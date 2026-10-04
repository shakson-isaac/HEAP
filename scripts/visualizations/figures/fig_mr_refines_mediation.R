#!/usr/bin/env Rscript

# ============================================================================
# fig_mr_refines_mediation.R  [figure_id: fig_mr_refines_mediation]
# ----------------------------------------------------------------------------
# PROMOTED (2026-07-11) out of fig_mr_supplement, where this was panels d and f
# of a six-panel landscape composite -- i.e. the MR section's thesis was buried
# as the fourth panel of a supplementary figure that gets scaled down in print.
#
# The claim: observational mediation nominates many proteins; MR supports almost
# none of them. Module 3 flags 1,466 of 2,686 tested proteins as observational
# exposome->protein->disease mediators. Requiring the protein->disease leg to be
# MR-causal (Tier 1 cis P->D) AND colocalized leaves 8 (either platform).
#
#   a  the funnel     2,686 tested -> 1,466 observational -> 8 causal+colocalized
#   b  who survives   the six, as modifiable exposure -> protein -> disease chains,
#                     with colocalization PP.H4, two-arm replication, and whether
#
# Panel b is the honest answer to "so what": the causal core is small, and half of
#
# Input : module3 primary_total (observational mediators, via heap_proportion_mediated)
#         mr_edges/summary/mr_tiered_edges.tsv + MRmotifs.tsv
#         coloc results via load_coloc_results()
# Output: figures/supplement/module5/fig_mr_refines_mediation.{pdf,png} + data tsv
#
# Authored at 6.5in = the supplement's \textwidth (scale 1.0).
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("cannot locate common/ helpers")
  for (f in c("figure_paths", "load_heap_results", "plot_theme",
              "label_helpers", "export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_mr_refines_mediation")
BS <- 7.5

# ------------------------------------------------- a: the attrition funnel ---
md  <- load_module3_results(covarType = "base", family = "lasso", mode = "primary_total")
pm  <- heap_proportion_mediated(md)
obs <- pm[predictor == "PXS_total" & NIE_q < 0.05 & pm_consistent == TRUE]
obs_prot <- unique(obs$protID)
n_tested <- uniqueN(pm[predictor == "PXS_total"]$protID)
n_obs    <- length(obs_prot)

sd <- heap_resolve_output(file.path("mr_edges", "summary"), must_exist = TRUE)
# ARM SCOPE: an edge touching the protein is evaluated within its pQTL platform;
# Tier 1 = supported on EITHER platform, Tier 1+ = on both (see macros/numbers.tex).
ed <- fread(file.path(sd, "mr_tiered_edges.tsv"))
mr_causal <- unique(ed[edge_dir == "Pcis_to_D" &
                       mr_tier_final %in% c("Tier1", "Tier1plus")]$src_id)
co         <- load_coloc_results()
coloc_prot <- unique(co[edge_dir == "Pcis_to_D" & colocalized == TRUE]$protID)
n_coloc    <- length(intersect(obs_prot, intersect(mr_causal, coloc_prot)))

message(sprintf("MR refines M3: tested=%d  observational=%d  MR-causal+colocalized=%d (%.1f%%)",
                n_tested, n_obs, n_coloc, 100 * n_coloc / n_tested))

fn <- data.table(
  stage = factor(c("Proteins tested\nfor mediation",
                   "Module-3 observational\nmediators",
                   "MR-causal & colocalized\nprotein → disease"),
                 levels = rev(c("Proteins tested\nfor mediation",
                                "Module-3 observational\nmediators",
                                "MR-causal & colocalized\nprotein → disease"))),
  n   = c(n_tested, n_obs, n_coloc),
  col = c("#BDC3C7", "#9E77B0", "#1A6B30"))
fn[, pct := fifelse(100 * n / n_tested < 1,
                    sprintf("%.1f%% of tested", 100 * n / n_tested),
                    sprintf("%.0f%% of tested", 100 * n / n_tested))]

pa <- ggplot(fn, aes(n, stage, fill = col)) +
  geom_col(width = .62) +
  geom_text(aes(label = comma(n)), hjust = -0.15, size = 2.4,
            fontface = "bold", colour = "grey15") +
  geom_text(aes(label = pct), hjust = -0.15, vjust = 2.6, size = 1.8, colour = "grey45") +
  scale_fill_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, .34))) +
  labs(x = "Proteins", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 1, lineheight = .95),
        plot.margin = margin(10, 4, 2, 2))

# ------------------------------------------- b: who actually survives --------
core_edges <- unique(ed[edge_dir == "Pcis_to_D" & mr_tier_final %in% c("Tier1", "Tier1plus"),
                        .(Protein = src_id, Disease = tgt_id)])
# ARM SCOPE: colocalization is a protein-anchored result, so take it from whichever
# platform supports it -- keeping this at arm=="UKB" made panel (b) show six proteins
# while panel (a) counted eight.
cc   <- co[edge_dir == "Pcis_to_D" & colocalized == TRUE,
           .(PPH4 = max(`PP.H4`)), by = .(Protein = protID, Disease = target)]
core <- merge(core_edges, cc, by = c("Protein", "Disease"))
core <- core[Protein %in% obs_prot]

# Three states, not two. A binary "replicates in deCODE" flag marked ICAM1 and SOST
# as replicated when they are colocalized in deCODE ONLY -- the opposite of what the
# reader would take from it.
armk <- function(a) co[arm == a & edge_dir == "Pcis_to_D" & colocalized == TRUE,
                       paste(protID, target)]
kU <- armk("UKB"); kD <- armk("DECODE")
core[, kk := paste(Protein, Disease)]
core[, arm_status := fcase(kk %in% kU & kk %in% kD, "both platforms",
                           kk %in% kU,              "UKB Olink only",
                           default =                "deCODE SomaScan only")]
core[, arm_status := factor(arm_status,
       levels = c("UKB Olink only", "deCODE SomaScan only", "both platforms"))]
core[, replicated := arm_status == "both platforms"]
core[, kk := NULL]

# The exposure leg comes from the TIER-1 MEDIATOR table -- the same rule the main
# text, Fig 4b, Fig 4e and S14 use -- not from the largest raw beta_EP.
#
# The previous pick_exp() took the exposure with the biggest |beta_EP| out of
# MRmotifs.tsv with no significance or tier requirement whatsoever (53 candidates
# for ASGR1 x lipoprotein disorder alone), then dropped "sexual|age_first_had" by
# name. It therefore printed "Health score -> ASGR1 -> Lipoprotein disorder" while
# Fig 4e printed "TV time -> ASGR1 -> Lipoprotein disorder" for the same pair --
# two shipped artifacts naming different exposures for one triad. This is the same
# failure NUM-3 recorded for Fig 4e; the ad-hoc name filter also excluded
# age_first_had_sexual_intercourse, which IS one of ASGR1's three genuine Tier-1
# exposures, while admitting 52 that are not.
#
# Only 3 of the 8 colocalized causal proteins have a Tier-1 E->P->D chain at all
# (ASGR1, ADM, FURIN). PCSK9, ALCAM, PRSS8, ICAM1 and SOST are Tier-1 and colocalized on the
# protein->disease leg but have no Tier-1 exposure->protein edge, so they get NO
# exposure rather than an invented one. That gap is the finding, not a blemish:
# a colocalized causal protein is not automatically a lifestyle mediator.
tri  <- fread(file.path(heap_path(), "docs", "manuscript_stats", "module5", "mr_triad_motifs.tsv"))
mcol <- names(tri)[grepl("motif", names(tri), ignore.case = TRUE)][1]
med  <- tri[get(mcol) == "A Mediator (E->P->D)", .(Protein, Exposure, Disease)]

# Where a protein-disease pair has several Tier-1 exposures, prefer the one Fig 4e
# draws, so the two figures cannot disagree.
FIG4E <- c(ASGR1 = "time_spent_watching_television_tv", ADM = "pack_years_of_smoking",
           FURIN = "time_spent_watching_television_tv")
pick_exp <- function(P, D) {
  cand <- med[Protein == P & Disease == D]
  if (!nrow(cand)) return(data.table(Exposure = NA_character_))
  pref <- if (P %in% names(FIG4E)) cand[grepl(FIG4E[[P]], Exposure, fixed = TRUE)] else cand[0]
  if (nrow(pref)) return(pref[1, .(Exposure)])
  cand[order(Exposure)][1, .(Exposure)]
}
core <- core[, cbind(.SD, pick_exp(Protein, Disease)), by = .(Protein, Disease)]

# Drug annotations removed 2026-08-29. They were a hardcoded map in this script, not
# derived from DrugBank or any other source, and no shipped figure or table in the
# paper uses druggability data -- the only live drug claim is the Discussion sentence
# on PCSK9 and ASGR1, which cites two trial papers directly.

pretty_dz <- function(x) {
  x0 <- sub("^finngen_R12_", "", x)
  m <- c(E4_LIPOPROT = "Lipoprotein disorder", I9_AF = "Atrial fibrillation",
         T2D = "Type 2 diabetes", T2D_WIDE = "Type 2 diabetes",
         I9_HYPTENSESS = "Hypertension", M13_OSTEOPOROSIS = "Osteoporosis",
         M13_ARTHROSIS_KNEE = "Knee osteoarthritis", E4_OBESITY = "Obesity",
         E4_OBESITYCAL = "Obesity (caloric)", E4_OBESITYNAS = "Obesity (unspecified)")[x0]
  fb <- gsub("_", " ", tolower(sub("^[A-Z][0-9]+_", "", x0)))
  substr(fb, 1, 1) <- toupper(substr(fb, 1, 1))
  fifelse(is.na(m), fb, m)
}
core[, dz  := pretty_dz(Disease)]
core[, exp_lab := heap_exposure_label(Exposure)]
core[, chain := fifelse(is.na(Exposure),
                        sprintf("%s  ->  %s", Protein, dz),
                        sprintf("%s  ->  %s  ->  %s", exp_lab, Protein, dz))]
core[, has_chain := !is.na(Exposure)]
stopifnot(!any(duplicated(core$chain)))   # a duplicate here means an un-filtered arm
setorder(core, -PPH4)
core[, chain := factor(chain, levels = rev(chain))]

message(sprintf("  causal core: %d chains | %d replicate in deCODE",
                nrow(core), sum(core$replicated)))

pb <- ggplot(core, aes(PPH4, chain)) +
  geom_segment(aes(x = 0.8, xend = PPH4, yend = chain), colour = "grey85", linewidth = .4) +
  geom_point(aes(shape = arm_status, fill = arm_status), size = 1.8, colour = "#1A6B30") +
  geom_text(aes(label = sprintf("%.3f", PPH4)), hjust = -0.35, size = 1.8, colour = "grey30") +
  scale_shape_manual(values = c(`UKB Olink only` = 1, `deCODE SomaScan only` = 24,
                                `both platforms` = 21), name = NULL,
                     drop = FALSE) +
  scale_fill_manual(values = c(`UKB Olink only` = NA, `deCODE SomaScan only` = "#1A6B30",
                               `both platforms` = "#1A6B30"), name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(0.8, 1.06), breaks = c(0.8, 0.9, 1.0),
                     expand = expansion(mult = c(.02, 0))) +
  labs(x = "Colocalization PP.H4 (cis-pQTL vs disease)", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 1.5),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2),
        plot.margin = margin(10, 4, 2, 2))

p <- (pa / pb) + plot_layout(heights = c(1, 1.45)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  fn[,   .(panel = "a", key = gsub("\n", " ", as.character(stage)), value = as.numeric(n))],
  core[, .(panel = "b", key = as.character(chain), value = PPH4)]), use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module5",
                 formats = c("pdf", "png"), width = 6.5, height = 4.6, website = TRUE)

message("fig_mr_refines_mediation: done.")
