#!/usr/bin/env Rscript

# ============================================================================
# fig_mr_coloc.R  [figure_id: fig_mr_coloc]
# ----------------------------------------------------------------------------
# CONSOLIDATED colocalization figure (2026-07-11). Merges:
#   fig_mr_coloc              the ASGR1 exemplar locuszoom (this file, before)
#   fig_mr_coloc_summary      the PP.H4 landscape / yield / colocalized-hit list
#   fig_mr_supplement panel a the PP.H4-vs-PP.H3 colocalization GATE
#
# All three were cited by one generic sentence and told one story between them:
# a cis-pQTL that merely sits NEAR a disease signal is not evidence of causation.
# coloc separates ONE shared causal variant (PP.H4) from TWO distinct variants in
# LD (PP.H3). Most candidates fail that test -- that is the point of the gate.
#
#   a  the gate       PP.H4 vs PP.H3 over every cis-pQTL x outcome test. The PP.H3
#                     corner is LD-confounded, NOT causal.
#   b  the yield      colocalized / tested by arm and edge type. This panel exists to
#                     make the SCOPE explicit, because the two source figures quoted
#                     different denominators without ever saying so:
#                        18/65 = ALL cis-pQTL tests (UKB P->D 18, UKB P->E 31, deCODE P->D 16)
#                        12/34 = cis P->D only, both arms  (the "12 shared vs 15 LD-confounded")
#                     Both were correct; neither stated its scope.
#   c  the survivors  the colocalized hits (PP.H4 >= 0.8), ranked
#   d  what one looks like   ASGR1 cis-pQTL vs the FinnGen lipoprotein-disorder locus,
#                     sharing lead SNP rs55714927 (PP.H4 = 0.998) -- the existence proof
#                     that panels a-c are summarising something real. The arm is NAMED
#                     on the facet strip: panel (c) distinguishes arms, so this one must
#                     too (see the PRV-2 note at panel d).
#
# Input : support/coloc/web/<ARM>__<PROT>__<TARGET>_plot_table.tsv (rsID-matched
#         windows; NOT the legacy coloc_index.tsv cache -- see panel d)
#         load_coloc_results() for the systematic landscape
# Output: figures/supplement/module5/fig_mr_coloc.{pdf,png} + data tsv
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
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_mr_coloc")
BS  <- 7.5
THR <- 0.8   # PP.H4 >= 0.8 = colocalized ; PP.H3 >= 0.8 = LD-confounded

# ------------------------------------------------- systematic landscape ------
co <- load_coloc_results()
co[, status3 := fcase(`PP.H4` >= THR, "Colocalized (shared variant)",
                      `PP.H3` >= THR, "LD-confounded (distinct variants)",
                      default        = "Inconclusive")]
co[, status3 := factor(status3, levels = c("Colocalized (shared variant)",
                                           "LD-confounded (distinct variants)",
                                           "Inconclusive"))]
EDGE <- c(Pcis_to_D = "cis-pQTL to disease", Pcis_to_E = "cis-pQTL to exposure")
co[, edge := factor(EDGE[edge_dir], levels = unname(EDGE))]

pd <- co[edge_dir == "Pcis_to_D"]
message(sprintf("coloc: %d/%d colocalized overall | cis P->D only: %d colocalized, %d LD-confounded of %d",
                co[`PP.H4` >= THR, .N], nrow(co),
                pd[`PP.H4` >= THR, .N], pd[`PP.H3` >= THR, .N], nrow(pd)))

PAL <- c(`Colocalized (shared variant)`      = "#1A6B30",
         `LD-confounded (distinct variants)` = "#D95F0E",
         `Inconclusive`                      = "#B0B0B0")

# ---- a: the gate ------------------------------------------------------------
lab_a <- co[`PP.H4` >= THR][order(-`PP.H4`)][seq_len(min(5, .N))]
pa <- ggplot(co, aes(`PP.H3`, `PP.H4`)) +
  geom_hline(yintercept = THR, linetype = "dashed", colour = "grey60", linewidth = .3) +
  geom_vline(xintercept = THR, linetype = "dashed", colour = "grey60", linewidth = .3) +
  geom_point(aes(colour = status3, shape = arm), size = 1.3, alpha = .9) +
  geom_text_repel(data = lab_a, aes(label = protID), size = 1.7, colour = "grey20",
                  segment.size = .18, segment.colour = "grey65", min.segment.length = 0,
                  box.padding = .35, max.overlaps = Inf, seed = 2) +
  scale_colour_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = c(UKB = 16, DECODE = 17), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .5, 1)) +
  scale_y_continuous(limits = c(0, 1.12), breaks = c(0, .5, 1)) +
  labs(x = "PP.H3  (distinct variants in LD)",
       y = "PP.H4\n(shared causal variant)") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom", legend.box = "vertical",
        legend.key.size = unit(5, "pt"),
        legend.text = element_text(size = BS - 2.5),
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2)) +
  guides(colour = guide_legend(nrow = 3, order = 1),
         shape  = guide_legend(nrow = 1, order = 2))

# ---- b: the yield, with scope made explicit ---------------------------------
yl <- co[, .(tested = .N, colocalized = sum(`PP.H4` >= THR)), by = .(arm, edge)]
yl[, notcol := tested - colocalized]
yl[, grp := paste0(arm, "\n", edge)]
yl[, lab := sprintf("%d / %d", colocalized, tested)]
ylb <- melt(yl, id.vars = "grp", measure.vars = c("colocalized", "notcol"),
            variable.name = "k", value.name = "n")
ylb[, k := factor(fifelse(k == "colocalized", "Colocalized", "Not colocalized"),
                  levels = c("Colocalized", "Not colocalized"))]

pb <- ggplot(ylb, aes(n, grp, fill = k)) +
  geom_col(width = .6) +
  geom_text(data = yl, aes(x = tested, y = grp, label = lab), inherit.aes = FALSE,
            hjust = -0.12, size = 1.9, fontface = "bold", colour = "grey15") +
  scale_fill_manual(values = c(Colocalized = "#1A6B30", `Not colocalized` = "#D5D8DC"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, .30))) +
  labs(x = "loci tested", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 2, lineheight = .9),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 1.5),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2))

# ---- c: the colocalized hits ------------------------------------------------
DZMAP <- c(E4_LIPOPROT = "Lipoprotein disorder", I9_AF = "Atrial fibrillation",
           T2D = "Type 2 diabetes", T2D_WIDE = "Type 2 diabetes (wide)",
           I9_HYPTENSESS = "Hypertension", M13_OSTEOPOROSIS = "Osteoporosis",
           M13_ARTHROSIS_KNEE = "Knee osteoarthritis", E4_OBESITY = "Obesity",
           E4_OBESITYCAL = "Obesity (caloric)", E4_OBESITYNAS = "Obesity (unspecified)")
pretty_dz <- function(x) {
  x0 <- sub("^finngen_R12_", "", x)
  m  <- DZMAP[x0]
  fb <- gsub("_", " ", tolower(sub("^[A-Z][0-9]+_", "", x0)))
  substr(fb, 1, 1) <- toupper(substr(fb, 1, 1))
  fifelse(is.na(m), fb, m)
}
hit <- co[`PP.H4` >= THR][order(`PP.H4`)]
hit[, tgt := fifelse(edge_dir == "Pcis_to_E",
                     heap_exposure_label(target), pretty_dz(target))]
hit[, tgt := heap_americanize(substr(tgt, 1, 30))]
hit[, lab := sprintf("%s -> %s (%s)", protID, tgt, arm)]
hit[, lab := factor(lab, levels = lab)]
pc <- ggplot(hit, aes(`PP.H4`, lab)) +
  geom_segment(aes(x = THR, xend = `PP.H4`, yend = lab),
               colour = "grey85", linewidth = .35) +
  geom_point(aes(shape = arm), size = 1.3, colour = "#1A6B30") +
  scale_shape_manual(values = c(UKB = 16, DECODE = 17), name = NULL) +
  scale_x_continuous(limits = c(THR, 1.004), breaks = c(.8, .9, 1)) +
  labs(x = "coloc PP.H4", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 3),
        legend.position = "none",
        plot.margin = margin(10, 4, 2, 2))

# ---- d: the exemplar locus --------------------------------------------------
# Read the ARM-KEYED export under support/coloc/web/, not the legacy cached
# <lead>_<prot>_<dz>_plot_table.tsv at the coloc root. Two reasons, both found
# 2026-08-29 (register PRV-2):
#
#  1. THE LEGACY CACHE IS THE WRONG ARM. coloc_index.tsv carries nsnps=4853 /
#     PP.H4=0.9981874, which is the DECODE row of coloc_results.tsv; the UKB row
#     is 4171 / 0.9984287. So this panel was drawing the deCODE SomaScan pQTL
#     while panel (c) labels arms explicitly and the main text discusses UKB
#     Olink. Every p_trait1 differed by orders of magnitude (lead 6e-55 vs
#     5.5e-94) for exactly that reason -- two platforms, not two runs.
#  2. THE LEGACY WINDOW IS TRUNCATED. LocusZoom.R filtered a GRCh37 LD panel
#     with GRCh38 bounds, so the window lost the variants beyond the build
#     offset: 2,508 variants over 6.77-7.68 Mb where the correct window holds
#     2,815 over 6.68-7.68 Mb. export_coloc_web.R matches on rsID instead.
#
# The LD colouring itself was never wrong -- every variant shared between the
# two tables has an identical r2 (max |delta| 5e-05), and the lead variant is
# rs55714927 in both. Only the arm and the left edge of the window change.
#
# ARM is settable; the panel names whichever it drew, so the two can never again
# be confused for one another.
coloc_dir <- heap_resolve_output(file.path("support", "coloc"), must_exist = FALSE)
ARM       <- Sys.getenv("HEAP_COLOC_ARM", unset = "UKB")
web_dir   <- file.path(coloc_dir, "web")

# Pick the strongest colocalized cis->disease locus IN THIS ARM that has an export.
cand <- co[arm == ARM & edge_dir == "Pcis_to_D" & `PP.H4` >= THR][order(-`PP.H4`)]
cand[, stem := sprintf("%s__%s__%s", arm, protID, target)]
cand <- cand[file.exists(file.path(web_dir, paste0(stem, "_plot_table.tsv")))]
if (!nrow(cand))
  stop("no colocalized ", ARM, " cis->disease locus with a web export under ", web_dir,
       "\nRun: Rscript scripts/support/coloc/export_coloc_web.R")
# Prefer ASGR1: it is the exemplar the legend and the main text describe, and it
# is a Tier-1 MEDIATOR. Ranking on PP.H4 alone hands the panel to PCSK9 (1.0000),
# which is a biomarker -- a different claim from the one this figure supports.
PREFER <- Sys.getenv("HEAP_COLOC_PROT", unset = "ASGR1")
sel    <- Sys.getenv("HEAP_COLOC_LOCUS",
                     unset = if (PREFER %in% cand$protID) cand[protID == PREFER][1]$stem
                             else cand[1]$stem)
# Resolve the row index OUTSIDE the data.table frame: inside `cand[...]` the
# columns are in scope, so `cand$stem == sel` with `sel` named `stem` would
# compare the column against itself and silently match every row (which is how
# an ASGR1 window first came back carrying PCSK9's protein name and PP.H4).
i_row  <- which(cand$stem == sel)[1L]
if (is.na(i_row))
  stop("locus not among the colocalized ", ARM, " exports: ", sel)
stem   <- sel
row    <- cand[i_row]

pt    <- fread(file.path(web_dir, paste0(stem, "_plot_table.tsv")))
mt_fp <- file.path(web_dir, paste0(stem, "_meta.tsv"))
mt    <- if (file.exists(mt_fp)) fread(mt_fp) else data.table()
lead  <- if (nrow(mt) && "lead" %in% names(mt)) mt$lead[1] else row$lead_snp
prot  <- row$protID; dzl <- row$target; pph4 <- row$`PP.H4`
ARMLAB <- c(UKB = "UK Biobank Olink", DECODE = "deCODE SomaScan")[ARM]
if (is.na(ARMLAB)) ARMLAB <- ARM
message(sprintf("coloc panel d: %s | arm %s (%s) | %d variants %.3f-%.3f Mb | lead %s | PP.H4 %.4f",
                stem, ARM, ARMLAB, nrow(pt), min(pt$pos)/1e6, max(pt$pos)/1e6, lead, pph4))

pt[, posMb := pos / 1e6]
pt[, ldb := cut(fifelse(is.na(r2), 0, r2), breaks = c(-Inf, .2, .4, .6, .8, Inf),
                labels = c("[0,0.2)", "[0.2,0.4)", "[0.4,0.6)", "[0.6,0.8)", "[0.8,1]"))]
LDPAL <- c(`[0,0.2)` = "#3B6EA5", `[0.2,0.4)` = "#7FBEE0", `[0.4,0.6)` = "#7FBF7F",
           `[0.6,0.8)` = "#F5A623", `[0.8,1]` = "#D62728")

T1 <- sprintf("%s cis-pQTL (%s)", prot, ARMLAB)
T2 <- sprintf("%s (FinnGen)", pretty_dz(dzl))
long <- rbindlist(list(
  pt[, .(posMb, np = -log10(pmax(p_trait1, 1e-300)), ldb, snp, trait = T1)],
  pt[, .(posMb, np = -log10(pmax(p_trait2, 1e-300)), ldb, snp, trait = T2)]))
long[, trait := factor(trait, levels = c(T1, T2))]
ldp <- long[snp == lead]

pdz <- ggplot(long, aes(posMb, np)) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
             colour = "grey60", linewidth = .25) +
  geom_point(aes(colour = ldb), size = .45, alpha = .85) +
  geom_point(data = ldp, shape = 23, size = 1.5, fill = "#6A0DAD",
             colour = "grey15", stroke = .3) +
  geom_text(data = ldp[trait == T1], aes(label = lead), vjust = -1.2, size = 1.7,
            fontface = "italic", colour = "grey20") +
  facet_wrap(~ trait, scales = "free_y", nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.16))) +
  scale_colour_manual(values = LDPAL, name = expression(LD~r^2)) +
  labs(x = sprintf("Chromosome %s (Mb)", pt$chr[1]),
       y = expression(-log[10](italic(p)))) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = BS - 1, face = "bold"),
        legend.position = "right", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2.5),
        legend.title = element_text(size = BS - 1.5),
        plot.margin = margin(10, 4, 2, 2))
if (is.finite(pph4))
  pdz <- pdz + geom_text(
    data = data.table(posMb = min(long$posMb), np = max(long[trait == T1]$np),
                      trait = factor(T1, levels = c(T1, T2))),
    aes(x = posMb, y = np, label = sprintf("coloc PP.H4 = %.3f", pph4)),
    hjust = 0, vjust = 2.4, size = 2.0, colour = "grey25", inherit.aes = FALSE)

p <- ((pa | pb) / pc / pdz) +
  plot_layout(heights = c(1.15, 1.0, 0.95)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- co[, .(arm, protID, target, edge_dir, lead_snp, `PP.H3`, `PP.H4`,
              status = as.character(status3))]
heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module5",
                 formats = c("pdf", "png"), width = 6.5, height = 8.0, website = TRUE)

message("fig_mr_coloc: done (exemplar locus ", stem, ").")
