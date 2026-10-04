#!/usr/bin/env Rscript

# ============================================================================
# fig_mr_motif_overview.R  [figure_id: fig_mr_motif_overview]
# ----------------------------------------------------------------------------
# MR triad motifs under two evidence bars (UKB pQTL arm).
#
# WHY THIS WAS REBUILT (2026-07-11). The MR section carried TWO different motif
# definitions and never said so:
#   main Fig 4b  (build_mr_panelb_folded.R) recomputes every motif with ALL SIX
#                edges evaluated at TIER 1  -> mediator = 7 triads / 3 proteins
#   MRmotifs.tsv (build_mr_tables.R:382)    defines motifs on nominal SIGNIFICANCE
#                (padj) -> mediator = 84 triads / 25 proteins
# The old supplement figure plotted the MRmotifs numbers while the citing sentence
# quoted the Fig 4b ones, so the reader could not check the claim against the
# figure. They differ ONLY in the NEGATED edges: a triad whose disease->protein
# edge is significant but sub-Tier-1 is DISQUALIFIED as a mediator by MRmotifs,
# but counts as a clean mediator in Fig 4b. That single choice is FURIN
# (pack-years smoking / TV time -> FURIN -> hypertension), and it is the whole
# 7-vs-5 gap.
#
# CANONICAL RULE (author's call, 2026-07-11): TIER 1 THROUGHOUT -- one evidence
# bar for all six edges, as Fig 4b already does. A merely-nominal reverse edge
# should not veto a Tier-1 causal chain. This figure now applies that rule and
# reproduces the main text EXACTLY: mediator 7 / 3, disease-liability 15,127 / 498.
#
# NB the two motif sets are NOT nested (FURIN is a mediator at Tier 1 but not at
# significance; motif C actually GROWS, 4,829 -> 4,905, because its negation
# loosens). So this is shown as two bars per motif, NOT a stack -- a stacked
# "share reaching Tier 1" would be mathematically wrong here.
#
# The result is the section's thesis, quantified: raising the bar to Tier 1
# destroys the mediator motif (84 -> 7) but barely touches disease-liability
# (17,999 -> 15,127). Causal mediation does not survive rigor; reporting does.
#
# Input : mr_edges/summary/MRmotifs.tsv        via load_mr_table()
#         mr_edges/summary/mr_tiered_edges.tsv (per-edge tiers)
# Output: figures/supplement/module5/fig_mr_motif_overview.{pdf,png} + data tsv
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

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_mr_motif_overview")
ARM <- Sys.getenv("HEAP_MR_ARM", unset = "UKB")
BS  <- 7.5

# ---------------------------------------------------------------- data ------
# Read the counts the summariser already computed rather than re-deriving the motif
# rule here. Three copies of that rule existed (summarize_mr_triads.R,
# build_mr_panelb_folded.R and this file); a fourth divergence is how Fig 4b and
# this panel could silently disagree.
cf <- file.path(heap_path(), "docs", "manuscript_stats", "module5", "mr_motif_counts.tsv")
if (!file.exists(cf))
  stop("missing ", cf, "\nRun: Rscript scripts/analysis_summaries/summarize_mr_triads.R")
cnt <- fread(cf)

NAME <- c(`A Mediator (E->P->D)`       = "A  Mediator (E->P->D)",
          `B Biomarker`                = "B  Biomarker",
          `C Exposure-marker`          = "C  Exposure-marker",
          `D Reverse (P->E)`           = "D  Reverse (P->E)",
          `E Disease-liability (D->P)` = "E  Disease-liability (D->P)")
cnt[, lab := NAME[motif]]

# ARM SCOPE: an edge touching the protein is evaluated within its pQTL platform, so a
# motif is Tier 1 if EITHER platform supports it and Tier 1+ if BOTH do. Splitting the
# Tier-1 bar that way shows what each platform contributes -- and shows that the
# mediator motif is UKB-only, because deCODE has no Tier-1 Pcis->E edges and too few
# Tier-1 E->P edges to close a triad.
# NB stacking is invalid on a log axis -- the segments would not sum to the bar, and
# the total label lands mid-bar. Each quantity therefore gets its own dodged bar, so
# every length reads correctly against the scale.
dt <- rbindlist(list(
  cnt[, .(lab, bar = "Any significant edge", seg = "any significant edge",
          n = nominal_triads, prot = nominal_proteins)],
  cnt[, .(lab, bar = "Tier 1", seg = "Tier 1 (either platform)",
          n = tier1_triads, prot = tier1_proteins)],
  cnt[, .(lab, bar = "Tier 1", seg = "UKB Olink",
          n = ukb_triads, prot = NA_integer_)],
  cnt[, .(lab, bar = "Tier 1", seg = "deCODE SomaScan",
          n = decode_triads, prot = NA_integer_)],
  cnt[, .(lab, bar = "Tier 1", seg = "both platforms (Tier 1+)",
          n = tier1plus_triads, prot = NA_integer_)]))
dt <- dt[n > 0]
SEGL <- c("any significant edge", "Tier 1 (either platform)", "UKB Olink",
          "deCODE SomaScan", "both platforms (Tier 1+)")
dt[, lab := factor(lab, levels = rev(unname(NAME)))]
dt[, bar := factor(bar, levels = c("Any significant edge", "Tier 1"))]
dt[, seg := factor(seg, levels = SEGL)]
dt[, mlab := fifelse(is.na(prot), comma(n), sprintf("%s  (%d prot.)", comma(n), prot))]

PAL <- c(`any significant edge` = "#BDC3C7", `Tier 1 (either platform)` = "#1A5276",
         `UKB Olink` = "#2E86C1", `deCODE SomaScan` = "#E8745C",
         `both platforms (Tier 1+)` = "#1A6B30")

p <- ggplot(dt, aes(n, lab, fill = seg)) +
  geom_col(position = position_dodge2(preserve = "single", reverse = TRUE), width = .78) +
  geom_text(aes(label = mlab),
            position = position_dodge2(width = .78, preserve = "single", reverse = TRUE),
            hjust = -0.08, size = 1.5, colour = "grey25") +
  facet_wrap(~ bar, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = PAL, name = NULL) +
  scale_x_log10(labels = comma, breaks = c(1, 100, 10000),
                expand = expansion(mult = c(0, .55))) +
  labs(x = "Exposure-protein-disease triads carrying the motif (log scale)", y = NULL) +
  theme_heap(base_size = BS) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = BS - 1, face = "bold"),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 1.5),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(8, 4, 2, 2))

heap_emit_figure(p, figure_id, data = dt, category = "supplement", subdir = "module5",
                 formats = c("pdf", "png"), width = 6.5, height = 5.0, website = TRUE)

message("fig_mr_motif_overview: done (both platforms).")
