#!/usr/bin/env Rscript

# ============================================================================
# fig_mr_attrition.R  [figure_id: fig_mr_attrition]
# ----------------------------------------------------------------------------
# CONSOLIDATED Mendelian-randomization attrition. Merges the three figures that
# Supplementary Note 6 already narrated in a SINGLE sentence:
#   fig_mr_tier_funnel      tier-by-tier attrition per instrument class
#   fig_mr_hitrate_effsize  per-edge hit rate + evidence strength (|z|)
#   fig_mr_hit_retention    share of hits lost to heterogeneity / directional pleiotropy
#
# One claim: "most nominally significant edges fall away as the evidence filters
# are applied, and only a minority survive to the highest tiers and replicate
# across the two pQTL arms."
#   a  the funnel      44,996 pQTL-instrumented edges tested -> 2 replicated Tier 1+
#   b  where the hits are   hit rate and |z| by edge direction (UKB arm)
#   c  what kills them  % of hits retained after the sensitivity filters, both arms
#
# NB fig_mr_tier_funnel and fig_mr_hit_retention were ORPHANED: their only
# main-text home (results_m5_mr.tex:36) is commented out, so Note 6's references
# to them rendered as "Fig. ??". Folding them into this figure -- which the Fig 4
# legend cites -- repairs that.
#
# Input : mr_tables/ via load_mr_table() / load_mr_table_both()
#         ("tier_funnel", "sensitivity_by_edgedir", "mr_sensitivity_long")
# Output: figures/supplement/module5/fig_mr_attrition.{pdf,png} + data tsv
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

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_mr_attrition")
BS <- 7.5

PAL_ARM <- c(UKB = "#3FA7A7", DECODE = "#E8745C")
TEAL    <- "#3FA7A7"

EDGE_LEVELS <- c("E_to_P", "Pcis_to_D", "Ptrans_to_D", "E_to_D",
                 "Pcis_to_E", "Ptrans_to_E", "D_to_P", "D_to_E")
EDGE_PRETTY <- c(E_to_P = "E→P", Pcis_to_D = "Pcis→D", Ptrans_to_D = "Ptrans→D",
                 E_to_D = "E→D", Pcis_to_E = "Pcis→E", Ptrans_to_E = "Ptrans→E",
                 D_to_P = "D→P", D_to_E = "D→E")
fac <- function(x) factor(EDGE_PRETTY[as.character(x)], levels = rev(EDGE_PRETTY[EDGE_LEVELS]))

# ============================================================================
# Rebuilt 2026-08-29 (author review of Supp Fig 26). Three things were wrong:
#
#  1. Panel (a) collapsed the eight directed edge types into three LANES while
#     (b)-(d) resolved all eight, so the figure changed its unit halfway down.
#  2. Nothing showed the gate sequence that actually produces a tier. The old
#     "retained after sensitivity filters" panel used `sens_pass`, which is the
#     STRICT clean-only rule (het and Egger both unflagged) -- but the tiering
#     does not use it. build_mr_tables.R:245 gates on `robust_pass`, which is
#     `clean OR rescued_presso OR rescued_median`. Attrition was therefore shown
#     under a filter the pipeline never applies, making the pruning look far
#     harsher than it is: 2,728 hits are flagged, but 2,544 of them are RECOVERED.
#  3. The recovery itself was invisible, so a reader could not see that MR-PRESSO
#     and the weighted median are doing most of the work.
#
# The gate order below is read off build_mr_tables.R:272-292, not off the prose:
#   significant (BH q<0.05, within edge_dir)
#     -> instrument sufficiency  (nsnp >= 3, OR cis, which is exempt)
#     -> robustness              (clean OR rescued by MR-PRESSO / weighted median)
#     -> direction               (Steiger significant AND forward)
#     -> class tail              (trans -> Tier 2; cis -> Tier 1 pending coloc,
#                                 demoted to Tier 2 if not colocalized)
#     -> Tier 1+                 (replicated in the other pQTL arm)
# Tier 1 for a protein-exposure/disease edge is therefore a CIS statement: the
# trans directions top out at Tier 2 by construction, which panel (a) now shows
# as a hard stop rather than leaving the reader to infer it.
# ============================================================================

sl <- load_mr_table("mr_sensitivity_long", "UKB")
te <- load_mr_table("mr_tiered_edges", "UKB")[dataset == "UKB"]
K  <- c("edge_dir", "src_id", "tgt_id")
m  <- merge(sl[dataset == "UKB"], te[, ..K, with = FALSE][, .SD, .SDcols = K][
              , tmp := TRUE][, tmp := NULL], by = K, all.x = TRUE)
m  <- merge(sl[dataset == "UKB"], te[, c(K, "mr_tier_final"), with = FALSE],
            by = K, all.x = TRUE)
m  <- m[edge_dir %in% EDGE_LEVELS]

# the gates, in pipeline order
m[, instr_ok := mr_hit & (nsnp >= 3 | edge_class == "cis")]
m[, rob_ok   := instr_ok & robust_pass]
m[, dir_ok   := rob_ok & !direction_unresolved & steiger_ok]

# ---- a: the gate sequence as a flowchart -------------------------------------
# The review asked for the pruning to be DRAWN, not inferred from a bar chart.
# A CONSORT-style trunk with the rescue as an explicit excursion is the only
# layout in which "fails a sensitivity test, then gets recovered by MR-PRESSO"
# is visible as a path rather than as a difference between two bars.
NN <- c(tested = nrow(m), significant = sum(m$mr_hit, na.rm = TRUE),
        instrumented = sum(m$instr_ok, na.rm = TRUE), robust = sum(m$rob_ok, na.rm = TRUE),
        directional = sum(m$dir_ok, na.rm = TRUE),
        tier1 = sum(m$mr_tier_final %in% c("Tier1", "Tier1plus")),
        tier1plus = sum(m$mr_tier_final == "Tier1plus"))

STEP <- data.table(
  y = 7:1,
  lab = c("Directed edges tested",
          "Significant  (BH q < 0.05, within direction)",
          "Instrument-sufficient  (>= 3 SNPs, or cis)",
          "Robust  (clean, or rescued)",
          "Direction established  (Steiger sig. & forward)",
          "Tier 1 / Tier 1+",
          "Tier 1+  (replicated in deCODE)"),
  n = as.integer(NN))
STEP[, txt := sprintf("%s\nn = %s", lab, comma(n))]
STEP[, hi := y <= 2]

DROP <- data.table(
  y = c(6.5, 5.5, 4.5, 3.5, 2.5, 1.5),
  ndrop = as.integer(-diff(NN)),
  why = c("not significant", "too few instruments",
          "heterogeneous / pleiotropic\nand not rescued",
          sprintf("direction unresolved (%s)\nor reverse (%s)",
                  comma(m[rob_ok == TRUE & direction_unresolved == TRUE, .N]),
                  comma(m[rob_ok == TRUE & direction_unresolved == FALSE & steiger_ok == FALSE, .N])),
          "trans-only, or cis not\ncolocalized",
          "not replicated in the\nsecond pQTL arm"))
DROP[, txt := sprintf("- %s  %s", comma(ndrop), why)]

n_flagged <- m[instr_ok == TRUE & clean == FALSE, .N]
n_lost    <- m[instr_ok == TRUE & clean == FALSE & !rescued_presso & !rescued_median, .N]
BOXW <- 3.55; BOXH <- 0.42
# Aesthetics follow theme_heap so the flowchart reads as a panel of this figure
# rather than as pasted-in artwork: grey40 strokes at 0.4 (the panel border
# weight), grey95 box fill (the facet-strip fill used in panel b), the same teal
# and green as the bars, and a panel border so all four panels share one frame.
STROKE <- "grey40"; SW <- 0.4
pa <- ggplot() +
  geom_segment(data = STEP[y > 1], aes(x = 0, xend = 0, y = y - BOXH, yend = y - 1 + BOXH),
               arrow = arrow(length = unit(3.2, "pt"), type = "closed"),
               linewidth = SW, colour = STROKE) +
  geom_rect(data = STEP, aes(xmin = -BOXW/2, xmax = BOXW/2,
                             ymin = y - BOXH, ymax = y + BOXH, fill = hi),
            colour = STROKE, linewidth = SW) +
  geom_text(data = STEP, aes(0, y, label = txt), size = 1.9, colour = "grey15",
            lineheight = .95) +
  geom_segment(data = DROP, aes(x = 0, xend = BOXW/2 + 0.10, y = y, yend = y),
               linewidth = .3, colour = "grey65") +
  geom_text(data = DROP, aes(BOXW/2 + 0.18, y, label = txt), hjust = 0,
            size = 1.7, colour = "grey35", lineheight = .95) +
  # the rescue excursion: same amber the rescue bars in panel (c) use
  annotate("rect", xmin = -BOXW/2 - 2.30, xmax = -BOXW/2 - 0.22, ymin = 3.02, ymax = 5.22,
           fill = "#FCF4E6", colour = "#E8A33D", linewidth = SW) +
  annotate("text", x = -BOXW/2 - 1.26, y = 4.12, size = 1.62, colour = "grey15",
           lineheight = 1.22, vjust = 0.5,
           label = sprintf("flagged by Q or Egger: %s\n\nMR-PRESSO rescues %s\nweighted median rescues %s\neither: %s\n\nlost: %s",
                           comma(n_flagged),
                           comma(m[instr_ok == TRUE & clean == FALSE & rescued_presso == TRUE, .N]),
                           comma(m[instr_ok == TRUE & clean == FALSE & rescued_median == TRUE, .N]),
                           comma(m[instr_ok == TRUE & clean == FALSE & (rescued_presso | rescued_median), .N]),
                           comma(n_lost))) +
  annotate("segment", x = -BOXW/2 - 0.22, xend = -BOXW/2 - 0.02, y = 4.12, yend = 4.12,
           arrow = arrow(length = unit(3.2, "pt"), type = "closed"),
           linewidth = SW, colour = "#E8A33D") +
  scale_fill_manual(values = c(`FALSE` = "grey95", `TRUE` = "#DCEBE0"), guide = "none") +
  coord_cartesian(xlim = c(-BOXW/2 - 2.42, BOXW/2 + 2.30), ylim = c(0.45, 7.55),
                  clip = "off", expand = FALSE) +
  theme_void(base_size = BS) +
  theme(panel.border = element_rect(fill = NA, colour = STROKE, linewidth = SW),
        plot.margin = margin(8, 5, 4, 5))

# ---- b: the same sequence, per edge direction --------------------------------
STAGES <- c("Tested", "Significant", "Instrumented", "Robust", "Directional",
            "Tier 1/1+", "Tier 1+")
fn <- m[, .(Tested = .N,
            Significant  = sum(mr_hit,   na.rm = TRUE),
            Instrumented = sum(instr_ok, na.rm = TRUE),
            Robust       = sum(rob_ok,   na.rm = TRUE),
            Directional  = sum(dir_ok,   na.rm = TRUE),
            `Tier 1/1+` = sum(mr_tier_final %in% c("Tier1", "Tier1plus")),
            `Tier 1+`   = sum(mr_tier_final == "Tier1plus")), by = edge_dir]
fl <- melt(fn, id.vars = "edge_dir", measure.vars = STAGES,
           variable.name = "stage", value.name = "n")
fl[, edge_lab := factor(EDGE_PRETTY[as.character(edge_dir)], levels = EDGE_PRETTY[EDGE_LEVELS])]
fl[, stage := factor(stage, levels = STAGES)]
fl[, tier_stage := stage %in% c("Tier 1/1+", "Tier 1+")]

message(sprintf("MR gate sequence (UKB): %s tested -> %s significant -> %s instrumented -> %s robust -> %s directional -> %s Tier 1/1+ -> %s Tier 1+",
                comma(sum(fn$Tested)), comma(sum(fn$Significant)), comma(sum(fn$Instrumented)),
                comma(sum(fn$Robust)), comma(sum(fn$Directional)),
                comma(sum(fn[["Tier 1/1+"]])), comma(sum(fn[["Tier 1+"]]))))

pb <- ggplot(fl, aes(stage, pmax(n, 1), fill = tier_stage)) +
  geom_col(data = fl[n > 0], width = .74) +
  geom_text(aes(label = ifelse(n == 0, "0", comma(n))), vjust = -0.35,
            size = 1.35, colour = "grey25") +
  facet_wrap(~ edge_lab, nrow = 2) +
  scale_fill_manual(values = c(`FALSE` = TEAL, `TRUE` = "#1A6B30"), guide = "none") +
  scale_y_log10(expand = expansion(mult = c(0, .30)), labels = comma) +
  labs(x = NULL, y = "Directed edges (log)") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = BS - 3),
        strip.text = element_text(size = BS - 1.5, face = "bold"),
        plot.margin = margin(10, 4, 2, 2))

# ---- c: what the robustness gate does -- flag, then rescue ------------------
# The panel the review asked for: an edge that trips heterogeneity or the Egger
# intercept is NOT dropped. MR-PRESSO (outlier-corrected global test) and the
# weighted median are consulted, and either can carry it through.
rs <- m[instr_ok == TRUE, .(
          `Clean`                    = sum(clean, na.rm = TRUE),
          `Rescued: MR-PRESSO`       = sum(!clean &  rescued_presso & !rescued_median, na.rm = TRUE),
          `Rescued: weighted median` = sum(!clean & !rescued_presso &  rescued_median, na.rm = TRUE),
          `Rescued: both`            = sum(!clean &  rescued_presso &  rescued_median, na.rm = TRUE),
          `Lost (heterogeneous /\npleiotropic)` = sum(!clean & !rescued_presso & !rescued_median, na.rm = TRUE)),
        by = edge_dir]
RLEV <- c("Clean", "Rescued: MR-PRESSO", "Rescued: weighted median", "Rescued: both",
          "Lost (heterogeneous /\npleiotropic)")
rl <- melt(rs, id.vars = "edge_dir", variable.name = "outcome", value.name = "n")[n > 0]
rl[, outcome  := factor(outcome, levels = RLEV)]
rl[, edge_lab := factor(EDGE_PRETTY[as.character(edge_dir)], levels = rev(EDGE_PRETTY[EDGE_LEVELS]))]
RPAL <- c(`Clean` = "#B8D8D8", `Rescued: MR-PRESSO` = "#E8A33D",
          `Rescued: weighted median` = "#7B5EA7", `Rescued: both` = "#3FA7A7",
          `Lost (heterogeneous /\npleiotropic)` = "#C0392B")

n_flag <- rs[, sum(.SD), .SDcols = RLEV[-1]]
n_resc <- rs[, sum(.SD), .SDcols = RLEV[2:4]]
message(sprintf("robustness gate: %s flagged, %s rescued (%.0f%%), %s lost",
                comma(n_flag), comma(n_resc), 100 * n_resc / n_flag,
                comma(rs[[RLEV[5]]] |> sum())))

pc <- ggplot(rl, aes(n, edge_lab, fill = outcome)) +
  geom_col(width = .70, position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = RPAL, name = NULL) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, .02))) +
  labs(x = "Share of significant, instrument-sufficient edges", y = NULL) +
  theme_heap(base_size = BS) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(panel.grid = element_blank(),
        legend.position = "top", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2.5),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2))

# ---- d: why an edge did not reach Tier 1 ------------------------------------
RSN <- c(insufficient_instruments = "Too few instruments",
         heterogeneous_pleiotropic = "Heterogeneous / pleiotropic\n(not rescued)",
         direction_unresolved = "Direction unresolved\n(Steiger n.s.)",
         reverse_direction = "Reverse direction",
         trans_only = "Trans-only (Tier 2 ceiling)",
         cis_not_colocalized = "Cis, not colocalized")
lost <- m[mr_hit == TRUE & !(mr_tier_final %in% c("Tier1", "Tier1plus")) &
          tier_reason %in% names(RSN)]
lc <- lost[, .N, by = .(edge_dir, tier_reason)]
lc[, reason   := factor(RSN[tier_reason], levels = unname(RSN))]
lc[, edge_lab := factor(EDGE_PRETTY[as.character(edge_dir)], levels = rev(EDGE_PRETTY[EDGE_LEVELS]))]
CPAL <- setNames(c("#9E9E9E", "#C0392B", "#E8A33D", "#7B5EA7", "#3FA7A7", "#1F6F8B"), unname(RSN))

pd <- ggplot(lc, aes(N, edge_lab, fill = reason)) +
  geom_col(width = .70) +
  scale_fill_manual(values = CPAL, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, .04)), labels = comma) +
  labs(x = "Significant edges not reaching Tier 1", y = NULL) +
  theme_heap(base_size = BS) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(panel.grid = element_blank(),
        legend.position = "top", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2.5),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2))

# ------------------------------------------------------------- assemble ------
p <- (pa / pb / (pc | pd)) +
  plot_layout(heights = c(1.30, 1.02, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  fl[, .(panel = "a", key = paste("UKB", edge_dir, stage, sep = " | "), value = as.numeric(n))],
  rl[, .(panel = "b", key = paste("UKB", edge_dir, gsub("\n", " ", outcome), sep = " | "), value = as.numeric(n))],
  lc[, .(panel = "c", key = paste("UKB", edge_dir, tier_reason, sep = " | "), value = as.numeric(N))],
  STEP[, .(panel = "a", key = paste("UKB", lab, sep = " | "), value = as.numeric(n))]),
  use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module5",
                 formats = c("pdf", "png"), width = 6.5, height = 7.1, website = TRUE)

message("fig_mr_attrition: done.")
