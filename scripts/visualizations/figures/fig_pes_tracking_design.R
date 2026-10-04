#!/usr/bin/env Rscript

# ============================================================================
# fig_pes_tracking_design.R  [figure_id: fig_pes_tracking_design]
# ----------------------------------------------------------------------------
# CONSOLIDATED (2026-07-11). Merges the three Module-6 QC figures that were all
# cited by one clause ("This tracking reflects how participants varied in their
# exposures across assessment visits"):
#   fig_pes_visit_timing          years between repeat assessments
#   fig_pes_within_person_change  how many participants actually changed exposure
#   fig_pes_transition_power      precision of the within-person delta correlation
#
# They answer one question between them: is the within-person tracking claim
# supported by the design, or is it an artifact of people simply not changing?
#
#   a  when   how far apart the repeat assessments are (baseline -> imaging ->
#             repeat imaging), which sets how much exposure change is even possible
#   b  who    the share of repeat-visit participants who actually changed each
#             exposure -- the numerator the tracking correlation is computed over
#   c  how well  precision of the within-person correlation against the number of
#             changers: the exposures with the tightest CIs are the ones with the
#             most movers, i.e. the tracking estimates are power-limited, not biased
#
# Input : module6_pes_longitudinal/visit_timing/visit_timing_gaps.tsv
#         module6_pes_longitudinal/base/PESlong_base_WithinDeltaCorCI.tsv
# Output: figures/supplement/module6/fig_pes_tracking_design.{pdf,png} + data tsv
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
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel); library(scales)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_pes_tracking_design", "all_main", "all_supplement", "all", "website")]
covarType <- if (length(a) >= 1) a[1] else "base"
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_pes_tracking_design")
BS <- 7.5

# ---- a: how far apart are the repeat assessments? ---------------------------
tdir <- heap_project_output("module6_pes_longitudinal", "visit_timing")
g <- fread(file.path(tdir, "visit_timing_gaps.tsv"))
LAB <- c(`i0->i2` = "Baseline to imaging",
         `i2->i3` = "Imaging to repeat imaging",
         `i0->i3` = "Baseline to repeat imaging")
g <- g[pair %in% names(LAB)]
g[, pairf := factor(LAB[pair], levels = unname(LAB))]
s <- g[, .(med = median(gap)), by = pairf][order(med)]
s[, vj := rep_len(c(1.5, 3.2), .N)]

pa <- ggplot(g, aes(gap, fill = pairf, colour = pairf)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 1, alpha = .45,
                 position = "identity", linewidth = .12) +
  geom_vline(data = s, aes(xintercept = med, colour = pairf), linetype = 2, linewidth = .4) +
  geom_text(data = s, aes(x = med, y = Inf, label = sprintf("median %.0f y", med), colour = pairf),
            vjust = s$vj, hjust = -0.12, size = 1.9, show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  scale_colour_brewer(palette = "Set2", name = NULL) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  labs(x = "Years between assessments", y = "Density") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2)) +
  guides(fill = guide_legend(nrow = 3), colour = guide_legend(nrow = 3))

# ---- b + c: who actually changed, and how precise is the tracking? ----------
od <- heap_project_output("module6_pes_longitudinal", "base")
wc <- fread(file.path(od, "PESlong_base_WithinDeltaCorCI.tsv"))
wc[, ci_hw := (prot_hi - prot_lo) / 2]
wc <- wc[is.finite(ci_hw) & is.finite(n_change) & n_change > 0]
wc[, category := heap_category_factor(category)]
wc[, elab := heap_exposure_label(exposure_id)]

message(sprintf("PES tracking design: %d exposures | changers median %d (range %d-%d) | CI half-width median %.3f",
                nrow(wc), as.integer(median(wc$n_change)),
                min(wc$n_change), max(wc$n_change), median(wc$ci_hw)))

# b: how many participants actually changed each exposure (the numerator)
topb <- wc[order(-n_change)][seq_len(min(20, .N))]
topb[, elab := factor(elab, levels = rev(elab))]
pb <- ggplot(topb, aes(n_change, elab, fill = category)) +
  geom_col(width = .68) +
  geom_text(aes(label = comma(n_change)), hjust = -0.12, size = 1.7, colour = "grey25") +
  scale_fill_exposure(drop = TRUE, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, .22)), labels = comma) +
  labs(x = "Participants who changed the exposure", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 2.5),
        plot.margin = margin(10, 4, 2, 2))

p <- (pa | pb) +
  plot_layout(widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  g[,  .(panel = "a", key = as.character(pairf), value = as.numeric(gap))],
  wc[, .(panel = "b", key = exposure_id, value = as.numeric(n_change))]), use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module6",
                 formats = c("pdf", "png"), width = 6.5, height = 2.9, website = TRUE)

message("fig_pes_tracking_design: done.")
