#!/usr/bin/env Rscript
# ============================================================================
# fig_pes_imaging_tracking.R   [figure_id: fig_pes_imaging_tracking]
# ----------------------------------------------------------------------------
# Module 6 SUPPLEMENT. Does the PES track change between the TWO REPEAT IMAGING
# visits (instance 2 -> 3), a short ~2-year interval where BOTH timepoints are
# held-out follow-ups (no baseline anchor)? For each exposure with >=15 held-out
# changers between visits 2 and 3:
#   dcor(2->3) = cor( y_raw(3)-y_raw(2), pred_prot(3)-pred_prot(2) ).
# Plotted against the baseline-anchored within-person Delta-correlation used in
# main Fig. 3c (baseline -> each imaging visit, ~10 y). Exposures that track over
# the long interval also track imaging-to-imaging -- a stronger, fully held-out
# test of longitudinal tracking. Time gaps (median): 0->2 ~10 y, 2->3 ~2 y.
#   Rscript scripts/visualizations/build_figures.R --figure fig_pes_imaging_tracking
# ============================================================================
local({
  cand <- c(file.path(getwd(),"scripts","visualizations","common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cand[dir.exists(cand)][1]
  for (f in c("figure_paths.R","load_heap_results.R","plot_theme.R","label_helpers.R","export_helpers.R")) source(file.path(cm,f))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_pes_imaging_tracking")
od <- heap_project_output("module6_pes_longitudinal","base")

# baseline-anchored Delta-correlation (the panel c value)
wc <- fread(file.path(od, "PESlong_base_WithinDeltaCorCI.tsv"))[, .(exposure_id, dcor_base = dcor_prot, category)]

# imaging-to-imaging (2 -> 3) Delta-correlation, computed from the held-out scores
fs <- list.files(od, pattern = "_HoldoutScores\\.tsv$", full.names = TRUE)
one <- function(f) {
  x <- tryCatch(fread(f, select = c("eid","instance","exposure_id","exposure_type","y_raw","pred_prot")), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) return(NULL)
  v2 <- x[instance == 2, .(eid, y2 = y_raw, p2 = pred_prot)]
  v3 <- x[instance == 3, .(eid, y3 = y_raw, p3 = pred_prot)]
  d  <- merge(v2, v3, by = "eid"); d[, `:=`(dY = y3 - y2, dP = p3 - p2)]
  d  <- d[is.finite(dY) & is.finite(dP)]
  nch <- sum(round(d$dY) != 0)
  if (nrow(d) < 30 || sd(d$dY) == 0 || nch < 15) return(NULL)
  data.table(exposure_id = x$exposure_id[1], exposure_type = x$exposure_type[1],
             n_pairs_23 = nrow(d), n_change_23 = nch, dcor_23 = cor(d$dY, d$dP))
}
i23 <- rbindlist(lapply(fs, one), fill = TRUE)
m <- merge(wc, i23, by = "exposure_id")
m[, category := heap_category_factor(category)]
m[, exposure_label := heap_exposure_label(exposure_id)]
pct_pos <- round(100 * mean(m$dcor_23 > 0))

EX   <- c("current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days","alcohol_intake_frequency_f1558_0_0",
          "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0","usual_walking_pace_f924_0_0","processed_meat_intake_f1349_0_0")
ELAB <- c("Current smoking","Alcohol frequency","Vigorous activity","Walking pace","Processed meat")
exr <- m[exposure_id %in% EX]; exr[, elab := ELAB[match(exposure_id, EX)]]

p <- ggplot(m, aes(dcor_base, dcor_23)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey60", linewidth = 0.4) +
  geom_point(aes(colour = category), size = 1.9, alpha = 0.7) +
  geom_point(data = exr, shape = 21, fill = "white", colour = "grey15", size = 2.6, stroke = 1.0) +
  geom_text_repel(data = exr, aes(label = elab), colour = "grey15", size = 2.5, fontface = "bold",
                  box.padding = 0.5, min.segment.length = 0, segment.colour = "grey55", max.overlaps = 30, seed = 9) +
  scale_colour_exposure(drop = TRUE) +
  coord_cartesian(ylim = c(-0.2, 0.7)) +
  labs(title = "PES tracks change between the two imaging visits",
       subtitle = NULL,
       x = "baseline-anchored within-person tracking (~10 y)",
       y = "imaging-to-imaging tracking  (visit 2 -> 3, ~2 y)") +
  theme_heap(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 7.8, colour = "grey35"),
        axis.title = element_text(face = "bold"), legend.position = "right", legend.title = element_blank())

heap_emit_figure(p, figure_id,
                 data = m[, .(exposure_id, exposure_label, category, exposure_type, dcor_base, n_pairs_23, n_change_23, dcor_23)],
                 category = "supplement", formats = c("pdf","png"),
                 width = 8.2, height = 5.8, website = TRUE)
message("fig_pes_imaging_tracking: done (", nrow(m), " exposures; ", pct_pos, "% track positively 2->3).")
