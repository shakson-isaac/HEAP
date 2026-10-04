#!/usr/bin/env Rscript

# ============================================================================
# fig_dose_response.R  [figure_id: fig_dose_response]
# ----------------------------------------------------------------------------
# REBUILT 2026-08-29. The cited sentence is "effect sizes were dose-dependent.
# For instance, the mean absolute effect of alcohol on the proteome scaled
# monotonically with reported alcohol intake frequency" -- a claim about EFFECT
# SIZE growing with dose. The former panel (a) plotted binned protein LEVEL for
# six exemplar pairs, which is a different quantity, and the sentence asking for
# it (results_m2_association.tex:51) had already been commented out. That panel
# also drew LEP four times of six: its cast deduplicated on exposure but not on
# protein, and LEP is the top protein for four separate activity exposures.
#
# It is dropped. Alcohol is now the worked example rather than the whole figure,
# and the general claim is shown across every exposure that supports it.
#
#   a  the case the text names: mean |beta| over the replicated proteome rises
#      monotonically with reported alcohol intake frequency
#   b  the same relation for every other exposure with >= 4 replicated ordinal
#      levels, spanning exercise, diet and deprivation
#
# Dropping the old panel (a) also drops the ~9 GB loader dependency: everything
# here comes from the module-2 replicated term table.
#
# Input : module2 replicated terms (load_module2_replicated)
# Output: figures/supplement/module2/fig_dose_response.{pdf,png} + data tsv
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
  library(data.table); library(ggplot2); library(patchwork)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_dose_response", "all_main", "all_supplement", "all", "website")]
covarType  <- if (length(a) >= 1) a[1] else "base"
experiment <- if (length(a) >= 2) a[2] else "M2_base_main"
figure_id  <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_dose_response")
BS <- 7.5
MIN_LEVELS <- 4L          # an exposure needs this many replicated levels to show a trend

m   <- load_module2_replicated(covarType = covarType, experiment = experiment)
sig <- m[replicated == TRUE & is.finite(beta_train) & is.finite(se_train)]

# An ordinal exposure emits one term per response level, named <Eid><level>.
sig[, lvl := sub(paste0("^", Eid), "", ID), by = Eid]
ord <- sig[grepl("^[0-9]+$", lvl)][, lvl := as.integer(lvl)]

# mean |beta| across the replicated proteome at each level, with a pooled SE
# The interval is on the MEAN across proteins, so its width has to come from how
# much |beta| varies BETWEEN proteins, not from how precisely each beta was
# measured. The previous sqrt(sum(se^2))/n used the within-term standard errors and
# ran about 3.5x too narrow: at the top alcohol level the terms' own SEs give a
# half-width of 0.0015 while the spread of |beta| across the 585 proteins gives
# 0.0096. Proteins are also correlated, so even this is optimistic -- it is a lower
# bound on the uncertainty, not an exact one.
agg <- ord[, .(n = .N, mean_abs = mean(abs(beta_train)),
               se_mean = sd(abs(beta_train)) / sqrt(.N)), by = .(Eid, Category, lvl)]
agg[, n_levels := uniqueN(lvl), by = Eid]
agg <- agg[n_levels >= MIN_LEVELS]
setorder(agg, Eid, lvl)
agg[, rank := seq_len(.N), by = Eid]                 # 1..k in reported order
agg[, `:=`(lo = mean_abs - 1.96 * se_mean, hi = mean_abs + 1.96 * se_mean)]

message(sprintf("dose-response: %d exposures with >=%d replicated levels",
                uniqueN(agg$Eid), MIN_LEVELS))
for (e in unique(agg$Eid))
  message(sprintf("  %-52s levels=%d  terms=%5d", substr(e, 1, 52),
                  agg[Eid == e, .N], agg[Eid == e, sum(n)]))

# =============================== a: the alcohol case =========================
# Level names come from the generated codebook, not a map typed into this file.
# docs/manuscript_stats/exposure_level_labels.tsv is built by
# scripts/support/build_exposure_level_labels.R, which derives the suffix ->
# meaning mapping from the UKB codings and validates it against the alcohol
# labels this figure used to hardcode.
LBL_F <- file.path(heap_path(), "docs", "manuscript_stats", "exposure_level_labels.tsv")
if (!file.exists(LBL_F))
  stop("missing ", LBL_F, "\nRun: Rscript scripts/support/build_exposure_level_labels.R")
LBL <- fread(LBL_F)
agg[, fid := suppressWarnings(as.integer(sub("^.*_f([0-9]+)_.*$", "\\1", Eid)))]
agg <- merge(agg, LBL[, .(fid = field_id, lvl = suffix, level_label)],
             by = c("fid", "lvl"), all.x = TRUE)
agg[is.na(level_label), level_label := paste("level", lvl)]
setorder(agg, Eid, lvl)

FIELD <- "alcohol_intake_frequency_f1558_0_0"
alc <- agg[Eid == FIELD]
if (!nrow(alc)) stop("no replicated alcohol_intake_frequency level terms")
wrap2 <- function(x) gsub("(.{1,14})(\\s|$)", "\\1\n", x)   # keep tick labels narrow
alc[, label := factor(trimws(wrap2(level_label)), levels = trimws(wrap2(level_label)))]

pa <- ggplot(alc, aes(label, mean_abs)) +
  geom_line(aes(group = 1), colour = "grey70", linewidth = .4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .16, linewidth = .35, colour = "grey45") +
  geom_point(aes(colour = lvl), size = 2.1) +
  scale_colour_gradient(low = "#FDD0A2", high = "#A63603", guide = "none") +
  scale_y_continuous(limits = c(0, max(alc$hi) * 1.12), expand = expansion(mult = c(0, .02))) +
  labs(title = "Alcohol intake frequency", x = NULL,
       y = expression("Mean absolute effect size  "*group("|",beta,"|"))) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(size = BS, face = "bold", hjust = 0.5),
        axis.text.x = element_text(size = BS - 2, lineheight = .85),
        plot.margin = margin(10, 4, 2, 2))

# ================= b: the same relation for every other exposure =============
oth <- agg[Eid != FIELD]
oth[, facet := heap_exposure_label(Eid)]
setorder(oth, Category, facet, rank)
oth[, facet := factor(facet, levels = unique(facet))]
lvl_order <- unique(oth[order(rank)]$level_label)
oth[, xlab := factor(level_label, levels = lvl_order)]

pb <- ggplot(oth, aes(xlab, mean_abs, colour = Category)) +
  geom_line(aes(group = facet), linewidth = .4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .12, linewidth = .3, alpha = .7) +
  geom_point(size = 1.1) +
  facet_wrap(~ facet, nrow = 2, scales = "free") +
  scale_colour_exposure(drop = TRUE, name = NULL) +
  scale_x_discrete(drop = TRUE) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL,
       y = expression("Mean  "*group("|",beta,"|"))) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = BS - 2, lineheight = .9),
        axis.text.x = element_text(size = BS - 3.2, angle = 40, hjust = 1),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2))

p <- (pa / pb) + plot_layout(heights = c(0.72, 1.6)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- agg[, .(panel = fifelse(Eid == FIELD, "a", "b"),
               key = paste(Eid, lvl, sep = " | "), x = as.numeric(rank), value = mean_abs)]

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module2",
                 formats = c("pdf", "png"), width = 6.5, height = 6.2, website = TRUE)

message("fig_dose_response: done.")
