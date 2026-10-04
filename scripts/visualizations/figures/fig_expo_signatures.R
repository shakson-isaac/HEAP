#!/usr/bin/env Rscript
# ============================================================================
# fig_expo_signatures.R  [figure_id: fig_expo_signatures]
# ----------------------------------------------------------------------------
# Per-exposure SIGNATURE strip: ALL 13 exposure categories as COLUMNS (x), each
# column's top proteins plotted as jittered dots on a SHARED R2 axis (y, honest
# magnitude). Columns ordered by reach (high -> low, left -> right); the
# environmental/socioeconomic columns are (near-)empty, showing only lifestyle
# exposures shape the proteome. Only the top-reach columns get protein labels
# (keeps the crowded left readable). Combines the full exposome (13 columns) +
# which proteins (labeled dots) + magnitude (shared y) + multi-contribution
# (a protein recurs across columns).
#
# Input : module1_predictive_r2_score_partition/.../predictive_r2_exposure_categories_*
# Output: figures/main/module1/fig_expo_signatures.{pdf,png} (+ CELL cell.{png,pdf})
# ============================================================================
local({
  cm <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers"))
    source(file.path(cm, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(1)

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_expo_signatures","all_main","all_supplement","all","website")]
covarType <- if (length(a) >= 1) a[1] else "base"
method    <- if (length(a) >= 2) a[2] else "lasso"
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_expo_signatures")
N_TOP    <- as.integer(Sys.getenv("HEAP_N_TOP", "8"))      # dots per exposure column
N_LAB    <- as.integer(Sys.getenv("HEAP_N_LAB", "3"))      # labeled proteins per labeled column
N_LABCAT <- as.integer(Sys.getenv("HEAP_N_LABCAT", "4"))   # only top-reach columns get labels
THR   <- 0.005
CELL <- nzchar(Sys.getenv("HEAP_CELL"))
BS <- if (CELL) 7 else 11; TTL <- if (CELL) 8.5 else 13
PSZ <- if (CELL) 1.2 else 1.9; LBL <- if (CELL) 1.95 else 2.8
B_W <- 3.85; B_H <- 4.15

ec <- load_module1_predictive_r2(covarType = covarType, method = method, level = "exposure_categories")
if ("method" %in% names(ec)) ec <- ec[get("method") == "score_unique_drop"]
pp <- ec[, .(r2 = mean(r2)), by = .(omic, category)]
catsig <- pp[, .(n = sum(r2 > THR)), by = category][order(-n)]
cats <- catsig$category                                          # all 13, by reach
top <- pp[r2 > THR][order(category, -r2)][, head(.SD, N_TOP), by = category]
top[, rk := seq_len(.N), by = category]
catn <- setNames(catsig$n, as.character(catsig$category))
xlev <- as.character(cats)                                       # highest reach at LEFT
top[, category := factor(as.character(category), levels = xlev)]
labcats <- xlev[seq_len(min(N_LABCAT, length(xlev)))]
lab <- top[rk <= N_LAB & as.character(category) %in% labcats]
xlabs <- setNames(sprintf("%s (%d)", heap_category_pretty(cats), catn[as.character(cats)]), xlev)

# annotate the (near-)empty environmental/socioeconomic tail so the white space carries the message
empty_cats <- intersect(xlev, names(catn)[catn == 0])
ymax <- max(top$r2, na.rm = TRUE)
ann <- NULL
if (length(empty_cats) >= 2) {
  xi <- match(empty_cats, xlev); x0 <- min(xi); x1 <- max(xi)
  # The bracket and its "no exposure-responsive proteins" note were removed
  # 2026-08-30: an empty category already reads as empty, and the threshold now
  # lives in the caption rather than on the panel.
  ann <- list()
}

jit <- position_jitter(width = 0.17, height = 0, seed = 1)
p <- ggplot(top, aes(category, r2, colour = category)) +
  geom_point(position = jit, size = PSZ, alpha = 0.85, stroke = 0) +
  ann +
  scale_colour_exposure(guide = "none") +
  scale_x_discrete(limits = xlev, labels = xlabs, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(title = if (CELL) "Top protein signature per exposure" else "Each exposure's top protein signatures",
       subtitle = if (CELL) NULL else sprintf("Module 1 (%s/%s): top %d proteins per exposure (shared R2); (n) = proteins with R2>0.5%%",
                          covarType, method, N_TOP),
       x = NULL, y = expression("predictive"~R^2),
       caption = if (CELL) NULL else "all 13 exposures shown; column (n) = proteins it shapes (R2 > 0.5%); top exposures' proteins labeled") +
  theme_heap(base_size = BS)

if (requireNamespace("ggrepel", quietly = TRUE))
  p <- p + ggrepel::geom_text_repel(data = lab, aes(label = omic),
           position = jit, size = LBL, fontface = "bold", colour = "grey15", max.overlaps = Inf,
           min.segment.length = 0, box.padding = 0.32, point.padding = 0.1, force = 2,
           segment.colour = "grey70", segment.size = 0.22, direction = "both", seed = 1)

p <- if (CELL) {
  p + theme(plot.title = element_text(size = TTL, hjust = 0.5, face = "bold"),
            plot.title.position = "panel", axis.title = element_text(face = "plain"),
            axis.text.x = element_text(size = rel(0.78), face = "bold", angle = 35, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = rel(0.85)),
            plot.caption = element_text(size = rel(0.58), colour = "grey45", hjust = 0),
            panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(), plot.margin = margin(3, 5, 2, 3))
} else {
  p + theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, face = "bold"),
            panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(linewidth = 0.3),
            panel.grid.minor = element_blank())
}

if (CELL) {
  FIGDIR_CELL <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module1")
  dir.create(FIGDIR_CELL, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(FIGDIR_CELL, paste0(figure_id, "_cell.png")), p, width = B_W, height = B_H, dpi = 400, bg = "white")
  ggsave(file.path(FIGDIR_CELL, paste0(figure_id, "_cell.pdf")), p, width = B_W, height = B_H, bg = "white")
} else {
  heap_emit_figure(p, figure_id, data = top, category = "main", formats = c("pdf","png"), width = 7.5, height = 4, website = TRUE)
}
message("fig_expo_signatures: done (13 exposure columns; reach ", paste(sprintf("%s=%d", cats, catn[as.character(cats)]), collapse=", "), ").")
