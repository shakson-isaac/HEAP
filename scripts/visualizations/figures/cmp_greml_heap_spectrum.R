#!/usr/bin/env Rscript
# ============================================================================
# cmp_greml_heap_spectrum.R -> E vs G scatter, HEAP vs GREML, ACCENTUATING
# exposure-responsive vs non-responsive proteins. Same 2,051 GREML-converged
# proteins. NO clipping (full range); a fixed flagship set (top exposomic by
# HEAP) is labeled in BOTH panels so you can see where each lands.
#   HEAP  = realized predictive R2 (unique drop-one PGS / PXS)
#   GREML = variance ceiling (V/Vp)
# Exposure-responsive = exposomic component >= 1% of variance (same in both).
# Output: figures/exploratory/module1/cmp_greml_heap_spectrum.png
# ============================================================================
local({
  cm <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers"))
    source(file.path(cm, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

pa <- file.path(heap_project_output("population_architecture"), "base", "grm_cutoff_0p025")
d <- fread(file.path(pa, "concordance_greml_vs_heap_r2.tsv")); d[, r2 := as.numeric(r2)]
gV <- dcast(d, protein ~ component, value.var = "greml"); setnames(gV, c("G","E","GxE"), c("gG","gE","gGxE"))
hV <- dcast(d, protein ~ component, value.var = "r2");    setnames(hV, c("G","E","GxE"), c("hG","hE","hGxE"))
m <- merge(gV, hV, by = "protein")[!is.na(gG) & !is.na(gE) & !is.na(hG) & !is.na(hE)]
m[, `:=`(hG = pmax(0,hG), hE = pmax(0,hE), gG = pmax(0,gG), gE = pmax(0,gE))]

THR <- 0.01
m[, resp_H := hE >= THR]; m[, resp_G := gE >= THR]
nH <- m[, sum(resp_H)]; nG <- m[, sum(resp_G)]
rho_E <- cor(m$hE, m$gE, method = "spearman", use = "complete.obs")
CELL <- nzchar(Sys.getenv("HEAP_CELL"))
BS <- if (CELL) 7 else 12
P0 <- if (CELL) 0.35 else 0.5; P1 <- if (CELL) 0.7 else 1.1; PF <- if (CELL) 1.0 else 1.5
LBLS <- if (CELL) 2.2 else 3.1; ANNS <- if (CELL) 2.5 else 3.6
FLAG <- m[order(-hE)][seq_len(if (CELL) 6L else 12L), protein]   # labeled in BOTH panels

COL_R <- HEAP_PAL_COMPONENT[["Exposome"]]; COL_N <- "grey80"
need_repel <- requireNamespace("ggrepel", quietly = TRUE)

scat <- function(dt, xv, yv, rv, ttl, sub, xlab, ylab, ncount) {
  fl <- dt[protein %in% FLAG]
  xr <- max(dt[[xv]], na.rm = TRUE); yr <- max(dt[[yv]], na.rm = TRUE)
  alab <- if (CELL) sprintf("%d exposure-\nresponsive", ncount)
          else sprintf("%d exposure-responsive\n(>= 1%% of variance)", ncount)
  p <- ggplot(dt, aes(.data[[xv]], .data[[yv]])) +
    geom_point(data = dt[get(rv) == FALSE], colour = COL_N, size = P0, alpha = 0.16, stroke = 0) +
    geom_point(data = dt[get(rv) == TRUE],  colour = COL_R, size = P1, alpha = 0.6, stroke = 0) +
    geom_hline(yintercept = THR, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    geom_point(data = fl, aes(.data[[xv]], .data[[yv]]), colour = "grey10", size = PF, stroke = 0) +
    annotate("text", x = xr * 0.98, y = yr * 0.30, hjust = 1, vjust = 1, colour = COL_R,
             fontface = "bold", size = ANNS, lineheight = 0.9, label = alab) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
    labs(title = ttl, subtitle = if (CELL) NULL else sub, x = xlab, y = ylab) +
    theme_heap(base_size = BS) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = if (CELL) rel(1) else rel(1.05)),
          plot.subtitle = element_text(hjust = 0.5, size = rel(0.78), colour = "grey35"),
          axis.title = element_text(face = "plain"), panel.grid.minor = element_blank(),
          panel.grid.major = if (CELL) element_blank() else element_line(linewidth = 0.25, colour = "grey92"))
  if (need_repel)
    p <- p + ggrepel::geom_text_repel(data = fl, aes(.data[[xv]], .data[[yv]], label = protein),
             size = LBLS, fontface = "bold", colour = "grey10", max.overlaps = Inf,
             min.segment.length = 0, box.padding = if (CELL) 0.4 else 0.55, point.padding = 0.2,
             segment.colour = "grey55", segment.size = 0.3, seed = 1)
  p
}
pA <- scat(m, "hG", "hE", "resp_H", "HEAP Prediction R2", "unique out-of-fold R2",
           expression("PGS predictive"~R^2), expression("PXS predictive"~R^2), nH)
pB <- scat(m, "gG", "gE", "resp_G", "GREML Variance Components", "variance component V/Vp",
           expression("Genetics: SNP-"*h^2), expression("Exposome:"~sigma[E]^2), nG)

# CELL mode: compact 2-scatter panel for the Fig 1 composite (panel b)
if (CELL) {
  OUT <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module1")
  VERT <- nzchar(Sys.getenv("HEAP_VERT"))
  if (VERT) {
    # portrait: HEAP over GREML, stacked — for the horizontal 3-column Fig 1
    pcell <- (pA / pB)
    ggsave(file.path(OUT, "cmp_greml_heap_spectrum_vert_cell.png"), pcell, width = 2.7, height = 4.5, dpi = 400, bg = "white")
    ggsave(file.path(OUT, "cmp_greml_heap_spectrum_vert_cell.pdf"), pcell, width = 2.7, height = 4.5, bg = "white")
    cat("wrote cmp_greml_heap_spectrum_vert_cell (portrait, green exposome)\n"); quit(save = "no")
  }
  pcell <- (pA | pB) + patchwork::plot_annotation(
    title = "Exposure-responsive proteins across the genetic spectrum (HEAP & GREML)",
    theme = theme(plot.title = element_text(face = "bold", size = 8.5, hjust = 0.5)))
  ggsave(file.path(OUT, "cmp_greml_heap_spectrum_cell.png"), pcell, width = 6.2, height = 2.35, dpi = 400, bg = "white")
  ggsave(file.path(OUT, "cmp_greml_heap_spectrum_cell.pdf"), pcell, width = 6.2, height = 2.35, bg = "white")
  cat("wrote cmp_greml_heap_spectrum_cell (green exposome)\n"); quit(save = "no")
}

# (c) cross-method agreement on the exposure axis (flagships labeled)
flc <- m[protein %in% FLAG]
pC <- ggplot(m, aes(hE, gE)) +
  geom_point(aes(colour = resp_H), size = 0.8, alpha = 0.4, stroke = 0) +
  scale_colour_manual(values = c("TRUE" = COL_R, "FALSE" = COL_N), guide = "none") +
  geom_smooth(method = "lm", se = FALSE, colour = "grey25", linewidth = 0.5) +
  geom_point(data = flc, colour = "grey10", size = 1.5, stroke = 0) +
  annotate("text", x = 0, y = max(m$gE), hjust = -0.1, vjust = 1, size = 3.6, colour = "grey20",
           label = sprintf("Spearman rho = %.2f", rho_E)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(title = "Both methods agree on which proteins are exposure-responsive",
       x = "HEAP exposomic R2 (realized)", y = "GREML exposomic V/Vp (ceiling)") +
  theme_heap(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.05)),
        panel.grid.minor = element_blank())
if (need_repel)
  pC <- pC + ggrepel::geom_text_repel(data = flc, aes(hE, gE, label = protein), size = 3.0,
           fontface = "bold", colour = "grey10", max.overlaps = Inf, min.segment.length = 0,
           box.padding = 0.5, segment.colour = "grey55", segment.size = 0.3, seed = 1)

p <- (pA | pB) / pC + plot_layout(heights = c(1, 0.92)) +
  plot_annotation(
    title = "Exposure-responsive vs non-responsive proteins (HEAP and GREML)",
    subtitle = "Exposomic component (y) vs genetic component (x, context only); top-12 exposure-responsive proteins labeled in both panels.",
    theme = theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
                  plot.subtitle = element_text(size = 10, colour = "grey35", hjust = 0.5)))
out <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module1/cmp_greml_heap_spectrum.png")
ggsave(out, p, width = 11, height = 10, dpi = 200, bg = "white")
cat(sprintf("HEAP exp-responsive %d | GREML exp-responsive %d | E Spearman %.2f | flagships: %s\nwrote %s\n",
            nH, nG, rho_E, paste(FLAG, collapse = ", "), out))
