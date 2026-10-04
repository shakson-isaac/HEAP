#!/usr/bin/env Rscript
# ============================================================================
# build_module4_composite_v4.R  (Fig5 interventions -- restructured)
# ----------------------------------------------------------------------------
# Reading order: schematic -> heatmap -> scatter -> network.
#   a  triangulation schematic (convergent proteins)  -- fig_m4_schematic_compact
#   b  interventional concordance heatmap              -- fig_m4_panel_b_cell
#   c  MR-causal proteins scatter                      -- fig_m4_panel_d_cell
#   d  shared-language network (full width, bottom)    -- fig_m4_shared_network_cell
# Top row = a + b + c (common height, schematic scaled down); bottom = d full
# width. Each PNG white-trimmed + placed 1:1; letters at each image top-left.
# Render the cells first (HEAP_CELL=1: panel_b, panel_d, shared_network) +
# tectonic fig_m4_schematic_compact.tex.
# ============================================================================
suppressPackageStartupMessages({ library(png); library(grid) })
fd <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module4")
nm <- c(a = "fig_m4_schematic_compact", b = "fig_m4_panel_b_cell",
        c = "fig_m4_panel_d_cell",      d = "fig_m4_shared_network_cell")

trim_white <- function(im, thr = 0.992, pad = 2L) {
  d <- dim(im)
  w <- if (length(d) == 3) (im[,,1] >= thr) & (im[,,2] >= thr) & (im[,,3] >= thr) else im >= thr
  rk <- which(rowSums(!w) > 0); ck <- which(colSums(!w) > 0)
  if (!length(rk) || !length(ck)) return(im)
  r0 <- max(1, min(rk) - pad); r1 <- min(d[1], max(rk) + pad)
  c0 <- max(1, min(ck) - pad); c1 <- min(d[2], max(ck) + pad)
  if (length(d) == 3) im[r0:r1, c0:c1, , drop = FALSE] else im[r0:r1, c0:c1, drop = FALSE]
}
img <- lapply(nm, function(n) {
  p <- file.path(fd, paste0(n, ".png"))
  if (!file.exists(p)) stop("missing panel PNG: ", p, " (render with HEAP_CELL=1)")
  trim_white(readPNG(p))
})
AR <- function(im) dim(im)[2] / dim(im)[1]

# ---- geometry (inches) -----------------------------------------------------
W <- 9.7; mL <- 0.15; mR <- 0.12; mT <- 0.12; mB <- 0.10; gx <- 0.30; gy <- 0.38
usable <- W - mL - mR
arA <- AR(img$a); arB <- AR(img$b); arC <- AR(img$c); arD <- AR(img$d)
aSC <- 0.72      # square schematic: its height and width move together, so its
                 # share is capped; the heatmap gave back width to pay for it
# top row (a scaled + b + c) fills the width at a common reference height H1
H1 <- (usable - 2 * gx) / (arA * aSC + arB + arC)
aH <- aSC * H1; aW <- arA * aH
bW <- arB * H1; cW <- arC * H1
# bottom row: network, full width
nW <- usable;  nH <- nW / arD
Htot <- mT + H1 + gy + nH + mB
top  <- Htot - mT
y1   <- top - H1                                # top-row bottom
y2   <- mB                                      # network bottom

place_one <- function(im, L, x, y0, w, h) {
  pushViewport(viewport(x = unit(x + w/2, "in"), y = unit(y0 + h/2, "in"),
                        width = unit(w, "in"), height = unit(h, "in")))
  grid.raster(im, interpolate = TRUE); popViewport()
  grid.text(L, x = unit(x - 0.02, "in"), y = unit(y0 + h, "in"), just = c("right", "top"),
            gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "sans", col = "#111111"))
}
draw <- function() {
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  # top row: a schematic (top-aligned, shorter), b heatmap, c scatter
  place_one(img$a, "a", mL,               top - aH, aW, aH)
  place_one(img$b, "b", mL + aW + gx,     y1,       bW, H1)
  place_one(img$c, "c", mL + aW + gx + bW + gx, y1, cW, H1)
  # bottom: network full width
  place_one(img$d, "d", mL,               y2,       nW, nH)
}

out_png <- file.path(fd, "Fig5_intervention_composite_v4.png")
out_pdf <- file.path(fd, "Fig5_intervention_composite_v4.pdf")
png(out_png, width = W, height = Htot, units = "in", res = 400, type = "cairo", bg = "white"); draw(); invisible(dev.off())
pdf(out_pdf, width = W, height = Htot); draw(); invisible(dev.off())
cat(sprintf("layout %.2f x %.2f in | top H1=%.2f (a %.2fx%.2f b %.2fx%.2f c %.2fx%.2f) | net %.2fx%.2f\n",
            W, Htot, H1, aW, aH, bW, H1, cW, H1, nW, nH))
cat("wrote", out_png, "\n")
