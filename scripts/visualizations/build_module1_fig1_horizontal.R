#!/usr/bin/env Rscript
# ============================================================================
# build_module1_fig1_horizontal.R -> horizontal 3-COLUMN Fig 1 (exposure-first).
#   a  portrait variance-partition + exposome schematic
#   b  stacked HEAP / GREML exposure-responsive scatters
#   c  per-exposure signature strip (all 13 exposures)
# Panels are trimmed, scaled to a COMMON HEIGHT and laid side-by-side; column
# widths follow each panel's native aspect ratio, so the landscape strip (c) is
# the widest column (room for its labels) and the portrait schematic/scatters
# (a,b) are narrow. Set HEAP_FIGH to change the common panel height (inches).
# Reads *_cell/*.png; writes figures/exploratory/module1/fig_module1_fig1_horizontal.{png,pdf}
# ============================================================================
suppressPackageStartupMessages({ library(png); library(grid) })
fd <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module1")
src <- c(a = Sys.getenv("HEAP_PANELA", "fig1_panelA_S1.png"),
         b = "cmp_greml_heap_spectrum_vert_cell.png",
         c = "fig_expo_signatures_cell.png")
labmap <- c(a = "a", b = "b", c = "c")

trim_white <- function(im, thr = 0.992, pad = 2L) {
  d <- dim(im)
  w <- if (length(d) == 3) (im[,,1] >= thr) & (im[,,2] >= thr) & (im[,,3] >= thr) else im >= thr
  rk <- which(rowSums(!w) > 0); ck <- which(colSums(!w) > 0)
  if (!length(rk) || !length(ck)) return(im)
  r0 <- max(1, min(rk)-pad); r1 <- min(d[1], max(rk)+pad); c0 <- max(1, min(ck)-pad); c1 <- min(d[2], max(ck)+pad)
  if (length(d) == 3) im[r0:r1, c0:c1, , drop = FALSE] else im[r0:r1, c0:c1, drop = FALSE]
}
img <- lapply(src, function(f) trim_white(readPNG(file.path(fd, f))))
ar  <- sapply(img, function(im) dim(im)[2] / dim(im)[1])   # width / height

# Target a fixed TOTAL width (Nature double-column ~180 mm = 7.087 in) and derive
# the common panel height from the aspect ratios (override with HEAP_FIGW).
Wt <- as.numeric(Sys.getenv("HEAP_FIGW", "7.087"))         # total figure width (in)
mL <- 0.10; mR <- 0.06; mT <- 0.18; mB <- 0.06; gx <- 0.16
Hp   <- (Wt - mL - mR - gx*(length(ar)-1)) / sum(ar)       # common panel height (in)
colW <- ar * Hp
W    <- Wt
Htot <- mT + Hp + mB

cells <- list(); xleft <- mL
for (k in names(src)) {
  cells[[length(cells)+1]] <- list(im = img[[k]], x0 = xleft, y0 = mB, w = colW[k], h = Hp, L = labmap[[k]])
  xleft <- xleft + colW[k] + gx
}
place <- function(p) {
  pushViewport(viewport(x = unit(p$x0 + p$w/2, "in"), y = unit(p$y0 + p$h/2, "in"),
                        width = unit(p$w, "in"), height = unit(p$h, "in")))
  grid.raster(p$im, interpolate = TRUE); popViewport()
  grid.text(p$L, x = unit(p$x0 - 0.01, "in"), y = unit(p$y0 + p$h + 0.02, "in"), just = c("left","bottom"),
            gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "sans", col = "#111111"))
}
draw <- function() { grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA)); invisible(lapply(cells, place)) }
out_png <- file.path(fd, "fig_module1_fig1_horizontal.png")
png(out_png, width = W, height = Htot, units = "in", res = 360, type = "cairo", bg = "white"); draw(); invisible(dev.off())
pdf(file.path(fd, "fig_module1_fig1_horizontal.pdf"), width = W, height = Htot); draw(); invisible(dev.off())
cat(sprintf("horizontal Fig1 %.2f x %.2f in ~ %.0f x %.0f mm | colW(in): a=%.2f b=%.2f c=%.2f | AR: %.2f/%.2f/%.2f\n",
            W, Htot, W*25.4, Htot*25.4, colW["a"], colW["b"], colW["c"], ar["a"], ar["b"], ar["c"]))
cat("wrote", out_png, "\n")
