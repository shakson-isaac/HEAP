#!/usr/bin/env Rscript
# ============================================================================
# build_module6_composite_landscape.R
# ----------------------------------------------------------------------------
# Module-6 landscape main figure, approximate 2x2. Each panel PNG is placed 1:1
# but FIRST auto-trimmed of its baked white border (ggplot plot.margin + title/
# axis padding ~0.1-0.15in/side) so the content fills its cell -- only the
# intentional gutters/margins remain. The grid layout is derived from the
# TRIMMED aspect ratios so a fills its top-left cell (a/b letters align), c/d
# align, right edges align, and the bottom row comes out taller than the top.
#   a = compact schematic   b = reads the exposome
#   c = tracks change        d = disease relevance (scatter + ladders)
# ============================================================================
suppressPackageStartupMessages({ library(png); library(grid) })
fd <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module6")
nm <- c(a="fig_module6_schematic_compact", b="fig_m6_panel_b_cell",
        c="fig_m6_panel_c_cell", d="fig_m6_panel_d_cell")

# crop the uniform white border off a panel image (keep a 2px hair of padding)
trim_white <- function(im, thr = 0.992, pad = 2L) {
  d <- dim(im)
  w <- if (length(d) == 3) (im[,,1] >= thr) & (im[,,2] >= thr) & (im[,,3] >= thr) else im >= thr
  rk <- which(rowSums(!w) > 0); ck <- which(colSums(!w) > 0)
  if (!length(rk) || !length(ck)) return(im)
  r0 <- max(1, min(rk) - pad); r1 <- min(d[1], max(rk) + pad)
  c0 <- max(1, min(ck) - pad); c1 <- min(d[2], max(ck) + pad)
  if (length(d) == 3) im[r0:r1, c0:c1, , drop = FALSE] else im[r0:r1, c0:c1, drop = FALSE]
}
img <- lapply(nm, function(n) trim_white(readPNG(file.path(fd, paste0(n, ".png")))))
ar  <- sapply(img, function(im) dim(im)[2] / dim(im)[1])      # trimmed width / height

# ---- layout (inches) derived from trimmed ARs ------------------------------
# top row fills width with a+b at height H1 (a fills -> aH==H1==bH, letters align);
# bottom fills with c+d at height H2. ar_c+ar_d < ar_a+ar_b -> H2 > H1 (bottom taller).
W  <- 9.7
mL<-0.17; mR<-0.12; mT<-0.12; mB<-0.10; gx<-0.20; gy<-0.16
usable <- W - mL - mR
H1 <- (usable - gx) / (ar["a"] + ar["b"])
H2 <- (usable - gx) / (ar["c"] + ar["d"])
aW<-ar["a"]*H1; bW<-ar["b"]*H1; cW<-ar["c"]*H2; dW<-ar["d"]*H2
Htot <- mT + H1 + gy + H2 + mB

topY0 <- Htot - mT - H1
cells <- list(
  a=list(im=img$a, x0=mL,       y0=topY0, w=aW, h=H1, L="a"),
  b=list(im=img$b, x0=mL+aW+gx, y0=topY0, w=bW, h=H1, L="b"),
  c=list(im=img$c, x0=mL,       y0=mB,    w=cW, h=H2, L="c"),
  d=list(im=img$d, x0=mL+cW+gx, y0=mB,    w=dW, h=H2, L="d"))

place <- function(p) {
  pushViewport(viewport(x=unit(p$x0+p$w/2,"in"), y=unit(p$y0+p$h/2,"in"),
                        width=unit(p$w,"in"), height=unit(p$h,"in")))
  grid.raster(p$im, interpolate=TRUE); popViewport()
  # letter sits in the gutter/margin just LEFT of the panel's top-left corner
  grid.text(p$L, x=unit(p$x0-0.02,"in"), y=unit(p$y0+p$h,"in"), just=c("right","top"),
            gp=gpar(fontsize=11, fontface="bold", fontfamily="sans", col="#111111"))
}
draw <- function() { grid.newpage(); grid.rect(gp=gpar(fill="white", col=NA)); invisible(lapply(cells, place)) }

out_png<-file.path(fd,"Module6_composite_landscape.png")
out_pdf<-file.path(fd,"Module6_composite_landscape.pdf")
png(out_png, width=W, height=Htot, units="in", res=400, type="cairo", bg="white"); draw(); invisible(dev.off())
pdf(out_pdf, width=W, height=Htot); draw(); invisible(dev.off())
cat(sprintf("layout %.2f x %.2f in (ar %.2f) ~ %.0f x %.0f mm | top row %.2f, bottom %.2f\n",
            W, Htot, W/Htot, W*25.4, Htot*25.4, H1, H2))
cat(sprintf("trimmed cells: a %.2fx%.2f  b %.2fx%.2f  c %.2fx%.2f  d %.2fx%.2f\n", aW,H1, bW,H1, cW,H2, dW,H2))
cat("wrote", out_png, "and", out_pdf, "\n")
