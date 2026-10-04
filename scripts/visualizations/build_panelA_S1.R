#!/usr/bin/env Rscript
# Stack refined S1 (exposome Sankey) over the variance-partition block -> panel a.
suppressPackageStartupMessages({ library(png); library(grid) })
fd <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module1")
trim_white <- function(im, thr = 0.992, pad = 3L) {
  d <- dim(im)
  w <- if (length(d) == 3) (im[,,1] >= thr) & (im[,,2] >= thr) & (im[,,3] >= thr) else im >= thr
  rk <- which(rowSums(!w) > 0); ck <- which(colSums(!w) > 0)
  if (!length(rk) || !length(ck)) return(im)
  r0 <- max(1, min(rk)-pad); r1 <- min(d[1], max(rk)+pad); c0 <- max(1, min(ck)-pad); c1 <- min(d[2], max(ck)+pad)
  if (length(d) == 3) im[r0:r1, c0:c1, , drop = FALSE] else im[r0:r1, c0:c1, drop = FALSE]
}
top <- trim_white(readPNG(file.path(fd, "fig1_exp_S1.png")))
bot <- trim_white(readPNG(file.path(fd, "fig1_partition_block.png")))
art <- dim(top)[2] / dim(top)[1]; arb <- dim(bot)[2] / dim(bot)[1]
Wp <- 3.5                       # common panel width (in); fonts ~consistent (both ~6.5cm authored)
ht <- Wp / art; hb <- Wp / arb
gap <- 0.06
Htot <- ht + gap + hb
out <- file.path(fd, "fig1_panelA_S1.png")
png(out, width = Wp, height = Htot, units = "in", res = 400, type = "cairo", bg = "white")
grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
pushViewport(viewport(x = unit(Wp/2, "in"), y = unit(Htot - ht/2, "in"),
                      width = unit(Wp, "in"), height = unit(ht, "in")))
grid.raster(top, interpolate = TRUE); popViewport()
pushViewport(viewport(x = unit(Wp/2, "in"), y = unit(hb/2, "in"),
                      width = unit(Wp, "in"), height = unit(hb, "in")))
grid.raster(bot, interpolate = TRUE); popViewport()
invisible(dev.off())
cat(sprintf("panel a (S1 + partition) %.2f x %.2f in | AR %.2f -> %s\n", Wp, Htot, Wp/Htot, out))
