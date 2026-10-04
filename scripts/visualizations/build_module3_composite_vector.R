#!/usr/bin/env Rscript
# ============================================================================
# build_module3_composite_vector.R  -> editable vector Fig 3
# Lays the Fig-3 panel VECTOR PDFs onto one canvas with the SAME scale-aware layout
# as the raster composite (fig3): row1 [a schematic (small) | b scale grid],
# row2 [c pleiotropy | d forest]. White border auto-trimmed from each panel PNG and
# applied to the PDF via \includegraphics[trim,clip]. Compile with tectonic in the dir.
# ============================================================================
suppressPackageStartupMessages(library(png))
fd <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3")
rows <- list(
  list(L=c("a","b"), stem=c("fig_module3_schematic_compact","fig_mediation_scale_main_cell"),
       dpi=c(300,400), s=c(0.72,1)),   # schematic enlarged 0.55 -> 0.72 (author, 2026-08-30)
  list(L=c("c","d"), stem=c("fig_mediation_pleiotropy_cell","fig_mediation_forest_cell"),
       dpi=c(400,400), s=c(1,1)))
panel_info <- function(stem, dpi) {
  im <- readPNG(file.path(fd, paste0(stem,".png"))); d <- dim(im)
  w <- if (length(d)==3) (im[,,1]>=.992)&(im[,,2]>=.992)&(im[,,3]>=.992) else im>=.992
  rk <- which(rowSums(!w)>0); ck <- which(colSums(!w)>0); H<-d[1]; W<-d[2]
  tl<-(min(ck)-1)/W; tr<-(W-max(ck))/W; tt<-(min(rk)-1)/H; tb<-(H-max(rk))/H
  pW<-W/dpi; pH<-H/dpi
  list(stem=stem, pageW=pW, pageH=pH, tl=tl, tr=tr, tt=tt, tb=tb, cW=pW*(1-tl-tr), cH=pH*(1-tt-tb))
}
info <- lapply(rows, function(r) Map(panel_info, r$stem, r$dpi))
ars  <- lapply(info, function(ri) sapply(ri, function(i) i$cW/i$cH))
scl  <- lapply(rows, function(r) r$s)

W<-9.7; mL<-0.16; mR<-0.10; mT<-0.12; mB<-0.10; gx<-0.24; gy<-0.30; HMAX<-4.3
usable <- W - mL - mR
Hrow <- sapply(seq_along(rows), function(i) min((usable-(length(ars[[i]])-1)*gx)/sum(ars[[i]]*scl[[i]]), HMAX))
band <- sapply(seq_along(rows), function(i) max(scl[[i]])*Hrow[i])
roww <- sapply(seq_along(rows), function(i) sum(ars[[i]]*scl[[i]]*Hrow[i]) + (length(ars[[i]])-1)*gx)
Htot <- mT + sum(band) + (length(rows)-1)*gy + mB

trim_opt <- function(i) sprintf("trim=%.4fin %.4fin %.4fin %.4fin, clip",
                                i$tl*i$pageW, i$tb*i$pageH, i$tr*i$pageW, i$tt*i$pageH)
L <- c("\\documentclass[border=0pt]{standalone}",
       "\\usepackage{graphicx}\\usepackage{tikz}\\renewcommand{\\familydefault}{\\sfdefault}",
       "\\begin{document}", "\\begin{tikzpicture}[x=1in,y=1in]",
       sprintf("\\fill[white] (0,0) rectangle (%.4f,%.4f);", W, Htot))
yTop <- Htot - mT
for (i in seq_along(rows)) {
  H <- Hrow[i]; x0 <- mL + max(0,(usable-roww[i])/2)
  for (j in seq_along(rows[[i]]$stem)) {
    s <- scl[[i]][j]; ph <- s*H; w <- ars[[i]][j]*ph; y0 <- yTop - band[i]/2 - ph/2; ii <- info[[i]][[j]]
    L <- c(L,
      sprintf("\\node[anchor=south west,inner sep=0] at (%.4f,%.4f) {\\includegraphics[%s,width=%.4fin,height=%.4fin]{%s.pdf}};",
              x0, y0, trim_opt(ii), w, ph, ii$stem),
      sprintf("\\node[anchor=east,font=\\bfseries\\large] at (%.4f,%.4f) {%s};", x0-0.03, yTop-0.06, rows[[i]]$L[j]))
    x0 <- x0 + w + gx
  }
  yTop <- yTop - band[i] - gy
}
L <- c(L, "\\end{tikzpicture}", "\\end{document}")
writeLines(L, file.path(fd, "Module3_composite_vector.tex"))
cat(sprintf("wrote Module3_composite_vector.tex | canvas %.2f x %.2f in\n", W, Htot))
