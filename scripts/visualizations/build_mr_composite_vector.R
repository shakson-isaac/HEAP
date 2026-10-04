#!/usr/bin/env Rscript
# Emit a TikZ doc laying the panel VECTOR PDFs onto one 7.0in-wide canvas, 2-row
# layout: top = [tier ladder | DAG] (=panel a) | b (motifs); bottom = c | d | e.
# Panels scaled to a common row height (aspect preserved); white borders auto-
# trimmed via includegraphics[trim,clip]. Compile -> one editable vector PDF.
suppressPackageStartupMessages(library(png))
OUT   <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module5")
SCHEM <- file.path(Sys.getenv("HEAP_ROOT", "/n/groups/patel/shakson_ukb/HEAP"), "mr_schematics")
DPI <- 300
P <- list(
  fun=list(dir=OUT,   stem="_panel_tier_ladder"),
  a  =list(dir=SCHEM, stem="fig_mr_evidence_dag"),
  b  =list(dir=OUT,   stem="_panel_b_folded"),
  c  =list(dir=OUT,   stem="_panel_c"),
  d  =list(dir=OUT,   stem="_panel_d"),
  e  =list(dir=OUT,   stem="fig_mr_main_mediators_dag_cell"))
meas <- function(p){ im<-readPNG(file.path(p$dir,paste0(p$stem,".png"))); d<-dim(im)
  w <- if(length(d)==3) (im[,,1]>=.992)&(im[,,2]>=.992)&(im[,,3]>=.992) else im>=.992
  rk<-which(rowSums(!w)>0); ck<-which(colSums(!w)>0); H<-d[1]; W<-d[2]
  tl<-(min(ck)-1)/W; tr<-(W-max(ck))/W; tt<-(min(rk)-1)/H; tb<-(H-max(rk))/H
  pageW<-W/DPI; pageH<-H/DPI
  list(pageW=pageW,pageH=pageH,tl=tl,tr=tr,tt=tt,tb=tb, cW=pageW*(1-tl-tr), cH=pageH*(1-tt-tb)) }
I <- lapply(P, meas); ar <- sapply(I, function(i) i$cW/i$cH)
trim_opt <- function(i) sprintf("trim=%.4fin %.4fin %.4fin %.4fin, clip",
                                i$tl*i$pageW, i$tb*i$pageH, i$tr*i$pageW, i$tt*i$pageH)
inc <- function(nm,x,y,w,h) sprintf(
  "\\node[anchor=south west,inner sep=0] at (%.4f,%.4f) {\\includegraphics[%s,width=%.4fin,height=%.4fin]{%s}};",
  x,y,trim_opt(I[[nm]]),w,h, file.path(P[[nm]]$dir, paste0(P[[nm]]$stem,".pdf")))
ltr <- function(x,y,L) sprintf("\\node[anchor=north west,font=\\fontsize{8}{9}\\selectfont\\bfseries] at (%.4f,%.4f) {%s};",x,y,L)

# --- layout (2 rows: [a-ladder|DAG|b] / [c|d|e]) ------------------------------
W<-7.0; mL<-0.14; mR<-0.08; mT<-0.10; mB<-0.07; gx<-0.12; gf<-0.10; FUNH<-0.70
usable <- W - mL - mR
# top row: funnel (0.8H) | DAG | b  -> fill usable
H1 <- (usable - gf - gx)/(FUNH*ar["fun"] + ar["a"] + ar["b"])
funW<-ar["fun"]*FUNH*H1; aW<-ar["a"]*H1; bW<-ar["b"]*H1
# bottom row: c | d | e  (e authored taller -> gets more vertical room; c & d centred beside it)
EK  <- 1.06                                  # e row-height / (c,d) row-height
Hcd <- (usable - 2*gx)/(ar["c"]+ar["d"]+EK*ar["e"]); He <- EK*Hcd
cW<-ar["c"]*Hcd; dW<-ar["d"]*Hcd; eW<-ar["e"]*He
H2 <- He                                     # bottom-row band height (= e)
Htot <- mT + H1 + 0.14 + H2 + mB
y1 <- Htot - mT - H1; y2 <- mB

L <- c("\\documentclass[border=0pt]{standalone}",
       "\\usepackage{graphicx}\\usepackage{tikz}\\renewcommand{\\familydefault}{\\sfdefault}",
       "\\begin{document}","\\begin{tikzpicture}[x=1in,y=1in]",
       sprintf("\\fill[white] (0,0) rectangle (%.4f,%.4f);", W, Htot))
# row1: funnel (left, vertically centred in row) + DAG + b
L <- c(L, inc("fun", mL, y1 + (H1-FUNH*H1)/2, funW, FUNH*H1))
L <- c(L, inc("a", mL+funW+gf, y1, aW, H1))
L <- c(L, inc("b", mL+funW+gf+aW+gx, y1, bW, H1))
# row2: c | d (top-aligned with e so titles/letters line up) | e (fills band)
ycd <- y2 + (He-Hcd)
L <- c(L, inc("c", mL, ycd, cW, Hcd),
          inc("d", mL+cW+gx, ycd, dW, Hcd),
          inc("e", mL+cW+dW+2*gx, y2, eW, He))
# letters
L <- c(L, ltr(mL-0.02, y1+H1+0.05, "a"), ltr(mL+funW+gf+aW+gx-0.02, y1+H1+0.05, "b"),
          ltr(mL-0.02, ycd+Hcd+0.05, "c"), ltr(mL+cW+gx-0.02, ycd+Hcd+0.05, "d"),
          ltr(mL+cW+dW+2*gx-0.02, y2+He+0.05, "e"))
L <- c(L, "\\end{tikzpicture}","\\end{document}")
writeLines(L, file.path(OUT,"Module5_composite_vector.tex"))
cat(sprintf("wrote Module5_composite_vector.tex | canvas %.2f x %.2f in\n", W, Htot))
