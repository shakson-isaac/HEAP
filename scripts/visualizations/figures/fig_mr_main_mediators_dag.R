#!/usr/bin/env Rscript
# ============================================================================
# fig_mr_main_mediators.R — Panel D of the MR main figure.
# Fully-decomposed cis-anchored causal MEDIATORS: for one exemplar exposure per
# protein, the exposure→protein→disease chain with β [95% CI] on every leg, the
# exposure→disease TOTAL effect, and the % mediated. The protein→disease (causal)
# leg is shown for BOTH instrument arms (UKB Olink + deCODE SomaScan). ASGR1 hero.
# CELL render mode per docs/MULTIPANEL_FIGURE_GUIDE.md.
# ============================================================================
local({
  cand <- c(file.path(getwd(),"scripts","visualizations","common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset="fig_mr_main_mediators_dag")
CELL <- nzchar(Sys.getenv("HEAP_CELL"))
BS  <- if (CELL) 9 else 10
LBL <- if (CELL) 3.3 else 2.7
CIL <- LBL - 0.55
B_W <- if (CELL) 6.9 else 7.2; B_H <- if (CELL) 3.3 else 4.8
PAL_ARM <- c(UKB="#2C7FB8", DECODE="#D95F0E")

sd <- heap_resolve_output(file.path("mr_edges","summary"), must_exist=TRUE)
u <- fread(file.path(sd,"MRmotifs.tsv"),
           select=c("Protein","Exposure","Disease","motif_A_mediator",
                    "beta_EP","se_EP","beta_PDcis","se_PDcis","beta_ED","se_ED"))
dec <- fread(file.path(sd,"DECODE","MRmotifs.tsv"),
             select=c("Protein","Disease","beta_PDcis","se_PDcis"))
dec <- unique(dec, by=c("Protein","Disease"))
co <- load_coloc_results()[edge_dir=="Pcis_to_D"]
# Module-1 exposomic R² (unique "E" block) — shades the protein nodes.
m1 <- load_module1_predictive_r2(covarType="base",method="lasso",level="coarse",experiment="M1_base_lasso")
eR2tab <- m1[block=="E", .(eR2=mean(r2,na.rm=TRUE)), by=.(Protein=omic)][is.finite(eR2)]
ER2MAX <- max(eR2tab$eR2, na.rm=TRUE)

pretty_dz <- function(x){ x <- sub("^finngen_R12_","",x)
  c(E4_LIPOPROT="Lipoprotein disorder", I9_AF="Atrial fibrillation",
    T2D="Type 2 diabetes", I9_HYPTENSESS="Hypertension")[x] }
# THIS is the copy that reaches the manuscript: build_mr_composite_vector.R takes
# panel e from fig_mr_main_mediators_dag_cell, and sync_figures.R maps
# Fig4 = exploratory/module5/Module5_composite_vector. fig_mr_main_mediators.R
# holds an identical PICK -- keep the two in sync or delete one.
#
# One exemplar exposure per protein. EVERY ROW HERE MUST BE A TIER-1 MEDIATOR --
# i.e. present with motif "A Mediator (E->P->D)" in
#   HEAP/docs/manuscript_stats/module5/mr_triad_motifs.tsv
# which is the rule the main text, Fig 4b and S14 all use. There are exactly 6
# such triads (7 before the 2026-08-11 symmetric-Steiger fix, which retired
# ASGR1 x usual_walking_pace) and they span exactly these 3 proteins, so the protein cast is the
# complete Tier-1 mediator set; only the exposure is a choice. The options are:
#   ASGR1  -> TV time | walking pace | never-smoking | age first sex  (all vs LIPOPROT)
#   ADM    -> pack-years only                                        (vs I9_AF)
#   FURIN  -> pack-years | TV time                                   (vs HYPTENSESS)
#
# Do NOT pick from MRmotifs.tsv's precomputed motif_A_mediator flags: that is the
# NOMINAL rule (84 triads) and it is NOT nested with Tier-1. FURIN x TV time is
# Tier-1 but nominal-FALSE; the reverse also happens. Selecting on the nominal
# flag is what produced the previous cast (NUM-3): ADM x walking pace and
# FURIN x walking pace are both ABSENT from the Tier-1 table, and FURIN x walking
# pace failed the nominal rule too, so it was a mediator under neither definition.
# Both traced to usual_walking_pace, the frailty proxy, whose ADM/FURIN edges are
# tiered "E Disease-liability (D->P)" -- the reverse of the arrow drawn here.
PICK <- data.table(
  # Row order is top-to-bottom in the panel. ASGR1 leads as the two-arm hero;
  # FURIN then ADM, so the two TV-time triads sit together and pack-years closes.
  Protein=c("ASGR1","FURIN","ADM"), hero=c(TRUE,FALSE,FALSE),
  match=c("time_spent_watching_television_tv","time_spent_watching_television_tv","pack_years_of_smoking"),
  dz=c("finngen_R12_E4_LIPOPROT","finngen_R12_I9_HYPTENSESS","finngen_R12_I9_AF"),
  exp_lab=c("TV time","TV time","pack-years"))

med <- rbindlist(lapply(seq_len(nrow(PICK)), function(i){
  p <- PICK$Protein[i]
  r <- u[Protein==p & grepl(PICK$match[i], Exposure) & Disease==PICK$dz[i]][1]
  dr <- dec[Protein==p & Disease==PICK$dz[i]]
  ci <- function(b,s) c(b, b-1.96*s, b+1.96*s)
  ep <- ci(r$beta_EP,r$se_EP); pdU <- ci(r$beta_PDcis,r$se_PDcis); ed <- ci(r$beta_ED,r$se_ED)
  pdD <- if (nrow(dr)) ci(dr$beta_PDcis,dr$se_PDcis) else c(NA,NA,NA)
  cph <- co[protID==p & target==r$Disease]
  coloc <- if (nrow(cph)) max(cph$`PP.H4`, na.rm=TRUE) else NA_real_
  data.table(Protein=p, hero=PICK$hero[i], exp_lab=PICK$exp_lab[i],
             disease=pretty_dz(r$Disease),
             ep_b=ep[1],ep_lo=ep[2],ep_hi=ep[3], pdU_b=pdU[1],pdU_lo=pdU[2],pdU_hi=pdU[3],
             pdD_b=pdD[1],pdD_lo=pdD[2],pdD_hi=pdD[3], ed_b=ed[1],ed_lo=ed[2],ed_hi=ed[3],
             coloc=coloc, prop=100*(r$beta_EP*r$beta_PDcis)/r$beta_ED)
}))
med <- merge(med, eR2tab, by="Protein", all.x=TRUE, sort=FALSE)
med[, coloc_ok := is.finite(coloc) & coloc>=0.8]
med[, coloc_lab := fifelse(coloc_ok, sprintf("coloc PP.H4 %.3f  (shared variant)", coloc),
                    fifelse(is.finite(coloc), sprintf("coloc PP.H4 %.2f  (LD-confounded)", coloc),
                            "coloc not tested"))]
# ---- triad-DAG: E→D direct (grey) vs E→P→D mediated (green); mediated-path
#      width = % mediated; protein shaded by exposomic R² ------------------------
Ex <- 0.2; Dx <- 8.9; Px <- (Ex+Dx)/2              # symmetric triad: apex centred over the base
# NOTE: Ph is set PER ROW below and is a column from that point on -- do not reintroduce a scalar.
# base widened from 0-6.3 (2026-08-31): the leg betas had ~2.8 units to sit in
# between the node box and the disease name and needed ~2.4, so they crowded the
# protein. The wider base gives each side ~3.0.
# per-row apex height: single-arm (Olink-only) triads are flatter than the two-arm ASGR1 triad,
# then stack all three with even gaps so ADM/FURIN take less height and ASGR1 keeps its room
med[, Ph := fifelse(is.na(pdD_b), 1.02, 1.52)]
{ HBOT <- 0.55; GAP <- 0.58; ht <- med$Ph + 0.80; yv <- numeric(nrow(med))
  yv[nrow(med)] <- HBOT
  for (i in (nrow(med)-1):1) yv[i] <- yv[i+1] + ht[i+1] + GAP + HBOT
  med[, yc := yv] }
# VARIANTS (HEAP_E_VAR): the panel is placed at ~50% in the composite, so a
# label.size below ~0.1mm is clamped to the device hairline and stops getting
# thinner -- borderless is the only way to lose the box. The arrowheads were
# 0.03in (2.2pt) and landed under the protein label, so the legs did not read
# as directed edges at all.
#   1 = borderless labels, bulkier heads, same positions
#   2 = borderless + betas moved OFF the legs (arrows visible end to end)
#   3 = no white grounds at all; betas ride the leg ends, heaviest arrows
VAR <- Sys.getenv("HEAP_E_VAR", "2")          # 2 is the shipped layout
VTAG <- if (nzchar(Sys.getenv("HEAP_E_VARTAG"))) paste0("_v", VAR) else ""   # side-by-side compares only
med[, w_med := if (VAR == "3") 0.75 + prop/20 else 0.55 + prop/22]
AHD <- c("1"=0.055, "2"=0.085, "3"=0.070, "6"=0.085)[[VAR]]
ar  <- arrow(length=unit(if(CELL)AHD else 0.048,"in"), type="closed")
ar2 <- arrow(length=unit(if(CELL)AHD*0.75 else 0.04,"in"), type="closed")
LSZ <- NA_real_                                # borderless in every variant
LFIL<- if (VAR == "3") NA else "white"         # V3 drops the white ground too
SHORT <- VAR == "6"          # point estimate only: a short label can hug its edge
fmtci  <- function(b,lo,hi) if (SHORT) sprintf("β = %.2f", b) else sprintf("β = %.2f [%.2f, %.2f]", b, lo, hi)
fmtci2 <- function(b,lo,hi) sprintf("%.2f [%.2f, %.2f]", b, lo, hi)   # arm name carries the "β ="
# leg-label anchors, per variant (columns, so aes() stays plain)
LGx <- c(Ex+0.24, Px-0.62); RGx <- c(Px+0.62, Dx-0.24)   # the two leg spans
MIDL <- mean(LGx) - 0.18; MIDR <- mean(RGx) + 0.45      # leg midpoints, eased off the node
# APEX HEIGHT IS PER ROW (Ph above: 1.52 for the two-arm ASGR1 triad, 1.02 for the
# flatter single-arm ones), so the label anchor has to be per row too. A scalar
# anchor put every label at the same absolute height -- level with ASGR1's node,
# but 0.50 ABOVE FURIN's and ADM's, which is what made those two look unmoored.
# The clearance each label needs also follows its own leg: half a label times that
# leg's slope, plus half the text height.
LHW <- 1.36                                              # half a full-CI label, data units
med[, midy := (0.12 + Ph - 0.34)/2]                      # this row's leg midpoint
med[, clr  := LHW*((Ph - 0.46)/diff(LGx)) + 0.155 + 0.06]
med[, `:=`(lx_ep = switch(VAR, "1"=(Ex+Px)/2, "2"=MIDL, "3"=Ex+0.70, "6"=MIDL),
           ly_ep = yc + switch(VAR, "1"=Ph*0.48, "2"=midy+clr, "3"=Ph*0.26, "6"=midy+0.34),
           lx_pd = switch(VAR, "1"=(Px+Dx)/2, "2"=MIDR, "3"=Dx-0.70, "6"=mean(RGx)),
           ly_pdU= yc + switch(VAR, "1"=Ph*0.63, "2"=midy+clr, "3"=Ph*0.42, "6"=midy+0.34),
           ly_pdD= yc + switch(VAR, "1"=Ph*0.20, "2"=midy+clr+0.44, "3"=Ph*0.14, "6"=midy-0.34))]
# Labels stay HORIZONTAL (author, 2026-08-31: diagonal on-edge labels rejected).

p <- ggplot(med) +
  # direct E→D: thin grey reference arrow
  geom_segment(aes(x=Ex+0.30, xend=Dx-0.30, y=yc-0.10, yend=yc-0.10),
               arrow=ar2, colour="grey65", linewidth=0.8, lineend="round") +
  # mediated E→P→D: thin green, subtle width ∝ % mediated
  geom_segment(aes(x=Ex+0.24, y=yc+0.12, xend=Px-0.62, yend=yc+Ph-0.34, linewidth=w_med),
               arrow=ar, colour="#2C7A3F", lineend="round") +
  geom_segment(aes(x=Px+0.62, y=yc+Ph-0.34, xend=Dx-0.24, yend=yc+0.12, linewidth=w_med),
               arrow=ar, colour="#2C7A3F", lineend="round") +
  scale_linewidth_identity() +
  # nodes
  geom_text(aes(x=Ex-0.12, y=yc, label=exp_lab), hjust=1, size=LBL, fontface="italic", colour="grey25") +
  geom_label(aes(x=Px, y=yc+Ph, label=Protein, fontface=ifelse(hero,"bold","plain"), fill=eR2),
             size=LBL+0.7, label.size=NA, label.padding=unit(0.10,"lines"), colour="grey10") +
  geom_text(aes(x=Dx+0.12, y=yc, label=disease), hjust=0, size=LBL, colour="grey15") +
  # % mediated + coloc: single line above the node (clear of the node box)
  geom_text(aes(x=Px, y=yc+Ph+0.82, label=sprintf("%.0f%% mediated  ·  coloc %.3f", prop, coloc)),
            size=LBL, fontface="bold", colour="#1A6B30") +
  # leg betas: HORIZONTAL, boxed (white box lifts them off the arrows)
  geom_label(aes(x=lx_ep, y=ly_ep, label=fmtci(ep_b,ep_lo,ep_hi)), size=CIL-0.2,
             colour="grey40", fill=LFIL, label.size=LSZ, label.padding=unit(0.07,"lines"), label.r=unit(0.04,"lines")) +
  geom_label(aes(x=lx_pd, y=ly_pdU, label=paste0("Olink ", fmtci(pdU_b,pdU_lo,pdU_hi))), size=CIL-0.2,
             colour=PAL_ARM["UKB"], fill=LFIL, fontface="bold", label.size=LSZ, label.padding=unit(0.07,"lines"), label.r=unit(0.04,"lines")) +
  geom_label(data=med[!is.na(pdD_b)], aes(x=lx_pd, y=ly_pdD, label=paste0("SomaScan ", fmtci(pdD_b,pdD_lo,pdD_hi))), size=CIL-0.2,
             colour=PAL_ARM["DECODE"], fill=LFIL, fontface="bold", label.size=LSZ, label.padding=unit(0.07,"lines"), label.r=unit(0.04,"lines")) +
  geom_label(aes(x=Px, y=yc-0.12, label=paste0("total ", fmtci(ed_b,ed_lo,ed_hi))), size=CIL-0.2,
             colour="grey45", fill=LFIL, label.size=LSZ, label.padding=unit(0.06,"lines"), label.r=unit(0.04,"lines")) +
  scale_fill_gradient(low="#EAF4E5", high="#2C7A3F", name="PXS R²",
                      limits=c(0, ER2MAX), breaks=c(0,0.05,0.10,0.15)) +
  scale_x_continuous(limits=c(-1.7, 12.4), expand=c(0,0)) +
  scale_y_continuous(limits=c(min(med$yc)-0.42, max(med$yc+med$Ph)+0.94), expand=c(0,0)) +
  coord_cartesian(clip="off") +
  labs(title="Cis-anchored causal mediators") +
  theme_void(base_size=BS) +
  theme(plot.title=element_text(hjust=0.5, size=if(CELL)11.0 else 14, face="bold", margin=margin(b=1)),
        plot.margin=margin(2,3,2,4), legend.position="right", legend.margin=margin(0,0,0,2),
        legend.title=element_text(size=if(CELL)7.5 else 9), legend.text=element_text(size=if(CELL)6.8 else 8),
        legend.key.height=unit(if(CELL)26 else 26,"pt"), legend.key.width=unit(if(CELL)7 else 8,"pt"))

if (CELL) {
  ggsave(heap_figure_path(paste0(figure_id, VTAG, "_cell.pdf"), category="exploratory", subdir="module5"), p, width=B_W, height=B_H, device=cairo_pdf)
  ggsave(heap_figure_path(paste0(figure_id, VTAG, "_cell.png"), category="exploratory", subdir="module5"), p, width=B_W, height=B_H, dpi=300)
} else {
  heap_emit_figure(p, figure_id, data=med, category="exploratory", subdir="module5",
                   formats=c("pdf","png"), width=B_W, height=B_H, website=FALSE)
}
message("mediators: ", paste(med$Protein, paste0(round(med$prop),"%"), collapse=", "))
