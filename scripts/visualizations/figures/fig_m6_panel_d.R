#!/usr/bin/env Rscript
# ============================================================================
# fig_m6_panel_d.R  -> Module 6 panel d (Q3): is the exposure signal disease-relevant?
# ----------------------------------------------------------------------------
# Disease signal is near-universal (the proteome is a broad health readout), so
# "predicts disease" alone does not tell you whether a score is a trustworthy,
# actionable exposure biomarker. The two properties that DO are READING the
# exposure and TRACKING within-person change. So the scatter axes are those two:
#   x = how well the proteome reads the exposure (held-out incremental, beyond
#       covariates); y = how well the score tracks within-person change (panel c
#       Delta-correlation). Point SIZE = disease relevance (best held-out C-index
#       gain over a 15-disease panel). No magnitude cutoffs -- the three regions
#       are illustrative archetypes:
#     right + top  = reads + tracks  -> modifiable, monitorable biomarker
#     right + low  = reads, no track -> fixed / cumulative (e.g. pack-years)
#     left         = does not read   -> confounded / non-specific (frailty)
# RIGHT: held-out Cox C-index ladders (covariates o -> +PES * -> +self-report x)
#   for one exemplar per region, with the modifiability tag. Hero contrast:
#   pack-years (fixed) vs current smoking (modifiable), both predict COPD.
#   COPD rather than emphysema: the same contrast holds on both outcomes, but
#   COPD carries ~3x the events (2,146 vs 653) and, unlike emphysema, is inside
#   the 181-disease grid the Results cite, so every drawn point is lookupable.
# Reads support/module6_quadrant_scan.R + _ladders.R + panel-b/c CI caches.
# ============================================================================
local({
  cand <- c(file.path(getwd(),"scripts","visualizations","common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cand[dir.exists(cand)][1]
  for (f in c("figure_paths.R","load_heap_results.R","plot_theme.R","label_helpers.R","export_helpers.R")) source(file.path(cm,f))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })

od  <- heap_project_output("module6_pes_longitudinal","base")
mpd <- heap_project_output("module6_pes_longitudinal","multipes_disease")

# --- render profile: standalone (default) vs small composite CELL ----------
# CELL mode = the version that goes in the 2x2 main figure: scatter + ladders
# only (no heatmap), Nature-spec small fonts, rendered at the exact cell size.
CELL  <- nzchar(Sys.getenv("HEAP_CELL"))
BS    <- if (CELL) 7    else 10     # base_size
LBL   <- if (CELL) 1.95 else 2.8    # primary point labels (mm)
LBL2  <- if (CELL) 1.75 else 2.45   # secondary labels (mm)
TTL   <- if (CELL) 8    else 11     # sub-panel title (pt)
MAXSZ <- if (CELL) 5.0  else 7.5    # scatter max point area
PSZ   <- if (CELL) 2.9  else 4.0    # ladder PES point
PSZ2  <- if (CELL) 1.7  else 2.3    # ladder cov / E points
LGT   <- if (CELL) 6    else 7.5    # legend title
LGX   <- if (CELL) 5.5  else 7.2    # legend text
BPAD  <- if (CELL) 0.35 else 0.8    # repel box padding
D_W   <- 5.05; D_H <- 2.85         # cell size (in); WIDER than c; taller (grow bottom row for vertical spacing -- keep C_H==D_H so fonts stay 1:1)

# --- read (held-out incremental, covariate skill floored at no-skill) ---
hdr <- as.data.table(load_module6_pes_longitudinal("base","holdout"))
hdr[, r2n := suppressWarnings(as.numeric(r2))]; hdr[, aucn := suppressWarnings(as.numeric(auc))]
inc <- hdr[, {
  typ <- if (all(is.na(r2n))) "binary" else "continuous"
  mv  <- function(mod) if (typ=="continuous") mean(r2n[model==mod], na.rm=TRUE) else mean(aucn[model==mod], na.rm=TRUE)
  flo <- if (typ=="continuous") 0 else 0.5
  .(incr = mv("prot_plus_cov") - max(mv("cov_only"), flo))
}, by=exposure_id]
# --- tracks (within-person Delta-correlation) ---
wc  <- fread(file.path(od,"PESlong_base_WithinDeltaCorCI.tsv"))[, .(exposure_id, dcor=dcor_prot, tracks=prot_lo>0)]
# --- disease (best held-out PES C-index gain over a clean panel) ---
scanAll <- fread(file.path(mpd,"quadrant_scan.tsv"))
best <- scanAll[!disease %in% c("Obesity","Alcoholic liver disease") & events>=500 & is.finite(dC_pes)][order(exposure_id,-dC_pes)][, .SD[1], by=exposure_id]

CAT2 <- c(modifiable="#009E73", fixed="#E69F00", confounded="#D55E00", classifier="#999999")  # colourblind-safe (Okabe-Ito)
EX <- data.table(
  exposure_id = c("current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days","pack_years_of_smoking_f20161_0_0",
                  "alcohol_intake_frequency_f1558_0_0","number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
                  "coffee_intake_f1498_0_0","usual_walking_pace_f924_0_0"),
  disease = c("COPD","COPD","Alcohol-use disorder","Type-2 diabetes","Type-2 diabetes","Type-2 diabetes"),
  cat2  = c("modifiable","fixed","modifiable","modifiable","classifier","confounded"),
  short = c("Current smoking","Pack-years smoking","Alcohol intake","Vigorous activity","Coffee","Walking pace"),
  dz    = c("COPD","COPD","alcohol-use disorder","type-2 diabetes","type-2 diabetes","type-2 diabetes"),
  pt    = c("Current smoking → COPD","Pack-years → COPD","Alcohol → alcohol-use disorder","Vig. activity → T2D","Coffee → T2D","Walking pace → T2D"),
  interp= c("reads + tracks + predicts → modifiable, monitorable biomarker",
            "reads + predicts but does NOT track → fixed, cumulative damage",
            "reads + tracks + predicts → modifiable biomarker",
            "reads + tracks + predicts → modifiable exposure biology",
            "reads + tracks well but no disease signal → exposure classifier (diet)",
            "predicts disease but does NOT read the exposure → confounding / reverse causation"))
lad5  <- fread(file.path(mpd,"quadrant_ladders.tsv"))[, .(exposure_id, disease, C0, C1_PES, C2_E)]
ladX  <- scanAll[exposure_id=="pack_years_of_smoking_f20161_0_0" & disease=="COPD", .(exposure_id, disease, C0, C1_PES, C2_E)]
d <- merge(EX, rbind(lad5, ladX), by=c("exposure_id","disease"))
d <- merge(d, inc, by="exposure_id"); d <- merge(d, wc, by="exposure_id", all.x=TRUE)

# ---- bootstrap C-indices (B=400 out-of-bag) replace the single 70/30 split ----
# The ladder above comes from one split with seed 42 and carries no uncertainty,
# which is why the Fig 6 caption had to disclose that it disagrees with S18
# (NUM-2: 3 of 12 estimates fell outside the table's CI). S18 covers every pair
# plotted here, so the figure can simply use it: same cohort, same models, but
# out-of-bag resampled. Agreement with the split is r = 0.988, mean |diff| 0.005.
BOOTF <- file.path(mpd, "pes_disease_ladder.tsv")   # the B=400 bootstrap table that ships as S18
if (file.exists(BOOTF)) {
  bb <- fread(BOOTF); bn <- names(bb)
  DZRE <- c(Emphysema="j43_first_reported_emphysema", COPD="j44_first_reported_other_chronic",
            `Alcohol-use disorder`="f10_first_reported_mental", `Type-2 diabetes`="e11_first_reported_non_insulin",
            `Ischaemic heart`="i25_first_reported_chronic_ischaem", `Lipid disorder`="e78_first_reported_disorders",
            Depression="f32_first_reported_depressive")
  bb[, `:=`(ek=get(bn[2]), dzn=get(bn[3]))]
  hit <- rbindlist(lapply(seq_len(nrow(d)), function(i){
    r <- bb[ek==d$exposure_id[i] & grepl(DZRE[[d$disease[i]]], dzn)][1]
    if (!nrow(r) || is.na(r$ek)) return(NULL)
    data.table(exposure_id=d$exposure_id[i], disease=d$disease[i],
               bC0=r$C_cov, bE=r$C_covE, bP=r$C_covP, bPlo=r$C_covP_lo, bPhi=r$C_covP_hi)}))
  if (nrow(hit)) {
    d <- merge(d, hit, by=c("exposure_id","disease"), all.x=TRUE)
    d[is.finite(bC0), `:=`(C0=bC0, C2_E=bE, C1_PES=bP)]      # bootstrap point estimates
    message(sprintf("  ladder: bootstrap C-index for %d/%d exemplars", sum(is.finite(d$bP)), nrow(d)))
  }
}
if (!"bPlo" %in% names(d)) d[, `:=`(bPlo=NA_real_, bPhi=NA_real_)]
d[, tracks := ifelse(is.na(tracks), FALSE, tracks)]

# ----------------------------------------------------------------------------
# LEFT: reads x tracks scatter, point size = disease relevance
# ----------------------------------------------------------------------------
q <- merge(merge(best[, .(exposure_id, dC_pes)], inc, by="exposure_id"), wc, by="exposure_id")  # need both axes
q[, `:=`(incr = pmax(pmin(incr, 0.47), -0.04), dcorc = pmax(pmin(dcor, 0.82), -0.13), dC_sz = pmin(pmax(dC_pes,0), 0.20))]
# extra scatter labels (notable cloud points; colored + labeled here but not in the ladders/heatmap)
# Only Strenuous sport kept. The two extra smoking encodings (Smoking status: current,
# Current smoking: No) were redundant with the "Current smoking" exemplar and crowded the
# top-right corner -> dropped to the grey cloud so the labeled points stay easy to read.
EXTRA <- data.table(
  exposure_id=c("types_of_physical_activity_in_last_4_weeks_f6164_0_0.multi_Strenuous_sports"),
  ecat=c("modifiable"),
  eshort=c("Strenuous sport"))
labelmap <- rbind(d[, .(exposure_id, ecat=as.character(cat2), eshort=short)], EXTRA)
q <- merge(q, labelmap, by="exposure_id", all.x=TRUE)
qL <- ggplot(q, aes(incr, dcorc)) +
  geom_hline(yintercept=0, colour="grey88", linewidth=0.4) +
  geom_point(data=q[is.na(ecat)], aes(size=dC_sz), colour="grey72", alpha=0.40) +
  geom_point(data=q[!is.na(ecat)], aes(size=dC_sz, colour=ecat)) +
  ggrepel::geom_text_repel(data=q[!is.na(ecat)], aes(label=eshort, colour=ecat), fontface="bold", size=LBL,
                           box.padding=BPAD, point.padding=0.3, min.segment.length=0, seed=3, max.overlaps=20, show.legend=FALSE) +
  scale_colour_manual(NULL, values=CAT2, breaks=c("modifiable","fixed","confounded","classifier"),
      labels=c("modifiable biomarker (reads + tracks)","fixed / cumulative (reads, no track)",
               "confounded (doesn't read)","classifier (reads, no disease)")) +
  scale_size_area("disease gain", max_size=MAXSZ, breaks=c(0.03,0.10,0.18)) +
  guides(colour=guide_legend(order=1, ncol=2, override.aes=list(size=if (CELL) 2.2 else 3.4)), size=guide_legend(order=2)) +
  scale_x_continuous(if (CELL) expression("reads:  incremental "*R^2*"/AUC (beyond covariates)") else "reads the exposure  (beyond covariates) →", limits=c(-0.05,0.48), breaks=seq(0,0.4,0.2)) +
  scale_y_continuous(if (CELL) expression("tracks:  within-person "*Delta*"-correlation") else "tracks within-person change →", limits=c(-0.14,0.86)) +
  labs(title=if (CELL) "Reads × tracks" else "Which exposure signals to trust") +
  theme_heap(base_size=BS) +
  theme(plot.title=element_text(face="bold", size=TTL, hjust=if (CELL) 0.5 else 0), plot.title.position="panel",
        axis.title=element_text(face=if (CELL) "plain" else "bold"),
        legend.position="bottom", legend.box="vertical",
        legend.title=element_text(size=LGT), legend.text=element_text(size=LGX), legend.margin=margin(1,1,1,1),
        legend.spacing.y=unit(0.01,"cm"), legend.key.height=unit(if (CELL) 0.22 else 0.32,"cm"),
        panel.grid.major = if (CELL) element_blank() else element_line(colour="grey92"),
        plot.margin=margin(5,3,2,7))

# ----------------------------------------------------------------------------
# RIGHT: exemplar disease-connection heatmap (specificity overlay on the grid)
# Biomarkers (read+track) light up FEW, biologically-specific diseases; confounders
# (don't read) light up the whole row (frailty / reverse causation is non-specific).
# ----------------------------------------------------------------------------
DZORD <- c("Type-2 diabetes","Lipid disorder","Ischaemic heart disease","Hypertension","Heart failure","Atrial fibrillation",
           "Stroke","COPD","Emphysema","Asthma","Chronic kidney disease","Depression","Dementia","Lung cancer","Knee osteoarthritis")
HEX <- rbind(EX[, .(exposure_id, short, cat2=as.character(cat2))],
             data.table(exposure_id="nap_during_day_f1190_0_0", short="Daytime napping", cat2="confounded"))
HEX[, rorder := match(cat2, c("modifiable","fixed","classifier","confounded"))]; HEX <- HEX[order(rorder)]
hh <- merge(scanAll[disease %in% DZORD, .(exposure_id, disease, dC_pes)], HEX, by="exposure_id")
hh[, disease := factor(heap_americanize(disease), levels=heap_americanize(DZORD))]; hh[, short := factor(short, levels=rev(HEX$short))]
hh[, gain := pmin(pmax(dC_pes,0), 0.15)]
ylabcol <- setNames(CAT2[HEX$cat2], HEX$short)
hM <- ggplot(hh, aes(disease, short, fill=gain)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(data=hh[dC_pes>=0.05], aes(label=sprintf("%.2f", dC_pes)), size=2.2, colour="white") +
  scale_fill_gradient("held-out C-index gain from the PES", low="#F2F6FB", high="#08306B", limits=c(0,0.15), breaks=c(0,0.05,0.10,0.15)) +
  scale_x_discrete(position="top") +
  labs(title="Disease-connection profile: biomarkers are specific; confounders are broad", x=NULL, y=NULL) +
  theme_heap(base_size=10) +
  theme(axis.text.x.top=element_text(angle=45, hjust=0, size=7), axis.text.y=element_text(colour=ylabcol[levels(hh$short)], face="bold", size=8.5),
        panel.grid=element_blank(), plot.title=element_text(face="bold", size=10.5),
        legend.position="bottom", legend.key.height=unit(0.28,"cm"), legend.title=element_text(size=8))

# ----------------------------------------------------------------------------
# TOP-RIGHT: exemplar C-index ladders (covariates -> +PES -> +E): objective vs self-report
# ----------------------------------------------------------------------------
d[, cat2f := factor(cat2, levels=c("confounded","classifier","fixed","modifiable"))]
d2 <- d[order(cat2f, C1_PES-C0)]; d2[, y := .I]; d2[, mshape := ifelse(tracks, 16L, 15L)]
d2[, mtag := ifelse(cat2 %in% c("classifier","confounded"), "", ifelse(tracks, "●  modifiable", "■  fixed / cumulative"))]
d2[, lab_l := if (CELL) short else pt]   # cell: short exposure name; standalone: "exposure → disease"
# cell sub-line = the ARCHETYPE only (disease moved next to the ladder line, below)
d2[, sub_show := if (CELL) cat2 else mtag]
pts <- rbindlist(list(d2[, .(y, m="Cov", x=C0)], d2[, .(y, m="E", x=C2_E)], d2[, .(y, m="PES", x=C1_PES)]))
# cell row: exposure (bold) -> archetype (italic) on the left; disease at the RIGHT end of each ladder line.
YL <- if (CELL) 0.78 else 0.30   # exposure name
YM <- if (CELL) 0.50 else 0.07   # archetype tag (under the exposure)
YP <- if (CELL) 0.22 else 0.62   # ladder line + the disease label (right of the PES point)
# drawn marker key (open circle / x / filled circle); unicode bullets do NOT render in the PDF font
lkey <- data.frame(kx=c(0.51,0.74,0.97), ky=7.5, ks=c(21L,4L,16L), kl=c("covariates","+ self-report","+ PES"))  # spread across the width so labels never collide
keylayers <- if (CELL) list(
  geom_point(data=lkey, aes(x=kx, y=ky), shape=lkey$ks, colour="grey25", fill="white", size=PSZ2, stroke=0.8, inherit.aes=FALSE),
  geom_text(data=lkey, aes(x=kx+0.013, y=ky, label=kl), hjust=0, size=LBL2, colour="grey30", inherit.aes=FALSE)) else NULL
lR <- ggplot() + keylayers +
  geom_text(data=d2, aes(x=0.50, y=y+YL, label=lab_l, colour=cat2f), hjust=0, fontface="bold", size=LBL) +
  geom_text(data=if (CELL) d2 else d2[mtag!=""], aes(x=0.50, y=y+YM, label=sub_show, colour=cat2f), hjust=0, size=LBL2, fontface="italic") +
  geom_text(data=d2, aes(x=C1_PES+0.028, y=y+YP, label=dz), hjust=0, size=LBL2, fontface="italic", colour="grey35") +   # disease at the RIGHT end of each ladder line (gap clears the PES dot)
  geom_segment(data=d2, aes(x=C0, xend=C1_PES, y=y+YP, yend=y+YP, colour=cat2f), linewidth=if (CELL) 0.7 else 1.0, alpha=0.4) +
  # NB no error bars: the bootstrap CIs are ~0.024 C-index units, 4% of the panel
  # width against a marker that is already ~1.8% -- they are illegible at cell size
  # and add nothing. \Tref{pes_cindex} carries them per pair.
  geom_point(data=pts[m=="Cov"], aes(x, y+YP), shape=21, fill="white", colour="grey45", size=PSZ2, stroke=0.7) +
  geom_point(data=pts[m=="E"], aes(x, y+YP), shape=4, colour="grey35", size=PSZ2, stroke=1.0) +
  geom_point(data=d2, aes(C1_PES, y+YP, colour=cat2f, shape=factor(mshape)), size=PSZ) +
  scale_shape_manual(values=c(`16`=16,`15`=15), guide="none") + scale_colour_manual(values=CAT2, guide="none") +
  scale_x_continuous("held-out Cox C-index for disease", limits=c(0.50,1.10), breaks=seq(0.55,0.90,0.10)) +
  scale_y_continuous(NULL, breaks=NULL, limits=c(0.6, if (CELL) 7.7 else 7.0), expand=expansion(add=c(0.05,0.10))) +
  labs(title=if (CELL) "Add PES vs self-report" else "Objective vs self-report: covariates / + self-report / + PES",
       subtitle=NULL) +
  theme_heap(base_size=BS) + theme(panel.grid.major=element_blank(),
        plot.title=element_text(face="bold", size=if (CELL) TTL-0.5 else 9.8, hjust=if (CELL) 0.5 else 0),
        plot.title.position="panel", axis.title=element_text(face=if (CELL) "plain" else "bold"),
        plot.margin=margin(8,5,2,9))   # top: key/title room; left: gap from the scatter; right: disease labels

FIGDIR <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module6")
if (CELL) {
  # main-figure cell: scatter + ladders only (heatmap dropped -> supplement)
  pcell <- wrap_plots(qL, lR, widths=c(1.08, 1.06))   # give the ladder a bit more width (the decluttered scatter can spare it) so labels breathe
  ggsave(file.path(FIGDIR,"fig_m6_panel_d_cell.png"), pcell, width=D_W, height=D_H, dpi=400, bg="white")
  ggsave(file.path(FIGDIR,"fig_m6_panel_d_cell.pdf"), pcell, width=D_W, height=D_H, bg="white")
  message("panel d CELL (scatter + ladders) done")
} else {
  p <- wrap_plots(wrap_plots(qL, lR, widths=c(1.15, 1)), hM, ncol=1, heights=c(1.45, 1)) +
    plot_annotation(title="Disease relevance: which exposure signals are trustworthy, modifiable, and specific",
          theme=theme(plot.title=element_text(face="bold", size=12.5)))
  ggsave(file.path(FIGDIR,"fig_m6_panel_d.png"), p, width=13.5, height=9.2, dpi=140, bg="white")
  message("panel d v6 (grid + ladders + heatmap) done")
}
