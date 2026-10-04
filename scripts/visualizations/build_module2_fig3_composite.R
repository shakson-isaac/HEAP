#!/usr/bin/env Rscript
# ============================================================================
# build_module2_fig3_composite.R                       == Manuscript Fig 3 (Module 2) ==
# ----------------------------------------------------------------------------
# COMBINED Module-2 figure -- 4 panels, 3 rows (research-backed; see
# project_module2_figure_reporting_research): a full-width ExWAS banner + the
# convergence headline + breadth + known-biology.
#   a ExWAS Manhattan (full-width banner) -> one point per exposure x protein pair
#        (E-block joint F-test p_E_block, 22,240 replicated), UNSIGNED, categories
#        ordered by impact, lead protein/category (white-halo labels) + counts
#   b exposure -> program -> tissue (full-width tripartite flow) -> exemplar
#        exposures (one true direction each) routed through biological program
#        clusters to organ systems; exposure->program edge colour = direction,
#        program->tissue grey = curated routing (edges precomputed upstream)
#   c category footprint   -> breadth x median|b| bubble, exemplar protein/category
#   d recovers known biology (small) -> curated canonical exposure->protein pairs
# (health-behavior axis scatter -> supplement [covered by fig_health_behavior_arms];
#  effect-by-category, GxE hubs, alcohol dose-response, eVar scatter -> supplement.)
#
# Run: module load gcc/14.2.0 R/4.4.2; export HEAP_PATHS_FILE=.../00_paths.R
#      Rscript scripts/visualizations/build_module2_fig3_composite.R [covarType] [experiment]
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]; if (is.na(common)) stop("cannot locate common/")
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers","program_clusters"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(png); library(grid); library(ggrepel); library(ggnewscale) })

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
covarType  <- if (length(a) >= 1) a[1] else "base"
experiment <- if (length(a) >= 2) a[2] else "M2_base_main"

BS <- 7; TTL <- 8.5; LBL <- 1.95
theme_cell <- function(base = BS)
  theme_heap(base_size = base) +
  theme(panel.grid = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(face = "plain"),
        plot.title = element_text(size = TTL, hjust = 0.5, face = "bold"),
        plot.title.position = "panel", plot.subtitle = element_text(size = 6, colour = "grey35", hjust = 0.5),
        plot.margin = margin(3, 3, 2, 3), legend.key.size = unit(0.18, "cm"),
        legend.text = element_text(size = 5.2), legend.title = element_text(size = 6))
CELLDIR <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "main/module2/fig3_cells")
dir.create(CELLDIR, recursive = TRUE, showWarnings = FALSE)
emit_cell <- function(p, id, w, h) {
  ggsave(file.path(CELLDIR, paste0(id, ".png")), p, width = w, height = h, dpi = 400, bg = "white")
  # cairo_pdf: vector output that keeps Unicode glyphs (β, ↑, ↓); base pdf() drops them
  ggsave(file.path(CELLDIR, paste0(id, ".pdf")), p, width = w, height = h, bg = "white", device = cairo_pdf)
}

# ============================================================================
# DATA
# ============================================================================
m <- load_module2_replicated(covarType, experiment = experiment)
pthr <- attr(m, "pval_thresh"); if (is.null(pthr) || !is.finite(pthr)) pthr <- 0.05 / nrow(m)
m <- m[is.finite(p_train) & is.finite(beta_train)]
eff <- m[replicated == TRUE & is.finite(se_train)]
tr <- load_module2_results(covarType, "train", experiment = experiment)
te <- load_module2_results(covarType, "test",  experiment = experiment)

# ============================================================================
# a (LEAD) category footprint: breadth x effect bubble, exemplar protein labeled
# ============================================================================
eb <- merge(tr$statFblock[, .(ID, omicID, Category, p_tr = p_E_block)],
            te$statFblock[, .(ID, omicID, p_te = p_E_block)], by = c("ID","omicID"))
eb <- eb[is.finite(p_tr) & is.finite(p_te)]; thrE <- 0.05/nrow(eb); sigE <- eb[p_tr < thrE & p_te < thrE]
ALLOW <- c("MMP12","APOM","LEP","FABP4","GDF15","CXCL17","CEACAM5","SELENOP","KIT","IGFBP1","IGFBP2",
           "CRP","PIGR","RETN","HGF","WFDC2","CA4","CEACAM6","IL6","CDCP1","PLAUR","NPPB","CST5","LPL")
exq <- function(ct){ d <- eff[Category == ct][order(-abs(beta_train))]; pk <- d[omicID %in% ALLOW][1]
  if (is.na(pk$omicID)) pk <- d[1]; sprintf("%s %s", pk$omicID, ifelse(pk$beta_train > 0, "↑", "↓")) }
medb <- eff[, .(med_abs = median(abs(beta_train))), by = Category]
agg_a <- sigE[, .(n_assoc = .N, n_prot = uniqueN(omicID), n_exp = uniqueN(ID)), by = Category]
agg_a <- merge(agg_a, medb, by = "Category"); agg_a[, Cat := heap_category_factor(Category)]
agg_a[, lab := sprintf("%s\n%s", heap_category_pretty(Category), sapply(Category, exq))]
ta <- nrow(sigE); tp <- uniqueN(sigE$omicID); texp <- uniqueN(sigE$ID)
pa <- ggplot(agg_a, aes(n_prot, med_abs, colour = Cat, size = n_assoc)) +
  geom_point(alpha = 0.85) +
  geom_text_repel(aes(label = lab), size = 1.95, lineheight = 0.82, show.legend = FALSE,
                  max.overlaps = Inf, box.padding = 0.5, seed = 1, segment.colour = "grey70", segment.size = 0.2) +
  scale_colour_exposure(drop = TRUE, guide = "none") +
  scale_size_continuous(range = c(1.5, 6), name = "# associations", breaks = c(1000, 4000, 7000)) +
  scale_x_log10(expand = expansion(mult = c(0.10, 0.12))) +
  labs(title = "Exposure Category Median Effect on the Proteome",
       x = "# of significant proteins", y = expression("median effect  "*abs(beta))) +
  theme_cell() + theme(legend.position = "right", legend.key.height = unit(0.16,"cm"))

# ============================================================================
# b (HEADLINE) health-behavior axis: two opposing protein arms
#   summed active days/week (well-read in Module 6) vs pack-years smoking
# ============================================================================
ppf <- function(pat) eff[grepl(pat, ID), .(bb = mean(beta_train)), by = omicID]
acts <- ppf("summed_days_activity"); smks <- ppf("pack_years_of_smoking")
dsc <- merge(setnames(copy(acts),"bb","activity"), setnames(copy(smks),"bb","smoking"), by = "omicID")
r_axis <- cor(dsc$activity, dsc$smoking)
PROT_ARM <- c("ADIPOQ","PON1","PON3","APOD","APOF","APOM","PLTP","IGFBP1","LPL")
INFL_ARM <- c("LEP","FABP4","GDF15","MMP12","IL1RN","OSM","CFH","CST3")
dsc[, arm := fifelse(omicID %in% PROT_ARM, "protective", fifelse(omicID %in% INFL_ARM, "inflammatory", "other"))]
labd <- dsc[arm != "other"]; xrs <- range(dsc$activity); yrs <- range(dsc$smoking)
pscat <- ggplot(dsc, aes(activity, smoking)) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-Inf, ymax=0, fill="#117733", alpha=0.06) +
  annotate("rect", xmin=-Inf, xmax=0, ymin=0, ymax=Inf, fill="#CC3344", alpha=0.06) +
  geom_hline(yintercept=0, colour="grey80", linewidth=0.3) + geom_vline(xintercept=0, colour="grey80", linewidth=0.3) +
  geom_abline(slope=-1, intercept=0, linetype="22", colour="grey65", linewidth=0.3) +
  geom_point(colour="grey62", size=0.45, alpha=0.5) +
  geom_point(data=labd, aes(colour=arm), size=1.1) +
  ggrepel::geom_text_repel(data=labd, aes(label=omicID, colour=arm), size=LBL-0.15, fontface="italic",
                           show.legend=FALSE, max.overlaps=Inf, seed=1, box.padding=0.28, segment.size=0.2) +
  scale_colour_manual(values=c(protective="#117733", inflammatory="#B2182B"), guide="none") +
  annotate("text", x=xrs[2]*0.98, y=yrs[1]*0.84, hjust=1, label="↑ activity, ↓ smoking\n(protective)", colour="#117733", size=1.75, fontface="bold", lineheight=0.85) +
  annotate("text", x=xrs[1]*0.98, y=yrs[2]*0.95, hjust=0, label="↓ activity, ↑ smoking\n(inflammation)", colour="#B2182B", size=1.75, fontface="bold", lineheight=0.85) +
  labs(title="Health-behavior axis: two opposing arms",
       subtitle=sprintf("per protein: more active days/week (x) vs more smoking (y); r=%.2f, ~all opposite", r_axis),
       x="effect of more active days/week (β)", y="effect of more smoking (β)") +
  theme_cell() + theme(plot.subtitle = element_text(size = 5.4, colour = "grey35", hjust = 0))

# (effect-size-by-category panel demoted to supplement: fig_effectsize_category)

# ============================================================================
# c recovers known biology
# ============================================================================
CUR <- rbindlist(list(
  data.table(field="pack_years_of_smoking", protein=c("GDF15","CXCL17","CEACAM5","MMP12"), expo="Smoking", bio="lung/stress"),
  data.table(field="alcohol_intake_frequency", protein="APOM", expo="Alcohol", bio="HDL/liver"),
  data.table(field="vigorous_physical_activity|summed_days_activity", protein=c("LEP","FABP4"), expo="Activity", bio="adiposity"),
  data.table(field="oily_fish_intake", protein="SELENOP", expo="Oily fish", bio="selenium")), fill=TRUE)
# The pair is SELECTED on the training fit but PLOTTED on the held-out test fit:
# these are the effects the reader should judge, and their SEs are the honest ones
# (the test split is ~1/4 the size, so the intervals are ~2x the training bars).
sel <- rbindlist(lapply(seq_len(nrow(CUR)), function(i){ r <- CUR[i]; d <- eff[omicID==r$protein & grepl(r$field, ID)]
  if (!nrow(d)) return(NULL); d <- d[which.max(abs(beta_train))]
  data.table(expo=r$expo, bio=r$bio, protein=r$protein, Category=as.character(d$Category),
             beta=d$beta_test, se=d$se_test) }))
sel[, `:=`(lo = beta-1.96*se, hi = beta+1.96*se, dir = ifelse(beta>0,"↑","↓"))]
sel[, strip := factor(sprintf("%s · %s", expo, bio), levels = unique(sprintf("%s · %s", CUR$expo, CUR$bio)))]
setorder(sel, strip, beta); sel[, ylab := factor(protein, levels = unique(protein))]
# coloured by EXPOSURE CATEGORY (the palette carrying a and b), not by sign --
# the sign is already legible from which side of zero the point sits on.
pc <- ggplot(sel, aes(beta, ylab, colour = Category)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
  # Smoking's held-out CI is ~0.03 wide against an axis spanning ~1.6, so a 1.9pt
  # marker swallowed it whole. Smaller point, bar drawn over it, so a tight
  # interval still shows either side instead of vanishing. At 1.1 the marker was
  # 3.1pt -- exactly the smoking CI's width -- so it had to go below that.
  geom_point(size = 0.85) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.26, linewidth = 0.4) +
  facet_grid(strip ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = HEAP_ECAT_COLORS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(title = "Protein Effects by Exposure Category",
       x = expression(beta*" (held-out test set)"), y = NULL) +
  theme_cell() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", size = 5.6),
        axis.text.y = element_text(face = "italic", size = 6), panel.spacing.y = unit(2, "pt"),
        plot.title.position = "plot", plot.title = element_text(size = TTL, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 5.0, colour = "grey35", hjust = 0.5))

# (pathway-themes, tissue-themes & GxE-hub panels demoted to supplement:
#  fig_pathway_themes, fig_tissue_themes, fig_gxe_summary / fig_gxe_composite)

# ============================================================================
# a (ANCHOR, full-width banner) entire-proteome ExWAS Manhattan
#   ONE point per (exposure field, protein) PAIR via the E-block joint F-test
#   (p_E_block; matches the 22,240 headline). UNSIGNED: the F-test is
#   directionless (and smoking-type categories mix oppositely-coded fields, so a
#   single sign is ill-defined) -- direction lives in the axis (b) and known-bio
#   (d). x = exposures grouped by category, categories ORDERED BY IMPACT
#   (# replicated). Recognizable lead protein per category + per-category counts.
# ============================================================================
MCAP <- 150
fb <- merge(tr$statFblock[, .(Eid, omicID, Category, pE_tr = p_E_block)],
            te$statFblock[, .(Eid, omicID, pE_te = p_E_block)], by = c("Eid", "omicID"))
fb <- fb[is.finite(pE_tr) & is.finite(pE_te)]
thrM <- 0.05/nrow(fb); fb[, rep := pE_tr < thrM & pE_te < thrM]
fb[, Category := heap_category_factor(Category)]
fb[, mlp := pmin(-log10(pE_tr + 1e-300), MCAP)]
fb[, fnrep := sum(rep), by = Eid]
mcatord <- fb[, .(cn = sum(rep)), by = Category][cn > 0][order(-cn)]   # drop empty categories
fb <- fb[as.character(Category) %in% as.character(mcatord$Category)]
fb[, Category := factor(as.character(Category), levels = as.character(mcatord$Category))]
setorder(fb, Category, -fnrep, Eid); fb[, fidx := as.integer(factor(Eid, levels = unique(Eid)))]
mxax <- fb[, .(mid = mean(range(fidx)), right = max(fidx)+0.5, nrep = sum(rep)), by = Category]; setorder(mxax, mid)
mbonf <- -log10(thrM + 1e-300)
mlead <- rbindlist(lapply(as.character(mcatord$Category), function(ct){
  d0 <- eff[as.character(Category)==ct][order(-abs(beta_train))]; pk <- d0[omicID %in% ALLOW][1]$omicID
  if (is.na(pk)) pk <- d0[1]$omicID
  d <- fb[as.character(Category)==ct & omicID==pk & rep==TRUE]; if (!nrow(d)) return(NULL)
  d[which.min(pE_tr), .(Category=ct, fidx, mlp, lead=pk)] }))
mlead[, Category := factor(Category, levels = levels(fb$Category))]
pmanh <- ggplot(fb, aes(fidx, mlp)) +
  geom_vline(xintercept=head(mxax$right,-1), colour="grey92", linewidth=0.2) +
  geom_point(data=fb[rep==FALSE], colour="grey85", size=0.08, alpha=0.16, stroke=0) +
  geom_point(data=fb[rep==TRUE], aes(colour=Category), size=0.38, alpha=0.65, stroke=0) +
  geom_hline(yintercept=mbonf, linetype="dashed", colour="blue", linewidth=0.3) +
  annotate("text", x=1, y=mbonf, label="Bonferroni", vjust=-0.5, hjust=0, size=1.8, colour="blue") +
  ggrepel::geom_text_repel(data=mlead, aes(fidx, mlp, label=lead, colour=Category), size=2.0, fontface="italic",
                           seed=1, max.overlaps=Inf, box.padding=0.45, point.padding=0.2, force=2, force_pull=1.4,
                           direction="both", ylim=c(NA, MCAP*1.08), min.segment.length=0.25, segment.size=0.2, segment.colour="grey80",
                           bg.color="white", bg.r=0.14, show.legend=FALSE) +
  scale_colour_exposure(drop=FALSE, guide="none") +
  scale_x_continuous(breaks=mxax$mid, labels=heap_category_pretty(mxax$Category), expand=expansion(mult=0.02)) +
  scale_y_continuous(limits=c(0, MCAP*1.16), expand=expansion(mult=c(0,0))) +
  labs(title="Entire-proteome ExWAS",
       x=NULL, y=expression(-log[10]~"(P"[E]*")")) +
  theme_cell() +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=6), axis.ticks.x=element_blank(),
        panel.grid.major.x=element_blank(), plot.subtitle=element_text(size=5.0, colour="grey35", hjust=0))

# ============================================================================
# b (TRIPARTITE) exposure -> biological program -> tissue
#   exemplar exposures (one true direction each) -> program clusters -> organ
#   systems. exposure->program edges colour = direction; program->tissue edges
#   grey = curated routing. Edges precomputed (no analysis in plotter) by
#   scripts/analysis_summaries/module2_program_tissue_edges.R.
# ============================================================================
PTdir <- file.path(heap_project_output("module4_enrichment"), "program_tissue")
EP <- fread(file.path(PTdir, "exposure_program_edges.tsv"))
PT <- fread(file.path(PTdir, "program_tissue_edges.tsv"))
PROG <- HEAP_PROGRAM_LEVELS; TISS <- HEAP_TISSUE_LEVELS; EXord <- as.data.table(HEAP_TRIPARTITE_EXEMPLARS)
PT <- PT[clust %in% PROG & organ %in% TISS & n_exp >= 4]
EP <- EP[clust %in% PROG & lab %in% EXord$lab]
ny <- nrow(EXord); ey <- setNames(seq(ny, 1, length.out = ny), EXord$lab)
spr <- function(n) seq(ny, 1, length.out = n)
py <- setNames(spr(length(PROG)), PROG); ty <- setNames(spr(length(TISS)), TISS)
EP[, `:=`(x0 = 0, y0 = ey[lab], x1 = 1, y1 = py[clust])]
PT[, `:=`(x0 = 1, y0 = py[clust], x1 = 2, y1 = ty[organ])]
EP[, dirf := factor(dir, levels = c("up","down"))]; EP[, w := npath/max(npath)]; PT[, w := n_exp/max(n_exp)]
RED <- "#B2182B"; BLU <- "#2166AC"
plab <- c("Innate immune"="Innate\nimmune","Adaptive immune"="Adaptive\nimmune","ECM / proteoglycan"="ECM /\nproteoglycan",
  "Neuronal / synaptic"="Neuronal /\nsynaptic","Growth-factor / RTK"="Growth-factor\n/ RTK",
  "Glycan / lipid metabolism"="Glycan / lipid\nmetabolism","Vascular / hemostasis / RAAS"="Vascular /\nhemostasis / RAAS","Muscle"="Muscle")
enode <- data.table(x = 0, y = ey[EXord$lab], lab = EXord$lab, col = ifelse(EXord$grp == "Harmful", RED, BLU))
pnode <- data.table(x = 1, y = py[PROG], lab = plab[PROG]); tnode <- data.table(x = 2, y = ty[TISS], lab = TISS)
# width key (bottom strip): flared wedge, thin->thick = fewer->more (counts); one per edge type
wedge <- function(x0, x1, yc, h0, h1)
  data.table(x = c(x0, x0, x1, x1), y = c(yc - h0, yc + h0, yc + h1, yc - h1))
# ONE key, right of the tissue column: direction (colour) above width (counts),
# so the reader looks in a single place instead of a ggplot legend on the right
# and two hand-drawn wedges under the flow.
KX <- 2.46; KWD <- 0.30                       # key column x, glyph width
KYa <- ny*0.60; KYb <- KYa - 0.52             # the two width wedges
KYu <- KYa + 1.18; KYd <- KYa + 0.62          # the two colour swatches
wEP <- wedge(KX, KX+KWD, KYa, 0.012, 0.10); wPT <- wedge(KX, KX+KWD, KYb, 0.012, 0.10)
kcol <- data.table(x = KX, xe = KX+KWD, y = c(KYu, KYd), col = c(RED, BLU),
                   lab = c("increased", "decreased"))
ppath <- ggplot() +
  geom_curve(data = PT, aes(x0, y0, xend = x1, yend = y1, linewidth = w), colour = "grey80", alpha = .55, curvature = -0.13, lineend = "round") +
  geom_curve(data = EP, aes(x0, y0, xend = x1, yend = y1, linewidth = w, colour = dirf), alpha = .9, curvature = -0.13, lineend = "round") +
  geom_point(data = tnode, aes(x, y), size = 1.6, colour = "grey20") +
  geom_point(data = enode, aes(x, y), size = 1.6, colour = enode$col) +
  geom_text(data = enode, aes(x - 0.05, y, label = lab), hjust = 1, size = LBL, fontface = "bold", colour = enode$col) +
  geom_label(data = pnode, aes(x, y, label = lab), size = LBL - 0.25, fontface = "bold", label.size = NA,
             fill = "white", lineheight = 0.82, label.padding = unit(0.05, "lines")) +
  geom_text(data = tnode, aes(x + 0.05, y, label = lab), hjust = 0, size = LBL, fontface = "bold") +
  annotate("text", x = 0, y = ny + 0.85, label = "EXPOSURE", size = 1.95, fontface = "italic", colour = "grey45") +
  annotate("text", x = 1, y = ny + 0.85, label = "PROGRAM", size = 1.95, fontface = "italic", colour = "grey45") +
  annotate("text", x = 2, y = ny + 0.85, label = "TISSUE", size = 1.95, fontface = "italic", colour = "grey45") +
  geom_polygon(data = wEP, aes(x, y), fill = "grey25") +
  geom_polygon(data = wPT, aes(x, y), fill = "grey72") +
  geom_segment(data = kcol, aes(x = x, xend = xe, y = y, yend = y), colour = kcol$col, linewidth = 1.5, lineend = "round") +
  geom_text(data = kcol, aes(x = xe + 0.05, y = y, label = lab), hjust = 0, size = LBL - 0.05, colour = kcol$col) +
  annotate("text", x = KX+KWD+0.05, y = KYa, hjust = 0, size = LBL - 0.05, colour = "grey25",
           label = "exposure → program") +
  annotate("text", x = KX+KWD+0.05, y = KYb, hjust = 0, size = LBL - 0.05, colour = "grey45",
           label = "program → tissue") +
  annotate("text", x = KX, y = KYa + 0.21, hjust = 0, size = LBL - 0.45, colour = "grey55", label = "fewer") +
  annotate("text", x = KX+KWD, y = KYa + 0.21, hjust = 1, size = LBL - 0.45, colour = "grey55", label = "more") +
  scale_colour_manual(values = c(up = RED, down = BLU), guide = "none") +
  scale_linewidth(range = c(0.3, 2.4), guide = "none") +
  scale_x_continuous(limits = c(-0.92, 3.28)) + scale_y_continuous(limits = c(0.45, ny + 1.15)) +
  labs(title = "Pathway and Tissue Enrichment of Exposure Associations") +
  theme_void(base_size = BS) +
  theme(plot.title = element_text(size = TTL, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 5.0, colour = "grey35", hjust = 0.5),
        legend.position = "none",
        plot.margin = margin(3, 3, 2, 3))

# ============================================================================
# RENDER CELLS  (a ExWAS Manhattan / b pathway+tissue / c bubble | d known-biology)
#   axis scatter (-> supplement: fig_health_behavior_arms), effect-by-cat, GxE -> supp
# ============================================================================
FULLW<-8.94; cW<-5.40; dW<-3.30; rA<-2.80; rB<-3.70; rC<-3.00
emit_cell(pmanh,"a_exwas_manhattan",FULLW,rA)
emit_cell(ppath,"b_pathtissue",FULLW,rB)
emit_cell(pa,"c_category_bubble",cW,rC);  emit_cell(pc,"d_known_biology",dW,rC)

# ============================================================================
# COMPOSE
# ============================================================================
trim_white <- function(im, thr=0.992, pad=2L) {
  d<-dim(im); w<-if(length(d)==3)(im[,,1]>=thr)&(im[,,2]>=thr)&(im[,,3]>=thr) else im>=thr
  rk<-which(rowSums(!w)>0); ck<-which(colSums(!w)>0); if(!length(rk)||!length(ck)) return(im)
  r0<-max(1,min(rk)-pad);r1<-min(d[1],max(rk)+pad);c0<-max(1,min(ck)-pad);c1<-min(d[2],max(ck)+pad)
  if(length(d)==3) im[r0:r1,c0:c1,,drop=FALSE] else im[r0:r1,c0:c1,drop=FALSE]
}
rd <- function(id) trim_white(readPNG(file.path(CELLDIR, paste0(id,".png"))))
# layout (a,c,d,b content order): row1 Manhattan (full) / row2 bubble | known-biology /
# row3 pathway+tissue (full). Panels are RE-LETTERED by position: a=Manhattan,
# b=category bubble, c=known-biology, d=pathway/tissue. Cell file names keep their
# original content tags (a_exwas_manhattan etc.); only position + on-figure letter move.
W<-9.20; mL<-0.16; mR<-0.10; mT<-0.10; mB<-0.10; gx2<-0.24; gy<-0.30
usable<-W-mL-mR; cW<-5.40; dW<-usable-gx2-cW
rManh<-2.80; rPair<-3.00; rPath<-3.70; Htot<-mT+rManh+rPair+rPath+2*gy+mB
yManh<-Htot-mT-rManh; yPair<-yManh-gy-rPair; yPath<-yPair-gy-rPath
cells <- list(
  list(im=rd("a_exwas_manhattan"),x=mL,        y=yManh,w=usable,h=rManh, L="a"),
  list(im=rd("c_category_bubble"),x=mL,        y=yPair,w=cW,    h=rPair, L="b"),
  list(im=rd("d_known_biology"),  x=mL+cW+gx2, y=yPair,w=dW,    h=rPair, L="c"),
  list(im=rd("b_pathtissue"),     x=mL,        y=yPath,w=usable,h=rPath, L="d"))
place <- function(p) {
  pushViewport(viewport(x=unit(p$x+p$w/2,"in"), y=unit(p$y+p$h/2,"in"), width=unit(p$w,"in"), height=unit(p$h,"in")))
  grid.raster(p$im, interpolate=TRUE); popViewport()
  grid.text(p$L, x=unit(p$x-0.02,"in"), y=unit(p$y+p$h,"in"), just=c("right","top"),
            gp=gpar(fontsize=11, fontface="bold", fontfamily="sans", col="#111111"))
}
draw <- function(){ grid.newpage(); grid.rect(gp=gpar(fill="white",col=NA)); invisible(lapply(cells,place)) }
OUT <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "main/module2")
png(file.path(OUT,"fig_exposomic_composite.png"), width=W, height=Htot, units="in", res=400, type="cairo", bg="white"); draw(); invisible(dev.off())
cairo_pdf(file.path(OUT,"fig_exposomic_composite.pdf"), width=W, height=Htot); draw(); invisible(dev.off())
cat(sprintf("Module2 composite %.2f x %.2f in\n", W, Htot)); cat("wrote", file.path(OUT,"fig_exposomic_composite.{png,pdf}"), "\n")
