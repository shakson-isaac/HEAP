#!/usr/bin/env Rscript
# Panel c (final): top exposome-responsive proteins as category-stacked bars.
# Each protein (ordered by PXS R², shown in the row label) -> a horizontal bar
# whose length = # distinct exposures with a Tier-1 E→P MR edge, segmented and
# colored by exposure category. Shows which proteins the exposome drives and the
# category mix behind each. png + vector pdf for the composite.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
OUT <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module5")
for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers"))
  source(file.path("scripts/visualizations/common", paste0(f, ".R")))
sd <- heap_resolve_output(file.path("mr_edges","summary"), must_exist=TRUE)

co <- load_module1_predictive_r2(covarType="base",method="lasso",level="coarse",experiment="M1_base_lasso")
eR2 <- co[block=="E", .(eR2=mean(r2,na.rm=TRUE)), by=.(Protein=omic)][is.finite(eR2)]
# ARM SCOPE: an edge touching the protein is evaluated within its pQTL platform;
# Tier 1 = supported on EITHER platform, Tier 1+ = on both (see macros/numbers.tex).
ed <- fread(file.path(sd,"mr_tiered_edges.tsv"))[edge_dir=="E_to_P" &
                                                 mr_tier_final %in% c("Tier1","Tier1plus")]
cm <- unique(fread(file.path(sd,"MRmotifs.tsv"), select=c("Exposure","ExposureCategory")))
catmap <- setNames(cm$ExposureCategory, cm$Exposure)
# rank the top exposome-responsive proteins that actually carry E->P edges
captured <- unique(ed$tgt_id)
top <- eR2[Protein %in% captured][order(-eR2)][1:12]
ep <- ed[tgt_id %in% top$Protein, .(Protein=tgt_id, exposure=src_id)]
ep[, cat := catmap[exposure]]; ep <- ep[!is.na(cat)]
catord <- ep[, .(n=uniqueN(exposure)), by=cat][order(-n)]
CATS <- catord$cat
cnt <- ep[, .(n=uniqueN(exposure)), by=.(Protein,cat)]
ABBR <- c(Diet_Weekly="Diet", Smoking="Smoking", Deprivation_Indices="SES",
          Exercise_Freq="Exercise", Exercise_MET="Activity", Sexual_Factors="Sexual",
          Alcohol="Alcohol", Sun_Exposure="Sun", Sleep="Sleep")
ablab <- function(x) ifelse(x %in% names(ABBR), ABBR[x], gsub("_"," ",x))
plev <- top[order(eR2)]$Protein
top[, Protein := factor(Protein, levels=plev)]
cnt[, Protein := factor(Protein, levels=plev)][, cat := factor(cat, levels=CATS)]
ylab <- setNames(sprintf("%s  (%.2f)", top$Protein, top$eR2), as.character(top$Protein))

p <- ggplot(cnt, aes(n, Protein, fill=cat)) +
  geom_col(width=0.74, colour="white", linewidth=0.2) +
  scale_fill_manual(values=HEAP_ECAT_COLORS, labels=ablab, name=NULL) +
  scale_y_discrete(labels=ylab) +
  scale_x_continuous(expand=expansion(mult=c(0,0.04))) +
  labs(title="Exposome-responsive proteins\nwith MR E→P hits",
       x="# exposures with a Tier-1 E→P edge", y="protein  (PXS R²)") +
  theme_heap(base_size=8) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=7.5, face="bold"), axis.title=element_text(face="plain", size=7),
        plot.title=element_text(size=11.5, face="bold", hjust=0.5, lineheight=0.92),
        legend.position="top", legend.direction="horizontal", legend.justification="left",
        legend.key.size=unit(6,"pt"), legend.text=element_text(size=5.6), legend.title=element_blank(),
        legend.box.margin=margin(-2,0,-3,0), legend.spacing.x=unit(1.5,"pt"),
        legend.margin=margin(0,0,0,0), plot.margin=margin(2,2,2,2))

# c & d authored at the same height so their text scales consistently in the composite
ggsave(file.path(OUT,"_panel_c.png"), p, width=3.9, height=3.3, dpi=300, bg="white")
ggsave(file.path(OUT,"_panel_c.pdf"), p, width=3.9, height=3.3, bg="white", device=cairo_pdf)
cat("wrote _panel_c.png/.pdf\n")
