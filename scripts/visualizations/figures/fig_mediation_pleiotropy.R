#!/usr/bin/env Rscript
# fig_mediation_pleiotropy.R  [figure_id: fig_mediation_pleiotropy]
# The mediation spectrum: every mediator protein placed by disease-PLEIOTROPY (#
# diseases it mediates, x) vs its strongest mediated effect (y). The gradient runs
# from many DISEASE-SPECIFIC mediators (left) to a few PLEIOTROPIC shared reporters
# (right, e.g. LEP/ADM/GDF15). No MR -- pure mediation characterization.
# Reads disease_mediators.tsv (module3_disease_specificity.R).
local({
  cm <- c(file.path(getwd(),"scripts","visualizations","common"),
          "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cm[dir.exists(cm)][1]
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers")) source(file.path(cm,paste0(f,".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })
figure_id <- "fig_mediation_pleiotropy"; CELL <- nzchar(Sys.getenv("HEAP_CELL")); BS <- if (CELL) 7 else 11
dm <- fread(file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3/disease_mediators.tsv"))
p <- dm[, .(pleiotropy=pleiotropy[1], max_eff=max(abs(dom_NIE-1))*100, n_exp=uniqueN(dom_cat),
            mean_PM=mean(PM[is.finite(PM) & PM>0 & PM<=1], na.rm=TRUE)), by=protID]
p[, class := fcase(pleiotropy<=3,"disease-specific (<=3)", pleiotropy>=20,"pleiotropic hub (>=20)", default="intermediate")]
p[, class := factor(class, levels=c("disease-specific (<=3)","intermediate","pleiotropic hub (>=20)"))]
HUBS <- p[pleiotropy>=20][order(-max_eff)][1:8, protID]    # notable hubs (spread by effect)
SPEC <- p[pleiotropy<=3][order(-max_eff)][1:6, protID]
p[, lab := fifelse(protID %in% c(HUBS,SPEC), protID, "")]
COL <- c("disease-specific (<=3)"="#5E3C99","intermediate"="grey70","pleiotropic hub (>=20)"="#C51B7D")
n_spec <- p[pleiotropy<=3,.N]; n_hub <- p[pleiotropy>=20,.N]

pp <- ggplot(p, aes(pleiotropy, max_eff)) +
  geom_point(aes(colour=class, size=n_exp), alpha=0.75, stroke=0) +
  geom_text_repel(aes(label=lab), size=if(CELL)2.0 else 2.5, max.overlaps=Inf, force=3,
                  min.segment.length=0, segment.colour="grey65", segment.size=0.25,
                  box.padding=0.55, point.padding=0.25, colour="grey12") +
  scale_colour_manual(values=COL, name=NULL) +
  scale_size_continuous(range=c(if(CELL)0.6 else 1, if(CELL)3.5 else 5), breaks=c(1,3,6,9), name="# driving\nexposures") +
  scale_x_continuous(trans="log10") +
  annotate("text", x=1.05, y=max(p$max_eff)*0.98, hjust=0, vjust=1, size=if(CELL)2.2 else 2.9, colour="#5E3C99",
           fontface="bold", label=sprintf("%d disease-specific\nmediators", n_spec)) +
  # Twice relocated. Over the pink cloud it was unreadable; a white-filled box
  # lifted it off the points but then covered them and the COL6A3 label. It now
  # sits in the empty top-right corner as PLAIN text, mirroring the purple count
  # at top-left -- the two class totals read as a pair, left-to-right along the
  # pleiotropy axis they describe.
  annotate("text", x=max(p$pleiotropy)*0.98, y=max(p$max_eff)*0.98, hjust=1, vjust=1,
           size=if(CELL)2.2 else 2.9, colour="#C51B7D", fontface="bold", lineheight=0.9,
           label=sprintf("%d pleiotropic\nshared reporters", n_hub)) +
  labs(title=if(CELL) NULL else "The mediation spectrum: disease-specific intermediaries to pleiotropic shared reporters",
       subtitle=if(CELL) NULL else "each mediator protein: x = # diseases it exposomically mediates (pleiotropy), y = strongest exposomic mediated effect; size = # driving exposures",
       x="Disease pleiotropy  (# diseases mediated, log)",
       y=expression(paste("Strongest exposomic mediated effect  |", NIE[E], "|  (% per SD)"))) +
  theme_heap(base_size=BS) +
  theme(plot.subtitle=element_text(size=if(CELL)6 else 8, colour="grey35"),
        axis.title=element_text(face=if(CELL)"plain" else "bold"), panel.grid.minor=element_blank(),
        legend.position="right", legend.text=element_text(size=if(CELL)5.5 else 8.5),
        legend.title=element_text(size=if(CELL)6.5 else 9.5), legend.key.size=unit(if(CELL)0.3 else 0.42,"cm"))

if (CELL) { B<-file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3/fig_mediation_pleiotropy_cell")
  ggsave(paste0(B,".png"), pp, width=5.6, height=4.4, dpi=400, bg="white"); ggsave(paste0(B,".pdf"), pp, width=5.6, height=4.4, bg="white") } else
  heap_emit_figure(pp, figure_id, data=p, category="exploratory", formats=c("pdf","png"), width=8, height=5.6, website=FALSE)
message("fig_mediation_pleiotropy done")
