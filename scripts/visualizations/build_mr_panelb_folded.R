#!/usr/bin/env Rscript
# Panel b (folded): the A–E motif EVIDENCE TABLE merged with the triad-FREQUENCY
# bars. Left = the ✓/✗/○ signature over the 6 edges (the table from the old
# schematic, edge numbering matches panel a); right = how many E–P–D triads show
# each motif. One panel does the motif key + the frequencies.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
OUT <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module5")
source("scripts/visualizations/common/figure_paths.R")
sd <- heap_resolve_output(file.path("mr_edges","summary"), must_exist=TRUE)
# ARM SCOPE: an edge touching the protein is evaluated WITHIN its pQTL platform, so a
# motif is assembled within a panel and is Tier 1 when either platform supports it.
# Read the counts summarize_mr_triads.R already computed under that rule instead of
# re-deriving it -- pooling the two arms' edge sets here is a DIFFERENT rule that lets
# a deCODE reverse edge veto a UKB triad, which drops FURIN from the mediator bar while
# panel (e) still draws it.
cf <- file.path(heap_path(), "docs", "manuscript_stats", "module5", "mr_motif_counts.tsv")
if (!file.exists(cf))
  stop("missing ", cf, "\nRun: Rscript scripts/analysis_summaries/summarize_mr_triads.R")
cnt <- fread(cf)
KEY <- c(`A Mediator (E->P->D)`="A", `B Biomarker`="B", `C Exposure-marker`="C",
         `D Reverse (P->E)`="D", `E Disease-liability (D->P)`="E")
cnt[, motif := KEY[motif]]
tab <- cnt[, .(motif, n_triads = tier1_triads, n_prot = tier1_proteins)]
setorder(tab, motif)
tab[, name := c(A="mediator",B="biomarker",C="exposure-marker",D="protein→exposure",E="disease-liability")[motif]]
COL <- c(A="#1A6B30",B="#2C7FB8",C="#3FA66A",D="#9E77B0",E="#D95F0E")
LEV <- c("E","D","C","B","A")                       # A at top (ggplot y is bottom-up)
mlab <- setNames(sprintf("%s  %s", tab$motif, tab$name), tab$motif)
tab[, motif := factor(motif, levels=LEV)]

# ---- evidence signatures over edges 1..6 (2=yes ✓, 1=either ○, 0=no ✗) -----
MAT <- list(A=c(2,2,2,0,0,0), B=c(2,0,2,0,2,0), C=c(2,0,2,0,0,1),
            D=c(0,1,1,2,1,1), E=c(1,1,1,1,2,2))
ev <- rbindlist(lapply(names(MAT), function(m) data.table(motif=m, edge=1:6, s=MAT[[m]])))
ev[, state := c("0"="no","1"="either","2"="yes")[as.character(s)]]
ev[, motif := factor(motif, levels=LEV)]
ev[, col := COL[as.character(motif)]]
edge_names <- c("E→P","P→D","E→D","P→E","D→P","D→E")

mat <- ggplot(ev, aes(edge, motif)) +
  geom_point(data=ev[state=="either"], shape=1, colour="grey60", size=2.7, stroke=0.7) +
  geom_point(data=ev[state=="no"],     shape=4, colour="grey78", size=2.0, stroke=0.8) +
  geom_point(data=ev[state=="yes"],    aes(colour=col), shape=16, size=3.3) +
  scale_colour_identity() +
  scale_x_continuous(breaks=1:6, labels=sprintf("%d\n%s", 1:6, edge_names),
                     position="top", limits=c(0.5,6.5), expand=c(0,0)) +
  scale_y_discrete(labels=mlab, expand=expansion(add=0.42)) +
  labs(x=NULL, y=NULL) +
  theme_minimal(base_size=10) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.text.x.top=element_text(size=8.6, lineheight=0.85, colour="grey35"),
        axis.text.y=element_text(size=9.5, face="bold", hjust=0))

bar <- ggplot(tab, aes(n_triads, motif, fill=motif)) +
  geom_col(width=0.62) +
  geom_text(aes(label=sprintf("%s  (%d proteins)", formatC(n_triads,big.mark=",",format="d"), n_prot)),
            hjust=-0.06, size=3.2, colour="grey25") +
  scale_fill_manual(values=COL, guide="none") +
  scale_x_log10(expand=expansion(mult=c(0,2.1)), breaks=c(100,10000), labels=c("100","10k")) +
  scale_y_discrete(expand=expansion(add=0.55)) +
  labs(x="# E–P–D triads", y=NULL) +
  theme_minimal(base_size=10) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.line.x=element_line(colour="grey40", linewidth=0.3),
        axis.line.y=element_line(colour="grey40", linewidth=0.3),
        axis.ticks.x=element_line(colour="grey45", linewidth=0.3), axis.ticks.length=unit(2.2,"pt"),
        axis.text.y=element_blank(), axis.title.x=element_text(size=8))

pb <- (mat | bar) + plot_layout(widths=c(1.28,1.12)) +
  plot_annotation(title="Tier 1 MR Motif E–P–D Triads",
                  theme=theme(plot.title=element_text(size=10.5, face="bold", hjust=0.5),
                              plot.margin=margin(2,2,1,2)))
ggsave(file.path(OUT,"_panel_b_folded.png"), pb, width=7.0, height=2.46, dpi=300, bg="white")
ggsave(file.path(OUT,"_panel_b_folded.pdf"), pb, width=7.0, height=2.46, bg="white", device=cairo_pdf)
cat("wrote _panel_b_folded.png/.pdf\n")
