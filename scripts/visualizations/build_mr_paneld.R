#!/usr/bin/env Rscript
# Panel d (final): proteins are disease RESPONDERS, not drivers. Diverging bars,
# one per top protein: orange (left) = # diseases that causally raise it (D→P,
# responder); green (right) = # diseases it causally drives (P→D, driver). The
# left bars are long, the right bars near-zero. png + vector pdf for the composite.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
OUT <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module5")
source("scripts/visualizations/common/figure_paths.R"); source("scripts/visualizations/common/plot_theme.R")
sd <- heap_resolve_output(file.path("mr_edges","summary"), must_exist=TRUE)
# ARM SCOPE: an edge touching the protein is evaluated within its pQTL platform;
# Tier 1 = supported on EITHER platform, Tier 1+ = on both (see macros/numbers.tex).
ed <- fread(file.path(sd,"mr_tiered_edges.tsv"))[mr_tier_final %in% c("Tier1","Tier1plus")]
# COUNT CONDITIONS, NOT ENDPOINTS. FinnGen ships nested codings of the same
# condition -- obesity appears as E4_OBESITY / E4_OBESITYCAL / E4_OBESITYNAS and
# type 2 diabetes as T2D / T2D_WIDE -- so an endpoint count rewards a protein for
# how finely its condition happens to be subdivided. 77% of the 492 responder
# proteins shrink when the families are collapsed (median 3 endpoints -> 1
# condition), and it reaches the driver side too: PRSS8's "two diseases" are
# E4_OBESITY and E4_OBESITYCAL, which is obesity twice.
collapse_dz <- function(x) {
  x <- sub("^finngen_R12_", "", x)
  fcase(grepl("^E4_OBESITY", x),    "obesity",
        x %in% c("T2D", "T2D_WIDE"), "type 2 diabetes",
        grepl("^I9_HEARTFAIL", x),   "heart failure",
        default = x)
}
dp <- ed[edge_dir=="D_to_P", .(n_resp=uniqueN(collapse_dz(src_id))), by=.(Protein=tgt_id)]
pd <- ed[edge_dir %in% c("Pcis_to_D","Ptrans_to_D"),
         .(n_drive=uniqueN(collapse_dz(tgt_id))), by=.(Protein=src_id)]
b <- merge(dp, pd, by="Protein", all=TRUE); for(c in c("n_resp","n_drive")) b[is.na(get(c)),(c):=0]
n_resp_only <- b[n_resp>0 & n_drive==0,.N]; n_drive <- b[n_drive>0,.N]
# ALWAYS carry the drivers. Ranking on (n_resp + n_drive) alone drew only PRSS8 and
# PCSK9 of the six, leaving ASGR1, ADM and FURIN -- the three mediators the section
# is built around -- off the panel the text cites for the driver count. The cast is
# now the six drivers plus the strongest responders, so every protein named in the
# Results is visible here.
# Carry every protein the Results names. The drivers, because the text quotes
# their count; and GDF15 / LGALS9 / IL1RN / HGF, because the sentence after it
# calls them broad disease-responsive proteins -- under the old endpoint ranking
# three of those four were absent from the panel that sentence cites (GDF15 ranked
# 52nd, IL1RN 67th, HGF 13th). The remaining slots go to the broadest responders.
# Ranking on endpoints also left the cast to an alphabetical tie-break: 23 proteins
# tied at >=6 endpoints and the panel took whichever 8 came first in table order.
drivers  <- b[n_drive>0]$Protein
EXEMPLAR <- c("GDF15", "LGALS9", "IL1RN", "HGF")
named    <- intersect(EXEMPLAR, b$Protein)
fill     <- b[!Protein %in% c(drivers, named)][order(-n_resp)][1:5]$Protein
top      <- b[Protein %in% c(drivers, named, fill)]
m <- rbind(top[, .(Protein, n=-n_resp, dir="responder (disease→protein)")],
           top[, .(Protein, n= n_drive, dir="driver (protein→disease)")])
plev <- top[order(n_resp-n_drive)]$Protein
m[, Protein := factor(Protein, levels=plev)]

p <- ggplot(m, aes(n, Protein, fill=dir)) +
  geom_col(width=0.72) +
  geom_vline(xintercept=0, colour="grey40", linewidth=0.5) +
  annotate("text", x=-3.4, y=0.25, label="←  # diseases driving", size=2.3, colour="grey30", hjust=0.5) +
  annotate("text", x=0.25,  y=0.25, label="# diseases driven  →", size=2.3, colour="grey30", hjust=0) +
  scale_fill_manual(values=c("responder (disease→protein)"="#D95F0E","driver (protein→disease)"="#1A6B30"), name=NULL) +
  scale_x_continuous(labels=abs, breaks=seq(-6,2,2), limits=c(-7.2,3.4)) +
  scale_y_discrete(expand=expansion(add=c(1.35,0.6))) +
  labs(title="Tier-1 MR links per protein", x=NULL, y=NULL) +
  theme_heap(base_size=8) +
  theme(legend.position="top", legend.key.size=unit(7,"pt"), legend.text=element_text(size=6.3),
        legend.margin=margin(0,0,-2,0), panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=7.5,face="bold"), axis.title.x=element_text(size=6.5, face="plain"),
        plot.title=element_text(size=11.5,face="bold",hjust=0.5,lineheight=0.9),
        plot.margin=margin(2,3,2,2))
# height matched to panels c & e so all three scale by the same factor in the composite
ggsave(file.path(OUT,"_panel_d.png"), p, width=3.5, height=3.3, dpi=300, bg="white")
ggsave(file.path(OUT,"_panel_d.pdf"), p, width=3.5, height=3.3, bg="white", device=cairo_pdf)
cat(sprintf("wrote _panel_d.png/.pdf | responders=%d drivers=%d\n", n_resp_only, n_drive))
