#!/usr/bin/env Rscript
# Compact TIER FUNNEL — MR evidence tiers; box width funnels DOWN with stringency.
# Standalone shows the criteria; CELL mode (for the composite) drops the criteria
# prose (it lives in the figure legend) and renders at the exact cell size.
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })
OUT <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module5")
CELL <- nzchar(Sys.getenv("HEAP_CELL"))

L <- data.table(
  y    = 1:4,
  tier = c("Suggestive","Tier 2","Tier 1","Tier 1+"),
  adds = c("few instruments, fails sensitivity, or direction unresolved",
           "BH-significant; trans-only or cis not colocalized",
           "+ established direction (Steiger) + sensitivity-robust",
           "+ replicated across UKB Olink & deCODE SomaScan"),
  col  = c("grey78","#9ECAE1","#4292C6","#1A6B30"))
L[, w := 0.46 + 0.62*((4-y)/3)]
XLO <- -1.2; XHI <- 0.62      # cell x range; the tiles sit at x = 0, not its midpoint
xc <- max(L$w)/2 + 0.12

p <- ggplot(L) +
  geom_tile(aes(x=0, y=y, width=w, height=0.90, fill=col)) +
  geom_text(aes(x=0, y=y, label=tier), size=if(CELL)2.3 else 3.0, fontface="bold",
            colour=ifelse(L$tier %in% c("Tier 1","Tier 1+"),"white","grey20")) +
  geom_segment(aes(x=-0.95, xend=-0.95, y=0.62, yend=4.38),
               arrow=arrow(length=unit(if(CELL)0.06 else 0.09,"in"), type="closed"), colour="grey45", linewidth=0.5) +
  annotate("text", x=-1.06, y=2.5, label="increasing stringency", angle=90, size=if(CELL)3.0 else 2.5, colour="grey45") +
  annotate("segment", x=-0.78, xend=-0.78, y=2.60, yend=4.42, colour="#1A6B30", linewidth=0.6) +
  annotate("text", x=-0.70, y=3.45, label="main findings", angle=90, size=if(CELL)3.0 else 2.5, fontface="bold", colour="#1A6B30") +
  scale_fill_identity() +
  scale_y_continuous(limits=c(0.52, 4.48), expand=c(0,0)) +
  labs(title="MR evidence tiers") +
  coord_cartesian(clip="off") +
  theme_void(base_size=10) +
  theme(plot.title=element_text(hjust=if(CELL) (0-XLO)/(XHI-XLO) else 0.5, size=if(CELL)10.6 else 10.5, face="bold", margin=margin(b=1)),
        plot.margin=margin(2,3,2,3))

if (CELL) {
  p <- p + scale_x_continuous(limits=c(XLO, XHI), expand=c(0,0))   # no criteria column
  ggsave(file.path(OUT,"_panel_tier_ladder.png"), p, width=2.6, height=1.7, dpi=300, bg="white")
  ggsave(file.path(OUT,"_panel_tier_ladder.pdf"), p, width=2.6, height=1.7, bg="white")
} else {
  p <- p + geom_text(aes(x=xc, y=y, label=adds), hjust=0, size=2.5, colour="grey30") +
       annotate("text", x=xc, y=4.55, label="Tier 1+ : protein edges only (E↔D caps at Tier 1)", hjust=0, size=2.2,
                fontface="italic", colour="#1A6B30") +
       scale_x_continuous(limits=c(-1.2, 3.05))
  ggsave(file.path(OUT,"_panel_tier_ladder.png"), p, width=6.0, height=2.7, dpi=300, bg="white")
}
cat("wrote _panel_tier_ladder.png (cell=", CELL, ")\n")
