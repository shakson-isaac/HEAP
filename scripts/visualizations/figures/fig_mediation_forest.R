#!/usr/bin/env Rscript
# fig_mediation_forest.R  [fig_mediation_forest]  (panel d as a FOREST plot)
# Representative protein intermediaries across four exposure->disease axes, shown as
# a forest plot: mediated effect NIE HR/SD with 95% CI, per protein, grouped by axis,
# colored by exposure, disease labeled. Reads intermediaries_forest.tsv.
local({ cm <- c(file.path(getwd(),"scripts","visualizations","common"),"/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"); cm <- cm[dir.exists(cm)][1]
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers")) source(file.path(cm,paste0(f,".R"))) })
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
CELL <- nzchar(Sys.getenv("HEAP_CELL")); BS <- if (CELL) 7 else 11; NP <- as.integer(Sys.getenv("HEAP_FOR_N", unset="4"))
d <- fread(file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3/intermediaries_forest.tsv"))
CMB <- c("Endocrine/metabolic","Renal/GU","Circulatory")
AX <- list(
  list(lab="Smoking",           cat="Smoking",       sys="Respiratory"),
  list(lab="Alcohol",           cat="Alcohol",       sys="Digestive"),
  list(lab="Physical activity", cat="Exercise_Freq", sys=CMB),
  list(lab="Diet",              cat="Diet_Weekly",   sys=CMB))
rows <- rbindlist(lapply(seq_along(AX), function(i){ a<-AX[[i]]
  x <- d[dom_cat==a$cat & system %in% a$sys][order(-aa)][seq_len(min(NP,.N))]; x[, axis:=a$lab]; x }), fill=TRUE)
DZS <- c(
  "non insulin dependent diabetes mellitus"                       = "type-2 diabetes",
  "chronic renal failure"                                         = "chronic kidney disease",
  "acute renal failure"                                           = "acute kidney injury",
  "other diseases of liver"                                       = "liver disease",
  "obesity"                                                       = "obesity",
  "heart failure"                                                 = "heart failure",
  "other disorders of kidney and ureter not elsewhere classified" = "kidney disorder",
  "respiratory failure not elsewhere classified"                  = "respiratory failure",
  "other interstitial pulmonary diseases"                         = "interstitial lung disease",
  "bronchiectasis"                                                = "bronchiectasis",
  "other disorders of peritoneum"                                 = "peritoneal disorder",
  "pleural effusion not elsewhere classified"                     = "pleural effusion",
  "asthma"                                                        = "asthma",
  "other chronic obstructive pulmonary disease"                   = "COPD",
  "disorders of lipoprotein metabolism and other lipidaemias"     = "lipid disorder",
  "cholelithiasis"                                                = "gallstones",
  "gout"                                                          = "gout")
# an unmapped disease keeps its full UKB wording rather than being cut mid-phrase;
# truncation is what produced labels like "other disorders"
rows[, dzS := ifelse(disease %in% names(DZS), DZS[disease], disease)]
if (any(!rows$disease %in% names(DZS)))
  message("  NOTE unmapped disease label(s): ",
          paste(unique(rows$disease[!rows$disease %in% names(DZS)]), collapse = "; "))
rows[, lab := sprintf("%s — %s", protID, dzS)][, catf := heap_category_pretty(dom_cat)]
rows[, axis_f := factor(axis, levels=sapply(AX,`[[`,"lab"))]
rows <- rows[order(axis_f, NIE_HR)]; rows[, rk:=.I][, lab_f := factor(rk, levels=rk, labels=lab)]
pal <- HEAP_ECAT_COLORS; names(pal) <- heap_category_pretty(names(pal))
xmax <- max(rows$u95)*1.01

p <- ggplot(rows, aes(NIE_HR, lab_f, colour=catf)) +
  geom_vline(xintercept=1, linetype="dashed", colour="grey55", linewidth=0.3) +
  geom_errorbarh(aes(xmin=l95, xmax=u95), height=0.28, linewidth=0.5) +
  geom_point(size=if(CELL)1.6 else 2.3) +
  facet_grid(axis_f ~ ., scales="free_y", space="free_y", switch="y", labeller=labeller(axis_f=label_wrap_gen(15))) +
  scale_colour_manual(values=pal, guide="none") +
  scale_x_continuous(trans="log10", breaks=c(0.9,1.0,1.1,1.25)) +
  coord_cartesian(xlim=c(min(rows$l95)*0.99, xmax)) +
  labs(title=if(CELL) NULL else "Exposomic mediated effects of representative intermediaries (NIE, 95% CI)",
       subtitle=if(CELL) NULL else "exposomic mediated hazard ratio per SD through each protein, dominant exposure; dashed = no effect",
       x=expression(paste("Exposomic mediated effect  ", NIE[E], "  (HR / SD, log)")), y=NULL) +
  theme_heap(base_size=BS) +
  theme(plot.subtitle=element_text(size=if(CELL)6 else 8, colour="grey35"),
        strip.placement="outside", strip.text.y.left=element_text(angle=0, face="bold", size=if(CELL)5.8 else 8),
        axis.text.y=element_text(size=if(CELL)5.8 else 8), axis.title.x=element_text(face=if(CELL)"plain" else "bold"),
        panel.grid.major.y=element_blank(), panel.grid.minor=element_blank())

if (CELL) { B<-file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3/fig_mediation_forest_cell")
  ggsave(paste0(B,".png"), p, width=5.6, height=5.0, dpi=400, bg="white"); ggsave(paste0(B,".pdf"), p, width=5.6, height=5.0, bg="white") } else
  heap_emit_figure(p,"fig_mediation_forest",data=rows[,.(axis,protID,dom_cat,disease,NIE_HR,l95,u95,PM,n_cases)],category="exploratory",formats=c("pdf","png"),width=7.5,height=6,website=FALSE)
message("forest done")
