#!/usr/bin/env Rscript
# fig_mediation_scale_main.R  [figure_id: fig_mediation_scale_main]
# Curated 'scale of mediation' grid for the MAIN figure: representative lifestyle
# exposures x representative diseases (spanning organ systems); cell = # proteins
# with an FDR-significant mediated effect (NIE, n>=100), number printed. The full
# all-disease heatmap lives in the supplement (fig_mediation_scale).
# Reads cat_disease_full.tsv.
local({
  cm <- c(file.path(getwd(),"scripts","visualizations","common"),
          "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cm[dir.exists(cm)][1]
  for (f in c("figure_paths","load_heap_results","plot_theme","label_helpers","export_helpers")) source(file.path(cm,paste0(f,".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
figure_id <- "fig_mediation_scale_main"; CELL <- nzchar(Sys.getenv("HEAP_CELL")); BS <- if (CELL) 7 else 11
g <- fread(file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3/cat_disease_full.tsv"))
g[, dcode := sub("^age_([a-z][0-9]+)_.*","\\1", DZ_ID)]
EXP <- c("Exercise_Freq","Diet_Weekly","Alcohol","Smoking","Sleep")          # representative exposures
DZ  <- data.table(
  dcode=c("e11","e66","e78","n18","n17","i50","i25","k76","k80","j44","j45","m10","f10","f32"),
  lab  =c("type-2 diabetes","obesity","lipid disorder","chronic kidney dis.","acute kidney inj.",
          "heart failure","ischemic heart dis.","liver disease","gallstones","COPD","asthma",
          "gout","alcohol-use disord.","depression"),
  sys  =c("Metabolic","Metabolic","Metabolic","Renal","Renal","Circulatory","Circulatory",
          "Digestive","Digestive","Respiratory","Respiratory","MSK","Mental","Mental"))
SYSL <- c("Metabolic","Renal","Circulatory","Digestive","Respiratory","MSK","Mental")
grid <- CJ(category=EXP, dcode=DZ$dcode)
sub <- merge(grid, g[, .(category,dcode,n_prot,n_strong)], by=c("category","dcode"), all.x=TRUE)
sub[is.na(n_prot), n_prot:=0][is.na(n_strong), n_strong:=0]
sub <- merge(sub, DZ, by="dcode")
sub[, exp_f := factor(heap_category_pretty(category), levels=rev(heap_category_pretty(EXP)))]
sub[, dz_f := factor(lab, levels=DZ$lab)]
sub[, sys_f := factor(sys, levels=SYSL)]

p <- ggplot(sub, aes(dz_f, exp_f, fill=n_prot)) +
  geom_tile(colour="white", linewidth=0.9) +
  geom_text(aes(label=n_prot, colour=n_prot>=110), size=if(CELL)2.2 else 3.0, fontface="bold") +
  scale_fill_viridis_c(option="mako", direction=-1, trans="sqrt", breaks=c(0,25,100,250,500),
                       name="mediator\nproteins") +
  scale_colour_manual(values=c(`TRUE`="white",`FALSE`="grey20"), guide="none") +
  facet_grid(. ~ sys_f, scales="free_x", space="free_x") +
  labs(title=if(CELL) NULL else "The scale of mediation across representative exposures and diseases",
       subtitle=if(CELL) NULL else "cell = # proteins with an FDR-significant mediated effect (NIE, n>=100); grouped by organ system",
       x=NULL, y=NULL) +
  theme_heap(base_size=BS) +
  theme(plot.subtitle=element_text(size=if(CELL)6 else 8, colour="grey35"),
        axis.text.x=element_text(angle=40, hjust=1, size=if(CELL)5.8 else 8),
        axis.text.y=element_text(size=if(CELL)6.5 else 9.5, face="bold"),
        panel.grid=element_blank(), axis.ticks=element_blank(),
        panel.spacing.x=unit(if(CELL)2 else 3,"pt"),
        strip.placement="outside", strip.background=element_rect(fill="grey92", colour=NA),
        strip.text.x=element_text(size=if(CELL)5.8 else 8, face="bold"),
        legend.position="right", legend.text=element_text(size=if(CELL)5.5 else 8),
        legend.title=element_text(size=if(CELL)6.5 else 9), legend.key.size=unit(if(CELL)0.30 else 0.45,"cm"))

if (CELL) { B<-file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module3/fig_mediation_scale_main_cell")
  ggsave(paste0(B,".png"), p, width=6.6, height=3.0, dpi=400, bg="white"); ggsave(paste0(B,".pdf"), p, width=6.6, height=3.0, bg="white") } else
  heap_emit_figure(p, figure_id, data=sub[, .(category,dcode,disease=lab,system=sys,n_prot,n_strong)], category="main", formats=c("pdf","png"), width=9.6, height=3.8, website=TRUE)
message("fig_mediation_scale_main done")
