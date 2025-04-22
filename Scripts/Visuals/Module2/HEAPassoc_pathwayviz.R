#Libraries:
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(pbapply)
library(dplyr)
library(purrr)
library(plotly)
library(htmlwidgets)
library(qs)

#For heatmaps
library(ComplexHeatmap)
library(circlize)

#'*LOAD GSEA results*
setwd("/n/groups/patel/shakson_ukb/UK_Biobank/")
HEAPht <- qread("./Output/HEAPres/HEAPgsea.qs")


#Tissue Effect Dataframe:
# purrr::map2 iterate over the list created to generate one dataframe:
TissueEffects <- map2_dfr(HEAPht@HEAPtgsea, names(HEAPht@HEAPtgsea), ~ .x %>%
                        select(all_of(c("ID", "NES","p.adjust"))) %>%
                        mutate(cID = .y))

#This dataframe can be used to confirm all associations have FDR p-value < 0.05
TE_df <- TissueEffects %>%
  select(c("ID","p.adjust","cID")) %>%
  pivot_wider(names_from = "ID",
              values_from = "p.adjust",
              values_fill = 1) %>%
  column_to_rownames(var = "cID")

TE_df_dir <- TissueEffects %>%
  select(c("ID","NES","cID")) %>%
  pivot_wider(names_from = "ID",
              values_from = "NES",
              values_fill = 0) %>%
  column_to_rownames(var = "cID")


#Pathway Effects Dataframe:
PathEffects <- map2_dfr(HEAPht@HEAPpgsea, names(HEAPht@HEAPpgsea), ~ .x %>%
                            select(all_of(c("Description", "NES","p.adjust"))) %>%
                            mutate(cID = .y))

Path_dir <- PathEffects %>%
  select(c("Description","NES","cID")) %>%
  pivot_wider(names_from = "Description",
              values_from = "NES",
              values_fill = 0) %>%
  column_to_rownames(var = "cID")




## Enrichments of All Significant Exposures and Tissues (Supplementary Section)
ht <- Heatmap(as.matrix(TE_df_dir),
              name = "NES",
              col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              #c("#440154", "#21908C", "#FDE725")
              width = ncol(TE_df_dir)*unit(2, "mm"),
              height = nrow(TE_df_dir)*unit(2, "mm"),
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_column_dend = FALSE, 
              show_row_dend = FALSE,
              #row_split = 5,
              #row_km = 3,
              #column_gap = unit(8, "mm"),
              #row_gap = unit(20, "mm"),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_names_side = "left",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10), # Legend title font size
                labels_gp = gpar(fontsize = 10), # Legend labels font size
                legend_height = unit(2, "cm"), # Height of the legend
                legend_width = unit(1, "cm")   # Width of the legend
              ),
              border = T,
              rect_gp = gpar(col = "black", lwd = 0.5)#,
              #border = T,
             # rect_gp = gpar(col = "black", lwd = 0.5)
)


# Save plot:
filename <- paste0("./Figures/HEAP/SF2/AllEtissue_ht.png")
png(file=filename,
    height = 12, width = 12, units = "in", res = 500)
draw(ht)
dev.off()


ht <- Heatmap(as.matrix(Path_dir),
              name = "NES",
              col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              #c("#440154", "#21908C", "#FDE725")
              width = ncol(Path_dir)*unit(2, "mm"),
              height = nrow(Path_dir)*unit(2, "mm"),
              show_row_names = TRUE,
              show_column_names = TRUE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_column_dend = FALSE, 
              show_row_dend = FALSE,
              #row_split = 5,
              #row_km = 3,
              #column_gap = unit(8, "mm"),
              #row_gap = unit(20, "mm"),
              row_names_gp = gpar(fontsize = 6),
              column_names_gp = gpar(fontsize = 6),
              row_names_side = "left",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10), # Legend title font size
                labels_gp = gpar(fontsize = 10), # Legend labels font size
                legend_height = unit(2, "cm"), # Height of the legend
                legend_width = unit(1, "cm")   # Width of the legend
              ),
              border = T,
              rect_gp = gpar(col = "black", lwd = 0.5)
)

# Save plot:
filename <- paste0("./Figures/HEAP/SF2/AllEpath_ht.png")
png(file=filename,
    height = 12, width = 12, units = "in", res = 500)
draw(ht)
dev.off()




## Enrichments of Specific Exposures (Main Figure)
selectIDs <- c("number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
  "fresh_fruit_intake_f1309_0_0",
  "alcohol_intake_frequency_f1558_0_06",
  "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days")
TE_df_dirv2 <- TE_df_dir %>%
                  filter(rownames(TE_df_dir) %in% selectIDs) %>%
                  filter(!if_all(everything(), ~ . == 0))
Path_dirv2 <- Path_dir %>%
                  filter(rownames(Path_dir) %in% selectIDs) %>%
                  filter(!if_all(everything(), ~ . == 0))
rownames(TE_df_dirv2)
rownames(Path_dirv2)


#'*One heatmap:*
newdir <- merge(Path_dirv2, TE_df_dirv2, by = "row.names")
rownames(newdir) <- newdir$Row.names
newdir$Row.names <- NULL
newdir <- as.data.frame(t(newdir))
colnames(newdir) <- c("Alcohol","Smoking","Diet","Exercise")

#Remove rows with all 0s:
newdir <- newdir %>%
  filter(!if_all(everything(), ~ . == 0))

HEAPanno <- data.frame(
  type = c(rep("Pathway", each = ncol(Path_dirv2)),
           rep("Tissue", each = ncol(TE_df_dirv2))),
  id = c(colnames(Path_dirv2), colnames(TE_df_dirv2))
)

HEAPanno[match(rownames(newdir),HEAPanno$id),]$type

# Row annotation using the tissue group information
row_anno <- rowAnnotation(
  tissue_group = anno_simple(HEAPanno[match(rownames(newdir),HEAPanno$id),]$type, 
                             col = c("Pathway" = "lightblue", "Tissue" = "lightgreen"))
)

#order newdir:
newdir <- newdir[,c("Exercise","Diet","Alcohol","Smoking")]

#Get row split vector before name changes:
rsplitn <- factor(HEAPanno[match(rownames(newdir), HEAPanno$id),]$type,
                  levels = c("Tissue", "Pathway"))

#Shorten names:
rownames(newdir)[rownames(newdir) == "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell"] <- "Lymphoid-Non-Lymphoid Interactions"
rownames(newdir)[rownames(newdir) == "Diseases associated with glycosaminoglycan metabolism"] <- "Glycosaminoglycan Metabolism Diseases"
rownames(newdir)[rownames(newdir) == "Platelet activation, signaling and aggregation"] <- "Platelet activation, signaling, aggregation"
rownames(newdir)[rownames(newdir) == "Metabolism of Angiotensinogen to Angiotensins"] <- "Angiotensinogen to Angiotensins Metabolism"


#'*ComplexHeatmap w/ both Tissue/Pathways*
ht <- Heatmap(as.matrix(newdir),
              name = "NES",
              col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
              #c("#440154", "#21908C", "#FDE725")
              width = ncol(newdir)*unit(3.5, "mm"),
              height = nrow(newdir)*unit(3, "mm"),
              show_row_names = TRUE,
              show_column_names = TRUE, # Adjust width if necessary
              row_names_max_width = unit(12, "cm"),
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              show_column_dend = FALSE, 
              show_row_dend = FALSE,
              row_split = rsplitn,
              #row_split = 5,
              #row_km = 3,
              #column_gap = unit(8, "mm"),
              row_gap = unit(10, "mm"),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp = gpar(fontsize = 10),
              row_names_side = "left",
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10), # Legend title font size
                labels_gp = gpar(fontsize = 10), # Legend labels font size
                legend_height = unit(2, "cm"), # Height of the legend
                legend_width = unit(1, "cm")   # Width of the legend
              ),
              border = T,
              rect_gp = gpar(col = "black", lwd = 0.5))
draw(ht)

# Save plot:
filename <- paste0("./Figures/HEAP/F2/Etissuepath_v2.png")
png(file=filename,
    height = 8, width = 4, units = "in", res = 1000)
draw(ht)
dev.off()

filename <- paste0("./Figures/HEAP/F2/Etissuepath_v2.svg")
svg(file=filename,
    height = 8, width = 4)
draw(ht)
dev.off()

filename <- paste0("./Figures/HEAP/F2/Etissuepath_v2.pdf")
pdf(file=filename,
    height = 8, width = 4)
draw(ht)
dev.off()
