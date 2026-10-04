#!/usr/bin/env Rscript
# ============================================================================
# fig_category_biology.R  [figure_id: fig_category_biology]
# ----------------------------------------------------------------------------
# PER-CATEGORY biology of exposure-responsive proteins (Module 1 per-category
# GSEA, run_module1_category_gsea.R). For each exposure category, its TOP
# significant GO/Reactome terms as a faceted lollipop (term NAMES on the axis,
# x = -log10 p.adj), colored by biological THEME. The shared cell-surface
# signaling terms (Olink panel background) are KEPT and shown in grey, so each
# category's distinctive biology (smoking->leukocyte chemotaxis; sleep->
# gluconeogenesis/cholesterol/cortisol; diet->lipid transport; exercise->
# vascular/ECM/integrin) reads against that background. NOT diet/exercise-
# confounded like the aggregate exposomic ranking.
#
# Input : docs/manuscript_stats/module1_enrichment_bycategory/gsea_all_significant_bycategory.tsv
# Output: figures/supplement/module1/fig_category_biology.{pdf,png} (+ CELL cell.{png,pdf})
# ============================================================================
local({
  cm <- "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common"
  for (f in c("figure_paths","plot_theme","label_helpers","export_helpers"))
    source(file.path(cm, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_category_biology")
CELL <- nzchar(Sys.getenv("HEAP_CELL"))
N_TOP <- as.integer(Sys.getenv("HEAP_N_TOP", "6"))
BS <- if (CELL) 7 else 11; TTL <- if (CELL) 8.5 else 13
B_W <- if (CELL) 5.2 else 7.6; B_H <- if (CELL) 3.6 else 5.6

f <- file.path(heap_path(), "docs", "manuscript_stats", "module1_enrichment_bycategory",
               "gsea_all_significant_bycategory.tsv")
S <- fread(f)[NES > 0]

themes <- c(
  "Immune / inflammation"        = "inflammat|chemotax|leukocyte|neutrophil|cytokine|immune|opsonin|complement|interleukin|defense",
  "Vascular / angiogenesis"      = "vasculature|angiogen|blood vessel|wound",
  "ECM / adhesion"               = "extracellular matrix|collagen|encapsulating|integrin|adhesion|plasma membrane region|external side",
  "Lipid transport"              = "lipid transport|lipid localization|lipoprotein|sterol transport",
  "Glucose / sterol metabolism"  = "gluconeogen|hexose|monosacchar|glucose|sterol|cholesterol|carbohydrate|amino acid metabolic",
  "Endocrine (cortisol)"         = "glucocorticoid|corticosteroid|hormone",
  "Secretory / ER"               = "endoplasmic reticulum|secret|vesicle",
  "Cell-surface (background)" = "signaling receptor|molecular transducer|transmembrane signal|cell surface|signal transduction|receptor activity")
assign_theme <- function(desc) { d <- tolower(desc); for (t in names(themes)) if (grepl(themes[[t]], d)) return(t); "Other" }
S[, theme := vapply(Description, assign_theme, character(1))]

catord <- intersect(c("Exercise_Freq","Diet_Weekly","Alcohol","Smoking","Sleep","Exercise_MET"),
                    unique(S$category))
shrt <- function(x) { x <- gsub("_", " ", x); ifelse(nchar(x) > 42, paste0(substr(x, 1, 40), "..."), x) }
top <- S[order(category, p.adjust)][, head(.SD, N_TOP), by = category]
# theme assignment above matches on the ORIGINAL (official GO/Reactome) name; only the
# DISPLAY label is americanized, so a British source term cannot leak onto the figure.
top[, term := shrt(heap_americanize(Description))]
top[, lp := -log10(p.adjust)]
top[, category := factor(as.character(category), levels = catord)]
top[, catlab := factor(heap_category_pretty(as.character(category)), levels = heap_category_pretty(catord))]
top <- top[order(category, lp)]
top[, key := factor(paste(category, term, seq_len(.N)), levels = unique(paste(category, term, seq_len(.N)))), by = category]
top[, key := factor(paste0(as.integer(category), "|", term, "|", seq_len(.N)))]   # unique per row
top <- top[order(category, lp)]
top[, key := factor(key, levels = key)]

thlev <- c("Immune / inflammation","Vascular / angiogenesis","ECM / adhesion","Lipid transport",
           "Glucose / sterol metabolism","Endocrine (cortisol)","Secretory / ER","Other",
           "Cell-surface (background)")
thcol <- setNames(c("#B2182B","#D6604D","#3B7DB8","#2E9E48","#B35806","#762A83","#117733","#777777","#C9CDD2"), thlev)
top[, theme := factor(theme, levels = thlev)]

p <- ggplot(top, aes(lp, key, colour = theme)) +
  geom_segment(aes(x = 0, xend = lp, yend = key), linewidth = 0.5) +
  geom_point(size = if (CELL) 1.9 else 2.6, stroke = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_discrete(labels = function(k) sub("^[0-9]+\\|(.*)\\|[0-9]+$", "\\1", k)) +
  scale_colour_manual(values = thcol, drop = TRUE, name = NULL) +
  facet_wrap(~ catlab, scales = "free_y", ncol = 2) +
  labs(title = if (CELL) "Each exposure's distinctive protein biology" else "Per-category enriched protein biology",
       subtitle = NULL,
       x = expression(-log[10]~p[adj]), y = NULL) +
  theme_heap(base_size = BS) +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        strip.text = element_text(face = "bold")) +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE))

if (CELL) {
  p <- p + theme(plot.title = element_text(size = TTL, hjust = 0.5, face = "bold"),
                 axis.text.y = element_text(size = rel(0.66)), axis.title = element_text(size = rel(0.85)),
                 legend.key.size = unit(0.5, "lines"), legend.text = element_text(size = rel(0.6)),
                 plot.margin = margin(3, 4, 2, 3))
  OUT <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module1")
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(OUT, paste0(figure_id, "_cell.png")), p, width = B_W, height = B_H, dpi = 400, bg = "white")
  ggsave(file.path(OUT, paste0(figure_id, "_cell.pdf")), p, width = B_W, height = B_H, bg = "white")
} else {
  heap_emit_figure(p, figure_id, data = top, category = "supplement", formats = c("pdf","png"), width = 7.6, height = 5.6, website = TRUE)
}
message("fig_category_biology: done (", nrow(top), " rows, ", length(catord), " categories; names visualized, signaling kept in grey).")
