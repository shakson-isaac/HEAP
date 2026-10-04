#!/usr/bin/env Rscript

# ============================================================================
# fig_gxe_architecture.R  [figure_id: fig_gxe_architecture]
# ----------------------------------------------------------------------------
# CONSOLIDATED Module-2 polygenic-GxE architecture. Merges:
#   fig_gxe_summary    (component split / by category / top exposures / hub proteins)
#   fig_gxe_cis_loci   (genomic Manhattan of the replicated cis-GxE loci)
# and RETIRES fig_gxe_miami (a per-category Manhattan that restated the
# by-category panel below and shipped as a 16.7 MB rasterized PDF).
#
# The claim it serves: replicated polygenic GxE is "detectable but sparse,
# largely cis-dominated, concentrated at a small number of loci (e.g. FOLR3),
# and less extensive than the exposomic main effects". Each panel is one clause:
#   a  cis-dominated        cis vs trans vs joint-only
#   b  across the exposome  which categories carry the (few) interactions
#   c  hub proteins         a handful of proteins carry most of them (FOLR3)
#   d  concentrated loci    cis-GxE maps to the protein's OWN gene, so the
#                           interactions stack into a few genomic towers
#
# NB the covariate-sensitivity of these interactions is a SEPARATE figure
# (fig_gxe_spec_sensitivity): it is a caveat, not support -- adjusting for SES
# retains only ~29% of the replicated GxE hits.
#
# Input : module2/<experiment>/<covarType>/univar_assoc_*.rds (statFblock, train+test)
#         + Olink gene-coordinate map (olink_protein_map_3k_v1.tsv)
# Output: figures/supplement/module2/fig_gxe_architecture.{pdf,png} + data tsv
#
# Authored at 6.5in = the supplement's \textwidth (scale 1.0).
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("cannot locate common/ helpers")
  for (f in c("figure_paths", "load_heap_results", "plot_theme",
              "label_helpers", "export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_gxe_architecture", "all_main", "all_supplement", "all", "website")]
covarType  <- if (length(a) >= 1) a[1] else "base"
experiment <- if (length(a) >= 2) a[2] else "M2_base_main"
figure_id  <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_gxe_architecture")
BS <- 7.5

# ---------------------------------------------------------------- data ------
tr <- load_module2_results(covarType, "train", experiment = experiment)
te <- load_module2_results(covarType, "test",  experiment = experiment)
gx <- merge(tr$statFblock[, .(ID, omicID, Category, pj_tr = p_GxE_joint,
                              pc_tr = p_GcisxE, pt_tr = p_GtrxE)],
            te$statFblock[, .(ID, omicID, pj_te = p_GxE_joint,
                              pc_te = p_GcisxE, pt_te = p_GtrxE)], by = c("ID", "omicID"))
gx <- gx[is.finite(pj_tr) & is.finite(pj_te)]
thr <- 0.05 / nrow(gx)
sig <- gx[pj_tr < thr & pj_te < thr]
sig[, cis   := is.finite(pc_tr) & is.finite(pc_te) & pc_tr < thr & pc_te < thr]
sig[, trans := is.finite(pt_tr) & is.finite(pt_te) & pt_tr < thr & pt_te < thr]
sig[, component := factor(fifelse(cis, "cis", fifelse(trans, "trans", "joint-only")),
                          levels = HEAP_GXE_LEVELS)]
message(sprintf("GxE replicated: %d pairs | cis=%d trans=%d joint-only=%d | %d exposures, %d proteins",
                nrow(sig), sum(sig$cis), sum(sig$trans), sum(!sig$cis & !sig$trans),
                uniqueN(sig$ID), uniqueN(sig$omicID)))

# ---- a: component split -----------------------------------------------------
compA <- sig[, .N, by = component]
pa <- ggplot(compA, aes(N, component, fill = component)) +
  geom_col(width = .6) + scale_fill_gxe() +
  geom_text(aes(label = N), hjust = -0.25, size = 2.3, colour = "grey15") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "# associations", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), legend.position = "none",
        plot.margin = margin(10, 4, 2, 2))

# ---- b: by exposure category ------------------------------------------------
catB <- sig[, .N, by = .(Category, component)]
ordB <- sig[, .(tot = .N), by = Category][order(tot)]
catB[, Category := factor(heap_category_pretty(Category),
                          levels = heap_category_pretty(ordB$Category))]
pb <- ggplot(catB, aes(N, Category, fill = component)) +
  geom_col(width = .68) + scale_fill_gxe() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(x = "# associations", y = NULL, fill = NULL) +
  guides(fill = guide_legend(title = NULL, nrow = 1)) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 1),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        plot.margin = margin(10, 4, 2, 2))

# ---- c: hub proteins --------------------------------------------------------
hub <- sig[, .N, by = .(omicID, component)]
ordH <- sig[, .(tot = .N), by = omicID][order(-tot)][seq_len(min(12, .N))]
hub <- hub[omicID %in% ordH$omicID]
hub[, omicID := factor(omicID, levels = rev(ordH$omicID))]
pc <- ggplot(hub, aes(N, omicID, fill = component)) +
  geom_col(width = .68) + scale_fill_gxe() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "# associations", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(), legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        plot.margin = margin(10, 4, 2, 2))

# ---- d: cis-GxE genomic loci ------------------------------------------------
cisgxe <- sig[(cis)]
cisgxe[, Category := heap_category_factor(Category)]
cisgxe[, np := -log10(pj_tr + 1e-300)]

M <- fread(file.path(heap_omicspred(), "olink_protein_map_3k_v1.tsv"))
gmap <- unique(M[, .(omicID = HGNC.symbol, chr = as.character(chr), gstart = as.numeric(gene_start))])
gmap <- gmap[chr %in% as.character(1:22)]; gmap[, chr := as.integer(chr)]
gmap <- gmap[!is.na(gstart)][order(chr, gstart)][, .SD[1], by = omicID]
d <- merge(cisgxe, gmap, by = "omicID")

chrlen <- M[chr %in% as.character(1:22),
            .(len = max(as.numeric(gene_end), na.rm = TRUE)), by = .(chr = as.integer(chr))][order(chr)]
chrlen[, offset := cumsum(as.numeric(len)) - len]
d <- merge(d, chrlen[, .(chr, offset)], by = "chr")
d[, gpos := offset + gstart]
chrlen[, mid := offset + len / 2]

# one label per locus, at its top point
lab <- d[, .SD[which.max(np)], by = omicID]

pd <- ggplot(d, aes(gpos, np, colour = Category)) +
  geom_hline(yintercept = -log10(thr), linetype = "dashed",
             colour = "grey55", linewidth = .25) +
  geom_point(size = .7, alpha = .85) +
  geom_text_repel(data = lab, aes(label = omicID), size = 1.9, fontface = "italic",
                  colour = "grey20", segment.size = .2, min.segment.length = 0,
                  max.overlaps = Inf, box.padding = .25, seed = 3, show.legend = FALSE) +
  scale_colour_exposure() +
  scale_x_continuous(breaks = chrlen$mid, labels = chrlen$chr, expand = expansion(mult = .02)) +
  labs(x = "Chromosome", y = expression(-log[10]~P~"(joint GxE)"), colour = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.key.size = unit(7, "pt"),
        legend.text = element_text(size = BS - 2),
        axis.text.x = element_text(size = BS - 2),
        plot.margin = margin(10, 4, 2, 2))

# ------------------------------------------------------------- assemble ------
p <- ((pa | pb | pc) / pd) +
  plot_layout(heights = c(0.85, 1.05)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  compA[, .(panel = "a", key = as.character(component), value = as.numeric(N))],
  catB[,  .(panel = "b", key = paste(Category, component, sep = " | "), value = as.numeric(N))],
  hub[,   .(panel = "c", key = paste(omicID, component, sep = " | "), value = as.numeric(N))],
  d[,     .(panel = "d", key = paste(omicID, ID, sep = " | "), value = np)]), use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module2",
                 formats = c("pdf", "png"), width = 6.5, height = 6.4, website = TRUE)

message("fig_gxe_architecture: done.")
