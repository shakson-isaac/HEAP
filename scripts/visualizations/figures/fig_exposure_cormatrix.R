#!/usr/bin/env Rscript

# ============================================================================
# fig_exposure_cormatrix.R  [figure_id: fig_exposure_cormatrix]
# ----------------------------------------------------------------------------
# Supplementary QC: correlation matrix of the analysis exposures, from the
# canonical HEAP loader object (HEAP.rds). Reads exposures the same way the
# modules do — heap_filter_exposures() + as_pxs_baseline() — so the matrix
# reflects exactly the exposure set used in Modules 1/2/3.
#
# Input : heap_loader_rds (override with HEAP_LOADER_RDS during transition)
# Output: figures/supplement/fig_exposure_cormatrix.{pdf,png} + figures/data/...tsv
#
# NOTE: loads the ~9 GB HEAP.rds object; run on a node with >=12 GB RAM and not
# concurrently with another HEAP.rds-loading job.
#
# Run:
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#   HEAP_LOADER_RDS=/n/scratch/.../UKB_intermediate/HEAP.rds \
#     Rscript scripts/visualizations/figures/fig_exposure_cormatrix.R
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("cannot locate common/ helpers")
  for (f in c("figure_paths", "plot_theme", "label_helpers", "export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
  have_ggcorrplot <- requireNamespace("ggcorrplot", quietly = TRUE)
})

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_exposure_cormatrix")

# --- load canonical HEAP loader + apply the analysis exposure filter --------
if (!file.exists(heap_loader_rds))
  stop("HEAP loader object not found: ", heap_loader_rds,
       "\nSet HEAP_LOADER_RDS to the HEAP.rds path (run scripts/loaders/HEAP_loader.R first).")
message("Loading HEAP loader object (~9 GB): ", heap_loader_rds)
heap <- readRDS(heap_loader_rds)
heap <- heap_filter_exposures(heap)            # same exposure set as the modules
pxs  <- as_pxs_baseline(heap); rm(heap); gc()

# --- assemble a wide exposure matrix (eid x exposures) ----------------------
elist <- pxs$Elist
ex_wide <- Reduce(function(a, b) merge(a, b, by = "eid", all = TRUE),
                  lapply(elist, as.data.table))
ex_mat <- ex_wide[, setdiff(names(ex_wide), "eid"), with = FALSE]
# Build the SAME feature set the PXS trains on: one-hot-expand every categorical
# exposure into 0/1 dummies (Module 1 uses Matrix::sparse.model.matrix), so the
# matrix shows the FULL set of PXS INPUT FEATURES (continuous + ordinal + binary +
# per-level dummies), not just the continuous/ordinal subset (which was ~81 of ~169).
ex_df <- as.data.frame(ex_mat)
for (nm in names(ex_df)) {                          # coded categoricals -> clean factors
  if (is.character(ex_df[[nm]])) ex_df[[nm]] <- factor(ex_df[[nm]])
  if (is.factor(ex_df[[nm]]))    ex_df[[nm]] <- droplevels(ex_df[[nm]])
}
ok <- vapply(ex_df, function(x) sum(!is.na(x)) > 100 &&
               length(unique(x[!is.na(x)])) >= 2L, logical(1))
ex_df <- ex_df[, ok, drop = FALSE]
op <- options(na.action = "na.pass")
X <- model.matrix(~ . - 1, data = ex_df)            # one-hot; NA-preserving for pairwise cor
options(op)
base_exposure <- setNames(colnames(ex_df)[attr(X, "assign")], colnames(X))
message("PXS feature matrix: ", nrow(X), " participants x ", ncol(X),
        " features (one-hot expanded from ", ncol(ex_df), " base exposures)")

# Spearman = Pearson on ranks. Rank each feature once (NA-preserving), then
# Pearson pairwise. cor(method="spearman", use="pairwise.complete.obs") re-ranks
# every column PAIR and is prohibitively slow at 169 features x 500k rows; the
# rank-once approach is the standard fast equivalent.
Xr <- apply(X, 2L, rank, na.last = "keep")
cm <- cor(Xr, use = "pairwise.complete.obs")
cm[is.na(cm)] <- 0

# order features by broad category (via each feature's BASE exposure) so blocks show
ecat <- pxs$Eid_cat
if (!is.null(ecat)) {
  ecat <- as.data.table(ecat)
  setnames(ecat, names(ecat)[1:2], c("Eid", "Category"))
  fb     <- unname(base_exposure[colnames(cm)])           # feature -> base exposure id
  cat_of <- ecat$Category[match(fb, ecat$Eid)]
  o <- order(heap_broad_category(cat_of), cat_of, fb, colnames(cm), na.last = TRUE)
  cm <- cm[o, o]
}

# readable axis labels: "<base short label>[: <level>]" for one-hot dummies
relabel_feat <- function(feat) {
  base <- unname(base_exposure[feat])
  lvl  <- trimws(gsub("[._]+", " ", substring(feat, nchar(base) + 1L)))
  lab  <- heap_exposure_label(base)
  ifelse(nzchar(lvl), paste0(lab, ": ", lvl), lab)
}
dimnames(cm) <- list(relabel_feat(rownames(cm)), relabel_feat(colnames(cm)))

# --- plot (long-form tile heatmap; ggcorrplot if available) -----------------
if (have_ggcorrplot) {
  p <- ggcorrplot::ggcorrplot(cm, hc.order = FALSE, type = "full",
                              outline.color = NA) +
    scale_fill_gradient2(low = "#1B6CA8", mid = "white", high = "#E07B39",
                         midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
    labs(title = "Exposure correlation matrix",
         subtitle = NULL, x = NULL, y = NULL) +   # ggcorrplot defaults to Var1/Var2
    theme_heap() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
          axis.text.y = element_text(size = 4))
} else {
  dt <- as.data.table(as.table(cm)); setnames(dt, c("V1", "V2", "r"))
  dt[, V1 := factor(V1, levels = colnames(cm))]
  dt[, V2 := factor(V2, levels = colnames(cm))]
  p <- ggplot(dt, aes(V1, V2, fill = r)) + geom_tile() +
    scale_fill_gradient2(low = "#1B6CA8", mid = "white", high = "#E07B39",
                         midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Exposure correlation matrix",
         subtitle = NULL,
         x = NULL, y = NULL) +
    theme_heap() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
          axis.text.y = element_text(size = 4))
}

# emit; data = the correlation matrix (with rownames as a column)
cm_out <- data.table(exposure = rownames(cm), as.data.table(cm))
heap_emit_figure(p, figure_id, data = cm_out, category = "supplement",
                 formats = c("pdf", "png"), width = 9, height = 8, website = FALSE)
message("fig_exposure_cormatrix: done (", ncol(cm), " exposures).")
