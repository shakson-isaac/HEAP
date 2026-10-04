#!/usr/bin/env Rscript
# ============================================================================
# module6_pes_score_correlation.R  (support analysis)
# ----------------------------------------------------------------------------
# Builds the cross-exposure correlation matrix of proteome-only PES scores
# (pes_prot_z) from the Module 6 out-of-fold predictions (TrainOOF, instance 0).
# Writes a long-form correlation table + a hierarchical-cluster order for the
# fig_pes_score_correlation plotter (no plotting here).
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (!is.na(hit)) source(hit)
})
suppressPackageStartupMessages(library(data.table))

covarType  <- "base"
score_col  <- "pes_prot_z"     # proteome-only PES (z-scored)
in_dir     <- heap_project_output("module6_pes_longitudinal", covarType)
out_dir    <- heap_project_output("module6_pes_longitudinal", "pes_score_correlation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "_TrainOOF\\.tsv$", full.names = TRUE)
if (!length(files)) stop("No TrainOOF.tsv under ", in_dir)
get_exp <- function(f) sub(sprintf("^PESlong_%s_(.*)_TrainOOF\\.tsv$", covarType), "\\1", basename(f))

message_ts <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
message_ts("Reading ", length(files), " TrainOOF files (", score_col, ", instance 0)")
L <- lapply(files, function(f) {
  d <- tryCatch(fread(f, select = c("eid", "instance", score_col)), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  setnames(d, score_col, "score")
  d <- d[instance == 0 & is.finite(score)]
  if (!nrow(d)) return(NULL)
  d[, exposure_id := get_exp(f)]
  d[, .(eid, exposure_id, score)]
})
dt <- rbindlist(Filter(Negate(is.null), L))
message_ts("rows=", nrow(dt), " | exposures=", uniqueN(dt$exposure_id), " | eids=", uniqueN(dt$eid))

W <- dcast(dt, eid ~ exposure_id, value.var = "score")
M <- as.matrix(W[, -1L]); rownames(M) <- as.character(W$eid)
# pairwise so exposures with partly-different samples still correlate
C <- cor(M, use = "pairwise.complete.obs")
C[!is.finite(C)] <- NA_real_

# hierarchical cluster order (on 1 - |r| so strong +/- associations cluster)
d <- as.dist(1 - abs(replace(C, is.na(C), 0)))
hc <- hclust(d, method = "average")
ord <- colnames(C)[hc$order]

# write order + long-form correlations
fwrite(data.table(exposure_id = ord, cluster_rank = seq_along(ord)),
       file.path(out_dir, paste0("pes_score_correlation_order_", covarType, ".tsv")), sep = "\t")
long <- as.data.table(as.table(C))
setnames(long, c("exposure_i", "exposure_j", "r"))
fwrite(long, file.path(out_dir, paste0("pes_score_correlation_long_", covarType, ".tsv")), sep = "\t")
# also the wide matrix for convenience
wide <- data.table(exposure_id = rownames(C)); wide <- cbind(wide, as.data.table(C))
fwrite(wide, file.path(out_dir, paste0("pes_score_correlation_matrix_", covarType, ".tsv")), sep = "\t")

message_ts("Wrote correlation matrix (", ncol(C), " exposures) to ", out_dir)
