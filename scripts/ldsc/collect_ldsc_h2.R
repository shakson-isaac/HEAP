#!/usr/bin/env Rscript

# ============================================================================
# collect_ldsc_h2.R — parse LDSC --h2 logs into one tidy summary table
# ----------------------------------------------------------------------------
# Scans the LDSC h2 output directory (one <exposure>.log per exposure, written
# by scripts/ldsc/run_ldsc_h2.sh) and extracts, per exposure:
#   h2, h2_se            — SNP heritability on the observed scale (+ s.e.)
#   intercept, intercept_se — the LDSC intercept = MODEL-BASED genomic-inflation
#                          estimate (≈1 under no confounding; > 1 flags
#                          confounding/structure rather than polygenicity)
#   lambda_gc            — LDSC's own lambda_GC (compare to the QQ-based one)
#   mean_chi2            — mean test statistic
#   ratio, ratio_se      — (intercept-1)/(mean_chi2-1): fraction of inflation
#                          NOT due to polygenic signal
#   n_snp                — regression SNPs after HapMap3 merge
#
# Output: <ldsc_out_root>/ldsc_h2_summary.tsv (+ a copy under figures/data so a
# figure can plot intercept-vs-lambda_GC). Adds exposure category/label.
#
# Run:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R Rscript scripts/ldsc/collect_ldsc_h2.R
# Override the h2 dir with the LDSC_OUT_ROOT env var (default IGLOO HEAP output).
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

# --- locate workflow/00_paths.R (for heap_gwas + heap_config) ---------------
local({
  if (exists("heap_gwas", mode = "function")) return(invisible())
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("collect_ldsc_h2.R: cannot find workflow/00_paths.R ",
                       "(set HEAP_PATHS_FILE).")
  source(hit)
})

out_root <- Sys.getenv("LDSC_OUT_ROOT", unset = heap_gwas("ldsc"))
h2_dir   <- file.path(out_root, "h2")
if (!dir.exists(h2_dir))
  stop("LDSC h2 directory not found: ", h2_dir,
       "\nRun the LDSC array (slurm/ldsc/ldsc_h2_array.sh) first.")

logs <- list.files(h2_dir, pattern = "\\.log$", full.names = TRUE)
if (!length(logs))
  stop("No *.log files under ", h2_dir, " — run the LDSC h2 array first.")

# --- parse one ldsc.py --h2 log ---------------------------------------------
# Lines look like:
#   Total Observed scale h2: 0.0858 (0.0042)
#   Lambda GC: 1.1234
#   Mean Chi^2: 1.2345
#   Intercept: 1.0234 (0.0089)
#   Ratio: 0.0987 (0.0345)     [ or "Ratio < 0 (...)" or "Ratio: NA (...)" ]
#   After merging with regression SNP LD, 1175331 SNPs remain.
.num <- function(x) suppressWarnings(as.numeric(x))
grab_est_se <- function(L, key) {
  ln <- grep(key, L, value = TRUE, fixed = TRUE)
  if (!length(ln)) return(c(NA_real_, NA_real_))
  ln <- ln[length(ln)]
  m <- regmatches(ln, regexec("([-0-9.eE]+)[[:space:]]*\\(([-0-9.eE]+)\\)", ln))[[1]]
  if (length(m) == 3L) return(c(.num(m[2]), .num(m[3])))
  m1 <- regmatches(ln, regexec(":[[:space:]]*([-0-9.eE]+)", ln))[[1]]
  if (length(m1) == 2L) return(c(.num(m1[2]), NA_real_))
  c(NA_real_, NA_real_)
}
parse_log <- function(path) {
  L <- readLines(path, warn = FALSE)
  h2  <- grab_est_se(L, "Total Observed scale h2:")
  itc <- grab_est_se(L, "Intercept:")
  rat <- grab_est_se(L, "Ratio:")
  lam <- grab_est_se(L, "Lambda GC:")[1]
  mc  <- grab_est_se(L, "Mean Chi^2:")[1]
  nsnp_ln <- grep("SNPs remain", L, value = TRUE)
  nsnp <- if (length(nsnp_ln))
    .num(gsub("[^0-9]", "", regmatches(nsnp_ln[length(nsnp_ln)],
              regexec(",[[:space:]]*([0-9]+) SNPs remain", nsnp_ln[length(nsnp_ln)]))[[1]][2]))
    else NA_real_
  failed <- any(grepl("ERROR|Traceback", L))
  data.table(exposure = sub("\\.log$", "", basename(path)),
             h2 = h2[1], h2_se = h2[2],
             intercept = itc[1], intercept_se = itc[2],
             lambda_gc = lam, mean_chi2 = mc,
             ratio = rat[1], ratio_se = rat[2],
             n_snp = nsnp, failed = failed)
}

res <- rbindlist(lapply(logs, parse_log), fill = TRUE)

# --- annotate with exposure category + a readable label ----------------------
add_meta <- function(res) {
  f <- tryCatch(heap_config("exposure_sets", "analysis_exposures.tsv"),
                error = function(e) NA_character_)
  if (!is.na(f) && file.exists(f)) {
    ex <- fread(f)
    res <- merge(res, ex[, .(exposure = variable, category, variable_type)],
                 by = "exposure", all.x = TRUE)
  } else res[, c("category", "variable_type") := NA_character_]
  res[, label := {
    x <- gsub("_f[0-9]+_[0-9]+_[0-9]+$", "", exposure)
    x <- gsub("_+", " ", x); trimws(x)
  }]
  res[]
}
res <- add_meta(res)
res[, h2_z := h2 / h2_se]
setorder(res, -h2)

# --- write outputs -----------------------------------------------------------
summ_path <- file.path(out_root, "ldsc_h2_summary.tsv")
fwrite(res, summ_path, sep = "\t")
Sys.chmod(summ_path, mode = "0664")
message("Wrote ", nrow(res), " exposures -> ", summ_path)

# Also stage a copy under the figure-data tree (best-effort) for plotting.
fig_data <- tryCatch(heap_project_root("figures", "data"), error = function(e) NA_character_)
if (!is.na(fig_data)) {
  dir.create(fig_data, recursive = TRUE, showWarnings = FALSE)
  fwrite(res, file.path(fig_data, "ldsc_h2_summary.tsv"), sep = "\t")
}

# --- console digest ----------------------------------------------------------
ok <- res[!is.na(h2) & failed == FALSE]
message(sprintf("Parsed %d logs | %d with an h2 estimate | %d flagged failed",
                nrow(res), nrow(ok), sum(res$failed, na.rm = TRUE)))
if (nrow(ok)) {
  print(head(ok[, .(exposure = substr(label, 1, 34),
                    h2 = round(h2, 4), se = round(h2_se, 4),
                    intercept = round(intercept, 4),
                    lambda_gc = round(lambda_gc, 3),
                    ratio = round(ratio, 3))], 15))
}
