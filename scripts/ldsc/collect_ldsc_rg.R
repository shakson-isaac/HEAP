#!/usr/bin/env Rscript

# ============================================================================
# collect_ldsc_rg.R — parse LDSC --rg logs into one tidy pairwise table
# ----------------------------------------------------------------------------
# The rg stage (slurm/ldsc/ldsc_rg.sh) writes one log per "root" exposure
# (<root>.rg.log), each holding a "Summary of Genetic Correlation Results"
# table of that root vs every LATER exposure (upper triangle). This script
# scans all of them and extracts, per unordered exposure pair:
#   p1, p2   — the two exposure ids (file basenames, .sumstats.gz stripped)
#   rg       — genetic correlation
#   se       — standard error of rg
#   z, p     — rg z-score and p-value
#
# Output: <ldsc_out_root>/ldsc_rg_summary.tsv (+ a copy under figures/data so
# the genetic-correlation figure can read it). Each pair appears once.
#
# Run:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R Rscript scripts/ldsc/collect_ldsc_rg.R
# Override the rg dir with the LDSC_OUT_ROOT env var (default IGLOO HEAP output).
# ============================================================================

suppressPackageStartupMessages({ library(data.table) })

# --- locate workflow/00_paths.R (for heap_gwas + heap_project_root) ----------
local({
  if (exists("heap_gwas", mode = "function")) return(invisible())
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("collect_ldsc_rg.R: cannot find workflow/00_paths.R ",
                       "(set HEAP_PATHS_FILE).")
  source(hit)
})

out_root <- Sys.getenv("LDSC_OUT_ROOT", unset = heap_gwas("ldsc"))
rg_dir   <- file.path(out_root, "rg")
if (!dir.exists(rg_dir))
  stop("LDSC rg directory not found: ", rg_dir,
       "\nRun the rg stage (slurm/ldsc/ldsc_rg.sh) first.")

logs <- list.files(rg_dir, pattern = "\\.rg\\.log$", full.names = TRUE)
if (!length(logs))
  stop("No *.rg.log files under ", rg_dir, " — run the rg stage first.")

.base <- function(x) sub("\\.sumstats\\.gz$", "", basename(x))
.num  <- function(x) suppressWarnings(as.numeric(x))

# --- parse the "Summary of Genetic Correlation Results" table of one log -----
# The block looks like (whitespace-delimited, p1/p2 are full file paths):
#   Summary of Genetic Correlation Results
#   p1  p2  rg  se  z  p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se
#   <path>.sumstats.gz  <path>.sumstats.gz  0.12  0.05  2.4  0.016  ...
parse_rg_log <- function(path) {
  L <- readLines(path, warn = FALSE)
  hdr <- grep("Summary of Genetic Correlation Results", L, fixed = TRUE)
  if (!length(hdr)) return(NULL)
  colline <- hdr[1] + 1L
  if (colline > length(L)) return(NULL)
  cols <- strsplit(trimws(L[colline]), "[[:space:]]+")[[1]]
  rows <- list()
  i <- colline + 1L
  while (i <= length(L) && nzchar(trimws(L[i]))) {
    f <- strsplit(trimws(L[i]), "[[:space:]]+")[[1]]
    if (length(f) >= length(cols)) {
      v <- as.list(f[seq_along(cols)]); names(v) <- cols
      rows[[length(rows) + 1L]] <- v
    }
    i <- i + 1L
  }
  if (!length(rows)) return(NULL)
  dt <- rbindlist(rows, fill = TRUE)
  dt[, .(p1 = .base(p1), p2 = .base(p2),
         rg = .num(rg), se = .num(se), z = .num(z), p = .num(p))]
}

res <- rbindlist(lapply(logs, parse_rg_log), fill = TRUE)
if (!nrow(res)) stop("Parsed 0 rg rows from ", length(logs), " logs.")

# de-duplicate unordered pairs (keep the lower-se estimate if a pair recurs)
res <- res[p1 != p2]
res[, key := ifelse(p1 < p2, paste(p1, p2), paste(p2, p1))]
setorder(res, key, se)
res <- res[, .SD[1L], by = key]
res[, key := NULL]

setorder(res, -rg)
summ_path <- file.path(out_root, "ldsc_rg_summary.tsv")
fwrite(res, summ_path, sep = "\t")
Sys.chmod(summ_path, mode = "0664")
message("Wrote ", nrow(res), " exposure pairs -> ", summ_path)

# Also stage a copy under the figure-data tree (best-effort) for plotting.
fig_data <- tryCatch(heap_project_root("figures", "data"), error = function(e) NA_character_)
if (!is.na(fig_data)) {
  dir.create(fig_data, recursive = TRUE, showWarnings = FALSE)
  fwrite(res, file.path(fig_data, "ldsc_rg_summary.tsv"), sep = "\t")
}

ok <- res[is.finite(rg)]
sig <- ok[is.finite(p) & p < 0.05]
message(sprintf("Parsed %d logs | %d pairs with an rg | %d nominally significant (p<0.05)",
                length(logs), nrow(ok), nrow(sig)))
if (nrow(ok)) {
  print(head(ok[order(-abs(rg)),
               .(p1 = substr(p1, 1, 28), p2 = substr(p2, 1, 28),
                 rg = round(rg, 3), se = round(se, 3), p = signif(p, 2))], 15))
}
