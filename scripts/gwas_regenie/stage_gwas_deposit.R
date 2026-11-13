#!/usr/bin/env Rscript
# ============================================================================
# stage_gwas_deposit.R --- convert the REGENIE exposure GWAS into a
# submission-ready per-exposure deposit.
#
# NOT part of HEAP_Supplementary_Data.zip: 172 exposures x 7.78M variants is
# ~121 GB raw and ~51 GB compressed, against a 182 MB supplement archive. These
# ship as a Tier-3 deposit -- one file per trait, not one archive, which is also
# what the GWAS Catalog expects.
#
# BGZIP, NOT GZIP, AND AN INDEX BESIDE IT. The two formats are the same size to
# within 3% (measured: 37.9 vs 38.9 MB per million variants), but a bgzipped
# file with a .tbi can be queried by region, so someone who wants one locus
# reads a few KB instead of a 300 MB file.
#
# The files are published in the requester-pays bucket gs://heap-gwas (listed
# on heap.bio/downloads); the downloader's own Google Cloud project is billed:
#
#   gcloud storage cp gs://heap-gwas/<exposure>.tsv.bgz* . --billing-project=YOUR_PROJECT
#   tabix <exposure>.tsv.bgz 1:1-100000
#
# Column mapping (REGENIE -> GWAS-Catalog standard):
#   ID       -> variant_id          CHROM  -> chromosome
#   GENPOS   -> base_pair_location  ALLELE1-> effect_allele   (BETA is wrt ALLELE1)
#   ALLELE0  -> other_allele        A1FREQ -> effect_allele_frequency
#   BETA     -> beta                SE     -> standard_error
#   LOG10P   -> p_value  (= 10^-LOG10P; REGENIE reports -log10 p, not p)
#   N, INFO  carried through as n and info
#
# Covariates: the base specification (age, age^2, sex, age x sex, age^2 x sex,
# assessment centre, 20 genetic PCs) -- see prepare_gwas_exposures.R:147.
#
#   Run:  Rscript stage_gwas_deposit.R [--limit N] [--out DIR]
#         --limit is for a proof run; omit to convert all exposures.
# ============================================================================
suppressPackageStartupMessages(library(data.table))
setDTthreads(0)

local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
# htslib is not on the default PATH on O2; fall back to the module install.
htslib_bin <- function(tool) {
  hit <- Sys.which(tool)
  if (nzchar(hit)) return(unname(hit))
  cand <- Sys.glob(file.path("/n/app/htslib", "*", "bin", tool))
  if (length(cand)) return(sort(cand, decreasing = TRUE)[1])
  stop(sprintf("%s not found; module load htslib or add it to PATH", tool))
}
BGZIP <- htslib_bin("bgzip"); TABIX <- htslib_bin("tabix")

a     <- commandArgs(TRUE)
LIMIT <- if ("--limit" %in% a) as.integer(a[which(a == "--limit") + 1L]) else Inf
# --index converts ONE exposure, so the set can run as a SLURM array; --finalize
# assembles the manifest and README afterwards from the rows each task left
# behind. Running the whole set in one process also works and is what --limit
# does; it just takes 169 x 35 s in series.
IDX   <- if ("--index" %in% a) as.integer(a[which(a == "--index") + 1L]) else NA_integer_
FINAL <- "--finalize" %in% a
OUT   <- if ("--out" %in% a) a[which(a == "--out") + 1L] else
         heap_project_output("gwas_deposit")
PARTS <- file.path(OUT, ".manifest_parts")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(PARTS, recursive = TRUE, showWarnings = FALSE)

src <- list.files(heap_project_output("gwas"), pattern = "\\.regenie$",
                  recursive = TRUE, full.names = TRUE)
if (!length(src)) stop("no .regenie files found")
src <- sort(src)   # a stable order, so --index N means the same file every time

take <- if (!is.na(IDX)) {
  if (IDX < 1L || IDX > length(src))
    stop(sprintf("--index %d is outside 1..%d", IDX, length(src)))
  IDX
} else seq_len(min(length(src), LIMIT))

if (!FINAL)
  message(sprintf("%d exposure GWAS found; converting %s", length(src),
                  if (!is.na(IDX)) sprintf("index %d", IDX)
                  else if (is.finite(LIMIT)) paste("first", LIMIT) else "all"))

KEEP <- c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","LOG10P")
man <- vector("list", 0L)
for (i in if (FINAL) integer(0) else take) {
  f  <- src[i]; ex <- sub("\\.regenie$", "", basename(f))
  t0 <- Sys.time()
  d  <- fread(f, select = KEEP, showProgress = FALSE)
  setnames(d, c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE"),
              c("chromosome","base_pair_location","variant_id","other_allele",
                "effect_allele","effect_allele_frequency","info","n","beta","standard_error"))
  # REGENIE reports -log10(p); the Catalog wants p. Recomputing loses precision
  # below ~1e-300, so keep the exact -log10 value alongside.
  d[, p_value := 10^(-LOG10P)]
  setnames(d, "LOG10P", "neg_log10_p_value")
  setcolorder(d, c("variant_id","chromosome","base_pair_location","effect_allele",
                   "other_allele","effect_allele_frequency","beta","standard_error",
                   "p_value","neg_log10_p_value","n","info"))
  # tabix needs the file sorted by chromosome then position. REGENIE emits it
  # that way already, but a resort is cheap next to the 4 minutes of I/O and
  # turns an assumption into a guarantee.
  setorder(d, chromosome, base_pair_location)

  o   <- file.path(OUT, paste0(ex, ".tsv.bgz"))
  tmp <- file.path(tempdir(), paste0(ex, ".tsv"))
  fwrite(d, tmp, sep = "\t")
  if (system2(BGZIP, c("-f", "-@", "4", "-c", shQuote(tmp)), stdout = o) != 0L)
    stop("bgzip failed for ", ex)
  unlink(tmp)
  # -s/-b/-e are the 1-based columns of chromosome and position; -S 1 skips the
  # header row, which is not a comment line in the Catalog format.
  if (system2(TABIX, c("-f", "-s", "2", "-b", "3", "-e", "3", "-S", "1", shQuote(o))) != 0L)
    stop("tabix failed for ", ex)
  row <- data.table(
    exposure_id = ex, file = basename(o), n_variants = nrow(d),
    n_samples = as.integer(stats::median(d$n, na.rm = TRUE)),
    size_mb = round(file.info(o)$size / 1e6, 1),
    index_kb = round(file.info(paste0(o, ".tbi"))$size / 1e3, 1))
  fwrite(row, file.path(PARTS, paste0(ex, ".tsv")), sep = "\t")
  man[[length(man) + 1L]] <- row
  message(sprintf("  [%3d/%d] %-52s %8.1f MB  %.0fs", i, length(src),
                  substr(ex, 1, 50), file.info(o)$size / 1e6,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
parts <- list.files(PARTS, pattern = "\\.tsv$", full.names = TRUE)
M <- if (length(parts)) rbindlist(lapply(parts, fread)) else rbindlist(man)
if (!nrow(M)) stop("nothing converted and no manifest rows on disk")
setorder(M, exposure_id)
fwrite(M, file.path(OUT, "manifest.tsv"), sep = "\t")

if (!is.na(IDX) && !FINAL) {
  message("  wrote manifest row; run --finalize once the array drains")
  quit(save = "no")
}

writeLines(c(
  "HEAP exposure GWAS -- summary statistics",
  "========================================",
  "",
  sprintf("One bgzipped file per exposure (%d in this bundle), each with a tabix",
          nrow(M)),
  "index (.tbi) beside it. Columns follow the",
  "GWAS Catalog standard:",
  "",
  "  variant_id, chromosome, base_pair_location, effect_allele, other_allele,",
  "  effect_allele_frequency, beta, standard_error, p_value, neg_log10_p_value,",
  "  n, info",
  "",
  "effect_allele is REGENIE's ALLELE1; beta is the effect of that allele.",
  "neg_log10_p_value is REGENIE's native output and is exact; p_value is derived",
  "from it and underflows to 0 below about 1e-308, so use neg_log10_p_value for",
  "the strongest associations.",
  "",
  "Association model: REGENIE step 2, adjusted for the base covariate set --",
  "age, age^2, sex, age x sex, age^2 x sex, assessment centre, and 20 genetic",
  "principal components.",
  "",
  "manifest.tsv lists every exposure with its variant and sample counts.",
  "",
  "Reading one region without downloading the file:",
  "",
  "  tabix <url>/ever_smoked_f20160_0_0.tsv.bgz 19:44900000-44920000",
  "",
  "Reading a whole file:",
  "",
  "  curl -O <url>/ever_smoked_f20160_0_0.tsv.bgz",
  "",
  "R: data.table::fread(cmd = paste('gzip -dc', file))",
  "Python: pandas.read_csv(file, sep='\\t', compression='gzip')",
  "",
  "bgzip output is valid gzip, so any gzip reader works; the index only adds",
  "the option of reading part of it."),
  file.path(OUT, "README.txt"))
message(sprintf("\n%s\n  %d exposures, %.1f GB total",
                OUT, nrow(M), sum(M$size_mb) / 1000))
