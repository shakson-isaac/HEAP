#!/usr/bin/env Rscript
# ============================================================================
# stage_gwas_deposit.R --- convert the REGENIE exposure GWAS into a
# submission-ready per-exposure deposit.
#
# NOT part of HEAP_Supplementary_Data.zip: 169 exposures x 7.78M variants is
# ~122 GB raw and ~42 GB gzipped, against a 182 MB supplement archive. These
# ship as a Tier-3 deposit (GWAS Catalog + Zenodo), which is also what the
# GWAS Catalog expects -- one file per trait, not one archive.
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
a     <- commandArgs(TRUE)
LIMIT <- if ("--limit" %in% a) as.integer(a[which(a == "--limit") + 1L]) else Inf
OUT   <- if ("--out" %in% a) a[which(a == "--out") + 1L] else
         heap_project_output("gwas_deposit")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

src <- list.files(heap_project_output("gwas"), pattern = "\\.regenie$",
                  recursive = TRUE, full.names = TRUE)
if (!length(src)) stop("no .regenie files found")
message(sprintf("%d exposure GWAS found; converting %s",
                length(src), if (is.finite(LIMIT)) paste("first", LIMIT) else "all"))

KEEP <- c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","INFO","N","BETA","SE","LOG10P")
man <- vector("list", 0L)
for (i in seq_len(min(length(src), LIMIT))) {
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
  o <- file.path(OUT, paste0(ex, ".tsv.gz"))
  fwrite(d, o, sep = "\t", compress = "gzip")
  man[[length(man) + 1L]] <- data.table(
    exposure_id = ex, file = basename(o), n_variants = nrow(d),
    n_samples = as.integer(stats::median(d$n, na.rm = TRUE)),
    size_mb = round(file.info(o)$size / 1e6, 1))
  message(sprintf("  [%3d/%d] %-52s %8.1f MB  %.0fs", i, min(length(src), LIMIT),
                  substr(ex, 1, 50), file.info(o)$size / 1e6,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
M <- rbindlist(man)
fwrite(M, file.path(OUT, "manifest.tsv"), sep = "\t")

writeLines(c(
  "HEAP exposure GWAS -- summary statistics",
  "========================================",
  "",
  sprintf("One gzipped file per exposure (%d in this bundle). Columns follow the",
          nrow(M)),
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
  "manifest.tsv lists every exposure with its variant and sample counts."),
  file.path(OUT, "README.txt"))
message(sprintf("\n%s\n  %d exposures, %.1f GB total",
                OUT, nrow(M), sum(M$size_mb) / 1000))
