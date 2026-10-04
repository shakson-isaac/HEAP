#!/usr/bin/env Rscript

# ============================================================================
# summarize_gwas_loci.R
# ----------------------------------------------------------------------------
# Manuscript-ready catalogue of independent lead variants / loci for every
# completed exposure GWAS (REGENIE step 2). Two machine-readable tables + a
# citation-ready markdown summary, so the genetic architecture of each exposure
# is documented in one place and trivial to quote / drop into a supplement.
#
# Definitions:
#   * genome-wide significant : p < 5e-8   (-log10p >= 7.30103)
#   * suggestive              : p < 1e-5   (-log10p >= 5)
#   * independent lead locus  : greedy LD-free distance clump of the genome-wide
#       significant variants -- take the most significant variant, drop all within
#       +/- WINDOW bp on the same chromosome, repeat (heap_gwas_lead_variants).
#       LD-free, so two truly-LD-linked signals >WINDOW apart count separately;
#       annotate with FUMA / Open Targets for LD-aware gene mapping (see notes).
#
# Outputs (authoritative copies in IGLOO; committed copies under docs/):
#   <OUT>/gwas/gwas_lead_variants.tsv    one row per (exposure, lead locus)
#   <OUT>/gwas/gwas_locus_summary.tsv    one row per exposure (counts + lambda + top)
#   docs/manuscript_stats/gwas_loci/{SUMMARY.md, *.tsv}   (best-effort; repo)
#   figures/data/gwas_locus_summary.tsv  (for plotting / website)
#
# Reading every *.regenie (~7.8M variants, ~50 s) is expensive, so per-exposure
# results are CACHED to <OUT>/gwas/gwas_loci_cache/<exposure>.rds and only missing
# exposures are scanned -- so this is incremental and RE-RUNNABLE as more GWAS
# finish. Env knobs:
#   HEAP_GWAS_LOCI_REFRESH=1   recompute all (ignore cache)
#   HEAP_GWAS_LOCI_MAX=<n>     cap exposures scanned this run
#   GW_WINDOW_KB=<kb>          clump half-window (default 500)
#
# Run:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/analysis_summaries/summarize_gwas_loci.R
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("summarize_gwas_loci.R: cannot locate common/ helpers")
  for (f in c("figure_paths", "load_heap_results", "label_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({ library(data.table) })

GW_THR   <- -log10(5e-8)
SUG_THR  <- -log10(1e-5)
WINDOW   <- 1000 * suppressWarnings(as.numeric(Sys.getenv("GW_WINDOW_KB", unset = "500")))
if (!is.finite(WINDOW) || WINDOW <= 0) WINDOW <- 5e5

completed <- list_exposure_gwas(completed_only = TRUE)
if (!length(completed))
  stop("No completed exposure GWAS under gwas/regenie_step2/.")

out_gwas   <- heap_gwas()                       # IGLOO .../output/gwas
cache_dir  <- file.path(out_gwas, "gwas_loci_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
refresh    <- nzchar(Sys.getenv("HEAP_GWAS_LOCI_REFRESH", unset = ""))

# --- per-exposure scan: lead variants + QC scalars (cached) ------------------
scan_one <- function(ex) {
  g <- tryCatch(load_exposure_gwas(ex), error = function(e) NULL)
  if (is.null(g)) return(NULL)
  g <- g[is.finite(LOG10P) & chr %in% 1:22]
  leads <- heap_gwas_lead_variants(g, log10p_thresh = GW_THR, window = WINDOW)
  if (nrow(leads)) {
    setorder(leads, -LOG10P)
    leads[, locus := .I]
    leads <- leads[, .(exposure = ex, locus, chr, pos, rsid = ID,
                       effect_allele = ALLELE1, other_allele = ALLELE0,
                       eaf = A1FREQ, beta = BETA, se = SE,
                       log10p = LOG10P, p = 10^(-LOG10P), n = N,
                       locus_label = sprintf("chr%d:%.2fMb", chr, pos / 1e6))]
  } else {
    leads <- data.table(exposure = character(0), locus = integer(0),
                        chr = integer(0), pos = integer(0), rsid = character(0),
                        effect_allele = character(0), other_allele = character(0),
                        eaf = numeric(0), beta = numeric(0), se = numeric(0),
                        log10p = numeric(0), p = numeric(0), n = numeric(0),
                        locus_label = character(0))
  }
  top <- if (nrow(leads)) leads[1] else NULL
  qc <- data.table(
    exposure     = ex,
    n_variants   = nrow(g),
    lambda_gc    = heap_gwas_lambda(g),
    mean_chi2    = if ("CHISQ" %in% names(g)) mean(g$CHISQ, na.rm = TRUE) else NA_real_,
    n_suggestive = sum(g$LOG10P >= SUG_THR),
    n_gwsig      = sum(g$LOG10P >= GW_THR),
    n_lead       = nrow(leads),
    max_log10p   = max(g$LOG10P, na.rm = TRUE),
    top_rsid     = if (!is.null(top)) top$rsid  else NA_character_,
    top_locus    = if (!is.null(top)) top$locus_label else NA_character_,
    top_p        = if (!is.null(top)) top$p     else NA_real_)
  list(qc = qc, leads = leads)
}

# Which exposures still need scanning?
heap_safe_name_local <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)
cache_path <- function(ex) file.path(cache_dir, paste0(heap_safe_name_local(ex), ".rds"))
n_cached <- sum(file.exists(vapply(completed, cache_path, character(1))))
todo <- if (refresh) completed else
  completed[!file.exists(vapply(completed, cache_path, character(1)))]
qmax <- suppressWarnings(as.integer(Sys.getenv("HEAP_GWAS_LOCI_MAX", unset = "")))
if (!is.na(qmax) && qmax > 0L && length(todo) > qmax) todo <- todo[seq_len(qmax)]

message(sprintf("summarize_gwas_loci: %d completed | %d cached | %d to scan | window=%dkb",
                length(completed), if (refresh) 0L else n_cached,
                length(todo), as.integer(WINDOW / 1000)))

if (length(todo)) {
  allowed <- tryCatch(length(parallel::mcaffinity()), error = function(e) 1L)
  ncores  <- max(1L, min(6L, allowed, length(todo)))
  old_thr <- getDTthreads(); setDTthreads(max(1L, floor(allowed / ncores)))
  message(sprintf("  scanning with %d worker(s) (allowed cpus=%d)", ncores, allowed))
  chunks <- split(todo, ceiling(seq_along(todo) / ncores))
  for (ci in seq_along(chunks)) {
    ch <- chunks[[ci]]
    res <- if (ncores > 1L && requireNamespace("parallel", quietly = TRUE))
      parallel::mclapply(ch, function(ex) { message("  [scan] ", ex); scan_one(ex) },
                         mc.cores = ncores, mc.preschedule = FALSE)
    else lapply(ch, function(ex) { message("  [scan] ", ex); scan_one(ex) })
    for (i in seq_along(ch)) if (!is.null(res[[i]])) saveRDS(res[[i]], cache_path(ch[i]))
    message(sprintf("  cached %d/%d exposures (chunk %d/%d)",
                    sum(file.exists(vapply(completed, cache_path, character(1)))),
                    length(completed), ci, length(chunks)))
  }
  setDTthreads(old_thr)
}

# --- aggregate everything that is cached -------------------------------------
have <- completed[file.exists(vapply(completed, cache_path, character(1)))]
if (!length(have)) stop("No cached per-exposure results — nothing to aggregate.")
objs <- lapply(have, function(ex) readRDS(cache_path(ex)))
summ  <- rbindlist(lapply(objs, `[[`, "qc"),    fill = TRUE)
leads <- rbindlist(lapply(objs, `[[`, "leads"), fill = TRUE)

# per-lead-SNP MR instrument strength: F ≈ Wald chi-square = (beta/se)^2. These
# exposure GWAS feed two-sample MR (Module 5), so the per-locus F documents
# instrument strength; lead SNPs are p<5e-8 by construction so all are "strong"
# (F > ~29), and min_f per exposure is the weakest lead instrument.
if (nrow(leads)) leads[, fstat := (beta / se)^2]
fstats <- if (nrow(leads)) {
  leads[, .(min_f = min(fstat), median_f = as.numeric(median(fstat))), by = exposure]
} else {
  data.table(exposure = character(0), min_f = numeric(0), median_f = numeric(0))
}
summ <- merge(summ, fstats, by = "exposure", all.x = TRUE)

# annotate with exposure category + readable label
ex_tab <- tryCatch(heap_exposure_table(), error = function(e) NULL)
if (!is.null(ex_tab)) {
  lut <- ex_tab[, .(exposure = variable, category, variable_type)]
  summ  <- merge(summ,  lut, by = "exposure", all.x = TRUE)
  leads <- merge(leads, lut[, .(exposure, category)], by = "exposure", all.x = TRUE)
}
summ[, label := heap_exposure_label(exposure)]
leads[, label := heap_exposure_label(exposure)]
setorder(summ, -n_lead, -max_log10p)
setcolorder(summ, c("exposure", "label", "category", "variable_type",
                    "n_variants", "lambda_gc", "mean_chi2", "n_suggestive",
                    "n_gwsig", "n_lead", "min_f", "median_f",
                    "max_log10p", "top_rsid", "top_locus", "top_p"))
setorder(leads, exposure, locus)
setcolorder(leads, c("exposure", "label", "category", "locus", "locus_label",
                     "chr", "pos", "rsid", "effect_allele", "other_allele",
                     "eaf", "beta", "se", "p", "log10p", "fstat", "n"))

# --- write authoritative TSVs (IGLOO) + figures/data + committed docs ---------
w <- function(dt, path) { fwrite(dt, path, sep = "\t"); Sys.chmod(path, "0664"); path }
p_summ  <- w(summ,  file.path(out_gwas, "gwas_locus_summary.tsv"))
p_leads <- w(leads, file.path(out_gwas, "gwas_lead_variants.tsv"))
message("Wrote ", nrow(summ), " exposures -> ", p_summ)
message("Wrote ", nrow(leads), " lead loci -> ", p_leads)

fig_data <- tryCatch(heap_project_root("figures", "data"), error = function(e) NA)
if (!is.na(fig_data)) { dir.create(fig_data, recursive = TRUE, showWarnings = FALSE)
  fwrite(summ, file.path(fig_data, "gwas_locus_summary.tsv"), sep = "\t") }

# --- citation-ready markdown summary (best-effort to repo docs) --------------
n_exp  <- nrow(summ); n_inst <- sum(summ$n_lead > 0)
tot_loci <- sum(summ$n_lead); med_lam <- median(summ$lambda_gc, na.rm = TRUE)
top_by_loci <- head(summ[order(-n_lead)], 15)
md <- c(
  "# Exposure GWAS — lead-variant / loci summary",
  "",
  sprintf("_Generated %s by `scripts/analysis_summaries/summarize_gwas_loci.R`._",
          format(Sys.Date())),
  "",
  "## Headline",
  sprintf("- **%d** completed exposure GWAS summarised (REGENIE step 2).", n_exp),
  sprintf("- **%d / %d** have at least one genome-wide-significant locus (p < 5e-8).",
          n_inst, n_exp),
  sprintf("- **%d** independent lead loci in total (greedy %d-kb distance clump).",
          tot_loci, as.integer(WINDOW / 1000)),
  sprintf("- Median genomic-inflation lambda_GC = **%.3f** (full set).", med_lam),
  "",
  "## Most-instrumented exposures (top 15 by # independent loci)",
  "",
  "| Exposure | Category | # loci | # GW-sig SNPs | lambda_GC | top locus | top SNP | top p |",
  "|---|---|--:|--:|--:|---|---|--:|",
  apply(top_by_loci, 1, function(r) sprintf(
    "| %s | %s | %s | %s | %s | %s | %s | %s |",
    r[["label"]], r[["category"]], r[["n_lead"]], r[["n_gwsig"]],
    formatC(as.numeric(r[["lambda_gc"]]), format = "f", digits = 3),
    r[["top_locus"]], r[["top_rsid"]],
    formatC(as.numeric(r[["top_p"]]), format = "e", digits = 1))),
  "",
  "## Files",
  "- `gwas_lead_variants.tsv` — one row per (exposure, lead locus): rsID, chr:pos,",
  "  alleles, EAF, beta, se, p, N. The supplementary lead-SNP table.",
  "- `gwas_locus_summary.tsv` — one row per exposure: variant counts, lambda_GC,",
  "  mean chi^2, # loci, peak signal, top locus.",
  "",
  "## Notes / caveats",
  "- Loci are **LD-free distance clumps**, not LD-aware. For manuscript locus→gene",
  "  mapping run FUMA / Open Targets Genetics on `gwas_lead_variants.tsv` (no gene",
  "  position panel is staged in IGLOO yet).",
  "- `lambda_GC` here is over all ~7.8M variants; LDSC reports its own lambda_GC on",
  "  the ~1M HapMap3 SNPs plus the **intercept** (confounding vs polygenicity) —",
  "  see scripts/ldsc. Effect allele = REGENIE ALLELE1.")
docs <- tryCatch(heap_path("docs", "manuscript_stats", "gwas_loci"),
                 error = function(e) NA_character_)
if (!is.na(docs)) {
  ok <- tryCatch({ dir.create(docs, recursive = TRUE, showWarnings = FALSE)
    writeLines(md, file.path(docs, "SUMMARY.md"))
    fwrite(summ,  file.path(docs, "gwas_locus_summary.tsv"), sep = "\t")
    fwrite(leads, file.path(docs, "gwas_lead_variants.tsv"), sep = "\t"); TRUE },
    error = function(e) FALSE)
  if (ok) message("Wrote manuscript summary -> ", file.path(docs, "SUMMARY.md"))
  else message("NOTE: could not write to ", docs, " (repo not group-writable?) — ",
               "authoritative TSVs are under ", out_gwas)
}

# --- console digest ----------------------------------------------------------
message(sprintf("\nDONE: %d exposures | %d instrumented | %d total lead loci | median lambda_GC=%.3f",
                n_exp, n_inst, tot_loci, med_lam))
print(head(summ[, .(label = substr(label, 1, 30), category,
                    n_lead, n_gwsig, lambda_gc = round(lambda_gc, 3),
                    top = top_locus, top_rsid)], 15))
