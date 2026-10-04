#!/usr/bin/env Rscript
# ============================================================================
# support/coloc/run_coloc_shortlist.R
# ----------------------------------------------------------------------------
# Targeted colocalization (coloc.abf) for the cis-pQTL Tier-1 SHORTLIST — the
# "shortlist -> confirm" step. For each cis-pQTL P->D / P->E edge that reached
# Tier 1 (mr_tiered_edges.tsv, lane==pQTL & edge_class==cis & mr_tier==Tier1),
# colocalize the protein cis-pQTL signal with the outcome GWAS (FinnGen disease
# for P->D, UKB exposure REGENIE for P->E) over the protein's cis window, and
# record PP.H4 so build/compare can promote (coloc-confirmed) or demote
# (cis & not-colocalized -> Tier 2) per mr_schematics/mr_hit_flowchart.tex.
#
# coloc.abf (single-causal-variant ABF) needs NO LD reference. Generalises the
# legacy ModuleMR/COLOC/runColoc.R (deCODE x FinnGen only) to both arms + both
# outcome types. Window = lead cis instrument (clump) +/- COLOC_WINDOW_KB.
#
# Output -> support/coloc/:
#   <locus>_coloc_summary.tsv, <locus>_plot_table.tsv  (per candidate)
#   coloc_shortlist_results.tsv  (aggregate: keyed for the tiering join)
#   coloc_index.tsv              (locus -> PP.H4, for fig_mr_coloc)
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R CPUS=4 \
#     Rscript scripts/support/coloc/run_coloc_shortlist.R
#   (HEAP_COLOC_N=3 to test on the first 3 candidates)
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("Could not locate workflow/00_paths.R (set HEAP_PATHS_FILE).")
  source(hit)
})
suppressPackageStartupMessages({ library(data.table); library(coloc); library(arrow); library(parallel) })
nthreads <- suppressWarnings(as.integer(Sys.getenv("CPUS", "")))
if (is.na(nthreads) || nthreads < 1) nthreads <- max(1L, detectCores() - 1L)
setDTthreads(1L)

WINDOW_KB <- as.numeric(Sys.getenv("COLOC_WINDOW_KB", "500"))
PP4_CONFIRM <- as.numeric(Sys.getenv("COLOC_PP4", "0.8"))
P1 <- 1e-4; P2 <- 1e-4; P12 <- 1e-5
ts <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

summ_dir   <- heap_project_output("mr_edges", "summary")
out_dir    <- heap_project_output("support", "coloc"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
inst_ukb   <- heap_project_output("mr", "protein_inst")
inst_dec   <- heap_project_output("mr", "protein_inst_decode")
UKB_PQTL   <- "/n/groups/patel/IGLOO/UKB/pQTL"
DEC_PQTL   <- "/n/groups/patel/IGLOO/DECODE/pQTL/final_somascan_smp"
DEC_MAP    <- fread("/n/groups/patel/IGLOO/DECODE/pQTLmetadata/somascan_protein_map.tsv")
FINNGEN    <- "/n/groups/patel/IGLOO/FinnGen/SummaryStats"
FG_MAN     <- fread("/n/groups/patel/IGLOO/FinnGen/finngen_R12_manifest.tsv")
EXPO_GWAS  <- heap_gwas("regenie_step2")

# ---- shortlist -------------------------------------------------------------
tiered <- fread(file.path(summ_dir, "mr_tiered_edges.tsv"))
sl <- tiered[lane == "pQTL" & edge_class == "cis" & mr_tier == "Tier1",
             .(dataset, edge_dir, protein = src_id, outcome = tgt_id)]
sl <- unique(sl)
ncap <- suppressWarnings(as.integer(Sys.getenv("HEAP_COLOC_N", "")))
if (!is.na(ncap)) sl <- sl[1:min(ncap, nrow(sl))]
ts("coloc shortlist: ", nrow(sl), " cis Tier-1 candidates (", uniqueN(sl$protein), " proteins)")

# incremental: skip candidates already in coloc_shortlist_results.tsv (non-failed);
# set HEAP_COLOC_FORCE=1 to recompute all. Keeps prior results in the union output.
prev_fp <- file.path(out_dir, "coloc_shortlist_results.tsv")
prev <- NULL
if (file.exists(prev_fp) && Sys.getenv("HEAP_COLOC_FORCE", "") == "") {
  prev <- fread(prev_fp)
  donek <- prev[coloc_status != "failed", paste(dataset, edge_dir, protein, outcome)]
  sl <- sl[!(paste(dataset, edge_dir, protein, outcome) %in% donek)]
  ts("incremental: ", length(donek), " already coloc'd; ", nrow(sl), " NEW cis Tier-1 to run")
}
if (!nrow(sl)) { message("No new cis Tier-1 candidates to coloc — nothing to do."); quit(save = "no", status = 0) }

# ---- helpers ---------------------------------------------------------------
norm_chr <- function(x) gsub("^chr", "", as.character(x))
to_maf <- function(x) { v <- suppressWarnings(as.numeric(x)); ifelse(is.finite(v), pmin(v, 1 - v), NA_real_) }
is_palindromic <- function(a1, a2) { a1 <- toupper(a1); a2 <- toupper(a2)
  (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C") }

lead_window <- function(protein, arm) {                    # cis lead from clump -> chr/lo/hi
  f <- file.path(if (arm == "UKB") inst_ukb else inst_dec,
                 paste0("protein_", protein, "_cis_clumped.tsv"))
  if (!file.exists(f)) return(NULL)
  cl <- fread(f)
  if (!nrow(cl)) return(NULL)
  cl <- cl[order(pval.exposure)][1]
  list(chr = norm_chr(cl$chr.exposure), pos = as.integer(cl$pos.exposure),
       lo = max(1L, as.integer(cl$pos.exposure - WINDOW_KB * 1000)),
       hi = as.integer(cl$pos.exposure + WINDOW_KB * 1000),
       lead = if ("SNP" %in% names(cl)) cl$SNP[1] else NA_character_)
}

read_ukb_pqtl <- function(protein, w) {                    # parquet, cis-window pushdown
  hits <- list.files(UKB_PQTL, pattern = paste0("^", protein, "_.*\\.parquet$"), full.names = TRUE)
  if (!length(hits)) return(NULL)
  d <- tryCatch(as.data.table(
        arrow::open_dataset(hits[1]) |>
          dplyr::filter(CHROM == as.integer(w$chr), GENPOS >= w$lo, GENPOS <= w$hi) |>
          dplyr::collect()), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  pv <- if ("PVAL" %in% names(d)) d$PVAL else 10^(-as.numeric(d$LOG10P))
  data.table(SNP = d$rsid, chr = norm_chr(d$CHROM), pos = as.integer(d$GENPOS),
             beta = as.numeric(d$BETA), se = as.numeric(d$SE), pval = as.numeric(pv),
             effect_allele = toupper(d$ALLELE1), other_allele = toupper(d$ALLELE0),
             eaf = as.numeric(d$A1FREQ), N = as.numeric(d$N))
}
read_dec_pqtl <- function(protein, w) {
  mp <- DEC_MAP[EntrezGeneSymbol == protein]
  if (!nrow(mp)) return(NULL)
  f <- file.path(DEC_PQTL, mp$file[1]); if (!file.exists(f)) return(NULL)
  d <- fread(f, showProgress = FALSE)
  d[, chr := norm_chr(Chrom)]; d <- d[chr == w$chr & Pos >= w$lo & Pos <= w$hi]
  if (!nrow(d)) return(NULL)
  snp <- if ("rsids" %in% names(d)) tstrsplit(d$rsids, ",", fixed = TRUE, keep = 1)[[1]] else d$Name
  data.table(SNP = snp, chr = d$chr, pos = as.integer(d$Pos),
             beta = as.numeric(d$Beta), se = as.numeric(d$SE), pval = as.numeric(d$Pval),
             effect_allele = toupper(d$effectAllele), other_allele = toupper(d$otherAllele),
             eaf = as.numeric(d$ImpMAF), N = as.numeric(gsub(",", "", as.character(d$N))))
}
read_finngen <- function(id, w) {
  f <- file.path(FINNGEN, paste0(id, ".gz")); if (!file.exists(f)) return(NULL)
  d <- fread(f, showProgress = FALSE)
  setnames(d, "#chrom", "chrom", skip_absent = TRUE)
  d[, chr := norm_chr(chrom)]; d <- d[chr == w$chr & pos >= w$lo & pos <= w$hi]
  if (!nrow(d)) return(NULL)
  snp <- if ("rsids" %in% names(d)) d$rsids else paste0(d$chr, ":", d$pos, "_", d$ref, "_", d$alt)
  ph <- sub("^finngen_R12_", "", id); man <- FG_MAN[phenocode == ph]
  nca <- if (nrow(man)) as.numeric(man$num_cases[1]) else NA_real_
  nco <- if (nrow(man)) as.numeric(man$num_controls[1]) else NA_real_
  list(d = data.table(SNP = snp, chr = d$chr, pos = as.integer(d$pos),
         beta = as.numeric(d$beta), se = as.numeric(d$sebeta), pval = as.numeric(d$pval),
         effect_allele = toupper(d$alt), other_allele = toupper(d$ref), eaf = as.numeric(d$af_alt)),
       type = "cc", s = nca / (nca + nco), N = nca + nco)
}
read_exposure <- function(id, w) {
  f <- file.path(EXPO_GWAS, id, paste0(id, ".regenie")); if (!file.exists(f)) return(NULL)
  d <- fread(f, showProgress = FALSE)
  d[, chr := norm_chr(CHROM)]; d <- d[chr == w$chr & GENPOS >= w$lo & GENPOS <= w$hi]
  if (!nrow(d)) return(NULL)
  pv <- if ("PVAL" %in% names(d)) d$PVAL else 10^(-as.numeric(d$LOG10P))
  list(d = data.table(SNP = d$ID, chr = d$chr, pos = as.integer(d$GENPOS),
         beta = as.numeric(d$BETA), se = as.numeric(d$SE), pval = as.numeric(pv),
         effect_allele = toupper(d$ALLELE1), other_allele = toupper(d$ALLELE0), eaf = as.numeric(d$A1FREQ)),
       type = "quant", s = NA_real_, N = stats::median(as.numeric(d$N), na.rm = TRUE))
}

harmonize <- function(d1, d2) {
  d1 <- d1[!is.na(SNP) & SNP != "" & SNP != "."][!duplicated(SNP)]
  d2 <- d2[!is.na(SNP) & SNP != "" & SNP != "."][!duplicated(SNP)]
  m <- merge(d1, d2, by = "SNP", suffixes = c(".1", ".2"))
  if (!nrow(m)) return(NULL)
  same <- m$effect_allele.1 == m$effect_allele.2 & m$other_allele.1 == m$other_allele.2
  swap <- m$effect_allele.1 == m$other_allele.2  & m$other_allele.1 == m$effect_allele.2
  m <- m[same | swap]; if (!nrow(m)) return(NULL)
  sw <- swap[same | swap]; if (any(sw)) m[sw, beta.2 := -beta.2]
  pal <- is_palindromic(m$effect_allele.1, m$other_allele.1)
  e1 <- suppressWarnings(as.numeric(m$eaf.1))
  m[!(pal & is.finite(e1) & e1 > 0.42 & e1 < 0.58)]
}
ds <- function(m, k, type, s = NULL, N = NULL) {
  out <- list(snp = m$SNP, beta = as.numeric(m[[paste0("beta.", k)]]),
              varbeta = as.numeric(m[[paste0("se.", k)]])^2,
              MAF = to_maf(m[[paste0("eaf.", k)]]), type = type, N = N)
  if (type == "cc") out$s <- s
  out
}

# ---- run one candidate -----------------------------------------------------
run_one <- function(r) {
  arm <- r$dataset; protein <- r$protein; outcome <- r$outcome
  w <- lead_window(protein, arm)
  base <- list(dataset = arm, edge_dir = r$edge_dir, protein = protein, outcome = outcome,
               nsnps = NA_integer_, PP.H3 = NA_real_, PP.H4 = NA_real_,
               coloc_status = "failed", note = "")
  if (is.null(w)) { base$note <- "no_cis_clump"; return(base) }
  pq <- if (arm == "UKB") read_ukb_pqtl(protein, w) else read_dec_pqtl(protein, w)
  if (is.null(pq) || !nrow(pq)) { base$note <- "no_pqtl"; return(base) }
  is_dis <- startsWith(outcome, "finngen_R12_")
  oc <- if (is_dis) read_finngen(outcome, w) else read_exposure(outcome, w)
  if (is.null(oc) || !nrow(oc$d)) { base$note <- "no_outcome"; return(base) }
  m <- harmonize(pq, oc$d)
  if (is.null(m) || nrow(m) < 50) { base$note <- "too_few_snps"; base$nsnps <- if (is.null(m)) 0L else nrow(m); return(base) }
  N1 <- stats::median(pq$N, na.rm = TRUE)
  d1 <- ds(m, "1", "quant", N = N1)
  d2 <- ds(m, "2", oc$type, s = oc$s, N = oc$N)
  res <- tryCatch(coloc::coloc.abf(d1, d2, p1 = P1, p2 = P2, p12 = P12), error = function(e) NULL)
  if (is.null(res)) { base$note <- "coloc_error"; base$nsnps <- nrow(m); return(base) }
  S <- res$summary
  locus <- paste(arm, protein, outcome, sep = "__")
  fwrite(data.table(locus = locus, lead = w$lead, chr = w$chr, pos = w$pos, window_kb = WINDOW_KB,
                    nsnps = nrow(m), N_pqtl = N1, N_outcome = oc$N, s = oc$s,
                    PP.H0 = S["PP.H0.abf"], PP.H1 = S["PP.H1.abf"], PP.H2 = S["PP.H2.abf"],
                    PP.H3 = S["PP.H3.abf"], PP.H4 = S["PP.H4.abf"]),
         file.path(out_dir, paste0(locus, "_coloc_summary.tsv")), sep = "\t")
  fwrite(data.table(snp = m$SNP, chr = m$chr.1, pos = m$pos.1,
                    p_trait1 = m$pval.1, p_trait2 = m$pval.2, r2 = NA_real_),
         file.path(out_dir, paste0(locus, "_plot_table.tsv")), sep = "\t")
  pp4 <- unname(S["PP.H4.abf"]); pp3 <- unname(S["PP.H3.abf"])
  base$nsnps <- nrow(m); base$PP.H3 <- pp3; base$PP.H4 <- pp4
  base$coloc_status <- if (pp4 >= PP4_CONFIRM) "confirmed" else if (pp3 >= PP4_CONFIRM) "distinct" else "ambiguous"
  base
}

results <- rbindlist(mclapply(seq_len(nrow(sl)), function(i) {
  r <- sl[i]
  ts(sprintf("[%d/%d] %s %s -> %s", i, nrow(sl), r$dataset, r$protein, r$outcome))
  out <- tryCatch(run_one(r), error = function(e) {
    list(dataset = r$dataset, edge_dir = r$edge_dir, protein = r$protein, outcome = r$outcome,
         nsnps = NA_integer_, PP.H3 = NA_real_, PP.H4 = NA_real_, coloc_status = "failed",
         note = conditionMessage(e)) })
  as.data.table(out)
}, mc.cores = nthreads), fill = TRUE)

# union with prior results (incremental) so the table accumulates across runs
if (!is.null(prev))
  results <- unique(rbind(prev, results, fill = TRUE),
                    by = c("dataset", "edge_dir", "protein", "outcome"))
fwrite(results, file.path(out_dir, "coloc_shortlist_results.tsv"), sep = "\t")
# refresh coloc_index for fig_mr_coloc (locus -> PP.H4), merge with any staged loci
idx <- results[coloc_status != "failed", .(locus = paste(dataset, protein, outcome, sep = "__"),
                                           chr = NA, PP.H4)]
fwrite(idx, file.path(out_dir, "coloc_index_shortlist.tsv"), sep = "\t")

cat("\n================ COLOC SHORTLIST RESULTS ================\n")
print(results[order(-PP.H4), .(dataset, protein, outcome, nsnps, PP.H3 = round(PP.H3,3),
                               PP.H4 = round(PP.H4,3), coloc_status, note)])
cat(sprintf("\nconfirmed (PP.H4>=%.2f): %d / %d\n", PP4_CONFIRM,
            sum(results$coloc_status == "confirmed"), nrow(results)))
cat("Wrote: ", file.path(out_dir, "coloc_shortlist_results.tsv"), "\nDONE.\n", sep = "")
