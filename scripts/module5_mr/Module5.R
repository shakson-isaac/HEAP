#!/usr/bin/env Rscript

local({
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "workflow", "00_paths.R"),
    file.path(getwd(), "..", "..", "..", "workflow", "00_paths.R")
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit)) source(hit)
})

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(TwoSampleMR)
  library(ieugwasr)
  library(genetics.binaRies)
})

# Shared per-edge sensitivity analyses (Steiger / single-SNP / leave-one-out / MR-PRESSO)
local({
  cand <- c(
    if (exists("heap_script", mode = "function")) heap_script("module5_mr", "mr_sensitivity.R") else "",
    file.path(Sys.getenv("HEAP_ROOT", ""), "scripts", "module5_mr", "mr_sensitivity.R"),
    file.path(getwd(), "scripts", "module5_mr", "mr_sensitivity.R")
  )
  cand <- cand[nzchar(cand)]
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit)) stop("mr_sensitivity.R not found (looked in HEAP_ROOT/scripts/module5_mr)")
  source(hit)
})

# parquet support
if (!requireNamespace("arrow", quietly = TRUE)) {
  stop("Package 'arrow' is required (protein GWAS are parquet). Load/install arrow in this R environment.")
}

# ---------------------------
# Args: edge_type idx [n_chunks] [make_plots]
# edge_type: EP ED PD PE DE DP
# idx: 1..n_chunks
# ---------------------------
# WARM_ONLY (env HEAP_MR_WARM_ONLY): source this file for its CFG + instrument
# functions only (cache pre-warming) — skip arg parsing and the edge run.
WARM_ONLY  <- nzchar(Sys.getenv("HEAP_MR_WARM_ONLY"))
make_plots <- FALSE

if (!WARM_ONLY) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) stop("Usage: Rscript Module5.R <EP|ED|PD|PE|DE|DP> <chunk_idx> [n_chunks] [make_plots]")

  edge_type  <- toupper(args[1])
  idx        <- as.integer(args[2])
  n_chunks   <- as.integer(ifelse(length(args) >= 3, args[3], 5000L))
  make_plots <- as.integer(ifelse(length(args) >= 4, args[4], 0L)) == 1L

  valid_types <- c("EP","ED","PD","PE","DE","DP")
  if (!edge_type %in% valid_types) stop("edge_type must be one of: ", paste(valid_types, collapse=", "))
}

# ---------------------------
# Config (EDIT PATHS AS NEEDED)
# ---------------------------
CFG <- list(
  # HEAP-specific MR outputs (IGLOO-rooted canonical). HEAP_MR_EDGES_DIR override
  # lets the registry-driven incremental run point at a delta edge-list dir
  # (mr_edges/registry/delta) so only missing edges are chunked + run.
  edges_dir   = Sys.getenv("HEAP_MR_EDGES_DIR", unset = heap_project_output("mr_edges", "global_edges")),
  # Per-experiment edge-output base. The launcher exports HEAP_MR_OUTDIR from the
  # manifest's output_path column (e.g. .../mr_edges/MR_UKB_primary) so outputs are
  # experiment-keyed; falls back to the legacy shared dir for direct invocation.
  outdirbase  = Sys.getenv("HEAP_MR_OUTDIR", unset = heap_project_output("mr_edges", "edges")),

  # Instrument cache directories (IGLOO-rooted; reused across runs)
  inst_dir_exposure = heap_project_output("mr", "clumps"),
  inst_dir_protein  = heap_project_output("mr", "protein_inst"),
  inst_dir_disease  = heap_project_output("mr", "disease_inst"),

  ld_bfile  = heap_ldref_prefix("EUR"),   # IGLOO/LDref/EUR (legacy fallback)
  plink_bin = genetics.binaRies::get_plink_binary(),

  p_thr         = 5e-8,
  F_min         = 10,
  clump_r2      = 0.001,
  clump_kb      = 10000,
  cis_window_bp = 1e6,

  # Canonical HEAP GWAS exposure summary stats (IGLOO-rooted).
  # Written by slurm/gwas_regenie/gwas_regenie_exposures_*.sh after step2 completes.
  # Override via HEAP_REGENIE_STEP2_DIR env var for transition period.
  exposure_dir = Sys.getenv(
    "HEAP_REGENIE_STEP2_DIR",
    unset = heap_gwas("regenie_step2")
  ),

  # Shared IGLOO genetics resources (read-only; already canonical)
  protein_dir  = igloo_path("UKB", "pQTL"),
  disease_dir  = igloo_path("FinnGen", "SummaryStats"),
  finngen_manifest = igloo_path("FinnGen", "finngen_R12_manifest.tsv"),

  omicpred_map = heap_omicspred_or_legacy("UKB_Olink_multi_ancestry_models_val_results_portal.csv"),
  pQTLloc_file = igloo_path("UKB", "pQTLmetadata", "olink_protein_map_3k_v1.tsv")
)

dir.create(CFG$outdirbase, recursive=TRUE, showWarnings=FALSE)
dir.create(CFG$inst_dir_protein, recursive=TRUE, showWarnings=FALSE)
dir.create(CFG$inst_dir_disease, recursive=TRUE, showWarnings=FALSE)

# ---------------------------
# Helpers
# ---------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

set_if_present <- function(DT, from, to) {
  if (from %in% names(DT) && !(to %in% names(DT))) setnames(DT, from, to)
  DT
}

# Normalize instrument tables -> TwoSampleMR exposure schema
ensure_tsmr_exposure <- function(DT, label) {
  DT <- as.data.table(DT)
  
  if (!"SNP" %in% names(DT)) {
    if ("rsid" %in% names(DT)) setnames(DT, "rsid", "SNP")
    if ("ID"   %in% names(DT)) setnames(DT, "ID",   "SNP")
  }
  
  DT <- set_if_present(DT, "beta", "beta.exposure")
  DT <- set_if_present(DT, "se",   "se.exposure")
  DT <- set_if_present(DT, "pval", "pval.exposure")
  
  DT <- set_if_present(DT, "BETA", "beta.exposure")
  DT <- set_if_present(DT, "SE",   "se.exposure")
  DT <- set_if_present(DT, "PVAL", "pval.exposure")
  
  DT <- set_if_present(DT, "effect_allele", "effect_allele.exposure")
  DT <- set_if_present(DT, "other_allele",  "other_allele.exposure")
  DT <- set_if_present(DT, "ALLELE1", "effect_allele.exposure")
  DT <- set_if_present(DT, "ALLELE0", "other_allele.exposure")
  
  DT <- set_if_present(DT, "eaf", "eaf.exposure")
  DT <- set_if_present(DT, "A1FREQ", "eaf.exposure")
  
  DT <- set_if_present(DT, "samplesize", "samplesize.exposure")
  DT <- set_if_present(DT, "N", "samplesize.exposure")
  
  DT <- set_if_present(DT, "chr", "chr.exposure")
  DT <- set_if_present(DT, "CHROM", "chr.exposure")
  
  DT <- set_if_present(DT, "pos", "pos.exposure")
  DT <- set_if_present(DT, "GENPOS", "pos.exposure")
  
  if (!"exposure" %in% names(DT)) DT[, exposure := label]
  if (!"id.exposure" %in% names(DT)) DT[, id.exposure := label]
  
  req <- c("SNP","beta.exposure","se.exposure","pval.exposure","effect_allele.exposure","other_allele.exposure")
  miss <- setdiff(req, names(DT))
  if (length(miss)) stop("Missing exposure cols: ", paste(miss, collapse=", "), "\nHave: ", paste(names(DT), collapse=", "))
  DT
}

# Normalize outcome tables -> TwoSampleMR outcome schema
ensure_tsmr_outcome <- function(DT, label) {
  DT <- as.data.table(DT)
  
  if (!"SNP" %in% names(DT)) {
    if ("rsid" %in% names(DT)) setnames(DT, "rsid", "SNP")
    if ("ID"   %in% names(DT)) setnames(DT, "ID",   "SNP")
  }
  
  DT <- set_if_present(DT, "beta", "beta.outcome")
  DT <- set_if_present(DT, "se",   "se.outcome")
  DT <- set_if_present(DT, "pval", "pval.outcome")
  
  DT <- set_if_present(DT, "BETA", "beta.outcome")
  DT <- set_if_present(DT, "SE",   "se.outcome")
  DT <- set_if_present(DT, "PVAL", "pval.outcome")
  
  DT <- set_if_present(DT, "effect_allele", "effect_allele.outcome")
  DT <- set_if_present(DT, "other_allele",  "other_allele.outcome")
  DT <- set_if_present(DT, "ALLELE1", "effect_allele.outcome")
  DT <- set_if_present(DT, "ALLELE0", "other_allele.outcome")
  
  DT <- set_if_present(DT, "eaf", "eaf.outcome")
  DT <- set_if_present(DT, "A1FREQ", "eaf.outcome")
  
  DT <- set_if_present(DT, "samplesize", "samplesize.outcome")
  DT <- set_if_present(DT, "N", "samplesize.outcome")
  
  DT <- set_if_present(DT, "chr", "chr.outcome")
  DT <- set_if_present(DT, "CHROM", "chr.outcome")
  
  DT <- set_if_present(DT, "pos", "pos.outcome")
  DT <- set_if_present(DT, "GENPOS", "pos.outcome")
  
  if (!"outcome" %in% names(DT)) DT[, outcome := label]
  if (!"id.outcome" %in% names(DT)) DT[, id.outcome := label]
  
  req <- c("SNP","beta.outcome","se.outcome","pval.outcome","effect_allele.outcome","other_allele.outcome")
  miss <- setdiff(req, names(DT))
  if (length(miss)) stop("Missing outcome cols: ", paste(miss, collapse=", "), "\nHave: ", paste(names(DT), collapse=", "))
  DT
}

select_and_clump <- function(exp_tsmr) {
  inst <- exp_tsmr[pval.exposure <= CFG$p_thr]
  inst[, F_stat := (beta.exposure / se.exposure)^2]
  inst <- inst[F_stat > CFG$F_min]
  if (nrow(inst) < 1) stop("No instruments pass thresholds.")
  clump_data(inst, clump_r2=CFG$clump_r2, clump_kb=CFG$clump_kb,
             bfile=CFG$ld_bfile, plink_bin=CFG$plink_bin)
}

pick_preferred <- function(mr_dt) {
  if (is.null(mr_dt) || nrow(mr_dt) == 0) return(NULL)
  pref <- c("Inverse variance weighted", "Wald ratio", "MR Egger", "Weighted median", "Weighted mode")
  for (m in pref) {
    r <- mr_dt[method == m]
    if (nrow(r) >= 1) return(r[1])
  }
  mr_dt[1]
}

run_mr_edge <- function(inst_exp, out_tsmr, outdir, prefix) {
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  
  lockfile <- file.path(outdir, ".LOCK")
  if (file.exists(lockfile)) return(invisible(NULL))
  file.create(lockfile)
  on.exit(unlink(lockfile), add=TRUE)
  
  dat <- harmonise_data(as.data.frame(inst_exp), as.data.frame(out_tsmr), action=2)
  dat <- dat[!is.na(dat$beta.exposure) & !is.na(dat$beta.outcome), ]
  
  if (nrow(dat) == 0) {
    sumrow <- data.table(prefix=prefix, nsnp=0L, method=NA_character_, b=NA_real_, se=NA_real_, pval=NA_real_)
    fwrite(sumrow, file.path(outdir, paste0(prefix, "_summary.tsv")), sep="\t")
    return(invisible(sumrow))
  }
  
  if (nrow(dat) == 1) {
    mr_res <- mr(dat, method_list=c("mr_wald_ratio"))
  } else if (nrow(dat) == 2) {
    mr_res <- mr(dat, method_list=c("mr_ivw","mr_egger_regression"))
  } else {
    mr_res <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  }
  
  mr_dt <- as.data.table(mr_res)
  pref <- pick_preferred(mr_dt)
  
  sumrow <- data.table(
    prefix=prefix,
    nsnp=nrow(dat),
    method=pref$method %||% NA_character_,
    b=pref$b %||% NA_real_,
    se=pref$se %||% NA_real_,
    pval=pref$pval %||% NA_real_
  )
  
  fwrite(mr_dt, file.path(outdir, paste0(prefix, "_mr_methods.tsv")), sep="\t")
  fwrite(as.data.table(dat), file.path(outdir, paste0(prefix, "_harmonised.tsv")), sep="\t")
  fwrite(sumrow, file.path(outdir, paste0(prefix, "_summary.tsv")), sep="\t")
  
  het <- tryCatch(mr_heterogeneity(dat), error=function(e) NULL)
  ple <- tryCatch(mr_pleiotropy_test(dat), error=function(e) NULL)
  if (!is.null(het)) fwrite(as.data.table(het), file.path(outdir, paste0(prefix, "_heterogeneity.tsv")), sep="\t")
  if (!is.null(ple)) fwrite(as.data.table(ple), file.path(outdir, paste0(prefix, "_pleiotropy.tsv")), sep="\t")

  # Extra sensitivity analyses (Steiger / single-SNP / leave-one-out / MR-PRESSO)
  run_mr_sensitivity(dat, outdir, prefix)
  
  if (make_plots && nrow(dat) >= 3) {
    p_sc <- mr_scatter_plot(mr_res, dat)
    ggsave(file.path(outdir, paste0(prefix, "_scatter.png")), p_sc[[1]], width=6, height=5, dpi=300)
  }
  
  invisible(sumrow)
}

# ---------------------------
# Readers -> generic schema (SNP,beta,se,pval,effect_allele,other_allele,eaf,chr,pos)
# ---------------------------
read_exposure_gwas_generic <- function(ExID) {
  # Canonical simplified name <exposure>.regenie; fall back to the legacy
  # regenie_step2_<exp>_<exp>.regenie for outputs produced before the rename.
  f <- file.path(CFG$exposure_dir, ExID, paste0(ExID, ".regenie"))
  if (!file.exists(f)) {
    legacy <- file.path(CFG$exposure_dir, ExID,
                        paste0("regenie_step2_", ExID, "_", ExID, ".regenie"))
    if (file.exists(legacy)) f <- legacy else stop("Exposure GWAS not found: ", f)
  }
  DT <- fread(f, showProgress=FALSE)
  if (!"PVAL" %in% names(DT) && "LOG10P" %in% names(DT)) DT[, PVAL := 10^(-1 * LOG10P)]
  
  set_if_present(DT, "ID", "SNP")
  set_if_present(DT, "BETA", "beta")
  set_if_present(DT, "SE", "se")
  set_if_present(DT, "PVAL", "pval")
  set_if_present(DT, "ALLELE1", "effect_allele")
  set_if_present(DT, "ALLELE0", "other_allele")
  set_if_present(DT, "A1FREQ", "eaf")
  set_if_present(DT, "N", "samplesize")
  set_if_present(DT, "CHROM", "chr")
  set_if_present(DT, "GENPOS", "pos")
  
  DT
}

read_protein_gwas_generic <- function(protID) {
  om <- fread(CFG$omicpred_map, showProgress=FALSE)
  OlinkID <- om[Gene %in% protID, unique(Olink_ID)]
  if (length(OlinkID) != 1) stop("Protein->Olink mapping not unique: ", protID)
  
  protFiles <- list.files(CFG$protein_dir, full.names=TRUE)
  pf <- protFiles[grepl(OlinkID, protFiles)]
  if (length(pf) != 1) stop("Protein GWAS file not unique: ", protID)
  
  DT <- as.data.table(arrow::read_parquet(pf))
  if (!"PVAL" %in% names(DT) && "LOG10P" %in% names(DT)) DT[, PVAL := 10^(-1 * LOG10P)]
  
  set_if_present(DT, "rsid", "SNP")
  set_if_present(DT, "BETA", "beta")
  set_if_present(DT, "SE", "se")
  set_if_present(DT, "PVAL", "pval")
  set_if_present(DT, "ALLELE1", "effect_allele")
  set_if_present(DT, "ALLELE0", "other_allele")
  set_if_present(DT, "A1FREQ", "eaf")
  set_if_present(DT, "N", "samplesize")
  set_if_present(DT, "CHROM", "chr")
  set_if_present(DT, "GENPOS", "pos")
  
  DT
}

read_disease_gwas_generic <- function(dzID) {
  f <- file.path(CFG$disease_dir, paste0(dzID, ".gz"))
  DT <- fread(f, showProgress=FALSE)
  
  if (!"rsid" %in% names(DT)) {
    if ("rsids" %in% names(DT)) DT[, rsid := tstrsplit(rsids, ",", fixed=TRUE, keep=1)]
    if (!"rsid" %in% names(DT)) DT[, rsid := paste0(get("#chrom"), ":", pos, "_", ref, "_", alt)]
  }
  
  setnames(DT, "rsid", "SNP")
  set_if_present(DT, "beta", "beta")
  set_if_present(DT, "sebeta", "se")
  set_if_present(DT, "pval", "pval")
  set_if_present(DT, "alt", "effect_allele")
  set_if_present(DT, "ref", "other_allele")
  set_if_present(DT, "af_alt", "eaf")
  set_if_present(DT, "#chrom", "chr")
  set_if_present(DT, "pos", "pos")

  # FinnGen case/control N from the manifest (phenocode = id minus finngen_R12_),
  # enabling Steiger directionality on the binary disease outcome (.gz files carry
  # no per-SNP N).
  .pheno <- sub("^finngen_R12_", "", dzID)
  .man <- tryCatch(fread(CFG$finngen_manifest, showProgress=FALSE), error=function(e) NULL)
  if (!is.null(.man) && "phenocode" %in% names(.man)) {
    .mrow <- .man[phenocode == .pheno]
    if (nrow(.mrow) >= 1L) {
      .nca <- as.numeric(.mrow$num_cases[1]); .nco <- as.numeric(.mrow$num_controls[1])
      DT[, ncase := .nca][, ncontrol := .nco][, samplesize := .nca + .nco]
    }
  }

  DT
}

get_protein_gene_coords <- function(protID) {
  loc <- fread(CFG$pQTLloc_file, showProgress=FALSE)
  row <- loc[HGNC.symbol %in% protID]
  if (nrow(row) != 1) stop("Gene coords not unique: ", protID)
  list(
    chr = as.character(row$chr[1]),
    start = min(row$gene_start[1], row$gene_end[1]),
    end   = max(row$gene_start[1], row$gene_end[1])
  )
}

# ---------------------------
# Instruments (disk-cached)
# ---------------------------
get_exposure_instruments <- function(ExID) {
  f <- file.path(CFG$inst_dir_exposure, paste0("exposure_", ExID, "_clumped.csv"))
  if (!file.exists(f)) {
    # fallback: compute if missing (optional)
    g <- read_exposure_gwas_generic(ExID)
    exp <- ensure_tsmr_exposure(g, ExID)
    inst <- select_and_clump(exp)
    fwrite(inst, f)
  }
  inst <- fread(f, showProgress=FALSE)
  ensure_tsmr_exposure(inst, ExID)
}

get_disease_instruments <- function(dzID) {
  cache_file <- file.path(CFG$inst_dir_disease, paste0("disease_", dzID, "_clumped.tsv"))
  if (file.exists(cache_file)) {
    inst <- fread(cache_file, showProgress=FALSE)
    return(ensure_tsmr_exposure(inst, dzID))
  }
  g <- read_disease_gwas_generic(dzID)
  exp <- ensure_tsmr_exposure(g, dzID)
  inst <- select_and_clump(exp)
  inst <- ensure_tsmr_exposure(inst, dzID)
  fwrite(inst, cache_file, sep="\t")
  inst
}

get_protein_instruments <- function(protID, which=c("cis","trans")) {
  which <- match.arg(which)
  cache_file <- file.path(CFG$inst_dir_protein, paste0("protein_", protID, "_", which, "_clumped.tsv"))
  
  if (file.exists(cache_file)) {
    inst <- fread(cache_file, showProgress=FALSE)
    return(ensure_tsmr_exposure(inst, protID))
  }
  
  g <- read_protein_gwas_generic(protID)
  coords <- get_protein_gene_coords(protID)
  g[, chr := as.character(chr)]
  
  #'*cis- location identifier!* 
  # Below OK for UKB/Decode Ensembl coords as 
  # olink_protein_map_3k_v1.tsv and BioMart has proper start/end directionality 
  # start < end coord
  cis_lo <- coords$start - CFG$cis_window_bp
  cis_hi <- coords$end   + CFG$cis_window_bp
  
  #'Extra Careful Version: No need 
  #gene_lo <- min(coords$start, coords$end)
  #gene_hi <- max(coords$start, coords$end)
  #cis_lo <- max(1, gene_lo - CFG$cis_window_bp)
  #cis_hi <- gene_hi + CFG$cis_window_bp
  
  
  if (which == "cis") {
    sub <- g[chr == coords$chr & pos > cis_lo & pos < cis_hi]
  } else {
    sub <- g[!(chr == coords$chr & pos > cis_lo & pos < cis_hi)]
  }
  if (nrow(sub) < 1) stop("No variants left for ", protID, " [", which, "]")
  
  exp <- ensure_tsmr_exposure(sub, protID)
  inst <- select_and_clump(exp)
  inst <- ensure_tsmr_exposure(inst, protID)
  fwrite(inst, cache_file, sep="\t")
  inst
}

# ---------------------------
# Memoised outcomes (per job)
# ---------------------------
.memo <- new.env(parent=emptyenv())
get_cached <- function(key) if (exists(key, envir=.memo, inherits=FALSE)) get(key, envir=.memo) else NULL
set_cached <- function(key, val) { assign(key, val, envir=.memo); invisible(val) }

# Exposure outcomes use a SINGLE-SLOT cache (not the full .memo): exposure GWAS
# are large genome-wide tables (~750 MB each), and a PE/DE chunk spans ~50
# exposures, so memoising them ALL blows memory (60-100 GB -> OOM, the same
# anti-pattern get_prot_outcome avoids above). The PE/DE loops sort their edges
# by exposure, so a one-deep slot loads each exposure exactly once while its
# consecutive edges run, then evicts it when the exposure changes.
.exp_slot <- new.env(parent=emptyenv())
.exp_slot$id <- NULL; .exp_slot$val <- NULL
get_exp_outcome <- function(ExID) {
  if (!is.null(.exp_slot$id) && identical(.exp_slot$id, ExID)) return(.exp_slot$val)
  .exp_slot$id <- NULL; .exp_slot$val <- NULL  # drop prior large table before loading next
  gc(FALSE)
  g <- read_exposure_gwas_generic(ExID)
  out <- ensure_tsmr_outcome(g, ExID)
  .exp_slot$id <- ExID; .exp_slot$val <- out
  out
}

get_prot_outcome <- function(protID) {
  # NOT memoised: protein outcomes are full genome-wide tables (UKB parquet up to
  # ~1.2 GB) that rarely repeat within a chunk, so caching them all blows memory.
  # Read fresh; the prior table becomes GC-eligible once the caller reassigns.
  ensure_tsmr_outcome(read_protein_gwas_generic(protID), protID)
}

# Disease outcomes ALSO use a SINGLE-SLOT cache: like exposures, disease GWAS are
# large genome-wide tables, and a PD/ED chunk spans ~many diseases, so memoising
# them ALL blows memory (PD chunks OOM'd at 48G). The PD/ED loops sort their edges
# by disease, so a one-deep slot loads each disease once then evicts on change.
.dz_slot <- new.env(parent=emptyenv())
.dz_slot$id <- NULL; .dz_slot$val <- NULL
get_dz_outcome <- function(dzID) {
  if (!is.null(.dz_slot$id) && identical(.dz_slot$id, dzID)) return(.dz_slot$val)
  .dz_slot$id <- NULL; .dz_slot$val <- NULL  # drop prior large table before loading next
  gc(FALSE)
  g <- read_disease_gwas_generic(dzID)
  out <- ensure_tsmr_outcome(g, dzID)
  .dz_slot$id <- dzID; .dz_slot$val <- out
  out
}

# ---------------------------
# Load edges + chunk  (WARM_ONLY skips this: pre-warm only builds instrument caches)
# ---------------------------
if (!WARM_ONLY) {

edge_file <- file.path(CFG$edges_dir, paste0("edges_", edge_type, ".tsv"))
if (!file.exists(edge_file)) stop("Missing edge file: ", edge_file)

edges <- fread(edge_file, showProgress=FALSE)
edges[, chunk := ((.I - 1L) %% n_chunks) + 1L]
todo <- edges[chunk == idx]
if (nrow(todo) == 0) stop("No edges for ", edge_type, " in chunk ", idx, "/", n_chunks)

# PE/DE load a large exposure GWAS as the OUTCOME, PD/ED load a large disease
# GWAS as the OUTCOME; order this chunk's edges by that outcome so the single-slot
# get_exp_outcome / get_dz_outcome cache loads each one once (never holding more
# than one) instead of caching all of them and OOM-ing.
if (edge_type %in% c("PE", "DE")) setorder(todo, Exposure)
if (edge_type %in% c("PD", "ED")) setorder(todo, Disease)

# ---------------------------
# Logging
# ---------------------------
logdir <- file.path(CFG$outdirbase, "logs", edge_type)
dir.create(logdir, recursive=TRUE, showWarnings=FALSE)
logfile <- file.path(logdir, paste0("chunk_", idx, ".log"))
loglines <- character(0)

safe_edge <- function(expr, tag) {
  tryCatch({ eval(expr); paste(tag, "SUCCESS") },
           error=function(e) paste(tag, "FAILED:", e$message))
}

# ---------------------------
# Run edges
# ---------------------------
if (edge_type == "EP") {
  for (i in seq_len(nrow(todo))) {
    ExID <- todo$Exposure[i]; protID <- todo$Protein[i]
    outdir <- file.path(CFG$outdirbase, "E_to_P", ExID, protID)
    sumfile <- file.path(outdir, "E_to_P_summary.tsv")
    if (file.exists(sumfile)) { loglines <- c(loglines, paste("[E->P]", ExID, protID, "SKIP")); next }
    
    tag <- paste("[E->P]", ExID, "->", protID)
    loglines <- c(loglines, safe_edge(quote({
      inst_E <- get_exposure_instruments(ExID)
      out_P  <- get_prot_outcome(protID)
      run_mr_edge(inst_E, out_P, outdir, "E_to_P")
    }), tag))
  }
}

if (edge_type == "ED") {
  for (i in seq_len(nrow(todo))) {
    ExID <- todo$Exposure[i]; dzID <- todo$Disease[i]
    outdir <- file.path(CFG$outdirbase, "E_to_D", ExID, dzID)
    sumfile <- file.path(outdir, "E_to_D_summary.tsv")
    if (file.exists(sumfile)) { loglines <- c(loglines, paste("[E->D]", ExID, dzID, "SKIP")); next }
    
    tag <- paste("[E->D]", ExID, "->", dzID)
    loglines <- c(loglines, safe_edge(quote({
      inst_E <- get_exposure_instruments(ExID)
      out_D  <- get_dz_outcome(dzID)
      run_mr_edge(inst_E, out_D, outdir, "E_to_D")
    }), tag))
  }
}

# PD: protein -> disease (cis + trans)
if (edge_type == "PD") {
  for (i in seq_len(nrow(todo))) {
    protID <- todo$Protein[i]; dzID <- todo$Disease[i]
    out_D <- get_dz_outcome(dzID)
    
    outdir_cis <- file.path(CFG$outdirbase, "Pcis_to_D", protID, dzID)
    sum_cis <- file.path(outdir_cis, "Pcis_to_D_summary.tsv")
    if (!file.exists(sum_cis)) {
      tag <- paste("[Pcis->D]", protID, "->", dzID)
      loglines <- c(loglines, safe_edge(quote({
        inst_Pcis <- get_protein_instruments(protID, "cis")
        run_mr_edge(inst_Pcis, out_D, outdir_cis, "Pcis_to_D")
      }), tag))
    } else loglines <- c(loglines, paste("[Pcis->D]", protID, dzID, "SKIP"))
    
    outdir_tr <- file.path(CFG$outdirbase, "Ptrans_to_D", protID, dzID)
    sum_tr <- file.path(outdir_tr, "Ptrans_to_D_summary.tsv")
    if (!file.exists(sum_tr)) {
      tag <- paste("[Ptrans->D]", protID, "->", dzID)
      loglines <- c(loglines, safe_edge(quote({
        inst_Ptr <- get_protein_instruments(protID, "trans")
        run_mr_edge(inst_Ptr, out_D, outdir_tr, "Ptrans_to_D")
      }), tag))
    } else loglines <- c(loglines, paste("[Ptrans->D]", protID, dzID, "SKIP"))
  }
}

# PE: protein -> exposure (cis + trans)
if (edge_type == "PE") {
  for (i in seq_len(nrow(todo))) {
    protID <- todo$Protein[i]; ExID <- todo$Exposure[i]
    out_E <- get_exp_outcome(ExID)
    
    outdir_cis <- file.path(CFG$outdirbase, "Pcis_to_E", protID, ExID)
    sum_cis <- file.path(outdir_cis, "Pcis_to_E_summary.tsv")
    if (!file.exists(sum_cis)) {
      tag <- paste("[Pcis->E]", protID, "->", ExID)
      loglines <- c(loglines, safe_edge(quote({
        inst_Pcis <- get_protein_instruments(protID, "cis")
        run_mr_edge(inst_Pcis, out_E, outdir_cis, "Pcis_to_E")
      }), tag))
    } else loglines <- c(loglines, paste("[Pcis->E]", protID, ExID, "SKIP"))
    
    outdir_tr <- file.path(CFG$outdirbase, "Ptrans_to_E", protID, ExID)
    sum_tr <- file.path(outdir_tr, "Ptrans_to_E_summary.tsv")
    if (!file.exists(sum_tr)) {
      tag <- paste("[Ptrans->E]", protID, "->", ExID)
      loglines <- c(loglines, safe_edge(quote({
        inst_Ptr <- get_protein_instruments(protID, "trans")
        run_mr_edge(inst_Ptr, out_E, outdir_tr, "Ptrans_to_E")
      }), tag))
    } else loglines <- c(loglines, paste("[Ptrans->E]", protID, ExID, "SKIP"))
  }
}

# DE: disease -> exposure
if (edge_type == "DE") {
  for (i in seq_len(nrow(todo))) {
    dzID <- todo$Disease[i]; ExID <- todo$Exposure[i]
    outdir <- file.path(CFG$outdirbase, "D_to_E", dzID, ExID)
    sumfile <- file.path(outdir, "D_to_E_summary.tsv")
    if (file.exists(sumfile)) { loglines <- c(loglines, paste("[D->E]", dzID, ExID, "SKIP")); next }
    
    tag <- paste("[D->E]", dzID, "->", ExID)
    loglines <- c(loglines, safe_edge(quote({
      inst_D <- get_disease_instruments(dzID)
      out_E  <- get_exp_outcome(ExID)
      run_mr_edge(inst_D, out_E, outdir, "D_to_E")
    }), tag))
  }
}

# DP: disease -> protein
if (edge_type == "DP") {
  for (i in seq_len(nrow(todo))) {
    dzID <- todo$Disease[i]; protID <- todo$Protein[i]
    outdir <- file.path(CFG$outdirbase, "D_to_P", dzID, protID)
    sumfile <- file.path(outdir, "D_to_P_summary.tsv")
    if (file.exists(sumfile)) { loglines <- c(loglines, paste("[D->P]", dzID, protID, "SKIP")); next }
    
    tag <- paste("[D->P]", dzID, "->", protID)
    loglines <- c(loglines, safe_edge(quote({
      inst_D <- get_disease_instruments(dzID)
      out_P  <- get_prot_outcome(protID)
      run_mr_edge(inst_D, out_P, outdir, "D_to_P")
    }), tag))
  }
}

writeLines(loglines, logfile)
cat("DONE Stage B: ", edge_type, " chunk ", idx, "/", n_chunks, "\nLog: ", logfile, "\n", sep="")

}  # end if(!WARM_ONLY)
