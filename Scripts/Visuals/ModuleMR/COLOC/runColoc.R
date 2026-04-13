#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

# ============================================================
# Minimal coloc runner anchored on a lead SNP (± window_kb)
#
# Required:
#   --pqtl      DECODE SomaScan pQTL file (one protein)
#   --gwas      FinnGen GWAS file (one endpoint)
#   --snp       rsID anchor
#   --out       output prefix
#
# Optional:
#   --window_kb 500
#   --s         FinnGen case fraction (cases/(cases+controls)) [recommended for cc]
#   --n_gwas    FinnGen total sample size N (REQUIRED if GWAS file lacks N)
#   --p1 --p2 --p12 coloc priors
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  w <- which(args == flag)
  if (length(w) == 0) return(default)
  if (w == length(args)) stop("Missing value for ", flag)
  args[w + 1]
}

pqtl_path  <- get_arg("--pqtl")
gwas_path  <- get_arg("--gwas")
lead_snp   <- get_arg("--snp")
window_kb  <- suppressWarnings(as.numeric(get_arg("--window_kb", "500")))
out_prefix <- get_arg("--out", paste0("coloc_", lead_snp))

s_casefrac <- suppressWarnings(as.numeric(get_arg("--s", NA_character_)))
if (!is.finite(s_casefrac)) s_casefrac <- NA_real_

n_gwas_override <- suppressWarnings(as.numeric(get_arg("--n_gwas", NA_character_)))
if (!is.finite(n_gwas_override)) n_gwas_override <- NA_real_

p1  <- suppressWarnings(as.numeric(get_arg("--p1", "1e-4")))
p2  <- suppressWarnings(as.numeric(get_arg("--p2", "1e-4")))
p12 <- suppressWarnings(as.numeric(get_arg("--p12", "1e-5")))

if (is.null(pqtl_path) || is.null(gwas_path) || is.null(lead_snp) || is.null(out_prefix)) {
  stop(
    "Usage:\n",
    "  Rscript runColoc.R --pqtl <DECODE.gz> --gwas <FinnGen.gz> --snp rs55714927 --window_kb 500 --s 0.10 --n_gwas 300000 --out /path/prefix\n"
  )
}
if (!file.exists(pqtl_path)) stop("pQTL file not found: ", pqtl_path)
if (!file.exists(gwas_path)) stop("GWAS file not found: ", gwas_path)
if (!is.finite(window_kb) || window_kb <= 0) stop("--window_kb must be positive")

set_if_present <- function(DT, from, to) {
  if (from %in% names(DT) && !(to %in% names(DT))) setnames(DT, from, to)
  DT
}

read_decode_pqtl <- function(path) {
  DT <- fread(path, showProgress = FALSE)
  
  if ("rsids" %in% names(DT)) {
    DT[, SNP := tstrsplit(rsids, ",", fixed=TRUE, keep=1)]
    DT[is.na(SNP) | SNP == "" | SNP == ".", SNP := NA_character_]
  }
  if (!"SNP" %in% names(DT) || all(is.na(DT$SNP))) {
    if ("Name" %in% names(DT) && !"SNP" %in% names(DT)) setnames(DT, "Name", "SNP")
  }
  
  set_if_present(DT, "Chrom", "chr")
  if ("chr" %in% names(DT)) DT[, chr := gsub("^chr", "", chr)]
  set_if_present(DT, "Pos", "pos")
  
  set_if_present(DT, "Beta", "beta")
  set_if_present(DT, "SE",   "se")
  set_if_present(DT, "Pval", "pval")
  
  set_if_present(DT, "effectAllele", "effect_allele")
  set_if_present(DT, "otherAllele",  "other_allele")
  
  set_if_present(DT, "ImpMAF", "eaf")
  
  set_if_present(DT, "N", "N")
  if ("N" %in% names(DT)) {
    DT[, N := gsub(",", "", as.character(N))]
    DT[, N := suppressWarnings(as.numeric(N))]
  }
  
  DT[, chr := as.character(chr)]
  DT[, pos := as.integer(pos)]
  DT
}

read_finngen_gwas <- function(path) {
  DT <- fread(path, showProgress = FALSE)
  
  if (!"rsid" %in% names(DT)) {
    if ("rsids" %in% names(DT)) DT[, rsid := tstrsplit(rsids, ",", fixed=TRUE, keep=1)]
    if (!"rsid" %in% names(DT)) DT[, rsid := paste0(get("#chrom"), ":", pos, "_", ref, "_", alt)]
  }
  setnames(DT, "rsid", "SNP")
  
  set_if_present(DT, "beta",   "beta")
  set_if_present(DT, "sebeta", "se")
  set_if_present(DT, "pval",   "pval")
  set_if_present(DT, "alt",    "effect_allele")
  set_if_present(DT, "ref",    "other_allele")
  set_if_present(DT, "af_alt", "eaf")
  set_if_present(DT, "#chrom", "chr")
  set_if_present(DT, "pos",    "pos")
  
  # some releases have n_total; R12 may not
  set_if_present(DT, "n_total", "N")
  if ("N" %in% names(DT)) {
    DT[, N := gsub(",", "", as.character(N))]
    DT[, N := suppressWarnings(as.numeric(N))]
  }
  
  DT[, chr := gsub("^chr", "", as.character(chr))]
  DT[, pos := as.integer(pos)]
  DT
}

infer_scalar_N <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

is_palindromic <- function(a1, a2) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  (a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C")
}

to_maf <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  if (!is.finite(v)) return(NA_real_)
  min(v, 1 - v)
}

dedup_by_snp <- function(DT) {
  DT <- as.data.table(DT)
  if ("pval" %in% names(DT)) DT[, pval := suppressWarnings(as.numeric(pval))]
  DT[, beta_num := suppressWarnings(as.numeric(beta))]
  DT[, se_num   := suppressWarnings(as.numeric(se))]
  DT[, absb := abs(beta_num)]
  DT[, p_rank := if ("pval" %in% names(DT)) pval else NA_real_]
  DT[, hasN := if ("N" %in% names(DT)) as.integer(is.finite(suppressWarnings(as.numeric(N)))) else 0L]
  setorder(DT, SNP, -hasN, p_rank, -absb, se_num, na.last = TRUE)
  out <- DT[, .SD[1], by = SNP]
  drop_cols <- intersect(c("beta_num","se_num","absb","p_rank","hasN"), names(out))
  if (length(drop_cols)) out[, (drop_cols) := NULL]
  out
}

harmonize_two <- function(d1, d2) {
  d1 <- dedup_by_snp(d1)
  d2 <- dedup_by_snp(d2)
  
  m <- merge(d1, d2, by="SNP", suffixes=c(".1",".2"), allow.cartesian = FALSE)
  if (nrow(m) == 0) return(NULL)
  
  m[, `:=`(
    effect_allele.1 = toupper(effect_allele.1),
    other_allele.1  = toupper(other_allele.1),
    effect_allele.2 = toupper(effect_allele.2),
    other_allele.2  = toupper(other_allele.2)
  )]
  
  same <- (m$effect_allele.1 == m$effect_allele.2) & (m$other_allele.1 == m$other_allele.2)
  swap <- (m$effect_allele.1 == m$other_allele.2)  & (m$other_allele.1 == m$effect_allele.2)
  
  keep <- same | swap
  m <- m[keep]
  if (nrow(m) == 0) return(NULL)
  
  swap <- swap[keep]
  if (any(swap)) m[swap, beta.2 := -beta.2]
  
  pal <- is_palindromic(m$effect_allele.1, m$other_allele.1)
  if ("eaf.1" %in% names(m)) {
    e1 <- suppressWarnings(as.numeric(m$eaf.1))
    amb <- pal & is.finite(e1) & e1 > 0.42 & e1 < 0.58
    m <- m[!amb]
  } else {
    m <- m[!pal]
  }
  
  m
}

make_dataset <- function(m, which=c("1","2"), type=c("quant","cc"), s=NULL) {
  which <- match.arg(which)
  type  <- match.arg(type)
  
  beta <- suppressWarnings(as.numeric(m[[paste0("beta.", which)]]))
  se   <- suppressWarnings(as.numeric(m[[paste0("se.", which)]]))
  varb <- se^2
  
  eaf_col <- paste0("eaf.", which)
  maf <- if (eaf_col %in% names(m)) vapply(m[[eaf_col]], to_maf, numeric(1)) else rep(NA_real_, nrow(m))
  
  ds <- list(
    snp = m$SNP,
    beta = beta,
    varbeta = varb,
    MAF = maf,
    N = NA_real_,
    type = type
  )
  if (type == "cc") ds$s <- s
  ds
}

# ============================================================
# Main
# ============================================================
pqtl <- read_decode_pqtl(pqtl_path)
gwas <- read_finngen_gwas(gwas_path)

req_cols <- c("SNP","chr","pos","beta","se","effect_allele","other_allele")
if (length(setdiff(req_cols, names(pqtl)))) stop("pQTL missing cols: ", paste(setdiff(req_cols, names(pqtl)), collapse=", "))
if (length(setdiff(req_cols, names(gwas)))) stop("GWAS missing cols: ", paste(setdiff(req_cols, names(gwas)), collapse=", "))

lead_row <- pqtl[SNP == lead_snp]
if (nrow(lead_row) == 0) lead_row <- gwas[SNP == lead_snp]
if (nrow(lead_row) == 0) stop("Lead SNP not found in either file: ", lead_snp)

chr0 <- as.character(lead_row$chr[1])
pos0 <- as.integer(lead_row$pos[1])
lo <- max(1L, pos0 - as.integer(window_kb * 1000))
hi <- pos0 + as.integer(window_kb * 1000)

pqtl_reg <- pqtl[chr == chr0 & pos >= lo & pos <= hi]
gwas_reg <- gwas[chr == chr0 & pos >= lo & pos <= hi]

cat("Lead SNP:", lead_snp, "chr", chr0, "pos", pos0, "\n")
cat("Window:", window_kb, "kb =>", chr0, ":", lo, "-", hi, "\n")
cat("pQTL SNPs in region:", nrow(pqtl_reg), "\n")
cat("GWAS SNPs in region:", nrow(gwas_reg), "\n")

# Force dataset1 N from DECODE (scalar mode)
N_pqtl_scalar <- infer_scalar_N(pqtl_reg$N)
cat("Inferred pQTL N (scalar):", N_pqtl_scalar, "\n")
if (!is.finite(N_pqtl_scalar)) stop("Could not infer pQTL N from DECODE file.")

m <- harmonize_two(pqtl_reg, gwas_reg)
if (is.null(m) || nrow(m) < 50) {
  dir.create(dirname(out_prefix), recursive=TRUE, showWarnings=FALSE)
  fwrite(data.table(
    lead_snp=lead_snp, chr=chr0, pos=pos0, lo=lo, hi=hi,
    status="TOO_FEW_OVERLAP_SNPS", nsnps=ifelse(is.null(m), 0L, nrow(m)),
    PP.H0=NA_real_, PP.H1=NA_real_, PP.H2=NA_real_, PP.H3=NA_real_, PP.H4=NA_real_
  ), paste0(out_prefix, "_coloc_summary.tsv"), sep="\t")
  if (!is.null(m)) fwrite(m, paste0(out_prefix, "_harmonized_snps.tsv"), sep="\t")
  stop("Too few overlapping SNPs after harmonization.")
}

# Build datasets
d1 <- make_dataset(m, which="1", type="quant")
d1$N <- N_pqtl_scalar

disease_type <- if (is.finite(s_casefrac)) "cc" else "quant"
d2 <- if (disease_type == "cc") make_dataset(m, which="2", type="cc", s=s_casefrac) else make_dataset(m, which="2", type="quant")

# Force dataset2 N:
N_gwas_scalar <- infer_scalar_N(m$N.2)
if (!is.finite(N_gwas_scalar) && is.finite(n_gwas_override)) N_gwas_scalar <- n_gwas_override

if (!is.finite(N_gwas_scalar)) {
  stop("dataset 2 (GWAS): N is missing. Provide --n_gwas <total sample size> for this FinnGen endpoint.")
}
d2$N <- N_gwas_scalar

cat("DEBUG: d1$N =", d1$N, " is.finite =", is.finite(d1$N), "\n")
cat("DEBUG: d2 type =", disease_type, " s =", ifelse(is.finite(s_casefrac), s_casefrac, NA_real_), " N =", d2$N, "\n")

if (disease_type == "cc" && (!is.finite(d2$s) || d2$s <= 0 || d2$s >= 1)) {
  stop("Case fraction s must be in (0,1). Provide --s cases/(cases+controls).")
}

res <- coloc::coloc.abf(dataset1=d1, dataset2=d2, p1=p1, p2=p2, p12=p12)

s <- res$summary
out <- data.table(
  lead_snp=lead_snp,
  chr=chr0, pos=pos0, lo=lo, hi=hi,
  window_kb=window_kb,
  pqtl_file=pqtl_path,
  gwas_file=gwas_path,
  disease_type=disease_type,
  s_casefrac=ifelse(is.finite(s_casefrac), s_casefrac, NA_real_),
  nsnps=nrow(m),
  N_pqtl=N_pqtl_scalar,
  N_gwas=N_gwas_scalar,
  p1=p1, p2=p2, p12=p12,
  PP.H0=unname(s["PP.H0.abf"]),
  PP.H1=unname(s["PP.H1.abf"]),
  PP.H2=unname(s["PP.H2.abf"]),
  PP.H3=unname(s["PP.H3.abf"]),
  PP.H4=unname(s["PP.H4.abf"])
)

dir.create(dirname(out_prefix), recursive=TRUE, showWarnings=FALSE)
fwrite(out, paste0(out_prefix, "_coloc_summary.tsv"), sep="\t")
fwrite(m,   paste0(out_prefix, "_harmonized_snps.tsv"), sep="\t")

cat("\nCOLOC DONE\n")
print(out)
cat("\nWrote:\n  ",
    paste0(out_prefix, "_coloc_summary.tsv"), "\n  ",
    paste0(out_prefix, "_harmonized_snps.tsv"), "\n", sep="")