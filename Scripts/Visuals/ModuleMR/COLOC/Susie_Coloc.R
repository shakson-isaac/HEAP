#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(ggplot2)
})

# ============================================================
# SuSiE colocalization with LD from PLINK reference
# - Reads harmonized SNP table (from your coloc.abf run)
# - Intersects SNPs with LD reference via PLINK --write-snplist
# - Computes LD matrix (R) with PLINK
# - Runs coloc::runsusie() for both traits (returns coloc-compatible SuSiE fits)
# - Runs coloc::coloc.susie()
#
# Args:
#   --harm     <*_harmonized_snps.tsv>
#   --bfile    <PLINK LD reference prefix> (e.g., .../LDref/EUR)
#   --plink    <plink binary path> (default: plink)
#   --out      <output prefix>
#   --N1       trait1 sample size (pQTL)
#   --N2       trait2 sample size (GWAS)
#
# Optional:
#   --L        max signals (default 5)
#   --maf      MAF filter (default 0.01)
#   --s        case fraction for trait2 if case-control (FinnGen) (recommended)
#   --p1 --p2 --p12  coloc priors
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  w <- which(args == flag)
  if (length(w) == 0) return(default)
  if (w == length(args)) stop("Missing value for ", flag)
  args[w + 1]
}

harm_path <- get_arg("--harm")
bfile     <- get_arg("--bfile")
plink_bin <- get_arg("--plink", "plink")
out_pref  <- get_arg("--out", "susie_coloc")

N1 <- suppressWarnings(as.numeric(get_arg("--N1", NA_character_)))
N2 <- suppressWarnings(as.numeric(get_arg("--N2", NA_character_)))
if (!is.finite(N1) || !is.finite(N2)) stop("You must supply --N1 and --N2 as numeric sample sizes.")

L <- suppressWarnings(as.integer(get_arg("--L", "5")))
if (!is.finite(L) || L < 1) L <- 5L

maf_thr <- suppressWarnings(as.numeric(get_arg("--maf", "0.01")))
if (!is.finite(maf_thr) || maf_thr < 0) maf_thr <- 0.01

s_casefrac <- suppressWarnings(as.numeric(get_arg("--s", NA_character_)))
if (!is.finite(s_casefrac)) s_casefrac <- NA_real_

p1  <- suppressWarnings(as.numeric(get_arg("--p1", "1e-4")))
p2  <- suppressWarnings(as.numeric(get_arg("--p2", "1e-4")))
p12 <- suppressWarnings(as.numeric(get_arg("--p12", "1e-5")))

if (is.null(harm_path) || is.null(bfile)) {
  stop(
    "Usage:\n",
    "  Rscript Susie_Coloc.R --harm <*_harmonized_snps.tsv> --bfile <LDref/EUR> --plink <plink> --N1 35704 --N2 300000 --s 0.10 --L 2 --maf 0.01 --out /path/prefix\n"
  )
}
if (!file.exists(harm_path)) stop("Missing --harm file: ", harm_path)

dir.create(dirname(out_pref), recursive=TRUE, showWarnings=FALSE)

# ---------------------------
# Load harmonized SNPs
# ---------------------------
m <- fread(harm_path)

need <- c("SNP","beta.1","se.1","beta.2","se.2")
miss <- setdiff(need, names(m))
if (length(miss)) stop("harm file missing required columns: ", paste(miss, collapse=", "))

# chr/pos
if (!("chr.1" %in% names(m) && "pos.1" %in% names(m))) {
  if (!all(c("chr","pos") %in% names(m))) stop("harm file needs chr/pos columns (chr.1/pos.1 or chr/pos).")
  m[, `:=`(chr.1 = chr, pos.1 = pos)]
}

# MAF from eaf if present
if ("eaf.1" %in% names(m)) {
  m[, maf1 := pmin(as.numeric(eaf.1), 1 - as.numeric(eaf.1))]
} else {
  m[, maf1 := NA_real_]
}
if ("eaf.2" %in% names(m)) {
  m[, maf2 := pmin(as.numeric(eaf.2), 1 - as.numeric(eaf.2))]
} else {
  m[, maf2 := NA_real_]
}

# Choose a MAF column for filtering: prefer maf1, else maf2, else no filter
maf_use <- if (any(is.finite(m$maf1))) "maf1" else if (any(is.finite(m$maf2))) "maf2" else NA_character_
if (!is.na(maf_use)) {
  m <- m[is.finite(get(maf_use)) & get(maf_use) >= maf_thr]
}

# Keep finite stats
m <- m[is.finite(beta.1) & is.finite(se.1) & is.finite(beta.2) & is.finite(se.2)]
m <- m[!duplicated(SNP)]
if (nrow(m) < 200) stop("Too few SNPs after filtering for SuSiE: ", nrow(m))

chr0 <- as.character(m$chr.1[1])
lo <- min(m$pos.1, na.rm=TRUE)
hi <- max(m$pos.1, na.rm=TRUE)

cat("SNPs for SuSiE (pre-LDref intersect):", nrow(m), "\n")
cat("Region:", chr0, ":", lo, "-", hi, "\n")

# ---------------------------
# Intersect SNPs with LD reference (PLINK --write-snplist)
# ---------------------------
snplist_in <- tempfile(pattern="susie_snps_in_", fileext=".txt")
fwrite(m[, .(SNP)], snplist_in, col.names=FALSE)

tmp_pref1 <- tempfile(pattern="plink_snplist_")
cmd1 <- sprintf(
  "%s --bfile %s --extract %s --chr %s --from-bp %d --to-bp %d --write-snplist --out %s",
  shQuote(plink_bin), shQuote(bfile), shQuote(snplist_in),
  shQuote(chr0), as.integer(lo), as.integer(hi), shQuote(tmp_pref1)
)
cat("Running PLINK snplist intersect:\n", cmd1, "\n")
ret <- system(cmd1)
if (ret != 0) stop("PLINK snplist step failed with exit code ", ret)

kept_snps_file <- paste0(tmp_pref1, ".snplist")
if (!file.exists(kept_snps_file)) stop("PLINK did not create .snplist: ", kept_snps_file)

kept <- fread(kept_snps_file, header=FALSE)[[1]]
kept <- kept[!is.na(kept) & kept != ""]
cat("SNPs present in LDref:", length(kept), "\n")
if (length(kept) < 200) stop("Too few SNPs found in LD reference after intersect: ", length(kept))

# Filter and ORDER exactly as kept list
setkey(m, SNP)
m2 <- m[J(kept), nomatch=0]
cat("SNPs used for SuSiE (post-LDref intersect):", nrow(m2), "\n")
if (nrow(m2) < 200) stop("Too few SNPs remain after ordering/intersection.")

# ---------------------------
# Compute LD matrix (R) on kept SNPs
# ---------------------------
snplist_kept <- tempfile(pattern="susie_snps_kept_", fileext=".txt")
fwrite(data.table(SNP = m2$SNP), snplist_kept, col.names=FALSE)

ld_prefix <- tempfile(pattern="plink_ld_")
cmd2 <- sprintf(
  "%s --bfile %s --extract %s --chr %s --from-bp %d --to-bp %d --r square gz --out %s",
  shQuote(plink_bin), shQuote(bfile), shQuote(snplist_kept),
  shQuote(chr0), as.integer(lo), as.integer(hi), shQuote(ld_prefix)
)
cat("Running PLINK LD:\n", cmd2, "\n")
ret <- system(cmd2)
if (ret != 0) stop("PLINK LD step failed with exit code ", ret)

ld_file <- paste0(ld_prefix, ".ld.gz")
if (!file.exists(ld_file)) stop("PLINK did not produce LD file: ", ld_file)

R <- as.matrix(fread(cmd = paste("zcat", shQuote(ld_file)), header=FALSE))
# PLINK --r square gz has no SNP labels; assign them from the exact SNP order used in --extract
stopifnot(nrow(R) == nrow(m2))
colnames(R) <- m2$SNP
rownames(R) <- m2$SNP
storage.mode(R) <- "double"
diag(R) <- 1
R[lower.tri(R)] <- t(R)[lower.tri(R)]

if (nrow(R) != nrow(m2)) stop("LD matrix dimension mismatch: LD=", nrow(R), " data=", nrow(m2))

# Name LD matrix
# ---- Robust SNP/LD alignment ----
# sanitize SNP IDs (kills hidden whitespace mismatches)
m2[, SNP := trimws(as.character(SNP))]

# Ensure LD has dimnames; if it already does, sanitize them too
if (is.null(colnames(R))) {
  # safest: use current m2 order (but we'll still validate)
  colnames(R) <- m2$SNP
  rownames(R) <- m2$SNP
} else {
  colnames(R) <- trimws(as.character(colnames(R)))
  rownames(R) <- trimws(as.character(rownames(R)))
}

# Intersect to SNPs present in LD (coloc::runsusie requires LD names contain all snps)
snps_ld <- colnames(R)
snps_dt <- m2$SNP
common <- intersect(snps_dt, snps_ld)

cat("DEBUG: SNPs in m2:", length(snps_dt), " SNPs in LD:", length(snps_ld),
    " common:", length(common), "\n")

if (length(common) < 200) stop("Too few SNPs overlap between summary stats and LD after name matching.")

# Reorder m2 to match LD column order (critical!)
# Keep only SNPs in LD, in exactly the LD order
m2 <- m2[SNP %in% common]
setkey(m2, SNP)
m2 <- m2[J(snps_ld), nomatch=0]

# Now subset/reorder R to match m2 (double safety)
R <- R[m2$SNP, m2$SNP, drop=FALSE]

# Final sanity checks
stopifnot(all(m2$SNP == colnames(R)))
stopifnot(all(m2$SNP == rownames(R)))
cat("DEBUG: Alignment OK. Final SNP count:", nrow(m2), "\n")

# Make R symmetric PSD (helps stability)
R <- (R + t(R)) / 2
ev <- eigen(R, symmetric=TRUE)
ev$values[ev$values < 1e-8] <- 1e-8
R <- ev$vectors %*% (ev$values * t(ev$vectors))
diag(R) <- 1
R <- (R + t(R)) / 2
stopifnot(nrow(R) == nrow(m2))
colnames(R) <- m2$SNP
rownames(R) <- m2$SNP

# ---------------------------
# Build coloc datasets for runsusie
# ---------------------------
maf1 <- if (any(is.finite(m2$maf1))) as.numeric(m2$maf1) else if (any(is.finite(m2$maf2))) as.numeric(m2$maf2) else rep(NA_real_, nrow(m2))
maf2 <- if (any(is.finite(m2$maf2))) as.numeric(m2$maf2) else if (any(is.finite(m2$maf1))) as.numeric(m2$maf1) else rep(NA_real_, nrow(m2))

d1 <- list(
  snp = m2$SNP,
  beta = as.numeric(m2$beta.1),
  varbeta = (as.numeric(m2$se.1))^2,
  MAF = maf1,
  N = N1,
  type = "quant"
)

d2 <- list(
  snp = m2$SNP,
  beta = as.numeric(m2$beta.2),
  varbeta = (as.numeric(m2$se.2))^2,
  MAF = maf2,
  N = N2,
  type = if (is.finite(s_casefrac)) "cc" else "quant"
)
if (d2$type == "cc") {
  if (!(s_casefrac > 0 && s_casefrac < 1)) stop("--s must be in (0,1) for case-control trait2.")
  d2$s <- s_casefrac
}

# ---------------------------
# Run SuSiE via coloc wrapper + coloc.susie
# ---------------------------
stopifnot(nrow(R) == length(m2$SNP))
stopifnot(all(rownames(R) == m2$SNP))
stopifnot(all(colnames(R) == m2$SNP))

# Store LD inside the dataset lists (required by your coloc::runsusie version)
d1$LD <- R
d2$LD <- R

cat("Running runsusie trait1 (L=", L, ")...\n")
fit1 <- coloc::runsusie(d1, L = L)

cat("Running runsusie trait2 (L=", L, ")...\n")
fit2 <- coloc::runsusie(d2, L = L)

cat("Running coloc.susie...\n")
cs <- coloc::coloc.susie(fit1, fit2, p1=p1, p2=p2, p12=p12)

# ---------------------------
# Output
# ---------------------------
sum_file <- paste0(out_pref, "_susie_coloc_summary.tsv")
res_file <- paste0(out_pref, "_susie_coloc_results.tsv")
pip_file <- paste0(out_pref, "_pip_tracks.tsv")

fwrite(as.data.table(cs$summary), sum_file, sep="\t")
fwrite(as.data.table(cs$results), res_file, sep="\t")

# PIPs are in susie fits
pip_dt <- data.table(
  SNP = m2$SNP,
  chr = m2$chr.1,
  pos = m2$pos.1,
  pip_trait1 = fit1$pip,
  pip_trait2 = fit2$pip
)
fwrite(pip_dt, pip_file, sep="\t")

cat("Wrote:\n ", sum_file, "\n ", res_file, "\n ", pip_file, "\n", sep="")

# ---------------------------
# Quick plots
# ---------------------------
p_a <- ggplot(pip_dt, aes(x=pos, y=pip_trait1)) + geom_point(size=0.7) +
  labs(title="SuSiE PIP (Trait 1)", x=paste0("chr", chr0, " position"), y="PIP")
ggsave(paste0(out_pref, "_pip_trait1.png"), p_a, width=7, height=3.5, dpi=300)

p_b <- ggplot(pip_dt, aes(x=pos, y=pip_trait2)) + geom_point(size=0.7) +
  labs(title="SuSiE PIP (Trait 2)", x=paste0("chr", chr0, " position"), y="PIP")
ggsave(paste0(out_pref, "_pip_trait2.png"), p_b, width=7, height=3.5, dpi=300)

# Locuscompare-style scatter if p-values exist
if ("pval.1" %in% names(m2) && "pval.2" %in% names(m2)) {
  mp <- copy(m2)
  mp[, `:=`(
    logp1 = -log10(pmax(as.numeric(pval.1), 1e-300)),
    logp2 = -log10(pmax(as.numeric(pval.2), 1e-300))
  )]
  p_c <- ggplot(mp, aes(x=logp1, y=logp2)) + geom_point(size=0.7) +
    labs(title="LocusCompare-style plot", x="-log10(p) Trait 1", y="-log10(p) Trait 2")
  ggsave(paste0(out_pref, "_locuscompare.png"), p_c, width=5.5, height=5, dpi=300)
}

cat("DONE.\n")