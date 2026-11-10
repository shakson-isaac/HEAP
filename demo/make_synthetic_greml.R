#!/usr/bin/env Rscript
# ============================================================================
# demo/make_synthetic_greml.R
# ----------------------------------------------------------------------------
# Synthetic genotypes and proteins for the GREML (population architecture)
# demo. Run after make_synthetic_heap.R has written a HEAP object (exposures and
# covariates in the real structure) for this stage; this script adds
#   <igloo>/UKB/gwas/UKBallchr.{pgen,pvar,psam}   PLINK 2 genotypes, FID = IID = eid
#   <igloo>/UKB/gwas/ukb_nonimputed_snps.pvar     the SNP list the GRM is built on
#   planted_greml.tsv                             the planted variance shares
# and replaces the HEAP object's proteins with
#   y = g + e + gxe + covariate effects + noise
# where g comes from the SNPs, e from the exposures, and gxe has covariance equal
# to (genetic relationship) x (exposure similarity) -- the three kernels GREML
# fits. Shares of phenotypic variance: G 0.30, E 0.25, GxE 0.10, noise 0.35.
# 20 pairs of close relatives are planted so the 0.025 relatedness cutoff
# visibly removes 40 people.
#
# Usage: Rscript demo/make_synthetic_greml.R <stage_dir> <igloo_root> [n_snps] [n_proteins]
#   needs plink2 (PLINK2_BIN or on PATH)
# ============================================================================
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
dir <- a[1]; igloo <- a[2]
M <- if (length(a) >= 3) as.integer(a[3]) else 50000L
K <- if (length(a) >= 4) as.integer(a[4]) else 3L
plink2 <- Sys.getenv("PLINK2_BIN", unset = Sys.which("plink2"))
if (!nzchar(plink2)) stop("plink2 not found: set PLINK2_BIN or put plink2 on PATH")
set.seed(as.integer(Sys.getenv("HEAP_DEMO_SEED", "1")))

rds  <- file.path(dir, "HEAP_synthetic.rds")
heap <- readRDS(rds)
eid  <- heap$prot_baseline$eid; N <- length(eid)
if (N %% 4L) stop("n_participants must be a multiple of 4 (PLINK .bed packing)")
SH   <- c(G = 0.30, E = 0.25, GxE = 0.10)
R    <- 20L                        # draws averaged to build the GxE effect

# ---- exposures as the E kernel sees them: one numeric column per exposure ---
X <- Reduce(function(x, y) merge(x, y, by = "eid", all = TRUE), heap$E_baseline)
X <- as.matrix(X[match(eid, X$eid), setdiff(names(X), "eid")])
X <- scale(X); X[is.na(X)] <- 0; mE <- ncol(X)

# ---- genotypes, in blocks; accumulate g and the genetic factors of gxe -------
gw <- file.path(igloo, "UKB", "gwas"); dir.create(gw, recursive = TRUE, showWarnings = FALSE)
pre <- file.path(dir, "geno_tmp"); con <- file(paste0(pre, ".bed"), "wb")
writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
rel_a <- 1:20; rel_b <- 21:40       # eid 21..40 are close relatives of 1..20
B <- matrix(0, N, K); U <- array(0, c(N, R, K)); blk <- 5000L; p_all <- numeric(0)
for (s in seq(1L, M, by = blk)) {
  b <- min(blk, M - s + 1L); p <- runif(b, 0.1, 0.5); p_all <- c(p_all, p)
  G <- matrix(rbinom(N * b, 2L, rep(p, each = N)), N, b)
  share <- matrix(runif(20 * b) < 0.5, 20, b)          # relatives share ~half the genome
  G[rel_b, ][share] <- G[rel_a, ][share]
  code <- c(0L, 2L, 3L)[G + 1L]                         # PLINK .bed: A2 count 0/1/2
  arr  <- array(code, c(4L, N / 4L, b))
  writeBin(as.raw(arr[1, , ] + 4L * arr[2, , ] + 16L * arr[3, , ] + 64L * arr[4, , ]), con)
  Z <- sweep(sweep(G, 2, 2 * p), 2, sqrt(2 * p * (1 - p)), "/")
  B <- B + Z %*% matrix(rnorm(b * K, 0, 1 / sqrt(M)), b, K)
  for (k in seq_len(K)) U[, , k] <- U[, , k] + Z %*% matrix(rnorm(b * R, 0, 1 / sqrt(M)), b, R)
}
close(con)
ids <- sprintf("rs%07d", seq_len(M)); chr <- sort(rep(1:22, length.out = M))
pos <- ave(seq_len(M), chr, FUN = seq_along) * 1000L
fwrite(data.table(chr, ids, 0, pos, "A", "G"), paste0(pre, ".bim"),
       sep = "\t", col.names = FALSE)
fwrite(data.table(eid, eid, 0, 0, 0, -9), paste0(pre, ".fam"), sep = " ", col.names = FALSE)
st <- system2(plink2, c("--bfile", pre, "--make-pgen", "--out", file.path(gw, "UKBallchr")),
              stdout = FALSE, stderr = FALSE)
if (st != 0) stop("plink2 --make-pgen failed")
# The pipeline reads the first two .psam columns as FID and IID.
fwrite(data.table(`#FID` = eid, IID = eid, SEX = "NA"), file.path(gw, "UKBallchr.psam"), sep = "\t")
fwrite(data.table(`#CHROM` = chr, POS = pos, ID = ids, REF = "A", ALT = "G"),
       file.path(gw, "ukb_nonimputed_snps.pvar"), sep = "\t")
unlink(paste0(pre, c(".bed", ".bim", ".fam")))

# ---- proteins ------------------------------------------------------------------
cv <- heap$covars_baseline[match(eid, heap$covars_baseline$eid), ]
age <- cv$age_when_attended_assessment_centre_f21003_0_0
sex <- as.numeric(as.character(cv$sex_f31_0_0))
unit <- function(v) (v - mean(v)) / sd(v)
for (k in seq_len(K)) {
  g   <- unit(B[, k])
  e   <- unit(X %*% rnorm(mE, 0, 1 / sqrt(mE)))
  gxe <- unit(rowSums(U[, , k] * (X %*% matrix(rnorm(mE * R, 0, 1 / sqrt(mE)), mE, R))) / sqrt(R))
  y <- sqrt(SH["G"]) * g + sqrt(SH["E"]) * e + sqrt(SH["GxE"]) * gxe +
       sqrt(1 - sum(SH)) * rnorm(N) + 0.02 * (age - 55) + 0.2 * sex
  heap$prot_baseline[[paste0("protein_", k)]] <- as.numeric(y)
}
heap$prot_baseline <- heap$prot_baseline[, c("eid", paste0("protein_", seq_len(K)))]
heap$protIDs <- paste0("protein_", seq_len(K))
saveRDS(heap, rds)

fwrite(data.table(protein = paste0("protein_", seq_len(K)), planted_G = SH["G"],
                  planted_E = SH["E"], planted_GxE = SH["GxE"], n = N, n_snps = M,
                  related_pairs = length(rel_a)),
       file.path(dir, "planted_greml.tsv"), sep = "\t")
cat(sprintf("Synthetic GREML data: %d participants (%d related pairs), %d SNPs, %d proteins\n",
            N, length(rel_a), M, K))
