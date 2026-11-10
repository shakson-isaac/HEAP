#!/usr/bin/env Rscript
# ============================================================================
# demo/make_synthetic_gwas.R
# ----------------------------------------------------------------------------
# Synthetic GWAS summary statistics for the Mendelian randomization demo, in the
# file formats scripts/module5_mr/Module5.R (UK Biobank pQTL arm) and
# Module5_deCODE.R (deCODE arm) read, plus the LD reference panel they clump
# against. The unmodified MR code then selects instruments (p < 5e-8, F > 10),
# clumps them (r2 < 0.001 within 10 Mb), splits cis from trans pQTLs, runs
# IVW / Egger / weighted median / Steiger / heterogeneity / MR-PRESSO, and the
# downstream scripts apply the paper's tiering and triad motifs.
#
# Genetics: a reference panel of 10,000 individuals x ~900 SNPs on 3
# chromosomes (LD blocks). Every trait's summary statistics are drawn from that
# panel's LD: z = sqrt(N) * R b + MVN(0, R), with b the trait's causal effects.
#
# Four planted triads, exposure_k -> protein_k -> disease_k:
#   1  mediator        E1 -> P1 (0.4 SD) -> D1 (0.3 SD)
#   2  reverse         P2 -> E2 (0.5 SD): the protein causes the exposure
#   3  unresolved      E3 -> P3 with P3 almost entirely determined by E3, so the
#                      instruments explain the same variance in both (Steiger
#                      cannot establish a direction)
#   4  failed coloc    as triad 1, but the given colocalization result for
#                      P4 -> D4 is PP.H4 = 0.30 (below the 0.8 gate)
# plus 20 null exposures (exposure_5 ... exposure_24) tested against protein_1.
# Colocalization is not computed in the demo; coloc_results.tsv is given.
#
# Usage: Rscript demo/make_synthetic_gwas.R <igloo_root> <project_root>
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(arrow) })
a <- commandArgs(trailingOnly = TRUE); igloo <- a[1]; proj <- a[2]
out <- file.path(proj, "output")
set.seed(20261113)
mk <- function(...) { d <- file.path(...); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

N_REF <- 10000L; N_E <- 20000L; N_P_UKB <- 30000L; N_P_DEC <- 35000L
N_CASE <- 5000L; N_CTRL <- 50000L; N_D <- N_CASE + N_CTRL

# ---- LD blocks ---------------------------------------------------------------
genes <- data.table(protein = paste0("protein_", 1:4), chr = c(1L, 1L, 2L, 2L),
                    center = c(20e6, 80e6, 20e6, 80e6))
roles <- c(paste0("cis_P", 1:4),
           paste0("E1_h", 1:4), paste0("E3_h", 1:4), paste0("E4_h", 1:4),
           paste0("P2_t", 1:3),
           as.vector(outer(paste0("D", 1:4, "_h"), 1:3, paste0)),
           paste0("null_h", 1:6), paste0("noise_", 1:11))
slots <- CJ(chr = 1:3, pos = seq(2e6, 140e6, by = 3e6))
slots <- slots[!slots[, sapply(seq_len(.N), function(i)
  any(genes$chr == chr[i] & abs(genes$center - pos[i]) < 2.5e6))]]
small <- roles[!grepl("^cis_", roles)]
slots <- slots[sample(.N, length(small))]
B <- rbind(
  genes[, .(role = paste0("cis_P", 1:4), chr, start = center - 490e3, m = 60L, step = 980e3 / 59, rho = 0.95)],
  data.table(role = small, chr = slots$chr, start = slots$pos, m = 15L, step = 1e4, rho = 0.90))
setorder(B, chr, start); B[, block := .I]          # PLINK needs SNPs sorted by position

# ---- reference panel genotypes ----------------------------------------------
pairs <- list(c("A","C"), c("A","G"), c("C","A"), c("C","T"),
              c("G","A"), c("G","T"), c("T","C"), c("T","G"))     # no A/T or C/G
snp <- list(); G <- list(); R <- list()
for (b in seq_len(nrow(B))) {
  m <- B$m[b]; rho <- B$rho[b]; p <- runif(m, 0.10, 0.45)
  hap <- function() {
    z <- matrix(rnorm(N_REF * m), N_REF, m)
    for (j in 2:m) z[, j] <- rho * z[, j - 1] + sqrt(1 - rho^2) * z[, j]
    sweep(z, 2, qnorm(1 - p), ">") * 1L
  }
  g <- hap() + hap()
  al <- pairs[sample(length(pairs), m, replace = TRUE)]
  snp[[b]] <- data.table(block = b, chr = B$chr[b], pos = as.integer(round(B$start[b] + (seq_len(m) - 1) * B$step[b])),
                         alt = vapply(al, `[`, "", 1), ref = vapply(al, `[`, "", 2), af = colMeans(g) / 2)
  G[[b]] <- g; R[[b]] <- cor(g)
}
S <- rbindlist(snp); S[, `:=`(idx = .I, rsid = paste0("rs", 100000 + .I))]
M <- nrow(S)
lead <- function(r) { b <- B$block[B$role == r]; S$idx[S$block == b][ceiling(B$m[b] / 2)] }

# PLINK 1 binary panel (SNP-major; A1 = alt). Codes: 2 alt -> 00, 1 -> 10, 0 -> 11.
ld <- mk(igloo, "LDref")
fwrite(data.table(i = 1:N_REF, j = 1:N_REF, 0, 0, 0, -9), file.path(ld, "EUR.fam"), sep = " ", col.names = FALSE)
fwrite(S[, .(chr, rsid, 0, pos, alt, ref)], file.path(ld, "EUR.bim"), sep = "\t", col.names = FALSE)
con <- file(file.path(ld, "EUR.bed"), "wb"); writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
code <- c(3L, 2L, 0L)                                   # index by genotype + 1
for (b in seq_along(G)) for (j in seq_len(ncol(G[[b]]))) {
  cc <- matrix(code[G[[b]][, j] + 1L], nrow = 4)
  writeBin(as.raw(cc[1, ] + cc[2, ] * 4L + cc[3, ] * 16L + cc[4, ] * 64L), con)
}
close(con)

# ---- causal effects (standardized) -------------------------------------------
eff <- list(); z0 <- numeric(M)
put <- function(v, roles, size) { for (r in roles) v[lead(r)] <- v[lead(r)] + size; v }
hE <- sqrt(0.005); hCis <- sqrt(0.03); hD <- sqrt(0.003)
for (k in c(1, 3, 4)) eff[[paste0("exposure_", k)]] <- put(z0, paste0("E", k, "_h", 1:4), hE)
eff$protein_1 <- put(0.40 * eff$exposure_1, "cis_P1", hCis)
eff$protein_3 <- put(0.98 * eff$exposure_3, "cis_P3", sqrt(0.01))
eff$protein_4 <- put(0.40 * eff$exposure_4, "cis_P4", hCis)
eff$protein_2 <- put(put(z0, "cis_P2", hCis), paste0("P2_t", 1:3), sqrt(0.02))
eff$exposure_2 <- 0.50 * eff$protein_2
for (k in 1:4) {
  own <- put(z0, paste0("D", k, "_h", 1:3), hD)
  eff[[paste0("disease_", k)]] <- if (k == 2) own else own + 0.30 * eff[[paste0("protein_", k)]]
}
for (j in 1:20) eff[[paste0("exposure_", 4 + j)]] <- put(z0, paste0("null_h", sample(6, 3)), hE)

sim_z <- function(b, N) {                     # marginal z from LD: sqrt(N) R b + MVN(0, R)
  unlist(lapply(seq_len(nrow(B)), function(k) {
    i <- S[block == k, idx]; Rk <- R[[k]]
    sqrt(N) * drop(Rk %*% b[i]) + drop(crossprod(chol(Rk + diag(1e-6, nrow(Rk))), rnorm(length(i))))
  }))
}
stats <- function(b, N) {
  z <- sim_z(b, N); se <- 1 / sqrt(2 * S$af * (1 - S$af) * N)
  lp <- -(pnorm(-abs(z), log.p = TRUE) + log(2)) / log(10)
  data.table(S[, .(chr, pos, rsid, alt, ref, af)], N = N, beta = z * se, se = se,
             log10p = lp, p = pmax(10^-lp, 1e-300))
}

# ---- exposure GWAS: regenie step-2 format -----------------------------------
EX <- paste0("exposure_", 1:24)
for (e in EX) {
  d <- stats(eff[[e]], N_E); dir <- mk(out, "gwas", "regenie_step2", e)
  fwrite(d[, .(CHROM = chr, GENPOS = pos, ID = rsid, ALLELE0 = ref, ALLELE1 = alt,
               A1FREQ = af, N, BETA = beta, SE = se, LOG10P = log10p)],
         file.path(dir, paste0(e, ".regenie")), sep = " ")
}

# ---- protein pQTL: UK Biobank parquet + deCODE SomaScan TSV -----------------
PR <- paste0("protein_", 1:4); oid <- sprintf("OID9000%d", 1:4)
ukp <- mk(igloo, "UKB", "pQTL"); dcp <- mk(igloo, "DECODE", "pQTL", "final_somascan_smp")
for (k in 1:4) {
  d <- stats(eff[[PR[k]]], N_P_UKB)
  write_parquet(d[, .(CHROM = chr, GENPOS = pos, ALLELE0 = ref, ALLELE1 = alt, A1FREQ = af,
                      N, BETA = beta, SE = se, LOG10P = log10p, rsid)],
                file.path(ukp, sprintf("%s_P0000%d_%s_v1_Demo.parquet", PR[k], k, oid[k])))
  d <- stats(eff[[PR[k]]], N_P_DEC)
  fwrite(d[, .(Chrom = paste0("chr", chr), Pos = pos, Name = paste0("chr", chr, ":", pos),
               rsids = rsid, effectAllele = alt, otherAllele = ref, Beta = beta, SE = se,
               Pval = p, N, ImpMAF = af)],
         file.path(dcp, sprintf("demo_soma_%d_%s.txt", k, PR[k])), sep = "\t")
}
fwrite(data.table(Gene = PR, Olink_ID = oid),
       file.path(mk(igloo, "UKB", "OMICSPRED"), "UKB_Olink_multi_ancestry_models_val_results_portal.csv"))
fwrite(genes[, .(HGNC.symbol = protein, chr, gene_start = center - 5e3, gene_end = center + 5e3)],
       file.path(mk(igloo, "UKB", "pQTLmetadata"), "olink_protein_map_3k_v1.tsv"), sep = "\t")
fwrite(genes[, .(file = sprintf("demo_soma_%d_%s.txt", 1:4, protein), EntrezGeneSymbol = protein,
                 hgnc_symbol = protein, SomaID = sprintf("SL00000%d", 1:4), OlinkID = oid,
                 chromosome_name = chr, TSS = center)],
       file.path(mk(igloo, "DECODE", "pQTLmetadata"), "somascan_protein_map.tsv"), sep = "\t")

# ---- disease GWAS: FinnGen format + manifest --------------------------------
DZ <- paste0("disease_", 1:4); fg <- mk(igloo, "FinnGen", "SummaryStats")
for (dz in DZ) {
  d <- stats(eff[[dz]], N_D)
  fwrite(d[, .(`#chrom` = chr, pos, ref, alt, rsids = rsid, beta, sebeta = se, pval = p, af_alt = af)],
         file.path(fg, paste0(dz, ".gz")), sep = "\t")
}
fwrite(data.table(phenocode = DZ, num_cases = N_CASE, num_controls = N_CTRL),
       file.path(igloo, "FinnGen", "finngen_R12_manifest.tsv"), sep = "\t")

# ---- edge lists, triad inventory, given colocalization ----------------------
ge <- mk(out, "mr_edges", "global_edges"); invisible(mk(out, "mr", "clumps"))
EPd <- rbind(data.table(Exposure = EX[1:4], Protein = PR), data.table(Exposure = EX[5:24], Protein = PR[1]))
fwrite(EPd[, .(Exposure, ExposureCategory = "Demo", Protein)], file.path(ge, "edges_EP.tsv"), sep = "\t")
fwrite(data.table(Exposure = EX[1:4], ExposureCategory = "Demo", Protein = PR), file.path(ge, "edges_PE.tsv"), sep = "\t")
fwrite(data.table(Protein = PR, Disease = DZ), file.path(ge, "edges_PD.tsv"), sep = "\t")
fwrite(data.table(Protein = PR, Disease = DZ), file.path(ge, "edges_DP.tsv"), sep = "\t")
fwrite(data.table(Exposure = EX[1:4], Disease = DZ), file.path(ge, "edges_ED.tsv"), sep = "\t")
fwrite(data.table(Exposure = EX[1:4], Disease = DZ), file.path(ge, "edges_DE.tsv"), sep = "\t")
fwrite(data.table(Exposure = EX[1:4], ExposureCategory = "Demo", Protein = PR, Disease = DZ,
                  Disease_UKB = DZ, ICD10 = "X00"), file.path(ge, "mr_triads.tsv"), sep = "\t")
co <- data.table(protID = PR[-2], target = DZ[-2], PP.H3 = c(0.02, 0.02, 0.65),
                 PP.H4 = c(0.95, 0.95, 0.30), edge_dir = "Pcis_to_D")
fwrite(rbind(cbind(arm = "UKB", co), cbind(arm = "DECODE", co)),
       file.path(mk(out, "support", "coloc"), "coloc_results.tsv"), sep = "\t")

cat(sprintf("Synthetic GWAS: %d SNPs, LD panel n = %d; %d exposures, %d proteins x 2 pQTL platforms, %d diseases\n",
            M, N_REF, length(EX), length(PR), length(DZ)))
