#!/usr/bin/env Rscript
# ============================================================================
# support/coloc/run_coloc_locus.R
# ----------------------------------------------------------------------------
# Per-LOCUS coloc.abf runner for the cis-pQTL colocalization pipeline. Reads a
# protein cis-pQTL signal (deCODE SomaScan .txt.gz OR UKB-PPP parquet) and an
# outcome GWAS (FinnGen disease .gz for P->D, UKB REGENIE .regenie for P->E),
# extracts +/- window_kb around the cis lead SNP, harmonizes alleles
# (strand/palindrome-aware), and runs coloc::coloc.abf.
#
# Allele harmonization + coloc.abf logic are ported VERBATIM from the legacy
# engine (ModuleMR/COLOC/runColoc.R + support/coloc/run_coloc_shortlist.R):
#   - pQTL          -> type="quant", N = scalar (per-locus median)
#   - FinnGen disease -> type="cc",  s = cases/(cases+controls), N = cases+controls
#   - UKB exposure  -> type="quant", N = median(per-SNP N)
#
# Two invocation modes:
#   (A) MANIFEST mode (canonical, array-friendly):
#         Rscript run_coloc_locus.R --manifest <coloc_manifest.tsv> --row <i>
#       Reads row i (the `row` column, 1-based) and resolves all inputs from it.
#   (B) DIRECT mode (runColoc.R-compatible, ad-hoc):
#         Rscript run_coloc_locus.R --pqtl <file> --gwas <file> --snp <leadSNP> \
#           --window_kb 500 --arm UKB --protID ASGR1 --disease finngen_R12_E4_LIPOPROT \
#           [--s <casefrac>] [--n_gwas <N>] [--pqtl_format auto|decode|ukb_parquet] \
#           [--out <prefix>]
#
# Output (per locus, into output/support/coloc/per_locus/ by default):
#   <out>_coloc_summary.tsv   (lead_snp, chr, pos, nsnps, PP.H0..PP.H4, arm, protID, disease, ...)
#   <out>_harmonized_snps.tsv (the harmonized merged SNP table)
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/support/coloc/run_coloc_locus.R --manifest <tsv> --row <i>
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            file.path(getwd(), "workflow", "00_paths.R"),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("Could not locate workflow/00_paths.R (set HEAP_PATHS_FILE).")
  source(hit)
})
suppressPackageStartupMessages({ library(data.table); library(coloc) })
setDTthreads(as.integer(Sys.getenv("CPUS", "1")))

P1 <- 1e-4; P2 <- 1e-4; P12 <- 1e-5
MIN_SNPS <- as.integer(Sys.getenv("COLOC_MIN_SNPS", "50"))
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# ---- args ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  w <- which(args == flag)
  if (!length(w)) return(default)
  if (w[1] == length(args)) stop("Missing value for ", flag)
  args[w[1] + 1]
}
num_arg <- function(flag, default = NA_real_) {
  v <- suppressWarnings(as.numeric(get_arg(flag, NA_character_)))
  if (is.finite(v)) v else default
}

# ---- shared helpers (verbatim from run_coloc_shortlist.R) ------------------
norm_chr <- function(x) gsub("^chr", "", as.character(x))
to_maf <- function(x) { v <- suppressWarnings(as.numeric(x)); ifelse(is.finite(v), pmin(v, 1 - v), NA_real_) }
is_palindromic <- function(a1, a2) { a1 <- toupper(a1); a2 <- toupper(a2)
  (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C") }

DEC_MAP_FP <- "/n/groups/patel/IGLOO/DECODE/pQTLmetadata/somascan_protein_map.tsv"
FG_MAN_FP  <- "/n/groups/patel/IGLOO/FinnGen/finngen_R12_manifest.tsv"

# pQTL readers — return a uniform data.table windowed to +/- window around lead.
read_ukb_pqtl <- function(path, w) {
  if (!requireNamespace("arrow", quietly = TRUE))
    stop("arrow is required to read UKB-PPP parquet pQTL but is not installed.")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("dplyr is required for the parquet column pushdown but is not installed.")
  d <- tryCatch(as.data.table(
        arrow::open_dataset(path) |>
          dplyr::filter(CHROM == as.integer(w$chr), GENPOS >= w$lo, GENPOS <= w$hi) |>
          dplyr::collect()), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  pv <- if ("PVAL" %in% names(d)) d$PVAL else 10^(-as.numeric(d$LOG10P))
  data.table(SNP = d$rsid, chr = norm_chr(d$CHROM), pos = as.integer(d$GENPOS),
             beta = as.numeric(d$BETA), se = as.numeric(d$SE), pval = as.numeric(pv),
             effect_allele = toupper(d$ALLELE1), other_allele = toupper(d$ALLELE0),
             eaf = as.numeric(d$A1FREQ), N = as.numeric(d$N))
}
read_dec_pqtl <- function(path, w) {
  d <- fread(path, showProgress = FALSE)
  d[, chr := norm_chr(Chrom)]; d <- d[chr == w$chr & Pos >= w$lo & Pos <= w$hi]
  if (!nrow(d)) return(NULL)
  snp <- if ("rsids" %in% names(d)) tstrsplit(d$rsids, ",", fixed = TRUE, keep = 1)[[1]] else d$Name
  data.table(SNP = snp, chr = d$chr, pos = as.integer(d$Pos),
             beta = as.numeric(d$Beta), se = as.numeric(d$SE), pval = as.numeric(d$Pval),
             effect_allele = toupper(d$effectAllele), other_allele = toupper(d$otherAllele),
             eaf = as.numeric(d$ImpMAF), N = as.numeric(gsub(",", "", as.character(d$N))))
}
read_pqtl <- function(path, w, fmt) {
  if (fmt == "auto") fmt <- if (grepl("\\.parquet$", path)) "ukb_parquet" else "decode"
  if (fmt == "ukb_parquet") read_ukb_pqtl(path, w) else read_dec_pqtl(path, w)
}

# outcome readers — return list(d=DT, type, s, N)
read_finngen <- function(path, w, s_override = NA_real_, n_override = NA_real_) {
  d <- fread(path, showProgress = FALSE)
  setnames(d, "#chrom", "chrom", skip_absent = TRUE)
  d[, chr := norm_chr(chrom)]; d <- d[chr == w$chr & pos >= w$lo & pos <= w$hi]
  if (!nrow(d)) return(NULL)
  snp <- if ("rsids" %in% names(d)) d$rsids else paste0(d$chr, ":", d$pos, "_", d$ref, "_", d$alt)
  # s / N: prefer explicit overrides, else derive from manifest by phenocode
  s <- s_override; N <- n_override
  if ((!is.finite(s) || !is.finite(N)) && file.exists(FG_MAN_FP)) {
    man <- fread(FG_MAN_FP)
    ph  <- sub("\\.gz$", "", sub("^finngen_R12_", "", basename(path)))
    mr  <- man[phenocode == ph]
    if (nrow(mr)) {
      nca <- as.numeric(mr$num_cases[1]); nco <- as.numeric(mr$num_controls[1])
      if (!is.finite(s)) s <- nca / (nca + nco)
      if (!is.finite(N)) N <- nca + nco
    }
  }
  list(d = data.table(SNP = snp, chr = d$chr, pos = as.integer(d$pos),
         beta = as.numeric(d$beta), se = as.numeric(d$sebeta), pval = as.numeric(d$pval),
         effect_allele = toupper(d$alt), other_allele = toupper(d$ref), eaf = as.numeric(d$af_alt)),
       type = "cc", s = s, N = N)
}
read_exposure <- function(path, w) {
  d <- fread(path, showProgress = FALSE)
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
mk_ds <- function(m, k, type, s = NULL, N = NULL) {
  out <- list(snp = m$SNP, beta = as.numeric(m[[paste0("beta.", k)]]),
              varbeta = as.numeric(m[[paste0("se.", k)]])^2,
              MAF = to_maf(m[[paste0("eaf.", k)]]), type = type, N = N)
  if (type == "cc") out$s <- s
  out
}

# ---- resolve the locus parameters from manifest row or direct args ---------
manifest_fp <- get_arg("--manifest")
row_i       <- suppressWarnings(as.integer(get_arg("--row", NA_character_)))
out_dir_def <- heap_project_output("support", "coloc", "per_locus")

if (!is.null(manifest_fp)) {
  if (!file.exists(manifest_fp)) stop("manifest not found: ", manifest_fp)
  if (is.na(row_i)) stop("--manifest requires --row <i> (1-based `row` column).")
  man <- fread(manifest_fp)
  mr  <- if ("row" %in% names(man)) man[row == row_i] else man[row_i]
  if (!nrow(mr)) stop("row ", row_i, " not found in manifest.")
  mr <- mr[1]
  arm        <- mr$arm
  protID     <- mr$protID
  disease    <- mr$disease_or_exposure_id
  edge_dir   <- if ("edge_dir" %in% names(mr)) mr$edge_dir else NA_character_
  pqtl_path  <- mr$pqtl_path
  gwas_path  <- mr$gwas_path
  lead_snp   <- mr$lead_snp
  chr0       <- norm_chr(mr$chr)
  pos0       <- as.integer(mr$pos)
  window_kb  <- as.numeric(mr$window_kb)
  out_type   <- if ("outcome_type" %in% names(mr)) mr$outcome_type else if (startsWith(disease, "finngen_R12_")) "cc" else "quant"
  s_casefrac <- if ("s_casefrac" %in% names(mr)) as.numeric(mr$s_casefrac) else NA_real_
  n_gwas     <- if ("n_gwas" %in% names(mr)) as.numeric(mr$n_gwas) else NA_real_
  pqtl_fmt   <- "auto"
} else {
  pqtl_path  <- get_arg("--pqtl");  gwas_path <- get_arg("--gwas")
  lead_snp   <- get_arg("--snp")
  arm        <- get_arg("--arm", if (grepl("\\.parquet$", pqtl_path %||% "")) "UKB" else "DECODE")
  protID     <- get_arg("--protID", NA_character_)
  disease    <- get_arg("--disease", NA_character_)
  edge_dir   <- get_arg("--edge_dir", NA_character_)
  window_kb  <- num_arg("--window_kb", 500)
  s_casefrac <- num_arg("--s", NA_real_)
  n_gwas     <- num_arg("--n_gwas", NA_real_)
  pqtl_fmt   <- get_arg("--pqtl_format", "auto")
  out_type   <- if (!is.null(gwas_path) && grepl("finngen", basename(gwas_path %||% ""))) "cc" else
                if (is.finite(s_casefrac)) "cc" else "quant"
  chr0 <- NA_character_; pos0 <- NA_integer_   # derived from lead SNP below
  if (is.null(pqtl_path) || is.null(gwas_path) || is.null(lead_snp))
    stop("DIRECT mode requires --pqtl, --gwas, --snp (or use --manifest --row).")
}

out_prefix <- get_arg("--out", file.path(out_dir_def, paste(arm, protID, disease, sep = "__")))
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(pqtl_path)) stop("pQTL file not found: ", pqtl_path)
if (!file.exists(gwas_path)) stop("GWAS file not found: ", gwas_path)
if (!is.finite(window_kb) || window_kb <= 0) stop("window_kb must be positive.")

# If chr/pos not provided (direct mode), peek into pQTL/GWAS for the lead SNP.
if (is.na(chr0) || is.na(pos0)) {
  if (grepl("\\.parquet$", pqtl_path) && requireNamespace("arrow", quietly = TRUE) &&
      requireNamespace("dplyr", quietly = TRUE)) {
    lr <- tryCatch(as.data.table(arrow::open_dataset(pqtl_path) |>
            dplyr::filter(rsid == lead_snp) |>
            dplyr::select(CHROM, GENPOS, rsid) |> dplyr::collect()), error = function(e) NULL)
    if (!is.null(lr) && nrow(lr)) { chr0 <- norm_chr(lr$CHROM[1]); pos0 <- as.integer(lr$GENPOS[1]) }
  } else if (!grepl("\\.parquet$", pqtl_path)) {
    d <- fread(pqtl_path, showProgress = FALSE)
    if ("rsids" %in% names(d)) d[, .snp := tstrsplit(rsids, ",", fixed = TRUE, keep = 1)[[1]]]
    lr <- if (".snp" %in% names(d)) d[.snp == lead_snp] else d[get("Name") == lead_snp]
    if (nrow(lr)) { chr0 <- norm_chr(lr$Chrom[1]); pos0 <- as.integer(lr$Pos[1]) }
  }
  if (is.na(chr0) || is.na(pos0))
    stop("Could not resolve chr/pos for lead SNP ", lead_snp, " from pQTL; pass --manifest with a resolved row.")
}

lo <- max(1L, as.integer(pos0 - window_kb * 1000))
hi <- as.integer(pos0 + window_kb * 1000)
w  <- list(chr = chr0, pos = pos0, lo = lo, hi = hi, lead = lead_snp)

cat(sprintf("Locus: %s :: %s -> %s  (%s)\n", arm, protID, disease, edge_dir))
cat(sprintf("Lead SNP: %s  chr %s pos %d   window %g kb => %s:%d-%d\n",
            lead_snp, chr0, pos0, window_kb, chr0, lo, hi))

# ---- read + harmonize + coloc.abf ------------------------------------------
write_fail <- function(note, nsnps = 0L) {
  fwrite(data.table(lead_snp = lead_snp, chr = chr0, pos = pos0, window_kb = window_kb,
                    nsnps = nsnps, arm = arm, protID = protID, disease = disease,
                    edge_dir = edge_dir, outcome_type = out_type,
                    s_casefrac = s_casefrac, status = note,
                    PP.H0 = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_,
                    PP.H3 = NA_real_, PP.H4 = NA_real_),
         paste0(out_prefix, "_coloc_summary.tsv"), sep = "\t")
}

pq <- read_pqtl(pqtl_path, w, pqtl_fmt)
if (is.null(pq) || !nrow(pq)) { write_fail("no_pqtl_in_region"); stop("No pQTL SNPs in region.") }

is_dis <- (out_type == "cc")
oc <- if (is_dis) read_finngen(gwas_path, w, s_override = s_casefrac, n_override = n_gwas) else read_exposure(gwas_path, w)
if (is.null(oc) || !nrow(oc$d)) { write_fail("no_outcome_in_region"); stop("No outcome SNPs in region.") }

cat(sprintf("pQTL SNPs in region: %d   outcome SNPs in region: %d\n", nrow(pq), nrow(oc$d)))

m <- harmonize(pq, oc$d)
if (is.null(m) || nrow(m) < MIN_SNPS) {
  n <- if (is.null(m)) 0L else nrow(m)
  if (!is.null(m)) fwrite(m, paste0(out_prefix, "_harmonized_snps.tsv"), sep = "\t")
  write_fail("too_few_snps", n)
  stop("Too few overlapping SNPs after harmonization: ", n)
}

N1 <- stats::median(pq$N, na.rm = TRUE)
if (!is.finite(N1)) { write_fail("missing_pqtl_N", nrow(m)); stop("Could not infer pQTL N.") }

if (is_dis) {
  if (!is.finite(oc$s) || oc$s <= 0 || oc$s >= 1) { write_fail("bad_case_fraction", nrow(m)); stop("Case fraction s must be in (0,1).") }
  if (!is.finite(oc$N)) { write_fail("missing_gwas_N", nrow(m)); stop("FinnGen N missing; pass --n_gwas.") }
}

d1 <- mk_ds(m, "1", "quant", N = N1)
d2 <- mk_ds(m, "2", oc$type, s = oc$s, N = oc$N)

cat(sprintf("coloc.abf: nsnps=%d  d1(quant) N=%g   d2(%s) s=%s N=%g\n",
            nrow(m), N1, oc$type, ifelse(is.finite(oc$s), round(oc$s, 4), "NA"), oc$N))

res <- coloc::coloc.abf(dataset1 = d1, dataset2 = d2, p1 = P1, p2 = P2, p12 = P12)
S <- res$summary

out <- data.table(
  arm = arm, protID = protID, disease = disease, edge_dir = edge_dir,
  lead_snp = lead_snp, chr = chr0, pos = pos0, lo = lo, hi = hi, window_kb = window_kb,
  outcome_type = oc$type, s_casefrac = oc$s, N_pqtl = N1, N_outcome = oc$N,
  nsnps = nrow(m), p1 = P1, p2 = P2, p12 = P12,
  PP.H0 = unname(S["PP.H0.abf"]), PP.H1 = unname(S["PP.H1.abf"]),
  PP.H2 = unname(S["PP.H2.abf"]), PP.H3 = unname(S["PP.H3.abf"]),
  PP.H4 = unname(S["PP.H4.abf"]),
  pqtl_file = pqtl_path, gwas_file = gwas_path, status = "ok"
)

fwrite(out, paste0(out_prefix, "_coloc_summary.tsv"), sep = "\t")
fwrite(m,   paste0(out_prefix, "_harmonized_snps.tsv"), sep = "\t")

cat("\nCOLOC DONE\n")
print(out[, .(arm, protID, disease, nsnps, PP.H3 = round(PP.H3, 4), PP.H4 = round(PP.H4, 4))])
cat("\nWrote:\n  ", paste0(out_prefix, "_coloc_summary.tsv"),
    "\n  ", paste0(out_prefix, "_harmonized_snps.tsv"), "\n", sep = "")
