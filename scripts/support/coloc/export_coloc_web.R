#!/usr/bin/env Rscript
# ============================================================================
# export_coloc_web.R
# ----------------------------------------------------------------------------
# Headless per-locus export for the website's regional colocalization view.
#
# For every locus with a retained *_harmonized_snps.tsv, emit the two tables an
# interactive locuszoom needs and nothing else:
#
#   <locus>_plot_table.tsv   snp, chr, pos, p_trait1, p_trait2, r2
#   <locus>_genes.tsv        gene, start, end, strand
#
# The LD step is the reason this cannot be done in the browser: r2 to the lead
# variant comes from PLINK against the 1000G EUR reference panel, the same way
# LocusZoom.R (ModuleMR/COLOC/LocusZoom.R) computes it for the print figure. The
# logic here is lifted from that script deliberately -- same panel, same window,
# same lead -- so the site and the figure colour identically. What is dropped is
# the ggplot rendering, which for 70+ loci is the expensive part and which the
# site does not use.
#
# Genes come from EnsDb.Hsapiens.v86, matching the figure's gene track.
#
# Run:
#   module load gcc/14.2.0 R/4.4.2 plink/1.90b7.7_20241022
#   Rscript scripts/support/coloc/export_coloc_web.R
#   Rscript scripts/support/coloc/export_coloc_web.R --only ASGR1   # one protein
# ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(EnsDb.Hsapiens.v86)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  w <- which(args == flag)
  if (!length(w)) return(default)
  args[w + 1]
}

OUT   <- Sys.getenv("HEAP_OUTPUT", "/n/groups/patel/IGLOO/UKB/HEAP/output")
COLOC <- file.path(OUT, "support", "coloc")
DEST  <- file.path(COLOC, "web")
BFILE <- get_arg("--bfile", "/n/groups/patel/IGLOO/LDref/EUR")
PLINK <- get_arg("--plink", "plink")
ONLY  <- get_arg("--only", NULL)
WINDOW_KB <- as.numeric(get_arg("--window_kb", "500"))
MIN_SNPS  <- 200

dir.create(DEST, recursive = TRUE, showWarnings = FALSE)

# Harmonized tables live either flat in support/coloc (legacy naming,
# <lead>_<protein>_<disease>) or under per_locus/ (<arm>__<protein>__<disease>,
# written by the current runner). Take both.
harm <- c(
  list.files(COLOC, pattern = "_harmonized_snps\\.tsv$", full.names = TRUE),
  list.files(file.path(COLOC, "per_locus"), pattern = "_harmonized_snps\\.tsv$",
             full.names = TRUE)
)
if (!length(harm)) stop("no *_harmonized_snps.tsv found under ", COLOC)
if (!is.null(ONLY)) harm <- harm[grepl(ONLY, basename(harm), fixed = TRUE)]
cat(sprintf("Loci with a harmonized table: %d\n", length(harm)))

ensdb <- EnsDb.Hsapiens.v86

lead_from_summary <- function(stem, dir) {
  f <- file.path(dir, paste0(stem, "_coloc_summary.tsv"))
  if (!file.exists(f)) return(NA_character_)
  s <- fread(f, nrows = 1)
  if ("lead" %in% names(s)) as.character(s$lead[1]) else NA_character_
}

ok <- 0L; skipped <- character()

for (hf in harm) {
  stem <- sub("_harmonized_snps\\.tsv$", "", basename(hf))
  dir0 <- dirname(hf)
  res <- tryCatch({
    m <- fread(hf, showProgress = FALSE)
    setnames(m, "chr.1", "chr", skip_absent = TRUE)
    setnames(m, "pos.1", "pos", skip_absent = TRUE)
    m[, p1 := as.numeric(pval.1)]
    m[, p2 := as.numeric(pval.2)]

    lead <- lead_from_summary(stem, dir0)
    if (is.na(lead) || !nzchar(lead)) {
      # Fall back to the strongest pQTL variant, which is what the coloc
      # summary's `lead` column records anyway.
      lead <- m[which.min(p1), SNP][1]
    }
    if (!(lead %in% m$SNP)) stop("lead ", lead, " not in harmonized table")

    chr0 <- m[SNP == lead, chr][1]
    pos0 <- m[SNP == lead, pos][1]
    lo <- max(1L, pos0 - as.integer(WINDOW_KB * 1000))
    hi <- pos0 + as.integer(WINDOW_KB * 1000)

    m <- m[chr == chr0 & pos >= lo & pos <= hi]
    m <- m[is.finite(p1) & is.finite(p2)][!duplicated(SNP)]
    if (nrow(m) < MIN_SNPS) stop("only ", nrow(m), " SNPs in window")

    # ---- intersect with the LD panel -------------------------------------
    f_in <- tempfile(fileext = ".txt")
    fwrite(m[, .(SNP)], f_in, col.names = FALSE)
    p1pref <- tempfile(pattern = "snplist_")
    # Filter by variant ID ONLY -- no --chr/--from-bp/--to-bp.
    #
    # The LD panel is GRCh37 and the coloc data is GRCh38: rs55714927 is at
    # 17:7,080,316 in the panel and 17:7,176,997 here, and rs4550150 differs by
    # 1.76Mb. Passing GRCh38 bounds to a GRCh37 panel silently shifts the window
    # -- where the offset is small it still overlaps and looks like it worked,
    # where it is large (GFRA1) plink returns "No variants remaining".
    #
    # The extract list is already exactly the variants inside the window, so the
    # positional filter was redundant as well as wrong. Matching on rsID makes
    # the build irrelevant. NOTE: ModuleMR/COLOC/LocusZoom.R has the same
    # filters and therefore the same shifted window.
    cmd <- sprintf("%s --bfile %s --extract %s --write-snplist --out %s",
                   shQuote(PLINK), shQuote(BFILE), shQuote(f_in), shQuote(p1pref))
    if (system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) != 0)
      stop("plink snplist failed")
    kept_f <- paste0(p1pref, ".snplist")
    if (!file.exists(kept_f)) stop("no snplist produced")
    kept <- fread(kept_f, header = FALSE)[[1]]
    setkey(m, SNP)
    m2 <- m[J(kept), nomatch = 0]
    if (nrow(m2) < MIN_SNPS) stop("only ", nrow(m2), " SNPs overlap the LD panel")

    # The lead variant is not always in the 1000G EUR panel -- rs5112 (APOC1)
    # and rs17632542 (FURIN) are simply absent from it. Rather than drop a third
    # of the loci, anchor LD on the strongest pQTL variant that IS in the panel
    # and record that we did: r2 is then to a PROXY, and the site says so. An
    # unlabelled proxy would quietly misstate what the colours mean.
    anchor <- lead
    anchor_is_lead <- TRUE
    if (!(lead %in% m2$SNP)) {
      anchor <- m2[which.min(p1), SNP][1]
      anchor_is_lead <- FALSE
      if (is.na(anchor)) stop("no usable LD anchor in the panel")
      cat(sprintf("      lead %s absent from the LD panel; anchoring on %s\n",
                  lead, anchor))
    }

    # ---- LD matrix -> r2 against the lead --------------------------------
    f_keep <- tempfile(fileext = ".txt")
    fwrite(data.table(SNP = m2$SNP), f_keep, col.names = FALSE)
    ldpref <- tempfile(pattern = "ld_")
    cmd <- sprintf("%s --bfile %s --extract %s --r square gz --out %s",
                   shQuote(PLINK), shQuote(BFILE), shQuote(f_keep), shQuote(ldpref))
    if (system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE) != 0)
      stop("plink LD failed")
    ldf <- paste0(ldpref, ".ld.gz")
    if (!file.exists(ldf)) stop("no LD matrix produced")
    LD <- as.matrix(fread(cmd = paste("zcat", shQuote(ldf)), header = FALSE))
    storage.mode(LD) <- "double"
    if (nrow(LD) != nrow(m2)) stop("LD dim mismatch ", nrow(LD), " vs ", nrow(m2))
    m2[, SNP := trimws(as.character(SNP))]
    colnames(LD) <- m2$SNP
    j <- match(anchor, colnames(LD))
    if (is.na(j)) stop("LD anchor not in the LD matrix")
    m2[, r2 := as.numeric(LD[, j])^2]

    pt <- m2[, .(snp = SNP, chr = chr, pos = pos, p_trait1 = p1, p_trait2 = p2,
                 r2 = round(r2, 4))]
    fwrite(pt, file.path(DEST, paste0(stem, "_plot_table.tsv")), sep = "\t")

    fwrite(data.table(locus = stem, lead = lead, anchor = anchor,
                      anchor_is_lead = anchor_is_lead, n_variants = nrow(pt),
                      chr = chr0, lo = lo, hi = hi),
           file.path(DEST, paste0(stem, "_meta.tsv")), sep = "\t")

    # ---- gene track ------------------------------------------------------
    gr <- ensembldb::genes(
      ensdb,
      filter = AnnotationFilter::AnnotationFilterList(
        AnnotationFilter::SeqNameFilter(as.character(chr0)),
        AnnotationFilter::GeneStartFilter(hi, condition = "<="),
        AnnotationFilter::GeneEndFilter(lo, condition = ">=")))
    gd <- as.data.table(as.data.frame(gr))
    if (nrow(gd)) {
      gd <- gd[gene_biotype == "protein_coding" & nzchar(gene_name)]
      gd <- gd[, .(gene = gene_name, start, end,
                   strand = ifelse(strand == "-", "-", "+"))]
      gd <- gd[!duplicated(gene)][order(start)]
    } else {
      gd <- data.table(gene = character(), start = integer(),
                       end = integer(), strand = character())
    }
    fwrite(gd, file.path(DEST, paste0(stem, "_genes.tsv")), sep = "\t")

    cat(sprintf("  %-52s %5d SNPs  %3d genes  anchor %s%s\n",
                stem, nrow(pt), nrow(gd), anchor,
                if (anchor_is_lead) "" else sprintf(" (proxy for %s)", lead)))
    TRUE
  }, error = function(e) {
    cat(sprintf("  %-52s SKIP: %s\n", stem, conditionMessage(e)))
    FALSE
  })
  if (isTRUE(res)) ok <- ok + 1L else skipped <- c(skipped, stem)
}

cat(sprintf("\nExported %d locus/loci to %s\n", ok, DEST))
if (length(skipped)) cat(sprintf("Skipped %d: %s\n", length(skipped),
                                 paste(head(skipped, 5), collapse = ", ")))
