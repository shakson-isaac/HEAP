#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(locuszoomr)
  library(EnsDb.Hsapiens.v86)
  library(scales)
  library(grid)
})

# ============================================================
# LocusZoom-style plots (offline LD via PLINK) with:
# - Custom LD legend bar (small, ~1/5 height)
# - Plot titles passed as args
#
# Args:
#   --harm      <*_harmonized_snps.tsv>
#   --bfile     <LDref/EUR>
#   --plink     <plink binary> (default "plink")
#   --lead      <rsid> (default rs55714927)
#   --window_kb <kb> (default 500)
#   --gene_kb   <kb> (default 500)
#   --title1    <string> (default "Trait 1 (pQTL)")
#   --title2    <string> (default "Trait 2 (GWAS)")
#   --ld_title  <string> (default "LD (r^2)")
#   --out       <prefix>
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
lead_snp  <- get_arg("--lead", "rs55714927")
window_kb <- suppressWarnings(as.numeric(get_arg("--window_kb", "500")))
gene_kb   <- suppressWarnings(as.numeric(get_arg("--gene_kb", "500")))
title1    <- get_arg("--title1", "Trait 1 (pQTL)")
title2    <- get_arg("--title2", "Trait 2 (GWAS)")
ld_title  <- get_arg("--ld_title", "LD (r^2)")
out_pref  <- get_arg("--out", "rs55714927_locuszoomr_ld")

if (is.null(harm_path) || is.null(bfile)) stop("Missing --harm or --bfile")
if (!file.exists(harm_path)) stop("Missing --harm file: ", harm_path)
if (!is.finite(window_kb) || window_kb <= 0) stop("--window_kb must be positive")
if (!is.finite(gene_kb)   || gene_kb   <= 0) stop("--gene_kb must be positive")

dir.create(dirname(out_pref), recursive=TRUE, showWarnings=FALSE)

# ---------------------------
# LocusZoom-like LD bar palette (discrete)
# Labels are bin upper-bounds: 0, 0.2, 0.4, 0.6, 0.8, 1.0
# ---------------------------
LD_LEVELS <- c("0", "0.2", "0.4", "0.6", "0.8", "1.0")

LD_COLORS <- c(
  "0"   = "#3B4CC0",  # deep blue
  "0.2" = "#2C7FB8",  # blue
  "0.4" = "#41B6C4",  # cyan
  "0.6" = "#A6D96A",  # green
  "0.8" = "#FDAE61",  # orange
  "1.0" = "#D7191C"   # red
)

# ---------------------------
# Compact stacked LD legend bar (small)
# ---------------------------
make_ld_legend_small <- function(title = "LD (r^2)") {
  df <- data.frame(
    y = seq_along(LD_LEVELS),
    bin = factor(LD_LEVELS, levels = LD_LEVELS)
  )
  
  ggplot(df, aes(x = 1, y = y, fill = bin)) +
    geom_tile(width = 0.55, height = 0.85, color = "black", linewidth = 0.25) +
    scale_fill_manual(values = LD_COLORS, drop = FALSE) +
    scale_y_continuous(breaks = df$y, labels = LD_LEVELS, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = title) +
    theme_void(base_size = 10) +
    theme(
      plot.title  = element_text(hjust = 0.5, size = 10, margin = margin(b = 3)),
      axis.text.y = element_text(size = 8, margin = margin(l = 4)),
      plot.margin = margin(0, 0, 0, 0)
    ) +
    guides(fill = "none")
}

# ---------------------------
# Load harmonized table
# ---------------------------
m <- fread(harm_path)
need <- c("SNP","chr.1","pos.1","pval.1","pval.2")
miss <- setdiff(need, names(m))
if (length(miss)) stop("harm file missing required columns: ", paste(miss, collapse=", "))

m[, SNP := trimws(as.character(SNP))]
m[, chr := gsub("^chr","", as.character(chr.1))]
m[, pos := as.integer(pos.1)]
m[, p1  := as.numeric(pval.1)]
m[, p2  := as.numeric(pval.2)]

if (!(lead_snp %in% m$SNP)) stop("Lead SNP not present in harmonized file: ", lead_snp)

chr0 <- m[SNP==lead_snp, chr][1]
pos0 <- m[SNP==lead_snp, pos][1]
lo_bp <- max(1L, pos0 - as.integer(window_kb*1000))
hi_bp <- pos0 + as.integer(window_kb*1000)

m <- m[chr == chr0 & pos >= lo_bp & pos <= hi_bp]
m <- m[is.finite(p1) & is.finite(p2)]
m <- m[!duplicated(SNP)]
if (nrow(m) < 200) stop("Too few SNPs in window after filtering: ", nrow(m))

cat("Lead SNP:", lead_snp, " chr", chr0, " pos", pos0, "\n")
cat("Assoc window:", chr0, ":", lo_bp, "-", hi_bp, " (", window_kb, "kb)\n", sep="")
cat("SNPs in assoc window:", nrow(m), "\n")

# ---------------------------
# Intersect SNPs with LDref snps (PLINK --write-snplist)
# ---------------------------
snplist_in <- tempfile(fileext=".txt")
fwrite(m[, .(SNP)], snplist_in, col.names=FALSE)

tmp_pref1 <- tempfile(pattern="plink_snplist_")
cmd1 <- sprintf(
  "%s --bfile %s --extract %s --chr %s --from-bp %d --to-bp %d --write-snplist --out %s",
  shQuote(plink_bin), shQuote(bfile), shQuote(snplist_in),
  shQuote(chr0), as.integer(lo_bp), as.integer(hi_bp), shQuote(tmp_pref1)
)
cat("Running PLINK snplist intersect:\n", cmd1, "\n")
ret <- system(cmd1)
if (ret != 0) stop("PLINK snplist step failed with exit code ", ret)

kept_snps_file <- paste0(tmp_pref1, ".snplist")
if (!file.exists(kept_snps_file)) stop("No .snplist produced by PLINK")
kept <- fread(kept_snps_file, header=FALSE)[[1]]
kept <- kept[!is.na(kept) & kept != ""]
if (length(kept) < 200) stop("Too few SNPs overlap with LD reference: ", length(kept))

setkey(m, SNP)
m2 <- m[J(kept), nomatch=0]
if (!(lead_snp %in% m2$SNP)) stop("Lead SNP not in LDref overlap set: ", lead_snp)

# ---------------------------
# Compute LD matrix (PLINK --r square gz)
# ---------------------------
snplist_kept <- tempfile(fileext=".txt")
fwrite(data.table(SNP=m2$SNP), snplist_kept, col.names=FALSE)

ld_prefix <- tempfile(pattern="plink_ld_")
cmd2 <- sprintf(
  "%s --bfile %s --extract %s --chr %s --from-bp %d --to-bp %d --r square gz --out %s",
  shQuote(plink_bin), shQuote(bfile), shQuote(snplist_kept),
  shQuote(chr0), as.integer(lo_bp), as.integer(hi_bp), shQuote(ld_prefix)
)
cat("Running PLINK LD:\n", cmd2, "\n")
ret <- system(cmd2)
if (ret != 0) stop("PLINK LD step failed with exit code ", ret)

ld_file <- paste0(ld_prefix, ".ld.gz")
if (!file.exists(ld_file)) stop("Missing LD file: ", ld_file)

LD <- as.matrix(fread(cmd=paste("zcat", shQuote(ld_file)), header=FALSE))
storage.mode(LD) <- "double"
if (nrow(LD) != nrow(m2)) stop("LD dim mismatch: ", nrow(LD), " vs ", nrow(m2))

rownames(LD) <- trimws(as.character(m2$SNP))
colnames(LD) <- trimws(as.character(m2$SNP))
m2[, SNP := trimws(as.character(SNP))]

j <- match(lead_snp, colnames(LD))
if (is.na(j)) stop("Lead SNP not found among LD matrix colnames: ", lead_snp)

m2[, r2 := as.numeric(LD[, j])^2]

# Save plot table
plot_table <- m2[, .(snp=SNP, chr=chr, pos=pos, p_trait1=p1, p_trait2=p2, r2=r2)]
fwrite(plot_table, paste0(out_pref, "_plot_table.tsv"), sep="\t")

# ---------------------------
# Build locus objects for gene tracks
# ---------------------------
ensdb <- EnsDb.Hsapiens.v86
seqname_chr <- paste0("chr", chr0)
xrange_bp   <- c(lo_bp, hi_bp)

df1 <- data.frame(snp=m2$SNP, chrom=paste0("chr", m2$chr), pos=m2$pos, p=m2$p1)
df2 <- data.frame(snp=m2$SNP, chrom=paste0("chr", m2$chr), pos=m2$pos, p=m2$p2)

lz1 <- locuszoomr::locus(
  data      = df1,
  ens_db    = ensdb,
  seqname   = seqname_chr,
  xrange    = xrange_bp,
  index_snp = lead_snp,
  chrom     = "chrom",
  pos       = "pos",
  p         = "p",
  LD        = LD
)

lz2 <- locuszoomr::locus(
  data      = df2,
  ens_db    = ensdb,
  seqname   = seqname_chr,
  xrange    = xrange_bp,
  index_snp = lead_snp,
  chrom     = "chrom",
  pos       = "pos",
  p         = "p",
  LD        = LD
)

# ---------------------------
# Association plot with labeled LD bins (0,0.2,...,1.0) + red lead diamond
# ---------------------------
make_assoc_plot <- function(m2, pcol, title) {
  dt <- copy(m2)
  dt[, pos_mb := pos / 1e6]
  dt[, logp := -log10(pmax(get(pcol), 1e-300))]
  dt[, is_lead := (SNP == lead_snp)]
  
  gws <- -log10(5e-8)
  
  dt[, r2_clean := fifelse(!is.finite(r2), NA_real_, pmin(pmax(r2, 0), 1))]
  
  # Map to LocusZoom-style labeled bins by upper bound: 0,0.2,0.4,0.6,0.8,1.0
  dt[, ld_lab := fifelse(r2_clean == 0, "0",
                         fifelse(r2_clean <= 0.2, "0.2",
                                 fifelse(r2_clean <= 0.4, "0.4",
                                         fifelse(r2_clean <= 0.6, "0.6",
                                                 fifelse(r2_clean <= 0.8, "0.8", "1.0")))))]
  dt[, ld_lab := factor(ld_lab, levels = LD_LEVELS)]
  
  # Fade low LD so high LD pops
  alpha_pal <- c("0"=0.15, "0.2"=0.30, "0.4"=0.55, "0.6"=0.70, "0.8"=0.90, "1.0"=1.00)
  dt[, alpha_bin := alpha_pal[as.character(ld_lab)]]
  
  # draw low first, high last
  dt[, ord := as.integer(ld_lab)]
  setorder(dt, ord)
  
  ggplot(dt, aes(x=pos_mb, y=logp)) +
    geom_point(aes(color=ld_lab, alpha=alpha_bin), size=1.3, na.rm=TRUE) +
    scale_alpha_identity() +
    geom_hline(yintercept=gws, linetype=3) +
    
    geom_point(
      data = dt[is_lead == TRUE],
      aes(x=pos_mb, y=logp),
      shape=23, fill=LD_COLORS[["1.0"]], color="black",
      size=4.0, stroke=0.9,
      inherit.aes=FALSE
    ) +
    geom_text(
      data = dt[is_lead == TRUE],
      aes(x=pos_mb, y=logp, label=lead_snp),
      hjust=1.05, vjust=0.5, size=3.2,
      nudge_x = -0.01,
      inherit.aes=FALSE
    ) +
    
    scale_color_manual(values = LD_COLORS, breaks = LD_LEVELS, drop = FALSE) +
    labs(
      title=title,
      x=paste0("Chromosome ", chr0, " (Mb)"),
      y=expression(-log[10](p))
    ) +
    coord_cartesian(clip="off") +
    expand_limits(y = max(dt$logp, na.rm=TRUE) * 1.08) +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(hjust=0.5, size=16),
      legend.position = "none",  # we add custom bar legend separately
      plot.margin = margin(10, 10, 10, 10)
    )
}

p_assoc1 <- make_assoc_plot(m2, "p1", title1)
p_assoc2 <- make_assoc_plot(m2, "p2", title2)

# ---------------------------
# Gene tracks (protein-coding)
# ---------------------------
gene_lo <- max(1L, pos0 - as.integer(gene_kb*1000))
gene_hi <- pos0 + as.integer(gene_kb*1000)

lz1g <- lz1; lz1g$xrange <- c(gene_lo, gene_hi)
lz2g <- lz2; lz2g$xrange <- c(gene_lo, gene_hi)

p_gene1 <- locuszoomr::gg_genetracks(lz1g, filter_gene_biotype = "protein_coding")
p_gene2 <- locuszoomr::gg_genetracks(lz2g, filter_gene_biotype = "protein_coding")

# Stack assoc + genes
p1 <- p_assoc1 / p_gene1 + plot_layout(heights = c(4.2, 2.4))
p2 <- p_assoc2 / p_gene2 + plot_layout(heights = c(4.2, 2.4))

ggsave(paste0(out_pref, "_trait1.png"), p1, width=8, height=7.5, dpi=300)
ggsave(paste0(out_pref, "_trait2.png"), p2, width=8, height=7.5, dpi=300)

# ---------------------------
# Small LD bar legend (~1/5 height) via spacer stacking
# ---------------------------
ld_bar <- make_ld_legend_small(ld_title)

# Make the legend column: spacer / legend / spacer
# Heights ratio: 2 : 1 : 2  => legend is 1/5 of total
ld_col <- plot_spacer() / ld_bar / plot_spacer() + plot_layout(heights = c(2, 1, 2))

p_pair <- (p1 | p2 | ld_col) + plot_layout(widths = c(1, 1, 0.14))
ggsave(paste0(out_pref, "_paired.png"), p_pair, width=17, height=10, dpi=300)

cat("Wrote:\n  ",
    paste0(out_pref, "_plot_table.tsv"), "\n  ",
    paste0(out_pref, "_trait1.png"), "\n  ",
    paste0(out_pref, "_trait2.png"), "\n  ",
    paste0(out_pref, "_paired.png"), "\n", sep="")
cat("DONE.\n")










# #!/usr/bin/env Rscript
# 
# suppressPackageStartupMessages({
#   library(data.table)
#   library(ggplot2)
#   library(patchwork)
#   library(locuszoomr)
#   library(EnsDb.Hsapiens.v86)
#   library(scales)
# })
# 
# # ============================================================
# # LocusZoomR-style plots (offline LD via PLINK) + gene filtering
# # - LD (r^2) computed from PLINK --r square gz
# # - Association plot: binned LD + lead SNP diamond + GWS line
# # - Gene tracks: locuszoomr gg_genetracks() with:
# #     filter_gene_biotype = "protein_coding"
# #     filter_gene_name    = <whitelist>
# #
# # Args:
# #   --harm      <*_harmonized_snps.tsv>
# #   --bfile     <LDref/EUR>
# #   --plink     <plink binary> (default "plink")
# #   --lead      <rsid> (default rs55714927)
# #   --window_kb <kb> (default 500)
# #   --gene_kb   <kb> (default 500)  # recommended for gene display
# #   --out       <prefix>
# # ============================================================
# 
# args <- commandArgs(trailingOnly = TRUE)
# get_arg <- function(flag, default=NULL) {
#   w <- which(args == flag)
#   if (length(w) == 0) return(default)
#   if (w == length(args)) stop("Missing value for ", flag)
#   args[w + 1]
# }
# 
# harm_path <- get_arg("--harm")
# bfile     <- get_arg("--bfile")
# plink_bin <- get_arg("--plink", "plink")
# lead_snp  <- get_arg("--lead", "rs55714927")
# window_kb <- suppressWarnings(as.numeric(get_arg("--window_kb", "500")))
# gene_kb   <- suppressWarnings(as.numeric(get_arg("--gene_kb", "500")))
# out_pref  <- get_arg("--out", "rs55714927_locuszoomr_ld")
# 
# if (is.null(harm_path) || is.null(bfile)) stop("Missing --harm or --bfile")
# if (!file.exists(harm_path)) stop("Missing --harm file: ", harm_path)
# if (!is.finite(window_kb) || window_kb <= 0) stop("--window_kb must be positive")
# if (!is.finite(gene_kb)   || gene_kb   <= 0) stop("--gene_kb must be positive")
# 
# dir.create(dirname(out_pref), recursive=TRUE, showWarnings=FALSE)
# 
# # ---------------------------
# # Gene filter settings (EDIT)
# # ---------------------------
# whitelist <- c("ASGR1", "ASGR2") # , "DLG4", "ACAP1")
# gene_biotype <- "protein_coding"
# 
# # ---------------------------
# # Load harmonized table
# # ---------------------------
# m <- fread(harm_path)
# need <- c("SNP","chr.1","pos.1","pval.1","pval.2")
# miss <- setdiff(need, names(m))
# if (length(miss)) stop("harm file missing required columns: ", paste(miss, collapse=", "))
# 
# m[, SNP := trimws(as.character(SNP))]
# m[, chr := gsub("^chr","", as.character(chr.1))]
# m[, pos := as.integer(pos.1)]
# m[, p1  := as.numeric(pval.1)]
# m[, p2  := as.numeric(pval.2)]
# 
# if (!(lead_snp %in% m$SNP)) stop("Lead SNP not present in harmonized file: ", lead_snp)
# 
# chr0 <- m[SNP==lead_snp, chr][1]
# pos0 <- m[SNP==lead_snp, pos][1]
# lo_bp <- max(1L, pos0 - as.integer(window_kb*1000))
# hi_bp <- pos0 + as.integer(window_kb*1000)
# 
# m <- m[chr == chr0 & pos >= lo_bp & pos <= hi_bp]
# m <- m[is.finite(p1) & is.finite(p2)]
# m <- m[!duplicated(SNP)]
# if (nrow(m) < 200) stop("Too few SNPs in window after filtering: ", nrow(m))
# 
# cat("Lead SNP:", lead_snp, " chr", chr0, " pos", pos0, "\n")
# cat("Assoc window:", chr0, ":", lo_bp, "-", hi_bp, " (", window_kb, "kb)\n", sep="")
# cat("SNPs in assoc window:", nrow(m), "\n")
# 
# # ---------------------------
# # Intersect SNPs with LDref snps (PLINK --write-snplist)
# # ---------------------------
# snplist_in <- tempfile(fileext=".txt")
# fwrite(m[, .(SNP)], snplist_in, col.names=FALSE)
# 
# tmp_pref1 <- tempfile(pattern="plink_snplist_")
# cmd1 <- sprintf(
#   "%s --bfile %s --extract %s --chr %s --from-bp %d --to-bp %d --write-snplist --out %s",
#   shQuote(plink_bin), shQuote(bfile), shQuote(snplist_in),
#   shQuote(chr0), as.integer(lo_bp), as.integer(hi_bp), shQuote(tmp_pref1)
# )
# cat("Running PLINK snplist intersect:\n", cmd1, "\n")
# ret <- system(cmd1)
# if (ret != 0) stop("PLINK snplist step failed with exit code ", ret)
# 
# kept_snps_file <- paste0(tmp_pref1, ".snplist")
# if (!file.exists(kept_snps_file)) stop("No .snplist produced by PLINK")
# kept <- fread(kept_snps_file, header=FALSE)[[1]]
# kept <- kept[!is.na(kept) & kept != ""]
# if (length(kept) < 200) stop("Too few SNPs overlap with LD reference: ", length(kept))
# 
# setkey(m, SNP)
# m2 <- m[J(kept), nomatch=0]
# if (!(lead_snp %in% m2$SNP)) stop("Lead SNP not in LDref overlap set: ", lead_snp)
# 
# cat("SNPs used (post-LDref intersect):", nrow(m2), "\n")
# 
# # ---------------------------
# # Compute LD matrix (PLINK --r square gz)
# # ---------------------------
# snplist_kept <- tempfile(fileext=".txt")
# fwrite(data.table(SNP=m2$SNP), snplist_kept, col.names=FALSE)
# 
# ld_prefix <- tempfile(pattern="plink_ld_")
# cmd2 <- sprintf(
#   "%s --bfile %s --extract %s --chr %s --from-bp %d --to-bp %d --r square gz --out %s",
#   shQuote(plink_bin), shQuote(bfile), shQuote(snplist_kept),
#   shQuote(chr0), as.integer(lo_bp), as.integer(hi_bp), shQuote(ld_prefix)
# )
# cat("Running PLINK LD:\n", cmd2, "\n")
# ret <- system(cmd2)
# if (ret != 0) stop("PLINK LD step failed with exit code ", ret)
# 
# ld_file <- paste0(ld_prefix, ".ld.gz")
# if (!file.exists(ld_file)) stop("Missing LD file: ", ld_file)
# 
# LD <- as.matrix(fread(cmd=paste("zcat", shQuote(ld_file)), header=FALSE))
# storage.mode(LD) <- "double"
# if (nrow(LD) != nrow(m2)) stop("LD dim mismatch: ", nrow(LD), " vs ", nrow(m2))
# 
# rownames(LD) <- trimws(as.character(m2$SNP))
# colnames(LD) <- trimws(as.character(m2$SNP))
# m2[, SNP := trimws(as.character(SNP))]
# 
# j <- match(lead_snp, colnames(LD))
# if (is.na(j)) {
#   stop("Lead SNP not found among LD matrix colnames. lead=", lead_snp,
#        "\nExample LD colnames: ", paste(head(colnames(LD)), collapse=", "))
# }
# m2[, r2 := as.numeric(LD[, j])^2]
# 
# cat("DEBUG r2 finite:", sum(is.finite(m2$r2)), "/", nrow(m2), "\n")
# 
# # Save plot table
# plot_table <- m2[, .(snp=SNP, chr=chr, pos=pos, p_trait1=p1, p_trait2=p2, r2=r2)]
# fwrite(plot_table, paste0(out_pref, "_plot_table.tsv"), sep="\t")
# 
# # ---------------------------
# # Build locus objects (for gene tracks / locuszoomr plotting)
# # IMPORTANT: locuszoomr::locus expects seqname style matching EnsDb.
# # We will use "chr17" for the locus object because your earlier working locus used chr-prefix.
# # ---------------------------
# ensdb <- EnsDb.Hsapiens.v86
# seqname_chr <- paste0("chr", chr0)
# xrange_bp   <- c(lo_bp, hi_bp)
# 
# df1 <- data.frame(snp=m2$SNP, chrom=paste0("chr", m2$chr), pos=m2$pos, p=m2$p1)
# df2 <- data.frame(snp=m2$SNP, chrom=paste0("chr", m2$chr), pos=m2$pos, p=m2$p2)
# 
# lz1 <- locuszoomr::locus(
#   data      = df1,
#   ens_db    = ensdb,
#   seqname   = seqname_chr,
#   xrange    = xrange_bp,
#   index_snp = lead_snp,
#   chrom     = "chrom",
#   pos       = "pos",
#   p         = "p",
#   LD        = LD
# )
# 
# lz2 <- locuszoomr::locus(
#   data      = df2,
#   ens_db    = ensdb,
#   seqname   = seqname_chr,
#   xrange    = xrange_bp,
#   index_snp = lead_snp,
#   chrom     = "chrom",
#   pos       = "pos",
#   p         = "p",
#   LD        = LD
# )
# 
# # ---------------------------
# # Association plot (binned LD) with lead SNP highlight
# # ---------------------------
# make_assoc_plot <- function(m2, pcol, title) {
#   dt <- copy(m2)
#   dt[, pos_mb := pos / 1e6]
#   dt[, logp := -log10(pmax(get(pcol), 1e-300))]
#   dt[, is_lead := (SNP == lead_snp)]
#   
#   gws <- -log10(5e-8)
#   
#   dt[, r2_clean := fifelse(!is.finite(r2), NA_real_, pmin(pmax(r2, 0), 1))]
#   dt[, r2_bin := cut(
#     r2_clean,
#     breaks = c(-Inf, 0, 0.2, 0.4, 0.6, 0.8, Inf),
#     labels = c("0", "0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"),
#     right = TRUE
#   )]
#   dt[, r2_bin := factor(r2_bin,
#                         levels = c("0", "0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"))]
#   
#   # LocusZoom-ish palette (stronger separation)
#   pal <- c(
#     "0"       = "#BEBEBE",
#     "0–0.2"   = "#1F78B4",
#     "0.2–0.4" = "#33A02C",
#     "0.4–0.6" = "#A6D854",
#     "0.6–0.8" = "#FDBF6F",
#     "0.8–1.0" = "#E31A1C"
#   )
#   
#   # Fixed alpha per BIN (this is key)
#   alpha_pal <- c(
#     "0"       = 0.15,
#     "0–0.2"   = 0.30,
#     "0.2–0.4" = 0.55,
#     "0.4–0.6" = 0.70,
#     "0.6–0.8" = 0.90,
#     "0.8–1.0" = 1.00
#   )
#   dt[, alpha_bin := alpha_pal[as.character(r2_bin)]]
#   
#   # Draw order: low LD first, high LD last (on top)
#   dt[, r2_order := as.integer(r2_bin)]
#   setorder(dt, r2_order)  # ensures high LD points drawn last
#   
#   ggplot(dt, aes(x=pos_mb, y=logp)) +
#     geom_point(aes(color=r2_bin, alpha=alpha_bin), size=1.35, na.rm=TRUE) +
#     scale_alpha_identity() +
#     geom_hline(yintercept=gws, linetype=3) +
#     
#     # Lead SNP diamond (red)
#     geom_point(
#       data = dt[is_lead == TRUE],
#       aes(x=pos_mb, y=logp),
#       shape=23, fill="#E31A1C", color="black",
#       size=4.2, stroke=0.9,
#       inherit.aes=FALSE
#     ) +
#     # label to the LEFT
#     geom_text(
#       data = dt[is_lead == TRUE],
#       aes(x=pos_mb, y=logp, label=lead_snp),
#       hjust=1.05, vjust=0.5, size=3.4,
#       nudge_x = -0.01,
#       inherit.aes=FALSE
#     ) +
#     
#     scale_color_manual(
#       values = pal,
#       breaks = names(pal),   # force all legend entries + order
#       drop   = FALSE,
#       name   = expression(LD~(r^2))
#     ) +
#     labs(
#       title=title,
#       x=paste0("Chromosome ", chr0, " (Mb)"),
#       y=expression(-log[10](p))
#     ) +
#     coord_cartesian(clip="off") +
#     expand_limits(y = max(dt$logp, na.rm=TRUE) * 1.08) +
#     theme_classic(base_size = 15) +
#     theme(
#       plot.title = element_text(hjust=0.5, size=16),
#       legend.position = "right",
#       legend.key.height = unit(0.7, "lines"),
#       plot.margin = margin(10, 10, 10, 10)
#     )
# }
# p_assoc1 <- make_assoc_plot(m2, "p1", "Trait 1 (pQTL)")
# p_assoc2 <- make_assoc_plot(m2, "p2", "Trait 2 (GWAS)")
# 
# # ---------------------------
# # Gene tracks using locuszoomr with filters
# # ---------------------------
# gene_lo <- max(1L, pos0 - as.integer(gene_kb*1000))
# gene_hi <- pos0 + as.integer(gene_kb*1000)
# 
# lz1g <- lz1; lz1g$xrange <- c(gene_lo, gene_hi)
# lz2g <- lz2; lz2g$xrange <- c(gene_lo, gene_hi)
# 
# # NOTE: gg_genetracks supports filtering in locuszoomr >= versions that document it.
# # If your version errors on these args, remove them and keep only the xrange trick.
# p_gene1 <- locuszoomr::gg_genetracks(
#   lz1g,
#   filter_gene_biotype = gene_biotype #,
#   #filter_gene_name    = whitelist
# )
# 
# p_gene2 <- locuszoomr::gg_genetracks(
#   lz2g,
#   filter_gene_biotype = gene_biotype #,
#   #filter_gene_name    = whitelist
# )
# 
# # Stack assoc + genes
# p1 <- p_assoc1 / p_gene1 + plot_layout(heights=c(4.0, 6))
# p2 <- p_assoc2 / p_gene2 + plot_layout(heights=c(4.0, 6))
# 
# ggsave(paste0(out_pref, "_trait1.png"), p1, width=8, height=7.5, dpi=300)
# ggsave(paste0(out_pref, "_trait2.png"), p2, width=8, height=7.5, dpi=300)
# 
# p_pair <- (p1 | p2) + plot_layout(guides = "collect") & theme(legend.position="right")
# ggsave(paste0(out_pref, "_paired.png"), p_pair, width=16, height=10, dpi=300)
# 
# cat("Wrote:\n  ",
#     paste0(out_pref, "_plot_table.tsv"), "\n  ",
#     paste0(out_pref, "_trait1.png"), "\n  ",
#     paste0(out_pref, "_trait2.png"), "\n  ",
#     paste0(out_pref, "_paired.png"), "\n", sep="")
# cat("DONE.\n")
