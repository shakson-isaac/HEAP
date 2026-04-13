#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(locuszoomr)
  library(EnsDb.Hsapiens.v86)
  library(scales)
})

# ---------------------------
# Args
# ---------------------------
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
out_pref  <- get_arg("--out", "rs55714927_locuszoomr_ld")

if (is.null(harm_path) || is.null(bfile)) stop("Missing --harm or --bfile")
if (!file.exists(harm_path)) stop("Missing --harm file: ", harm_path)
if (!is.finite(window_kb) || window_kb <= 0) stop("--window_kb must be positive")

dir.create(dirname(out_pref), recursive=TRUE, showWarnings=FALSE)

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
cat("Window:", chr0, ":", lo_bp, "-", hi_bp, " (", window_kb, "kb)\n", sep="")
cat("SNPs in window (harm):", nrow(m), "\n")

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

cat("SNPs used (post-LDref intersect):", nrow(m2), "\n")

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

rownames(LD) <- m2$SNP
colnames(LD) <- m2$SNP

# r^2 to lead SNP (PLINK gives r)
#m2[, r2 := as.numeric((LD[, lead_snp])^2)[SNP]]
# ---- robust r^2 to lead SNP (positional) ----
m2[, SNP := trimws(as.character(SNP))]
colnames(LD) <- trimws(as.character(colnames(LD)))
rownames(LD) <- trimws(as.character(rownames(LD)))

j <- match(lead_snp, colnames(LD))
if (is.na(j)) {
  stop("Lead SNP not found among LD matrix colnames. lead=", lead_snp,
       "\nExample LD colnames: ", paste(head(colnames(LD)), collapse=", "))
}

r2_vec <- as.numeric(LD[, j])^2
m2[, r2 := r2_vec]   # SAME order as LD

cat("DEBUG lead col index:", j, "\n")
cat("DEBUG r2 summary:\n")
print(summary(m2$r2))
cat("DEBUG r2 finite:", sum(is.finite(m2$r2)), "/", nrow(m2), "\n")

# Save plot table
plot_table <- m2[, .(snp=SNP, chr=chr, pos=pos, p_trait1=p1, p_trait2=p2, r2=r2)]
fwrite(plot_table, paste0(out_pref, "_plot_table.tsv"), sep="\t")

# ---------------------------
# Build locus objects for gene tracks only
# ---------------------------
ensdb <- EnsDb.Hsapiens.v86
seqname_chr <- paste0("chr", chr0)
xrange_bp <- c(lo_bp, hi_bp)

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
# Build “LocusZoom-style” association plots with LD colouring (offline)
# ---------------------------
# make_assoc_plot <- function(m2, pcol, title) {
#   dt <- copy(m2)
#   dt[, pos_mb := pos / 1e6]
#   dt[, logp := -log10(pmax(get(pcol), 1e-300))]
#   lead_x <- dt[SNP == lead_snp, pos_mb][1]
#   
#   ggplot(dt, aes(x=pos_mb, y=logp, color=r2)) +
#     geom_point(size=1.6, alpha=0.95) +
#     geom_vline(xintercept=lead_x, linetype=2) +
#     labs(title=title,
#          x=paste0("Chromosome ", chr0, " (Mb)"),
#          y=expression(-log[10](p)),
#          color=expression(LD~r^2)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust=0.5),
#           legend.position = "left")
# }

# make_assoc_plot <- function(m2, pcol, title) {
#   dt <- copy(m2)
#   dt[, pos_mb := pos / 1e6]
#   dt[, logp := -log10(pmax(get(pcol), 1e-300))]
#   lead_x <- dt[SNP == lead_snp, pos_mb][1]
#   
#   ggplot(dt, aes(x=pos_mb, y=logp)) +
#     geom_point(aes(color=r2), size=1.6, alpha=0.95, na.rm=TRUE) +
#     geom_vline(xintercept=lead_x, linetype=2) +
#     scale_color_gradientn(
#       colours = c("grey80","skyblue3","turquoise3","gold","red3"),
#       limits = c(0,1),
#       oob = scales::squish,
#       name = expression(LD~r^2)
#     ) +
#     labs(title=title,
#          x=paste0("Chromosome ", chr0, " (Mb)"),
#          y=expression(-log[10](p))) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust=0.5),
#           legend.position = "left")
# }

# make_assoc_plot <- function(m2, pcol, title) {
#   dt <- copy(m2)
#   dt[, pos_mb := pos / 1e6]
#   dt[, logp := -log10(pmax(get(pcol), 1e-300))]
#   dt[, is_lead := (SNP == lead_snp)]
#   lead_x <- dt[is_lead == TRUE, pos_mb][1]
#   
#   # Genome-wide significance threshold line: p=5e-8
#   gws <- -log10(5e-8)
#   
#   ggplot(dt, aes(x=pos_mb, y=logp)) +
#     # Background points (low emphasis)
#     geom_point(aes(color=r2),
#                size=1.4,
#                alpha=0.55,
#                na.rm=TRUE) +
#     # Lead SNP as diamond on top
#     geom_point(data = dt[is_lead == TRUE],
#                aes(x=pos_mb, y=logp),
#                shape=23,  # filled diamond
#                fill="white",
#                color="black",
#                size=3.6,
#                stroke=1.0,
#                inherit.aes = FALSE) +
#     # Label lead SNP
#     geom_text(data = dt[is_lead == TRUE],
#               aes(x=pos_mb, y=logp, label=lead_snp),
#               vjust=-0.8, hjust=0.5, size=3.5,
#               inherit.aes = FALSE) +
#     # Vertical line at lead position
#     geom_vline(xintercept=lead_x, linetype=2) +
#     # Genome-wide significance line
#     geom_hline(yintercept=gws, linetype=3) +
#     # LD color scale (you already used this)
#     scale_color_gradientn(
#       colours = c("grey85","skyblue3","turquoise3","gold","red3"),
#       limits = c(0,1),
#       oob = scales::squish,
#       name = expression(LD~r^2)
#     ) +
#     labs(title=title,
#          x=paste0("Chromosome ", chr0, " (Mb)"),
#          y=expression(-log[10](p))) +
#     theme_bw() +
#     theme(
#       plot.title = element_text(hjust=0.5),
#       legend.position = "left"
#     )
# }

make_assoc_plot <- function(m2, pcol, title) {
  dt <- copy(m2)
  dt[, pos_mb := pos / 1e6]
  dt[, logp := -log10(pmax(get(pcol), 1e-300))]
  dt[, is_lead := (SNP == lead_snp)]
  lead_x <- dt[is_lead == TRUE, pos_mb][1]
  
  # Genome-wide significance: p = 5e-8
  gws <- -log10(5e-8)
  
  # Bin LD r^2 (classic LocusZoom)
  # include 0 as its own bin so "unlinked" points stand out
  dt[, r2_clean := fifelse(!is.finite(r2), NA_real_, pmin(pmax(r2, 0), 1))]
  dt[, r2_bin := cut(
    r2_clean,
    breaks = c(-Inf, 0, 0.2, 0.4, 0.6, 0.8, Inf),
    labels = c("0", "0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"),
    right = TRUE
  )]
  
  # Explicit ordering (so legend looks right)
  dt[, r2_bin := factor(r2_bin, levels = c("0", "0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0"))]
  
  ggplot(dt, aes(x=pos_mb, y=logp)) +
    geom_point(aes(color=r2_bin), size=1.4, alpha=0.9, na.rm=TRUE) +
    # lead SNP diamond + label
    geom_point(data = dt[is_lead == TRUE],
               aes(x=pos_mb, y=logp),
               shape=23, fill="white", color="black",
               size=3.6, stroke=1.0, inherit.aes=FALSE) +
    geom_text(data = dt[is_lead == TRUE],
              aes(x=pos_mb, y=logp, label=lead_snp),
              vjust=-0.8, hjust=0.5, size=3.5,
              inherit.aes=FALSE) +
    #geom_vline(xintercept=lead_x, linetype=2) +
    geom_hline(yintercept=gws, linetype=3) +
    scale_color_manual(values = c(
      "0"       = "grey80",
      "0–0.2"   = "steelblue3",
      "0.2–0.4" = "deepskyblue3",
      "0.4–0.6" = "turquoise3",
      "0.6–0.8" = "gold",
      "0.8–1.0" = "red3"
    )) +
    labs(
      title=title,
      x=paste0("Chromosome ", chr0, " (Mb)"),
      y=expression(-log[10](p)),
      color=expression(LD~(r^2))
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust=0.5),
      legend.position = "right"
    )
}

p_assoc1 <- make_assoc_plot(m2, "p1", "Trait 1 (pQTL)")
p_assoc2 <- make_assoc_plot(m2, "p2", "Trait 2 (GWAS)")

# Gene tracks from locuszoomr (offline via EnsDb)
# gg_genetracks returns a ggplot object
p_gene1 <- locuszoomr::gg_genetracks(lz1)
p_gene2 <- locuszoomr::gg_genetracks(lz2)

# smaller gene window (e.g., +/- 200kb) for cleaner labels
# gene_lo <- max(1L, pos0 - 200000L)
# gene_hi <- pos0 + 200000L
# lz1_genes <- lz1; lz1_genes$xrange <- c(gene_lo, gene_hi)
# lz2_genes <- lz2; lz2_genes$xrange <- c(gene_lo, gene_hi)
# 
# p_gene1 <- locuszoomr::gg_genetracks(lz1_genes)
# p_gene2 <- locuszoomr::gg_genetracks(lz2_genes)

# Stack assoc + genes
p1 <- p_assoc1 / p_gene1 + plot_layout(heights=c(3, 3))
p2 <- p_assoc2 / p_gene2 + plot_layout(heights=c(3, 3))

ggsave(paste0(out_pref, "_trait1.png"), p1, width=8, height=7, dpi=300)
ggsave(paste0(out_pref, "_trait2.png"), p2, width=8, height=7, dpi=300)

#p_pair <- p1 | p2
p_pair <- (p1 | p2) + plot_layout(guides = "collect")
ggsave(paste0(out_pref, "_paired.png"), p_pair, width=16, height=7, dpi=300)

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
# })
# 
# # ============================================================
# # locuszoomr LocusZoom-style plots for two traits with LD
# #
# # Inputs:
# #   --harm      harmonized SNP table from your coloc run (*_harmonized_snps.tsv)
# #               must contain: SNP, chr.1, pos.1, pval.1, pval.2
# #   --bfile     PLINK LD reference prefix (same used for MR clumping)
# #   --plink     plink binary path
# #   --lead      lead SNP rsID (default rs55714927)
# #   --window_kb window around lead (default 500)
# #   --out       output prefix (path ok)
# #
# # Outputs:
# #   <out>_trait1.png
# #   <out>_trait2.png
# #   <out>_paired.png
# #   <out>_plot_table.tsv
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
# out_pref  <- get_arg("--out", "rs55714927_ASGR1_E4_LIPOPROT_locuszoomr")
# 
# if (is.null(harm_path) || is.null(bfile)) {
#   stop(
#     "Usage:\n",
#     "  Rscript LocusZoomR_ASGR1_E4_LIPOPROT.R \\\n",
#     "    --harm /path/rs55714927_ASGR1_E4_LIPOPROT_harmonized_snps.tsv \\\n",
#     "    --bfile /n/groups/.../LDref/EUR \\\n",
#     "    --plink /path/to/plink \\\n",
#     "    --lead rs55714927 \\\n",
#     "    --window_kb 500 \\\n",
#     "    --out /path/rs55714927_ASGR1_E4_LIPOPROT\n"
#   )
# }
# if (!file.exists(harm_path)) stop("Missing --harm file: ", harm_path)
# if (!is.finite(window_kb) || window_kb <= 0) stop("--window_kb must be positive")
# 
# dir.create(dirname(out_pref), recursive=TRUE, showWarnings=FALSE)
# 
# # ---------------------------
# # Load harmonized table
# # ---------------------------
# m <- fread(harm_path)
# 
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
# 
# cat("Lead SNP:", lead_snp, "chr", chr0, "pos", pos0, "\n")
# cat("Window:", window_kb, "kb =>", chr0, ":", lo_bp, "-", hi_bp, "\n")
# cat("SNPs in window (harm):", nrow(m), "\n")
# if (nrow(m) < 200) stop("Too few SNPs in window after filtering: ", nrow(m))
# 
# # ---------------------------
# # Intersect SNPs with LD reference (PLINK --write-snplist)
# # ---------------------------
# snplist_in <- tempfile(pattern="lz_snps_in_", fileext=".txt")
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
# if (!file.exists(kept_snps_file)) stop("PLINK did not create .snplist: ", kept_snps_file)
# kept <- fread(kept_snps_file, header=FALSE)[[1]]
# kept <- kept[!is.na(kept) & kept != ""]
# cat("SNPs present in LDref:", length(kept), "\n")
# if (length(kept) < 200) stop("Too few SNPs overlap with LD reference: ", length(kept))
# 
# # Reorder to match PLINK's snplist order
# setkey(m, SNP)
# m2 <- m[J(kept), nomatch=0]
# cat("SNPs used (post-LDref intersect):", nrow(m2), "\n")
# if (!(lead_snp %in% m2$SNP)) stop("Lead SNP not in LD reference set: ", lead_snp)
# 
# # ---------------------------
# # Compute LD matrix (PLINK --r square gz), label it with SNP names
# # ---------------------------
# snplist_kept <- tempfile(pattern="lz_snps_kept_", fileext=".txt")
# fwrite(data.table(SNP = m2$SNP), snplist_kept, col.names=FALSE)
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
# if (!file.exists(ld_file)) stop("PLINK did not produce LD file: ", ld_file)
# 
# LD <- as.matrix(fread(cmd = paste("zcat", shQuote(ld_file)), header=FALSE))
# storage.mode(LD) <- "double"
# if (nrow(LD) != nrow(m2)) stop("LD dimension mismatch: LD=", nrow(LD), " m2=", nrow(m2))
# 
# # PLINK LD matrix has no row/col names -> attach from m2 order
# rownames(LD) <- m2$SNP
# colnames(LD) <- m2$SNP
# 
# # Save table used for plotting (nice for debugging/reuse)
# plot_table <- m2[, .(snp=SNP, chr=chr, pos=pos, p_trait1=p1, p_trait2=p2)]
# fwrite(plot_table, paste0(out_pref, "_plot_table.tsv"), sep="\t")
# 
# # ---------------------------
# # Build locuszoomr locus objects
# # Important: EnsDb uses "chr17" style seqname; we supply seqname accordingly.
# # ---------------------------
# ensdb <- EnsDb.Hsapiens.v86
# seqname_chr <- paste0("chr", chr0)
# xrange_bp <- c(lo_bp, hi_bp)
# 
# #df1 <- data.frame(snp=m2$SNP, chrom=m2$chr, pos=m2$pos, p=m2$p1)
# #df2 <- data.frame(snp=m2$SNP, chrom=m2$chr, pos=m2$pos, p=m2$p2)
# 
# df1 <- data.frame(snp=m2$SNP, chrom=paste0("chr", m2$chr), pos=m2$pos, p=m2$p1)
# df2 <- data.frame(snp=m2$SNP, chrom=paste0("chr", m2$chr), pos=m2$pos, p=m2$p2)
# 
# cat("DEBUG seqname:", seqname_chr, "\n")
# cat("DEBUG chrom unique:", paste(unique(df1$chrom), collapse=", "), "\n")
# cat("DEBUG pos range:", min(df1$pos), "-", max(df1$pos), "\n")
# cat("DEBUG nrows:", nrow(df1), "\n")
# 
# # locus() requires ens_db in your version
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
# # Plot (ggplot objects)
# # locus_ggplot is the ggplot-based LocusZoom plotter.
# # ---------------------------
# p1 <- locuszoomr::locus_ggplot(lz1) + ggtitle("Trait 1 (pQTL)")
# p2 <- locuszoomr::locus_ggplot(lz2) + ggtitle("Trait 2 (GWAS)")
# 
# ggsave(paste0(out_pref, "_trait1.png"), p1, width=7, height=5, dpi=300)
# ggsave(paste0(out_pref, "_trait2.png"), p2, width=7, height=5, dpi=300)
# 
# p_pair <- p1 | p2
# ggsave(paste0(out_pref, "_paired.png"), p_pair, width=14, height=5, dpi=300)
# 
# cat("Wrote:\n  ",
#     paste0(out_pref, "_plot_table.tsv"), "\n  ",
#     paste0(out_pref, "_trait1.png"), "\n  ",
#     paste0(out_pref, "_trait2.png"), "\n  ",
#     paste0(out_pref, "_paired.png"), "\n", sep="")
# cat("DONE:\n ")