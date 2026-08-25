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

# Required packages
library(data.table)
library(stringr)

# OmicsPred reference: prefer the shared IGLOO copy, fall back to legacy.
# base path is kept self-consistent (read + write cis_trans_snps under one root).
base_load_save_path <- if (dir.exists(heap_omicspred())) heap_omicspred() else
  legacy_ukb_path("Data", "OMICSPRED")
full_snp_file_path  <- file.path(base_load_save_path, "UKB_Olink_Multiancestry_scores")

# ensure output dirs exist
dir.create(file.path(base_load_save_path, "cis_trans_snps", "cis"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_load_save_path, "cis_trans_snps", "trans"), recursive = TRUE, showWarnings = FALSE)

# read Olink protein mapping file with gene locations
olink_map <- fread(file.path(base_load_save_path, "olink_protein_map_3k_v1.tsv"), sep = "\t")

# create gene_chr column with integer chromosome values where applicable
# keep non-numeric chr values as-is (e.g., "X", "Y", "MT")
olink_map[, gene_chr := ifelse(str_detect(chr, "^\\d+$"), as.integer(chr), chr)]

# list of individual files with SNP coefficients for plasma protein expression from Omicspred
ukb_protgs <- fread(file.path(base_load_save_path, "UKBprotGSfiles.txt"), header = FALSE)

# helper to read the 'pgs_name' from the header
read_pgs_name <- function(path) {
  con <- file(path, "r")
  on.exit(close(con))
  oid <- NA_character_
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (grepl("pgs_name", line, fixed = TRUE)) {
      # header format assumed like: "## pgs_name = <OID>" or "pgs_name=<OID>"
      parts <- strsplit(line, "=", fixed = TRUE)[[1]]
      if (length(parts) >= 2) {
        oid <- str_trim(parts[2])
      }
      break
    }
  }
  oid
}

for (fname in ukb_protgs[[1]]) {
  cat(fname, "\n")
  file_path <- file.path(full_snp_file_path, fname)
  
  # read Olink ID from file header
  oid <- read_pgs_name(file_path)
  if (is.na(oid) || oid == "") {
    cat("Could not find Olink ID (pgs_name) in header.\n")
    next
  }
  cat(sprintf("Olink ID found: %s\n", oid))
  
  # get gene location from olink_map
  oid_row <- olink_map[OlinkID == oid]
  if (nrow(oid_row) > 0) {
    gene_start <- as.numeric(oid_row$gene_start[1])
    gene_end   <- as.numeric(oid_row$gene_end[1])
    lower_bound <- min(gene_start, gene_end) - 1e6
    upper_bound <- max(gene_start, gene_end) + 1e6
    gene_chr    <- oid_row$gene_chr[1]  # can be integer or character
    cat(sprintf("Gene region for Olink ID %s: Chr %s: %g - %g\n",
                oid, as.character(gene_chr), lower_bound, upper_bound))
  } else {
    cat(sprintf("Olink ID %s not found in the mapping file.\n", oid))
    next
  }
  
  # read protein genetic score file, skipping first 10 header lines
  opgs <- fread(file_path, skip = 10, sep = "\t", check.names = TRUE)
  cat("opgs dim: ", paste(dim(opgs), collapse = " x "), "\n")
  
  # Ensure expected columns exist
  needed_cols <- c("chr_name", "chr_position", "rsid")
  if (!all(needed_cols %in% names(opgs))) {
    stop(sprintf("Missing required columns in %s. Found: %s",
                 fname, paste(names(opgs), collapse = ", ")))
  }
  
  # make types consistent for filtering
  # Convert chr_name to character for robust comparison (works for numeric + X/Y/MT)
  opgs[, chr_name_chr := as.character(chr_name)]
  gene_chr_chr <- as.character(gene_chr)
  
  # filter cis variants (same chr and within ±1Mb window)
  cis <- opgs[chr_name_chr == gene_chr_chr &
                chr_position >= lower_bound &
                chr_position <= upper_bound]
  
  cat("cis dim:   ", paste(dim(cis), collapse = " x "), "\n")
  
  # trans = everything not in cis by rsid
  trans <- opgs[!rsid %in% cis$rsid]
  cat("trans dim: ", paste(dim(trans), collapse = " x "), "\n\n")
  
  # write cis and trans to .txt files
  opgs_id <- sub("_model.txt$", "", fname)
  fwrite(cis[, !("chr_name_chr")],   # drop helper column
         file = file.path(base_load_save_path, "cis_trans_snps", "cis",
                          sprintf("%s_cis_variants.txt", opgs_id)),
         sep = "\t", quote = FALSE, na = "NA")
  
  fwrite(trans[, !("chr_name_chr")],
         file = file.path(base_load_save_path, "cis_trans_snps", "trans",
                          sprintf("%s_trans_variants.txt", opgs_id)),
         sep = "\t", quote = FALSE, na = "NA")
}
