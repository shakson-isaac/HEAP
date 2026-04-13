suppressPackageStartupMessages({
  library(data.table)
  library(pbapply)
})

CFG <- list(
  edges_dir  = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/global_edges",
  outdirbase = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges"
)

# -----------------------------
# Build expected file paths for one edge_type
# -----------------------------
expected_edge_paths <- function(edge_type, CFG, suffix = "mr_methods") {
  # suffix in {"mr_methods","heterogeneity","pleiotropy"}
  edge_type <- toupper(edge_type)
  
  edge_file <- file.path(CFG$edges_dir, paste0("edges_", edge_type, ".tsv"))
  edges <- fread(edge_file, showProgress = FALSE)
  
  make_paths <- function(edge_dir, src, tgt) {
    file.path(CFG$outdirbase, edge_dir, src, tgt, paste0(edge_dir, "_", suffix, ".tsv"))
  }
  
  if (edge_type == "EP") {
    paths <- make_paths("E_to_P", edges$Exposure, edges$Protein)
    meta  <- data.table(edge_type=edge_type, edge_dir="E_to_P", src_id=edges$Exposure, tgt_id=edges$Protein)
    
  } else if (edge_type == "ED") {
    paths <- make_paths("E_to_D", edges$Exposure, edges$Disease)
    meta  <- data.table(edge_type=edge_type, edge_dir="E_to_D", src_id=edges$Exposure, tgt_id=edges$Disease)
    
  } else if (edge_type == "PD") {
    paths_cis   <- make_paths("Pcis_to_D",   edges$Protein, edges$Disease)
    paths_trans <- make_paths("Ptrans_to_D", edges$Protein, edges$Disease)
    paths <- c(paths_cis, paths_trans)
    meta  <- rbind(
      data.table(edge_type=edge_type, edge_dir="Pcis_to_D",   src_id=edges$Protein, tgt_id=edges$Disease),
      data.table(edge_type=edge_type, edge_dir="Ptrans_to_D", src_id=edges$Protein, tgt_id=edges$Disease)
    )
    
  } else if (edge_type == "PE") {
    paths_cis   <- make_paths("Pcis_to_E",   edges$Protein, edges$Exposure)
    paths_trans <- make_paths("Ptrans_to_E", edges$Protein, edges$Exposure)
    paths <- c(paths_cis, paths_trans)
    meta  <- rbind(
      data.table(edge_type=edge_type, edge_dir="Pcis_to_E",   src_id=edges$Protein, tgt_id=edges$Exposure),
      data.table(edge_type=edge_type, edge_dir="Ptrans_to_E", src_id=edges$Protein, tgt_id=edges$Exposure)
    )
    
  } else if (edge_type == "DE") {
    paths <- make_paths("D_to_E", edges$Disease, edges$Exposure)
    meta  <- data.table(edge_type=edge_type, edge_dir="D_to_E", src_id=edges$Disease, tgt_id=edges$Exposure)
    
  } else if (edge_type == "DP") {
    paths <- make_paths("D_to_P", edges$Disease, edges$Protein)
    meta  <- data.table(edge_type=edge_type, edge_dir="D_to_P", src_id=edges$Disease, tgt_id=edges$Protein)
    
  } else stop("Unknown edge_type: ", edge_type)
  
  meta[, file := paths]
  meta
}

# Backwards-compatible alias for your existing code
expected_mr_method_paths <- function(edge_type, CFG) {
  expected_edge_paths(edge_type, CFG, suffix = "mr_methods")
}

# -----------------------------
# Safe fread that skips missing/0-byte files
# -----------------------------
safe_fread <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  sz <- file.info(path)$size
  if (is.na(sz) || sz == 0) return(NULL)  # skip 0B files
  # fread can still error on malformed files; catch
  tryCatch(
    fread(path, sep = "\t", showProgress = FALSE, ...),
    error = function(e) NULL
  )
}

# -----------------------------
# Load mr_methods for one edge_type (unchanged behavior)
# -----------------------------
load_mr_methods_one <- function(edge_type, CFG,
                                select_cols = c("method","nsnp","b","se","pval",
                                                "id.exposure","id.outcome","exposure","outcome")) {
  
  meta <- expected_edge_paths(edge_type, CFG, suffix = "mr_methods")
  
  keep <- file.exists(meta$file)
  meta <- meta[keep]
  if (nrow(meta) == 0) {
    message("No existing mr_methods for ", edge_type)
    return(data.table())
  }
  
  dt_list <- pbapply::pblapply(seq_len(nrow(meta)), function(i) {
    x <- safe_fread(meta$file[i], select = select_cols)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    cbind(meta[i, .(edge_type, edge_dir, src_id, tgt_id)], x)
  })
  
  dt_list <- Filter(Negate(is.null), dt_list)
  if (length(dt_list) == 0) return(data.table())
  rbindlist(dt_list, fill = TRUE)
}

# -----------------------------
# Load heterogeneity / pleiotropy for one edge_type
# -----------------------------
load_mr_sensitivity_one <- function(edge_type, CFG,
                                    suffix = c("heterogeneity","pleiotropy"),
                                    select_cols = NULL) {
  suffix <- match.arg(suffix)
  
  meta <- expected_edge_paths(edge_type, CFG, suffix = suffix)
  
  # keep only files that exist and are non-empty
  keep <- file.exists(meta$file) & (file.info(meta$file)$size > 0)
  keep[is.na(keep)] <- FALSE
  meta <- meta[keep]
  
  if (nrow(meta) == 0) {
    message("No existing ", suffix, " files for ", edge_type)
    return(data.table())
  }
  
  dt_list <- pbapply::pblapply(seq_len(nrow(meta)), function(i) {
    # select_cols NULL -> read all cols (safest because formats vary)
    x <- if (is.null(select_cols)) safe_fread(meta$file[i]) else safe_fread(meta$file[i], select = select_cols)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    cbind(meta[i, .(edge_type, edge_dir, src_id, tgt_id)], x)
  })
  
  dt_list <- Filter(Negate(is.null), dt_list)
  if (length(dt_list) == 0) return(data.table())
  rbindlist(dt_list, fill = TRUE)
}

# -----------------------------
# Convenience: load all 3 for one edge_type
# -----------------------------
load_mr_all_one <- function(edge_type, CFG,
                            mr_select = c("method","nsnp","b","se","pval",
                                          "id.exposure","id.outcome","exposure","outcome")) {
  list(
    mr  = load_mr_methods_one(edge_type, CFG, select_cols = mr_select),
    het = load_mr_sensitivity_one(edge_type, CFG, suffix = "heterogeneity"),
    ple = load_mr_sensitivity_one(edge_type, CFG, suffix = "pleiotropy")
  )
}

# -----------------------------
# Example usage
# -----------------------------
PD_all <- load_mr_all_one("PD", CFG)
EP_all <- load_mr_all_one("EP", CFG)
ED_all <- load_mr_all_one("ED", CFG)
DE_all <- load_mr_all_one("DE", CFG)
PE_all <- load_mr_all_one("PE", CFG)
DP_all <- load_mr_all_one("DP", CFG)

# Save the Sensitivity Results and Main Results together:
saveRDS(PD_all, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","PDres.rds"))
saveRDS(EP_all, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","EPres.rds"))

saveRDS(ED_all, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","EDres.rds"))
saveRDS(DE_all, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","DEres.rds"))

saveRDS(PE_all, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","PEres.rds"))
saveRDS(DP_all, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","DPres.rds"))


head(PD_all$het)
head(PD_all$ple)







# Your original objects (if you want to keep them)
mr_PD <- PD_all$mr
mr_EP <- EP_all$mr
mr_ED <- ED_all$mr
mr_DE <- DE_all$mr
mr_PE <- PE_all$mr
mr_DP <- DP_all$mr

het_PD <- PD_all$het
ple_PD <- PD_all$ple


# Example: load PD only (cis+trans)
mr_PD <- load_mr_methods_one("PD", CFG) #Took 8 minutes
mr_EP <- load_mr_methods_one("EP", CFG) #Took 8 minutes
mr_ED <- load_mr_methods_one("ED", CFG) #Took 8 minutes
mr_DE <- load_mr_methods_one("DE", CFG) #Took 8 minutes
mr_PE <- load_mr_methods_one("PE", CFG) #Took 8 minutes
mr_DP <- load_mr_methods_one("DP", CFG) #Took 8 minutes