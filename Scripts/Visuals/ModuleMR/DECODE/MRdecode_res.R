suppressPackageStartupMessages({
  library(data.table)
  library(pbapply)
})

CFG <- list(
  edges_dir  = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/global_edges",
  outdirbase = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE"
)

# Build expected mr_methods file paths for one edge_type
expected_mr_method_paths <- function(edge_type, CFG) {
  edge_type <- toupper(edge_type)
  edge_file <- file.path(CFG$edges_dir, paste0("edges_", edge_type, ".tsv"))
  edges <- fread(edge_file, showProgress = FALSE)
  
  if (edge_type == "EP") {
    # outdir: edges/E_to_P/ExID/protID/E_to_P_mr_methods.tsv
    paths <- file.path(CFG$outdirbase, "E_to_P", edges$Exposure, edges$Protein, "E_to_P_mr_methods.tsv")
    meta  <- data.table(edge_type=edge_type, edge_dir="E_to_P", src_id=edges$Exposure, tgt_id=edges$Protein)
    
  } else if (edge_type == "ED") {
    paths <- file.path(CFG$outdirbase, "E_to_D", edges$Exposure, edges$Disease, "E_to_D_mr_methods.tsv")
    meta  <- data.table(edge_type=edge_type, edge_dir="E_to_D", src_id=edges$Exposure, tgt_id=edges$Disease)
    
  } else if (edge_type == "PD") {
    paths_cis   <- file.path(CFG$outdirbase, "Pcis_to_D",   edges$Protein, edges$Disease, "Pcis_to_D_mr_methods.tsv")
    paths_trans <- file.path(CFG$outdirbase, "Ptrans_to_D", edges$Protein, edges$Disease, "Ptrans_to_D_mr_methods.tsv")
    paths <- c(paths_cis, paths_trans)
    meta  <- rbind(
      data.table(edge_type=edge_type, edge_dir="Pcis_to_D",   src_id=edges$Protein, tgt_id=edges$Disease),
      data.table(edge_type=edge_type, edge_dir="Ptrans_to_D", src_id=edges$Protein, tgt_id=edges$Disease)
    )
    
  } else if (edge_type == "PE") {
    paths_cis   <- file.path(CFG$outdirbase, "Pcis_to_E",   edges$Protein, edges$Exposure, "Pcis_to_E_mr_methods.tsv")
    paths_trans <- file.path(CFG$outdirbase, "Ptrans_to_E", edges$Protein, edges$Exposure, "Ptrans_to_E_mr_methods.tsv")
    paths <- c(paths_cis, paths_trans)
    meta  <- rbind(
      data.table(edge_type=edge_type, edge_dir="Pcis_to_E",   src_id=edges$Protein, tgt_id=edges$Exposure),
      data.table(edge_type=edge_type, edge_dir="Ptrans_to_E", src_id=edges$Protein, tgt_id=edges$Exposure)
    )
    
  } else if (edge_type == "DE") {
    paths <- file.path(CFG$outdirbase, "D_to_E", edges$Disease, edges$Exposure, "D_to_E_mr_methods.tsv")
    meta  <- data.table(edge_type=edge_type, edge_dir="D_to_E", src_id=edges$Disease, tgt_id=edges$Exposure)
    
  } else if (edge_type == "DP") {
    paths <- file.path(CFG$outdirbase, "D_to_P", edges$Disease, edges$Protein, "D_to_P_mr_methods.tsv")
    meta  <- data.table(edge_type=edge_type, edge_dir="D_to_P", src_id=edges$Disease, tgt_id=edges$Protein)
    
  } else stop("Unknown edge_type: ", edge_type)
  
  meta[, file := paths]
  meta
}

# Read all existing mr_methods files for one edge_type (with progress bar)
load_mr_methods_one <- function(edge_type, CFG,
                                select_cols = c("method","nsnp","b","se","pval","id.exposure","id.outcome","exposure","outcome")) {
  
  meta <- expected_mr_method_paths(edge_type, CFG)
  
  # Vectorized existence check (fast)
  keep <- file.exists(meta$file)
  meta <- meta[keep]
  
  if (nrow(meta) == 0) {
    message("No existing mr_methods for ", edge_type)
    return(data.table())
  }
  
  # Read with progress bar
  dt_list <- pbapply::pblapply(seq_len(nrow(meta)), function(i) {
    x <- fread(meta$file[i], sep="\t", showProgress=FALSE, select=select_cols)
    cbind(meta[i, .(edge_type, edge_dir, src_id, tgt_id)], x)
  })
  
  rbindlist(dt_list, fill=TRUE)
}

# Example: load PD only (cis+trans)
mr_PD <- load_mr_methods_one("PD", CFG) #Took 8 minutes
mr_EP <- load_mr_methods_one("EP", CFG) #Took 8 minutes
mr_ED <- load_mr_methods_one("ED", CFG) #Took 8 minutes
mr_DE <- load_mr_methods_one("DE", CFG) #Took 8 minutes
mr_PE <- load_mr_methods_one("PE", CFG) #Took 8 minutes
mr_DP <- load_mr_methods_one("DP", CFG) #Took 8 minutes

fwrite(mr_PD, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","DECODE","PDres.csv"))
fwrite(mr_EP, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","DECODE","EPres.csv"))
fwrite(mr_ED, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","DECODE","EDres.csv"))
fwrite(mr_DE, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","DECODE","DEres.csv"))
fwrite(mr_PE, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","DECODE","PEres.csv"))
fwrite(mr_DP, file = file.path("/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/","summary","DECODE","DPres.csv"))


# Or load everything
all_types <- c("EP","ED","PD","PE","DE","DP")
mr_all <- rbindlist(lapply(all_types, load_mr_methods_one, CFG=CFG), fill=TRUE)
