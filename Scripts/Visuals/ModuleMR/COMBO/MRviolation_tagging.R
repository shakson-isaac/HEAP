#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(pbapply)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# CONFIG
# ============================================================

CFG_UKB <- list(
  tag        = "UKB",
  edges_dir  = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/global_edges",
  outdirbase = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges"
)

CFG_DEC <- list(
  tag        = "DECODE",
  edges_dir  = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/global_edges",
  outdirbase = "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/edges_DECODE"
)

OUTDIR <- "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MRplots/SENSITIVITY_UKB_vs_DECODE"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

CFG_ANALYSIS <- list(
  adj_method = "BH",
  alpha_q    = 0.05,  # hit threshold on adjusted p-value
  alpha_sens = 0.05,  # heterogeneity Q_pval and Egger intercept p threshold
  keep_methods = c("Inverse variance weighted", "Wald ratio"),
  edge_types = c("EP","ED","PD","PE","DE","DP")
)

# ============================================================
# PATH BUILDER (works for BOTH UKB + DECODE layouts)
# ============================================================

expected_edge_paths <- function(edge_type, CFG, suffix = "mr_methods") {
  # suffix in {"mr_methods","heterogeneity","pleiotropy"}
  edge_type <- toupper(edge_type)
  
  edge_file <- file.path(CFG$edges_dir, paste0("edges_", edge_type, ".tsv"))
  edges <- fread(edge_file, showProgress = FALSE)
  
  make_paths <- function(edge_dir, src, tgt) {
    # NOTE: matches your UKB convention: <edge_dir>_<suffix>.tsv
    # and assumes DECODE is the same naming convention.
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

safe_fread <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  sz <- file.info(path)$size
  if (is.na(sz) || sz == 0) return(NULL)
  tryCatch(
    fread(path, sep = "\t", showProgress = FALSE, ...),
    error = function(e) NULL
  )
}

load_one_suffix <- function(edge_type, CFG,
                            suffix = c("mr_methods","heterogeneity","pleiotropy"),
                            select_cols = NULL) {
  suffix <- match.arg(suffix)
  meta <- expected_edge_paths(edge_type, CFG, suffix = suffix)
  
  keep <- file.exists(meta$file) & (file.info(meta$file)$size > 0)
  keep[is.na(keep)] <- FALSE
  meta <- meta[keep]
  
  if (nrow(meta) == 0) {
    message("[", CFG$tag, "] No existing ", suffix, " files for ", edge_type)
    return(data.table())
  }
  
  dt_list <- pbapply::pblapply(seq_len(nrow(meta)), function(i) {
    x <- if (is.null(select_cols)) safe_fread(meta$file[i]) else safe_fread(meta$file[i], select = select_cols)
    if (is.null(x) || nrow(x) == 0) return(NULL)
    cbind(meta[i, .(edge_type, edge_dir, src_id, tgt_id)], x)
  })
  
  dt_list <- Filter(Negate(is.null), dt_list)
  if (length(dt_list) == 0) return(data.table())
  out <- rbindlist(dt_list, fill = TRUE)
  out[, dataset := CFG$tag]
  out
}

load_bundle_all_edges <- function(CFG, edge_types,
                                  mr_select = c("method","nsnp","b","se","pval",
                                                "id.exposure","id.outcome","exposure","outcome")) {
  mr_all  <- rbindlist(lapply(edge_types, \(et) load_one_suffix(et, CFG, "mr_methods", select_cols = mr_select)), fill=TRUE)
  het_all <- rbindlist(lapply(edge_types, \(et) load_one_suffix(et, CFG, "heterogeneity", select_cols = NULL)), fill=TRUE)
  ple_all <- rbindlist(lapply(edge_types, \(et) load_one_suffix(et, CFG, "pleiotropy", select_cols = NULL)), fill=TRUE)
  list(mr = mr_all, het = het_all, ple = ple_all)
}

safe_padj <- function(p, method = "BH") {
  if (length(p) == 0) return(p)
  if (all(is.na(p))) return(rep(NA_real_, length(p)))
  p.adjust(p, method = method)
}

# ============================================================
# BUILD "MR HIT + SENSITIVITY" TABLE FOR ONE DATASET
# ============================================================

build_mr_sensitivity_table <- function(bundle, cfg_analysis) {
  mr_all  <- as.data.table(bundle$mr)
  het_all <- as.data.table(bundle$het)
  ple_all <- as.data.table(bundle$ple)
  
  if (nrow(mr_all) == 0) return(data.table())
  
  # Keep only main MR methods
  mr_main <- mr_all[method %in% cfg_analysis$keep_methods]
  
  # Adjust p-values within edge_dir (same as your UKB script)
  mr_main[, pval_adj := safe_padj(pval, method = cfg_analysis$adj_method), by = edge_dir]
  mr_main[, mr_hit := !is.na(pval_adj) & pval_adj < cfg_analysis$alpha_q]
  
  # --- heterogeneity (IVW Q_pval) ---
  het_ivw <- copy(het_all)
  if (nrow(het_ivw)) {
    if ("method" %in% names(het_ivw)) {
      het_ivw <- het_ivw[method == "Inverse variance weighted"]
    }
    # prefer Q_pval if present; otherwise try common alternates
    qcol <- c("Q_pval","Q_pval_IVW","Q_pval_ivw","Q_p","Qpval")[c("Q_pval","Q_pval_IVW","Q_pval_ivw","Q_p","Qpval") %in% names(het_ivw)][1]
    if (!is.na(qcol)) {
      het_ivw <- het_ivw[, .(dataset, edge_dir, src_id, tgt_id, het_pval = get(qcol))]
    } else {
      het_ivw <- data.table()
    }
  }
  
  # --- pleiotropy (Egger intercept p) ---
  ple_dt <- copy(ple_all)
  if (nrow(ple_dt)) {
    # Some MR packages store egger intercept test under method == "MR Egger"
    if ("method" %in% names(ple_dt)) {
      ple_dt <- ple_dt[method %in% c("MR Egger", "Egger", "MR-Egger", "Egger intercept")]
    }
    # typical column is pval; sometimes "pval" or "p.value"
    pcol <- c("pval","p.value","p_value","Pvalue","P")[c("pval","p.value","p_value","Pvalue","P") %in% names(ple_dt)][1]
    if (!is.na(pcol)) {
      ple_dt <- ple_dt[, .(dataset, edge_dir, src_id, tgt_id, egger_pval = get(pcol))]
    } else {
      ple_dt <- data.table()
    }
  }
  
  # Merge diagnostics onto MR rows
  setkeyv(mr_main, c("dataset","edge_dir","src_id","tgt_id"))
  mr_sens <- copy(mr_main)
  
  if (exists("het_ivw") && nrow(het_ivw)) {
    setkeyv(het_ivw, c("dataset","edge_dir","src_id","tgt_id"))
    mr_sens <- het_ivw[mr_sens]
  } else {
    mr_sens[, het_pval := NA_real_]
  }
  
  if (exists("ple_dt") && nrow(ple_dt)) {
    setkeyv(ple_dt, c("dataset","edge_dir","src_id","tgt_id"))
    mr_sens <- ple_dt[mr_sens]
  } else {
    mr_sens[, egger_pval := NA_real_]
  }
  
  # Flags + pass rule (missing=pass)
  mr_sens[, het_flag   := !is.na(het_pval)   & het_pval   < cfg_analysis$alpha_sens]
  mr_sens[, pleio_flag := !is.na(egger_pval) & egger_pval < cfg_analysis$alpha_sens]
  
  mr_sens[, sens_pass :=
            (is.na(het_pval)   | het_pval   >= cfg_analysis$alpha_sens) &
            (is.na(egger_pval) | egger_pval >= cfg_analysis$alpha_sens)
  ]
  
  mr_sens[, hit_after_sens := mr_hit & sens_pass]
  mr_sens
}

# ============================================================
# RUN: LOAD BOTH DATASETS
# ============================================================

message("Loading UKB bundle...")
UKB_bundle <- load_bundle_all_edges(CFG_UKB, CFG_ANALYSIS$edge_types)

message("Loading DECODE bundle...")
DEC_bundle <- load_bundle_all_edges(CFG_DEC, CFG_ANALYSIS$edge_types)

saveRDS(UKB_bundle, "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary/COMPARE_UKB_vs_DECODE/UKBmr.rds") 
saveRDS(DEC_bundle, "/n/groups/patel/shakson_ukb/UK_Biobank/Output/MR_edges/summary/COMPARE_UKB_vs_DECODE/DECODEmr.rds") 


message("Building sensitivity tables...")
UKB_mr_sens <- build_mr_sensitivity_table(UKB_bundle, CFG_ANALYSIS)
DEC_mr_sens <- build_mr_sensitivity_table(DEC_bundle, CFG_ANALYSIS)

mr_sens_all <- rbindlist(list(UKB_mr_sens, DEC_mr_sens), fill = TRUE)

# ============================================================
# SUMMARIES: VIOLATION RATES + HIT RETENTION
# ============================================================

# 1) overall per dataset
overall_summary <- mr_sens_all[, .(
  n_edges = .N,
  n_hits = sum(mr_hit, na.rm=TRUE),
  n_hits_after = sum(hit_after_sens, na.rm=TRUE),
  
  # diagnostics availability among hits
  n_hit_with_het = sum(mr_hit & !is.na(het_pval), na.rm=TRUE),
  n_hit_with_ple = sum(mr_hit & !is.na(egger_pval), na.rm=TRUE),
  
  # violations among hits (only where diagnostic exists)
  n_het_flag = sum(mr_hit & !is.na(het_pval)   & het_flag, na.rm=TRUE),
  n_ple_flag = sum(mr_hit & !is.na(egger_pval) & pleio_flag, na.rm=TRUE)
), by = dataset]

overall_summary[, hit_retention_pct := ifelse(n_hits > 0, 100*n_hits_after/n_hits, NA_real_)]
overall_summary[, het_flag_rate_pct := ifelse(n_hit_with_het > 0, 100*n_het_flag/n_hit_with_het, NA_real_)]
overall_summary[, ple_flag_rate_pct := ifelse(n_hit_with_ple > 0, 100*n_ple_flag/n_hit_with_ple, NA_real_)]

# 2) per dataset x edge_dir
edge_summary <- mr_sens_all[, .(
  n_edges = .N,
  n_hits = sum(mr_hit, na.rm=TRUE),
  n_hits_after = sum(hit_after_sens, na.rm=TRUE),
  
  n_hit_with_het = sum(mr_hit & !is.na(het_pval), na.rm=TRUE),
  n_hit_with_ple = sum(mr_hit & !is.na(egger_pval), na.rm=TRUE),
  
  n_het_flag = sum(mr_hit & !is.na(het_pval)   & het_flag, na.rm=TRUE),
  n_ple_flag = sum(mr_hit & !is.na(egger_pval) & pleio_flag, na.rm=TRUE)
), by = .(dataset, edge_dir)]

edge_summary[, hit_retention_pct := ifelse(n_hits > 0, 100*n_hits_after/n_hits, NA_real_)]
edge_summary[, het_flag_rate_pct := ifelse(n_hit_with_het > 0, 100*n_het_flag/n_hit_with_het, NA_real_)]
edge_summary[, ple_flag_rate_pct := ifelse(n_hit_with_ple > 0, 100*n_ple_flag/n_hit_with_ple, NA_real_)]

# ============================================================
# OUTPUT TABLES
# ============================================================

fwrite(overall_summary, file.path(OUTDIR, "Sensitivity_overall_summary_UKB_vs_DECODE.tsv"), sep="\t")
fwrite(edge_summary,    file.path(OUTDIR, "Sensitivity_by_edgeDir_summary_UKB_vs_DECODE.tsv"), sep="\t")
fwrite(mr_sens_all,     file.path(OUTDIR, "MR_sensitivity_table_ALL_edges_UKB_and_DECODE.tsv"), sep="\t")

# ============================================================
# PLOTS
# ============================================================

# Plot A: hit retention per edge_dir (UKB vs DECODE)
p_retention <- ggplot(edge_summary, aes(x = edge_dir, y = hit_retention_pct, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  theme_bw() +
  labs(
    x = NULL,
    y = "% MR hits retained after sensitivity",
    title = "MR hit retention after heterogeneity + pleiotropy filtering",
    subtitle = "Missing diagnostics = pass; hits defined by IVW/Wald with BH q<0.05 within edge_dir"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

ggsave(file.path(OUTDIR, "Hit_retention_after_sensitivity_by_edgeDir.png"),
       p_retention, width = 10, height = 5, dpi = 600)

# Plot B: violation rates among hits with diagnostic available
edge_long <- rbindlist(list(
  edge_summary[, .(dataset, edge_dir, diag="Heterogeneity", rate=het_flag_rate_pct)],
  edge_summary[, .(dataset, edge_dir, diag="Pleiotropy",    rate=ple_flag_rate_pct)]
), fill=TRUE)

p_viols <- ggplot(edge_long, aes(x = edge_dir, y = rate, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ diag, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(
    x = NULL,
    y = "% flagged among MR hits\n(with diagnostic available)",
    title = "Sensitivity violations among MR hits",
    subtitle = "Heterogeneity: IVW Cochran Q p<0.05; Pleiotropy: Egger intercept p<0.05"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())

ggsave(file.path(OUTDIR, "Violation_rates_among_hits_by_edgeDir.png"),
       p_viols, width = 10, height = 7, dpi = 600)

message("\nDone. Wrote outputs to:\n  ", OUTDIR, "\n", sep="")
print(overall_summary)