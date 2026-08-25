#!/usr/bin/env Rscript

# ============================================================================
# summarize_replicated_associations.R  —  Module 2 downstream analysis
# ----------------------------------------------------------------------------
# Canonical producer of the REPLICATED exposure->protein (E) and GxE association
# tables. "Replicated" = Bonferroni-significant in BOTH the train and test split
# of the Module 2 univariate associations.
#
# This replaces the legacy logic that was hidden inside a visualization script
# (analyHEAP_association.R building HEAPassoc@HEAPsig$...$sigBOTH, written out by
# Module2/HEAPassoc_table.R). Statistical analysis lives here, in the module;
# figures/tables downstream just read the output.
#
# Source of truth (per protein):
#   heap_project_output("module2", <covarType>)/univar_assoc_<idx>.rds
#     $train / $test = list(statE[[1]], statGxE[[2]], statR2[[3]], statFblock[[4]])
#   statE  cols : ID, Eid, Category, Estimate, Std. Error, t value, Pr(>|t|),
#                 R2, adj.R2, samplesize, omicID
#   statGxE cols: ID(=E:omic_GScis/GStrans), Eid, Category, E_id, E_term,
#                 G_component, Estimate, ..., Pr(>|t|), ..., omicID
#
# Output (canonical IGLOO HEAP module 2 area):
#   heap_project_output("module2", "ReplicatedEassoc.csv")     # read by Module5_load.R
#   heap_project_output("module2", "ReplicatedGxEassoc.csv")
#   heap_project_output("module2", "ReplicatedEassoc_<covarType>.csv")  (provenance copy)
#   heap_project_output("module2", "ReplicatedGxEassoc_<covarType>.csv")
#
# Replication rule (matches the legacy definition):
#   threshold = 0.05 / N   where N = number of associations in the train|test
#               outer join (standard Bonferroni over all tested associations)
#   keep rows with Pr(>|t|)_train < threshold AND Pr(>|t|)_test < threshold
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R \
#     Rscript scripts/module2_associations/summarize_replicated_associations.R Type5
# ============================================================================

local({
  candidates <- c(
    Sys.getenv("HEAP_PATHS_FILE", unset = ""),
    file.path(getwd(), "workflow", "00_paths.R"),
    file.path(getwd(), "..", "workflow", "00_paths.R"),
    "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R"
  )
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Could not locate workflow/00_paths.R (set HEAP_PATHS_FILE).")
  source(hit)
})

suppressPackageStartupMessages({ library(data.table) })

args      <- commandArgs(trailingOnly = TRUE)
covarType <- if (length(args) >= 1) args[1] else "base"   # canonical primary = base

# --- locate Module 2 outputs ------------------------------------------------
# IGLOO canonical, local staging fallback. After the 2026 covariate-restructure
# Module 2 writes experiment-nested dirs module2/<experiment>/<covarType>/;
# older runs were flat module2/<covarType>/. Try flat first, then nested
# (preferring the most-complete run), under both roots.
.has_univar <- function(d)
  length(list.files(d, pattern = "^univar_assoc_.*\\.rds$")) > 0L
.find_m2_dir <- function(root_fn) {
  flat <- root_fn("module2", covarType)
  if (dir.exists(flat) && .has_univar(flat)) return(flat)
  nested <- Sys.glob(file.path(root_fn("module2"), "*", covarType))
  nested <- nested[dir.exists(nested) & vapply(nested, .has_univar, logical(1))]
  if (length(nested) == 0L) return(NA_character_)
  # prefer the run with the most protein files (the complete one)
  nested[order(-vapply(nested,
    function(d) length(list.files(d, pattern = "^univar_assoc_.*\\.rds$")),
    integer(1)))][1]
}
m2_dir <- .find_m2_dir(heap_project_output)
if (is.na(m2_dir)) m2_dir <- .find_m2_dir(heap_output)
if (is.na(m2_dir))
  stop("No univar_assoc_*.rds for covarType=", covarType,
       " under module2/", covarType, " or module2/*/", covarType,
       "\nRun Module 2 for covarType=", covarType, " first.", call. = FALSE)

files <- list.files(m2_dir, pattern = "^univar_assoc_.*\\.rds$", full.names = TRUE)

message("Replicated-association summary | covarType=", covarType,
        " | ", length(files), " protein files\n  source: ", m2_dir)

# --- aggregate one component (1=statE, 2=statGxE) across proteins for a split -
aggregate_component <- function(files, comp_idx, split) {
  rbindlist(lapply(files, function(f) {
    obj <- readRDS(f)
    part <- obj[[split]]
    if (is.null(part) || length(part) < comp_idx || is.null(part[[comp_idx]]))
      return(NULL)
    as.data.table(part[[comp_idx]])
  }), fill = TRUE)
}

PCOL <- "Pr(>|t|)"
KEYS <- c("ID", "omicID")

# --- replication for a given component ---------------------------------------
build_replicated <- function(comp_idx, label) {
  train <- unique(aggregate_component(files, comp_idx, "train"))  # drop exact dup rows
  test  <- unique(aggregate_component(files, comp_idx, "test"))
  if (nrow(train) == 0L || nrow(test) == 0L)
    stop("Empty ", label, " table for ", covarType,
         " (train=", nrow(train), ", test=", nrow(test), ").", call. = FALSE)

  # If a key is duplicated after unique() (e.g. identical metrics under >1 model),
  # collapse to one row per key to keep the train/test join 1:1.
  train <- train[, .SD[1L], by = KEYS]
  test  <- test[,  .SD[1L], by = KEYS]

  nonkey_tr <- setdiff(names(train), KEYS)
  nonkey_te <- setdiff(names(test),  KEYS)
  setnames(train, nonkey_tr, paste0(nonkey_tr, "_train"))
  setnames(test,  nonkey_te, paste0(nonkey_te, "_test"))

  all_assoc <- merge(train, test, by = KEYS, all = TRUE)   # outer join = all tested
  all_assoc[, AssocID := paste0(omicID, ":", ID)]

  thr  <- 0.05 / nrow(all_assoc)                            # Bonferroni over all assoc
  p_tr <- paste0(PCOL, "_train"); p_te <- paste0(PCOL, "_test")
  sig  <- all_assoc[!is.na(get(p_tr)) & !is.na(get(p_te)) &
                    get(p_tr) < thr & get(p_te) < thr]

  message(sprintf("  %-5s: tested=%d  replicated(sigBOTH)=%d  Bonferroni p<%.3g",
                  label, nrow(all_assoc), nrow(sig), thr))
  list(sig = sig, threshold = thr, n_tested = nrow(all_assoc))
}

E   <- build_replicated(1L, "E")
GxE <- build_replicated(2L, "GxE")

# --- write canonical + provenance copies -------------------------------------
dir.create(heap_project_output("module2"), recursive = TRUE, showWarnings = FALSE)

write_one <- function(dt, canonical_name, provenance_name) {
  fwrite(dt, heap_project_output("module2", canonical_name), sep = ",")
  fwrite(dt, heap_project_output("module2", provenance_name), sep = ",")
}
write_one(E$sig,   "ReplicatedEassoc.csv",
          paste0("ReplicatedEassoc_",   covarType, ".csv"))
write_one(GxE$sig, "ReplicatedGxEassoc.csv",
          paste0("ReplicatedGxEassoc_", covarType, ".csv"))

# small machine-readable provenance record
meta <- data.table(
  covarType        = covarType,
  n_protein_files  = length(files),
  E_n_tested       = E$n_tested,   E_n_replicated   = nrow(E$sig),
  E_bonferroni_p   = E$threshold,
  GxE_n_tested     = GxE$n_tested, GxE_n_replicated = nrow(GxE$sig),
  GxE_bonferroni_p = GxE$threshold,
  generated_at     = as.character(Sys.time())
)
fwrite(meta, heap_project_output("module2", paste0("replicated_associations_meta_", covarType, ".csv")))

message("\nWrote (canonical, read by Module5_load.R):\n  ",
        heap_project_output("module2", "ReplicatedEassoc.csv"), "\n  ",
        heap_project_output("module2", "ReplicatedGxEassoc.csv"))
message("DONE.")
