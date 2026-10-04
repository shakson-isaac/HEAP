#!/usr/bin/env Rscript
# ============================================================================
# module6_pes_disease_ladder_annotate.R   (support: post-process the S18 ladder)
# ----------------------------------------------------------------------------
# Adds to multipes_disease/pes_disease_ladder.tsv the two things
# module6_pes_disease_bootstrap.R cannot record itself. Neither needs the
# bootstrap re-run, which takes an hour.
#
# 1. dz_id -- the disease MACHINE KEY. The bootstrap records `disease` as the
#    REGEX PATTERN it matched ("e11_first_reported_non_insulin"), not the column
#    that pattern resolved to ("age_e11_..._f130708_0_0"). Without it S18 ships
#    a readable name and an ICD chapter but nothing that joins to
#    \Tref{disease_list} (SUPPLEMENT_DATA_STANDARD.md section 3;
#    SUPP_TABLE_REVIEW.md E4). Resolved with the bootstrap's own rule:
#    grep(pattern, names(DZ_df), ignore.case=TRUE), first "^age_" hit.
#
# 2. label -- resynced from ladder_slate.tsv. The slate is authoritative for
#    labels because it carries annotations a generated label cannot reconstruct,
#    notably the "(broad control)" marks on the two walking-pace rows that drive
#    S18's Note column. When the slate was first generated those marks were lost
#    and Note shipped ALL-EMPTY, which section 3 forbids and which falsified the
#    legend clause describing it. Syncing here means a slate label fix never
#    requires re-bootstrapping.
#
# Idempotent. The HEAP.rds load (~45 s, ~40 GB) is SKIPPED when every row
# already carries a dz_id, so a label-only fix runs in seconds.
#
# Writes in place, leaving a .bak_annotate copy beside the file.
#   Run:  Rscript scripts/support/module6_pes_disease_ladder_annotate.R
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages(library(data.table))
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep="")

out_dir <- heap_project_output("module6_pes_longitudinal","multipes_disease")
ladder  <- file.path(out_dir, "pes_disease_ladder.tsv")
slate   <- file.path(out_dir, "ladder_slate.tsv")
if (!file.exists(ladder)) stop("missing ", ladder)
L <- fread(ladder)
n0 <- nrow(L)

# ---- 1. dz_id (only if needed) ---------------------------------------------
if (!"dz_id" %in% names(L) || anyNA(L$dz_id) || any(!nzchar(L$dz_id))) {
  msg("resolving dz_id -- loading HEAP.rds for names(DZ_df)")
  heap <- readRDS(heap_loader_rds); dznames <- names(heap$disease$DZ_df); rm(heap); gc()
  resolve <- function(pat) {
    ac <- grep(pat, dznames, value = TRUE, ignore.case = TRUE)
    ac <- ac[grepl("^age_", ac)]
    if (!length(ac)) return(NA_character_)
    ac[1]
  }
  pats <- unique(L$disease)
  LK <- data.table(disease = pats, dz_id = vapply(pats, resolve, character(1), USE.NAMES = FALSE),
                   n_matches = vapply(pats, function(p)
                     sum(grepl("^age_", grep(p, dznames, value=TRUE, ignore.case=TRUE))), integer(1)))
  if (anyNA(LK$dz_id)) stop("unresolved disease pattern(s): ", paste(LK[is.na(dz_id)]$disease, collapse=", "))
  if (nrow(LK[n_matches > 1])) {
    msg("WARNING: pattern(s) matching >1 age_ column; first hit used, as the bootstrap does:")
    print(LK[n_matches > 1])
  }
  if ("dz_id" %in% names(L)) L[, dz_id := NULL]
  L <- merge(L, LK[, .(disease, dz_id)], by = "disease", all.x = TRUE, sort = FALSE)
  lk_dir <- file.path(heap_root, "docs", "manuscript_stats", "module6")
  dir.create(lk_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(LK, file.path(lk_dir, "ladder_dz_lookup.tsv"), sep = "\t")
  msg("dz_id resolved for ", uniqueN(LK$disease), " distinct diseases")
} else {
  msg("dz_id already present on all ", nrow(L), " rows -- skipping the HEAP.rds load")
}

# ---- 2. label, resynced from the slate -------------------------------------
if (file.exists(slate)) {
  S <- fread(slate)[, .(exposure_id, disease, slate_label = label, slate_note = note)]
  L <- merge(L, S, by = c("exposure_id","disease"), all.x = TRUE, sort = FALSE)
  chg <- L[!is.na(slate_label) & slate_label != label]
  if (nrow(chg)) {
    msg("relabelling ", nrow(chg), " row(s) from the slate:")
    print(chg[, .(from = label, to = slate_label)])
  } else msg("labels already match the slate")
  miss <- L[is.na(slate_label)]
  if (nrow(miss)) msg("WARNING: ", nrow(miss), " ladder row(s) not in the slate; label left as-is")
  L[!is.na(slate_label), label := slate_label]
  if ("note" %in% names(L)) L[, note := NULL]
  setnames(L, "slate_note", "note")
  L[is.na(note), note := ""]
  L[, slate_label := NULL]
  msg("rows carrying a note: ", sum(nzchar(L$note)))
} else msg("WARNING: no ladder_slate.tsv; label sync skipped")

setcolorder(L, c("label","exposure_id","disease","dz_id"))
stopifnot(nrow(L) == n0, !anyNA(L$dz_id))
file.copy(ladder, paste0(ladder, ".bak_annotate"), overwrite = TRUE)
fwrite(L, ladder, sep = "\t")
msg("wrote ", ladder, " (", nrow(L), " rows)")
msg("rows carrying a '(broad control)' annotation: ", sum(grepl("broad control", L$label)))
