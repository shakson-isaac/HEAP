#!/usr/bin/env Rscript
# ============================================================================
# module6_visit_timing.R  (support analysis)
# ----------------------------------------------------------------------------
# Years between repeat assessments (baseline i0, imaging i2, repeat-imaging i3)
# for the proteomics repeat-measured participants. Reads per-instance assessment
# timing from HEAP.rds $covars_long and writes the pairwise gap distribution for
# fig_pes_visit_timing. (Heavy: loads HEAP.rds.)
# ============================================================================
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]; if (!is.na(hit)) source(hit)
})
suppressPackageStartupMessages(library(data.table))
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

out_dir <- heap_project_output("module6_pes_longitudinal", "visit_timing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

msg("Loading HEAP.rds covars_long")
heap <- readRDS(heap_loader_rds)
cv <- as.data.table(heap$covars_long)
rm(heap); invisible(gc())
msg("covars_long: ", nrow(cv), " rows; cols: ", paste(head(names(cv), 40), collapse = ", "))

# per-instance timing: prefer assessment_year (field 53), else per-visit age
yearcol <- intersect(c("assessment_year"), names(cv))
agecol  <- grep("age_when_attended_assessment_centre", names(cv), value = TRUE)[1]
if (length(yearcol)) {
  cv[, t := suppressWarnings(as.numeric(get(yearcol)))]; unit <- "assessment_year"
} else if (!is.na(agecol)) {
  cv[, t := suppressWarnings(as.numeric(get(agecol)))]; unit <- agecol
} else stop("No assessment_year or per-visit age column in covars_long")
msg("Using timing column: ", unit)

cv <- cv[is.finite(t) & instance %in% c(0, 2, 3), .(eid, instance, t)]
w <- dcast(cv, eid ~ instance, value.var = "t")
setnames(w, c("0", "2", "3"), c("i0", "i2", "i3"), skip_absent = TRUE)
gaps <- rbindlist(list(
  if ("i2" %in% names(w)) w[is.finite(i0) & is.finite(i2), .(eid, pair = "i0->i2", gap = i2 - i0)],
  if ("i3" %in% names(w)) w[is.finite(i0) & is.finite(i3), .(eid, pair = "i0->i3", gap = i3 - i0)],
  if (all(c("i2","i3") %in% names(w))) w[is.finite(i2) & is.finite(i3), .(eid, pair = "i2->i3", gap = i3 - i2)]
), fill = TRUE)
gaps <- gaps[is.finite(gap) & gap > 0]
gaps[, unit := unit]
fwrite(gaps, file.path(out_dir, "visit_timing_gaps.tsv"), sep = "\t")

s <- gaps[, .(n = .N, mean = round(mean(gap), 2), median = median(gap),
              q25 = quantile(gap, .25), q75 = quantile(gap, .75)), by = pair]
fwrite(s, file.path(out_dir, "visit_timing_summary.tsv"), sep = "\t")
msg("Wrote gaps to ", out_dir)
print(s)
