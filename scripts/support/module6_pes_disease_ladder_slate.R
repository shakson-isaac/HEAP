#!/usr/bin/env Rscript
# ============================================================================
# module6_pes_disease_ladder_slate.R   (support: the S18 ladder slate)
# ----------------------------------------------------------------------------
# Emits the exposure->disease slate that module6_pes_disease_bootstrap.R runs to
# produce multipes_disease/pes_disease_ladder.tsv (Supplementary Table S18).
#
# WHY THIS EXISTS. The previous slate, ladder_candidates.tsv, was a hand/screen
# selected list with no generator: rows 6-22 were simply the top 17 pairs by
# dC_beyondE out of a 165-exposure x 181-disease universe. Nothing recorded that
# rule, so the table read as an evidence locker while being a highlight reel
# (SUPPLEMENT_DATA_STANDARD.md section 4 -- a cap must be stated). It also
# covered only 2 of the 12 pairs Fig 6d draws a ladder for, so table and figure
# described the same analysis on almost disjoint slates (rule 1).
#
# THE RULE, which the S18 legend must state:
#   the union of
#     (A) every exposure->disease pair for which Fig 6d draws a C-index ladder
#         -- the curated slate in module6_quadrant_ladders.R plus the one pair
#         panel d pulls out of quadrant_scan.tsv (pack-years -> emphysema); and
#     (B) the 22 previously shipped screen-selected candidate pairs, retained so
#         the rebuild drops no row a reader may already have cited,
#   deduplicated by (exposure, ICD chapter).
#
# Disease patterns are regexes resolved by grep against names(DZ_df) with a
# leading "age_" filter, exactly as module6_quadrant_ladders.R and
# module6_pes_disease_bootstrap.R both do. One canonical pattern per ICD chapter:
# the old candidate file spelled the same column two ways (bare
# "e11_first_reported_non_insulin" and the fully resolved
# "age_e11_..._f130708_0_0"), which made one outcome look like two.
# Where the figure and the old slate disagreed on spelling, the figure wins.
#
# Writes multipes_disease/ladder_slate.tsv (cols: exposure_id, disease, label).
# ============================================================================
local({ cand <- c(Sys.getenv("HEAP_PATHS_FILE",""),"/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  h <- cand[nzchar(cand)&file.exists(cand)][1]; if(!is.na(h)) source(h) })
suppressPackageStartupMessages(library(data.table))
# heap_exposure_label() lives with the figure label helpers, not in 00_paths.R
local({ f <- file.path(Sys.getenv("HEAP_ROOT","/n/groups/patel/shakson_ukb/HEAP"),
                       "scripts","visualizations","common","label_helpers.R")
  if (file.exists(f)) source(f) else stop("missing label_helpers.R: ", f) })
msg <- function(...) cat(format(Sys.time(),"[%H:%M:%S] "),...,"\n",sep="")

out_dir <- heap_project_output("module6_pes_longitudinal","multipes_disease")

# ---- canonical regex + readable label per ICD chapter ----------------------
# labels match the dz_name map in HEAP_manuscript/supp/table_specs/pes_cindex.R
ICD <- fread(text='icd\tpattern\tlabel
E11\te11_first_reported_non_insulin\tType 2 diabetes
E66\te66_first_reported_obesity\tObesity
E78\te78_first_reported_disorders_of_lipoprotein\tLipid disorder
F10\tf10_first_reported_mental_and_behavioural_disorders_due_to_use_of_alcohol\tAlcohol use disorder
F17\tf17_first_reported_mental_and_behavioural_disorders_due_to_use_of_tobacco\tTobacco use disorder
F32\tf32_first_reported_depressive\tDepression
I25\ti25_first_reported_chronic_ischaemic\tIschaemic heart disease
J43\tj43_first_reported_emphysema\tEmphysema
J44\tj44_first_reported_other_chronic_obstructive\tCOPD
K70\tk70_first_reported_alcoholic_liver\tAlcoholic liver disease
K76\tk76_first_reported_other_diseases_of_liver\tOther diseases of liver', sep='\t', header=TRUE)

icd_of <- function(x) toupper(sub("^(?:age_)?([a-z][0-9]+)_.*$", "\\1", x))

# ---- (A) the pairs Fig 6d draws a ladder for ------------------------------
# 11 from module6_quadrant_ladders.R's curated slate + the pack-years -> emphysema
# row that fig_m6_panel_d.R pulls out of quadrant_scan.tsv (its `ladX`).
FIG <- fread(text='exposure_id\ticd
current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days\tJ44
current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days\tJ43
alcohol_intake_frequency_f1558_0_0\tF10
number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0\tE11
summed_days_activity_f22033_0_0\tE11
coffee_intake_f1498_0_0\tE11
oily_fish_intake_f1329_0_0\tI25
oily_fish_intake_f1329_0_0\tE78
nap_during_day_f1190_0_0\tF32
usual_walking_pace_f924_0_0\tE11
pm2_5_mean\tE11
pack_years_of_smoking_f20161_0_0\tJ43', sep='\t', header=TRUE)
FIG[, `:=`(src = "fig6d", prev_label = NA_character_)]

# ---- (B) the previously shipped candidate slate ---------------------------
prev_f <- file.path(out_dir, "ladder_candidates.tsv")
if (!file.exists(prev_f)) stop("missing previous slate: ", prev_f)
# keep the previous slate's label verbatim: it carries annotations the generated
# label cannot reconstruct -- notably the "(broad control)" marks on the two
# walking-pace rows, which drive the Note column of Supplementary Table S18.
# Regenerating those labels blindly silently empties that column (section 3
# forbids an all-empty column) and falsifies the legend clause describing it.
PREV <- fread(prev_f)[, .(exposure_id, icd = icd_of(disease), src = "prev",
                          prev_label = label)]

# ---- union, deduplicated by (exposure, ICD); figure rows win --------------
S <- unique(rbind(FIG, PREV), by = c("exposure_id","icd"))
# dedup keeps the Fig 6d row for a pair present in both slates, whose prev_label
# is NA -- so re-attach the previous label by key, or the annotation is lost for
# exactly the pairs both slates share (e.g. walking pace -> T2D, a broad control).
S[, prev_label := NULL]
S <- merge(S, PREV[, .(exposure_id, icd, prev_label)],
           by = c("exposure_id","icd"), all.x = TRUE)
stopifnot(!anyDuplicated(S[, .(exposure_id, icd)]))
unknown <- setdiff(S$icd, ICD$icd)
if (length(unknown)) stop("no canonical pattern for ICD: ", paste(unknown, collapse=", "))

S <- merge(S, ICD, by = "icd", sort = FALSE)
# The annotation travels in its OWN column, not inside the label. Reusing the
# previous label verbatim would keep "(broad control)" but also its exposure
# wording, so one exposure_id could display under two names in the same sheet
# ("Alcohol freq." and "Alcohol frequency"). Generate every label the same way,
# and carry the annotation separately -- which also stops the S18 transform
# having to parse free text to populate its Note column.
S[, note := fifelse(!is.na(prev_label) & grepl("broad control", prev_label, fixed = TRUE),
                    "Reverse-causation control", "")]
S[, label := paste0(heap_exposure_label(exposure_id), " -> ", label)]
setorder(S, icd, exposure_id)
OUT <- S[, .(exposure_id, disease = pattern, label, note)]

# every exposure must have a PES out-of-fold file or the bootstrap silently drops it
od <- heap_project_output("module6_pes_longitudinal","base")
miss <- unique(OUT$exposure_id)[!file.exists(file.path(od,
          paste0("PESlong_base_", unique(OUT$exposure_id), "_TrainOOF.tsv")))]
if (length(miss)) stop("no TrainOOF file for: ", paste(miss, collapse=", "))

fwrite(OUT, file.path(out_dir, "ladder_slate.tsv"), sep = "\t")
msg("wrote ", file.path(out_dir,"ladder_slate.tsv"), " (", nrow(OUT), " pairs; ",
    S[src=="fig6d", .N], " from Fig 6d, ", S[src=="prev", .N], " carried over)")
print(OUT[, .(exposure_id = substr(exposure_id,1,42), disease = substr(disease,1,34))])
