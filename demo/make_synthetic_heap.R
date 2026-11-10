#!/usr/bin/env Rscript
# ============================================================================
# demo/make_synthetic_heap.R
# ----------------------------------------------------------------------------
# Writes a SYNTHETIC analysis dataset with the same structure as the HEAP.rds
# object built by scripts/loaders/HEAP_loader.R, plus the per-protein genetic
# score files and OmicsPred ID map that Modules 1 and 2 read. Nothing here is
# derived from UK Biobank data: every value is simulated.
#
# Proteins, exposures and diseases have generic labels (protein_1, exposure_1,
# disease_1). The exposure set keeps the real structure -- 169 exposures in the
# real categories with the real types -- and the `base` covariates keep their
# real names. The modules read these names from config/, so this script also
# writes a demo copy of the two config files that define them; run_demo.sh runs
# the unmodified code against that copy:
#   <demo_config>/exposure_sets/analysis_exposures.tsv
#   <demo_config>/protein_sets/omicspred_proteins.txt
#
# Each protein is simulated as
#   y = 0.30 * cis score + 0.15 * trans score + 0.25 * x1 + 0.15 * x2
#       + 0.15 * x1 * cis score + covariate effects + N(0, 1) noise
# where x1 and x2 are two continuous exposures (written to planted_effects.tsv),
# so the demo output can be checked against a known answer.
#
# Usage:  Rscript demo/make_synthetic_heap.R <out_dir> <demo_config> [n_participants] [n_proteins]
#   reads the real config from $HEAP_REPO (the repository).
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(yaml) })

args   <- commandArgs(trailingOnly = TRUE)
out    <- args[1]
dcfg   <- args[2]
N      <- if (length(args) >= 3) as.integer(args[3]) else 2000L
NPROT  <- if (length(args) >= 4) as.integer(args[4]) else 6L

repo <- Sys.getenv("HEAP_REPO", unset = ".")
cfgp <- function(...) file.path(repo, "config", ...)
set.seed(20261109)

# ---- names: real structure, generic labels ----------------------------------
ae <- fread(cfgp("exposure_sets", "analysis_exposures.tsv"))[include == 1L]
ae[, `:=`(variable = paste0("exposure_", .I), field_id = NA_integer_, notes = "synthetic")]
base  <- read_yaml(cfgp("covariates", "covariate_sets.yml"))$covariate_sets$base$covariates
prots <- paste0("protein_", seq_len(NPROT))
pcol  <- prots

dir.create(file.path(dcfg, "exposure_sets"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(dcfg, "protein_sets"),  recursive = TRUE, showWarnings = FALSE)
fwrite(ae, file.path(dcfg, "exposure_sets", "analysis_exposures.tsv"), sep = "\t")
writeLines(prots, file.path(dcfg, "protein_sets", "omicspred_proteins.txt"))

eid <- seq_len(N)

# ---- covariates (base set) -------------------------------------------------
age <- round(runif(N, 40, 70))
sex <- rbinom(N, 1, 0.5)
cov <- data.frame(eid = eid)
for (v in base) {
  cov[[v]] <- switch(v,
    age_when_attended_assessment_centre_f21003_0_0 = age,
    sex_f31_0_0  = factor(sex),
    age2         = age^2,
    age_sex      = age * sex,
    age2_sex     = age^2 * sex,
    uk_biobank_assessment_centre_f54_0_0 = factor(sample(1:6, N, replace = TRUE)),
    rnorm(N))                                    # genetic principal components
}

# ---- exposures, by declared type -------------------------------------------
sim_exposure <- function(type) switch(type,
  binary     = rbinom(N, 1, 0.3),
  ordinal    = sample(0:4, N, replace = TRUE),
  continuous = rnorm(N))
E <- lapply(split(ae, ae$category), function(d) {
  df <- data.frame(eid = eid)
  for (i in seq_len(nrow(d))) df[[d$variable[i]]] <- sim_exposure(d$variable_type[i])
  df
})
# Missingness: 10% of participants lack one exposure value. (Module 1 fits all
# exposures jointly on complete cases, so missingness spread evenly over every
# field would leave almost no complete rows.)
miss_ids <- sample(eid, round(0.10 * N))
for (i in miss_ids) {
  v <- sample(ae$variable, 1); cat <- ae[variable == v, category]
  E[[cat]][i, v] <- NA
}

# Two continuous exposures carry the planted effects.
cont <- ae[variable_type == "continuous"]
x1n <- cont$variable[1]; x2n <- cont$variable[2]
getx <- function(v) { cat <- ae[variable == v, category]; E[[cat]][[v]] }
x1 <- getx(x1n); x2 <- getx(x2n)

# ---- genetic scores, OmicsPred map, proteins -------------------------------
op_id <- sprintf("OPGS%06d", seq_along(prots))
gdir  <- file.path(out, "igloo", "UKB")
dir.create(file.path(gdir, "ProtGScis"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(gdir, "ProtGStrans"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(gdir, "OMICSPRED"),   recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(Gene = prots, OMICSPRED_ID = op_id),
       file.path(gdir, "OMICSPRED", "UKB_Olink_multi_ancestry_models_val_results_portal.csv"))

m1 <- mean(x1, na.rm = TRUE); s1 <- sd(x1, na.rm = TRUE)
m2 <- mean(x2, na.rm = TRUE); s2 <- sd(x2, na.rm = TRUE)
zz <- function(v, m, s) { a <- (v - m) / s; a[is.na(a)] <- 0; a }
sim_protein <- function(cis, trans, xa, xb, age, sex)
  0.30 * cis + 0.15 * trans + 0.25 * zz(xa, m1, s1) + 0.15 * zz(xb, m2, s2) +
  0.15 * zz(xa, m1, s1) * cis + 0.02 * (age - 55) + 0.20 * sex + rnorm(length(cis))

prot <- data.frame(eid = eid); GS <- list()
for (j in seq_along(prots)) {
  cis <- rnorm(N); trans <- rnorm(N); GS[[j]] <- list(cis = cis, trans = trans)
  fwrite(data.table(`#IID` = eid, SCORE1_AVG = cis),
         file.path(gdir, "ProtGScis", paste0(op_id[j], ".sscore")), sep = "\t")
  fwrite(data.table(`#IID` = eid, SCORE1_AVG = trans),
         file.path(gdir, "ProtGStrans", paste0(op_id[j], ".sscore")), sep = "\t")
  y <- sim_protein(cis, trans, x1, x2, age, sex)
  y[runif(N) < 0.03] <- NA
  prot[[pcol[j]]] <- y
}

# ---- incident disease (Module 3) --------------------------------------------
# Disease 1: hazard rises 0.5 log-HR per SD of the FIRST protein, plus a small
# direct effect of x1, so only that protein mediates x1 -> disease 1.
# Disease 2: age only (no protein or exposure effect).
# Value = age at first report; NA = no event within 13 years of follow-up.
FU <- 13
dz_ids <- c("disease_1", "disease_2")
p1 <- prot[[pcol[1]]]; p1z <- (p1 - mean(p1, na.rm = TRUE)) / sd(p1, na.rm = TRUE); p1z[is.na(p1z)] <- 0
event_age <- function(lp, rate) { t <- rexp(N, rate * exp(lp)); ifelse(t < FU, age + t, NA) }
DZ_df <- data.frame(eid = eid,
  d1 = event_age(0.5 * p1z + 0.2 * zz(x1, m1, s1) + 0.03 * (age - 55), 0.017),
  d2 = event_age(0.03 * (age - 55), 0.017),
  recode_age_of_assessment_0_0 = age,
  recode_age_of_death_0_0 = NA_real_,
  age_of_removal_0_0 = NA_real_,
  age_of_lastfollowup = age + FU)
names(DZ_df)[2:3] <- dz_ids

# ---- repeat visits (Module 6) -----------------------------------------------
# 20% of participants return at instance 2 (+9 y); half of those also at
# instance 3 (+11 y). x1 changes between visits (visit-to-visit r = 0.7); every
# other exposure is carried forward; proteins are re-simulated from the same
# genetic scores and the new x1.
hold <- sort(sample(eid, round(0.20 * N)))
i3   <- sort(sample(hold, round(length(hold) / 2)))
visit <- function(ids, inst, dage) {
  k <- match(ids, eid)
  x1f <- 0.7 * x1[k] + sqrt(1 - 0.7^2) * s1 * rnorm(length(k)) + 0.3 * m1
  Ef <- lapply(E, function(df) { d <- df[k, , drop = FALSE]
    if (x1n %in% names(d)) d[[x1n]] <- x1f; cbind(d[1], instance = inst, d[-1]) })
  pf <- data.frame(eid = ids, instance = inst)
  for (j in seq_along(prots))
    pf[[pcol[j]]] <- sim_protein(GS[[j]]$cis[k], GS[[j]]$trans[k], x1f, x2[k],
                                 age[k] + dage, sex[k])
  cf <- cov[k, , drop = FALSE]
  a <- cf$age_when_attended_assessment_centre_f21003_0_0 + dage
  cf$age_when_attended_assessment_centre_f21003_0_0 <- a
  cf$age2 <- a^2; cf$age_sex <- a * sex[k]; cf$age2_sex <- a^2 * sex[k]
  list(E = Ef, prot = pf, cov = cbind(cf[1], instance = inst, cf[-1]))
}
V0 <- list(E = lapply(E, function(df) cbind(df[1], instance = 0L, df[-1])),
           prot = cbind(prot[1], instance = 0L, prot[-1]),
           cov  = cbind(cov[1], instance = 0L, cov[-1]))
V2 <- visit(hold, 2L, 9); V3 <- visit(i3, 3L, 11)
E_long <- setNames(lapply(names(E), function(n) rbind(V0$E[[n]], V2$E[[n]], V3$E[[n]])), names(E))
split_df <- data.frame(eid = eid, has_i0 = TRUE, has_i2 = eid %in% hold, has_i3 = eid %in% i3)
split_df$n_visits  <- 1L + split_df$has_i2 + split_df$has_i3
split_df$instances <- ifelse(split_df$has_i3, "0,2,3", ifelse(split_df$has_i2, "0,2", "0"))
split_df$analysis_set <- ifelse(split_df$has_i2 | split_df$has_i3,
                                "holdout_repeat_proteomics", "train_baseline_only")

# Module 6 reads its exposure list from this file (written by Module6_config.R in a full run).
m6 <- file.path(out, "igloo", "UKB", "HEAP", "output", "module6_pes_test")
dir.create(m6, recursive = TRUE, showWarnings = FALSE)
fwrite(ae[, .(exposure_id = variable,
              exposure_type = ifelse(variable_type == "binary", "binary", "continuous"))],
       file.path(m6, "exposure_specs.tsv"), sep = "\t")

# ---- assemble the HEAP object (same element names as HEAP_loader.R) --------
heap <- list(
  meta = list(object_version = "synthetic-demo", created_at = format(Sys.time()),
              instances = c(0L, 2L, 3L), canonical_instance = 0L, projID = NA_integer_,
              n_proteins = length(prots), n_eid_baseline = N,
              n_exposure_cats = length(E), exposure_cats = names(E),
              skip_disease = FALSE, prevalent_disease_cols = character(0)),
  protIDs       = prots,
  Elist_names   = names(E),
  Eid_cat       = data.frame(Eid = ae$variable, Category = ae$category),
  ordinalIDs    = ae[variable_type == "ordinal", variable],
  covars_list   = base,
  prot_long     = rbind(V0$prot, V2$prot, V3$prot),
  E_long        = E_long,
  covars_long   = rbind(V0$cov, V2$cov, V3$cov),
  split_df      = split_df,
  prot_baseline = prot,
  E_baseline    = E,
  covars_baseline = cov,
  disease       = list(DZ_df = DZ_df, DZ_ids = dz_ids)
)
saveRDS(heap, file.path(out, "HEAP_synthetic.rds"))

fwrite(data.table(protein = prots, planted_cis = 0.30, planted_trans = 0.15,
                  exposure_x1 = x1n, planted_x1 = 0.25,
                  exposure_x2 = x2n, planted_x2 = 0.15,
                  planted_x1_by_cis = 0.15,
                  disease_mediated = dz_ids[1], disease_null = dz_ids[2],
                  mediator_protein = prots[1], planted_protein_logHR = 0.5,
                  n_holdout = length(hold), visit_r_x1 = 0.7),
       file.path(out, "planted_effects.tsv"), sep = "\t")

cat(sprintf("Synthetic HEAP object: %d participants, %d proteins (%s), %d exposures in %d categories\n",
            N, length(prots), paste(prots, collapse = ", "), nrow(ae), length(E)))
cat("Planted exposure effects: ", x1n, " (0.25 SD), ", x2n, " (0.15 SD)\n", sep = "")
cat("Wrote ", normalizePath(out), "\n", sep = "")
