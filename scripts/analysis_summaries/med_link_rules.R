# ============================================================================
# med_link_rules.R --- THE definition of a significant mediated link.
# ----------------------------------------------------------------------------
# WHY THIS FILE EXISTS (NUM-6). The rule was implemented twice and the two
# copies disagreed. summarize_module3_mediation.R -- which writes the cited
# deposit -- required q<0.05 AND a sign-consistent proportion mediated, with BH
# over every link in the driver class: 22,270 significant at base.
# module3_spec_sensitivity.R -- which fed the robustness figure and table --
# required only q<0.05, with BH recomputed on the n_cases>=100 subset: 23,549.
# About 5% apart at every specification. Neither was wrong; they answered
# slightly different questions under the same name, and nothing surfaced it
# because the two artifacts had never shipped together.
#
# Resolved 2026-08-08 in favor of the DEPOSIT rule: it is stricter (a link whose
# proportion mediated is not sign-consistent cannot be read as mediation), and it
# is what the already-cited med_exposure_total.tsv uses, so harmonizing changes
# no published number.
#
# Everything that counts mediated links should call med_derive_links() rather
# than reimplementing this. summarize_module3_mediation.R still carries its own
# inline copy because it writes the cited base deposit and is deliberately left
# untouched; export_med_spec_deposit.R verifies byte-for-byte that the two agree,
# so that pair cannot drift silently.
#
# Expects one row per (protein, disease, predictor) with at least:
#   nie_logHR nde_logHR nie_se nie_p predictor_class driver_present
# ============================================================================
suppressPackageStartupMessages(library(data.table))

med_derive_links <- function(L) {
  stopifnot(is.data.table(L))
  L[, total_logHR := nie_logHR + nde_logHR]
  L[, total_HR    := exp(total_logHR)]
  # PM is defined only when the indirect effect points the same way as the total.
  # NOTE the clamp: the retired builder's header said it EXCLUDED PM>1; its code
  # clamped to 1 and kept the row, and the code is what produced the shipped
  # numbers. Following the comment instead shifts every proportion mediated.
  L[, prop_mediated := fifelse(sign(nie_logHR) == sign(total_logHR) & total_logHR != 0,
                               nie_logHR / total_logHR, NA_real_)]
  L[prop_mediated < 0, prop_mediated := NA_real_]
  L[, pm_clamped := !is.na(prop_mediated) & prop_mediated > 1]
  L[prop_mediated > 1, prop_mediated := 1]
  L[, q := NA_real_]
  L[!is.na(nie_p), q := p.adjust(nie_p, "BH"), by = predictor_class]
  # estimable separates a STRUCTURAL zero (no cis/trans variant, so no genetic
  # pathway -- the zero is the finding) from a cell where nothing was estimated.
  L[, estimable := !(is.na(nie_se) & driver_present %in% c(TRUE, "TRUE"))]
  unest <- which(!L$estimable)
  for (cc in c("nie_HR","nie_l95","nie_u95","nde_HR","nde_l95","nde_u95",
               "total_HR","prop_mediated","nie_logHR"))
    if (cc %in% names(L)) set(L, i = unest, j = cc, value = NA_real_)
  L[, sig := !is.na(q) & q < 0.05 & !is.na(prop_mediated)]
  L[]
}
