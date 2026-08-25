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

script_file <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_file) == 0L) getwd() else dirname(normalizePath(sub("^--file=", "", script_file[1L])))
source(file.path(script_dir, "common.R"))

plot_histogram <- function(values, main, xlab, out_path) {
  stats_values <- values[is.finite(values)]
  if (length(stats_values) == 0L) {
    return(invisible(NULL))
  }
  grDevices::pdf(out_path, width = 8, height = 6)
  hist(stats_values, breaks = 30, main = main, xlab = xlab, col = "grey80", border = "white")
  graphics::rug(stats_values)
  grDevices::dev.off()
}

plot_scatter <- function(x, y, protein_ids, out_path) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) {
    return(invisible(NULL))
  }
  grDevices::pdf(out_path, width = 8, height = 6)
  plot(
    x[keep], y[keep],
    xlab = "Population GxE variance proportion",
    ylab = "Predictive OOF delta-R2 for GxE",
    main = "Population GxE variance vs predictive OOF delta-R2",
    pch = 19,
    col = "steelblue"
  )
  abline(stats::lm(y[keep] ~ x[keep]), col = "firebrick", lwd = 2)
  graphics::text(x[keep], y[keep], labels = protein_ids[keep], pos = 3, cex = 0.6)
  grDevices::dev.off()
}

collect_summary_rows <- function(model_dir) {
  row_files <- Sys.glob(file.path(model_dir, "*_summary.tsv"))
  row_files <- sort(row_files)
  if (length(row_files) == 0L) {
    return(NULL)
  }
  do.call(rbind, lapply(row_files, read_tsv, header = TRUE))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: summarize_population_architecture.R <config.R> <run_id> <covar_spec> [--model=primary|sensitivity] [--center-exposures=true|false]", call. = FALSE)
}

config_path <- args[1L]
run_id <- args[2L]
covar_spec_name <- args[3L]
opts <- parse_optional_args(args[-(1L:3L)])
model_name <- get_opt(opts, "model", "primary")
center_exposures <- as_bool(get_opt(opts, "center_exposures", TRUE), default = TRUE)
exposure_mode <- if (center_exposures) "centered" else "uncentered"
cfg <- load_config(config_path)
min_protein_n <- as.integer(get_opt(opts, "min_protein_n", cfg$min_protein_n %||% 2000L))
paths <- resolve_run_paths(cfg, run_id, covar_spec_name, exposure_mode = exposure_mode)
model_dir_override <- get_opt(opts, "model_dir", NULL)
summary_dir_override <- get_opt(opts, "summary_dir", NULL)
plots_dir_override <- get_opt(opts, "plots_dir", NULL)
model_dir <- if (!is.null(model_dir_override) && nzchar(model_dir_override)) {
  ensure_dir(model_dir_override)
} else if (model_name == "primary") {
  paths$models_primary
} else {
  paths$models_sensitivity
}
summary_dir <- if (!is.null(summary_dir_override) && nzchar(summary_dir_override)) {
  ensure_dir(summary_dir_override)
} else {
  paths$summary
}
plots_dir <- if (!is.null(plots_dir_override) && nzchar(plots_dir_override)) {
  ensure_dir(plots_dir_override)
} else {
  paths$plots
}
summary_path <- file.path(summary_dir, paste0("per_protein_summary_", model_name, ".tsv"))
collected <- collect_summary_rows(model_dir)
summary_df <- if (!is.null(collected)) {
  collected
} else if (file.exists(summary_path)) {
  read_tsv(summary_path, header = TRUE)
} else {
  stopf("Per-protein summary not found and no per-protein row files were available for model %s.", model_name)
}
write_tsv(summary_df, summary_path)

aggregate_summary <- data.frame(
  metric = c(
    "n_proteins",
    "n_converged",
    "median_n",
    "median_variance_G",
    "median_variance_E",
    "median_variance_GxE",
    "mean_variance_Covars_fixed"
  ),
  value = c(
    nrow(summary_df),
    sum(summary_df$converged %in% c(TRUE, "TRUE")),
    stats::median(as.numeric(summary_df$n), na.rm = TRUE),
    stats::median(as.numeric(summary_df$variance_G), na.rm = TRUE),
    stats::median(as.numeric(summary_df$variance_E), na.rm = TRUE),
    stats::median(as.numeric(summary_df$variance_GxE), na.rm = TRUE),
    mean(as.numeric(summary_df$variance_Covars_fixed), na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write_tsv(aggregate_summary, file.path(summary_dir, paste0("aggregate_summary_", model_name, ".tsv")))

plot_histogram(
  as.numeric(summary_df$variance_G),
  main = "Distribution of variance_G across proteins",
  xlab = "variance_G",
  out_path = file.path(plots_dir, paste0("variance_G_", model_name, ".pdf"))
)
plot_histogram(
  as.numeric(summary_df$variance_E),
  main = "Distribution of variance_E across proteins",
  xlab = "variance_E",
  out_path = file.path(plots_dir, paste0("variance_E_", model_name, ".pdf"))
)
plot_histogram(
  as.numeric(summary_df$variance_GxE),
  main = "Distribution of variance_GxE across proteins",
  xlab = "variance_GxE",
  out_path = file.path(plots_dir, paste0("variance_GxE_", model_name, ".pdf"))
)
plot_scatter(
  as.numeric(summary_df$variance_GxE),
  as.numeric(summary_df$predictive_gxe_delta_r2),
  summary_df$protein,
  out_path = file.path(plots_dir, paste0("variance_GxE_vs_predictive_deltaR2_", model_name, ".pdf"))
)

report_lines <- c(
  paste("# Population Architecture Report -", model_name),
  "",
  paste("- Run ID:", run_id),
  paste("- Covariate specification:", covar_spec_name),
  paste("- Exposure kernel mode:", exposure_mode),
  paste("- Proteins analysed:", nrow(summary_df)),
  paste("- Converged fits:", sum(summary_df$converged %in% c(TRUE, "TRUE"))),
  paste("- Median sample size:", stats::median(as.numeric(summary_df$n), na.rm = TRUE)),
  "",
  "## Estimands",
  "",
  "- `variance_G`, `variance_E`, and `variance_GxE` are REML variance proportions from the multi-kernel model.",
  "- `variance_Covars_fixed` is the variance share of the fitted fixed covariate predictor and is not a heritability analogue.",
  "- `variance_CovarsxG_sensitivity` and `variance_CovarsxE_sensitivity` are only reported for the sensitivity model.",
  "",
  "## QC notes",
  "",
  paste("- Proteins with `n < ", min_protein_n, "` were skipped.", sep = ""),
  paste("- Non-converged or boundary-constrained fits retained in the summary: ", sum(!summary_df$converged %in% c(TRUE, "TRUE")), ".", sep = ""),
  "- Predictive OOF comparisons are descriptive only and target a different estimand than the population variance components."
)
writeLines(report_lines, con = file.path(summary_dir, paste0("report_", model_name, ".md")))
timestamp_msg("Summary outputs written to", summary_dir)
