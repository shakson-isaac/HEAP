
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
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

timestamp_msg <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " ")))
}

stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

parse_optional_args <- function(args) {
  out <- list()
  if (length(args) == 0L) {
    return(out)
  }

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }

    keyval <- substring(arg, 3L)
    if (grepl("=", keyval, fixed = TRUE)) {
      parts <- strsplit(keyval, "=", fixed = TRUE)[[1L]]
      key <- gsub("-", "_", parts[1L], fixed = TRUE)
      value <- paste(parts[-1L], collapse = "=")
      out[[key]] <- value
    } else {
      key <- gsub("-", "_", keyval, fixed = TRUE)
      out[[key]] <- TRUE
    }
  }

  out
}

get_opt <- function(opts, key, default = NULL) {
  if (key %in% names(opts)) {
    opts[[key]]
  } else {
    default
  }
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  if (is.logical(x)) {
    return(isTRUE(x))
  }
  val <- tolower(as.character(x)[1L])
  if (val %in% c("1", "true", "t", "yes", "y")) {
    return(TRUE)
  }
  if (val %in% c("0", "false", "f", "no", "n")) {
    return(FALSE)
  }
  default
}

load_project_config <- function(path, required = c(
  "project_root", "loader_rds", "genotype_pfile", "ld_pruned_pvar",
  "gcta_bin", "plink2_bin", "output_root", "covariate_specs"
)) {
  if (!file.exists(path)) {
    stopf("Config file not found: %s", path)
  }

  # parent = globalenv() (NOT baseenv): the config file references HEAP_PATHS and the
  # path helpers (igloo_path, heap_loader_rds, heap_config, load_covariate_set_*, ...)
  # which 00_paths.R and config_helpers.R source into the GLOBAL environment. With
  # parent = baseenv() the config body cannot see them -> "object 'HEAP_PATHS' not found".
  env <- new.env(parent = globalenv())
  sys.source(path, envir = env)
  if (!exists("cfg", envir = env, inherits = FALSE)) {
    stopf("Config file must define `cfg`: %s", path)
  }

  cfg <- get("cfg", envir = env, inherits = FALSE)
  missing <- required[!required %in% names(cfg)]
  if (length(missing) > 0L) {
    stopf("Missing config entries: %s", paste(missing, collapse = ", "))
  }

  cfg
}

load_config <- function(path) {
  load_project_config(path)
}

resolve_threads <- function(cfg) {
  env_threads <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
  cfg_threads <- suppressWarnings(as.integer(cfg$threads %||% 1L))
  candidates <- c(env_threads, cfg_threads, 1L)
  candidates <- candidates[is.finite(candidates) & !is.na(candidates) & candidates > 0L]
  as.integer(candidates[1L])
}

resolve_run_paths_from_root <- function(output_root, run_id, covar_spec, exposure_mode = "centered") {
  root <- file.path(output_root, run_id, covar_spec, exposure_mode)
  paths <- list(
    root = root,
    inputs = file.path(root, "inputs"),
    kernels = file.path(root, "kernels"),
    models = file.path(root, "models"),
    models_primary = file.path(root, "models", "primary"),
    models_sensitivity = file.path(root, "models", "sensitivity"),
    summary = file.path(root, "summary"),
    plots = file.path(root, "plots"),
    logs = file.path(root, "logs")
  )
  invisible(lapply(paths, ensure_dir))
  paths
}

resolve_run_paths <- function(cfg, run_id, covar_spec, exposure_mode = "centered") {
  runtime_output_root <- Sys.getenv("POPARCH_OUTPUT_ROOT", unset = "")
  output_root <- if (nzchar(runtime_output_root)) runtime_output_root else cfg$output_root
  resolve_run_paths_from_root(output_root, run_id, covar_spec, exposure_mode = exposure_mode)
}

define_pxs_class <- function() {
  if (!methods::isClass("PXSconstruct")) {
    methods::setClass(
      "PXSconstruct",
      slots = c(
        Elist = "list",
        Elist_names = "character",
        Eid_cat = "data.frame",
        ordinalIDs = "character",
        UKBprot_df = "data.frame",
        protIDs = "character",
        covars_df = "data.frame",
        covars_list = "character"
      )
    )
  }
}

as_pxs <- function(x) {
  if (is.list(x) && !isS4(x)) {
    return(x)
  }
  define_pxs_class()
  if (!isS4(x)) {
    stopf("PXS object must be an S4 PXSconstruct or a list.")
  }
  list(
    Elist = x@Elist,
    Elist_names = x@Elist_names,
    Eid_cat = x@Eid_cat,
    ordinalIDs = x@ordinalIDs,
    UKBprot_df = x@UKBprot_df,
    protIDs = x@protIDs,
    covars_df = x@covars_df,
    covars_list = x@covars_list
  )
}

read_pxs_loader <- function(path) {
  if (!file.exists(path)) {
    stopf("PXS loader RDS not found: %s", path)
  }
  define_pxs_class()
  readRDS(path)
}

normalize_loader <- function(loader) {
  if (isS4(loader)) return(as_pxs(loader))
  if (!is.null(loader$E_baseline)) {
    return(list(
      Elist       = loader$E_baseline,
      Elist_names = loader$Elist_names,
      Eid_cat     = loader$Eid_cat,
      ordinalIDs  = loader$ordinalIDs,
      UKBprot_df  = loader$prot_baseline,
      protIDs     = loader$protIDs,
      covars_df   = loader$covars_baseline,
      covars_list = loader$covars_list
    ))
  }
  loader
}

prot_clean <- function(prot_id) {
  gsub("-", "_", prot_id)
}

align_to_ids <- function(df, ids, id_col = "eid") {
  ids_df <- data.frame(ids, stringsAsFactors = FALSE)
  names(ids_df) <- id_col
  merge(ids_df, df, by = id_col, all.x = TRUE, sort = FALSE)
}

read_omicpred_map <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    stopf("OMICSPRED map not found: %s", path %||% "<missing>")
  }
  data.table::fread(path, data.table = FALSE)
}

resolve_omicpred_id <- function(protein_id, omicpred_map) {
  idx <- match(protein_id, omicpred_map$Gene)
  op_id <- omicpred_map$OMICSPRED_ID[idx]
  if (length(op_id) == 0L || is.na(op_id) || !nzchar(op_id)) {
    stopf("No OMICSPRED_ID found for protein %s.", protein_id)
  }
  op_id
}

read_sscore_or_zero <- function(fpath, out_col, eids, strict = FALSE) {
  make_zero <- function(eids_in) {
    out <- data.frame(eid = as.integer(eids_in), stringsAsFactors = FALSE)
    out[[out_col]] <- 0
    out
  }

  if (!file.exists(fpath)) {
    msg <- sprintf("Missing score file: %s", fpath)
    if (strict) {
      stopf(msg)
    }
    warning(msg, call. = FALSE)
    return(make_zero(eids))
  }

  dt <- data.table::fread(fpath, data.table = FALSE)
  id_col <- intersect(c("eid", "IID", "#IID", "id", "ID"), names(dt))[1L]
  if (is.na(id_col)) {
    id_col <- names(dt)[1L]
  }

  score_col <- intersect(c("SCORE1_AVG", "SCORE1_SUM", "SCORE1", "score", "SCORE"), names(dt))[1L]
  if (is.na(score_col)) {
    score_col <- names(dt)[ncol(dt)]
  }

  out <- data.frame(
    eid = as.integer(dt[[id_col]]),
    stringsAsFactors = FALSE
  )
  out[[out_col]] <- as.numeric(dt[[score_col]])
  out <- out[out$eid %in% eids, , drop = FALSE]
  out
}

extract_protgs <- function(protein_id, eids, omicpred_map, which = c("cis", "trans"),
                           cfg, strict = FALSE) {
  which <- match.arg(which)
  score_dir <- if (which == "cis") cfg$gs_cis_dir %||% NULL else cfg$gs_tr_dir %||% NULL
  if (is.null(score_dir) || !nzchar(score_dir)) {
    stopf("Config is missing %s for protein-specific prep.", if (which == "cis") "gs_cis_dir" else "gs_tr_dir")
  }

  op_id <- resolve_omicpred_id(protein_id, omicpred_map)
  fpath <- file.path(score_dir, paste0(op_id, ".sscore"))
  out_col <- paste0(prot_clean(protein_id), if (which == "cis") "_GScis" else "_GStrans")
  read_sscore_or_zero(fpath, out_col = out_col, eids = eids, strict = strict)
}

build_gs_struct <- function(protein_id, pxs_loader, omicpred_map, cfg, strict = FALSE) {
  pxs <- as_pxs(pxs_loader)
  protein_name <- prot_clean(protein_id)
  if (!protein_name %in% names(pxs$UKBprot_df)) {
    stopf("Protein %s is not present in the loader phenotype table.", protein_name)
  }

  omic_orig <- as.data.frame(pxs$UKBprot_df[, c("eid", protein_name), drop = FALSE])
  eids <- omic_orig$eid

  cis <- extract_protgs(protein_id, eids, omicpred_map, which = "cis", cfg = cfg, strict = strict)
  trans <- extract_protgs(protein_id, eids, omicpred_map, which = "trans", cfg = cfg, strict = strict)
  omic_pgs <- merge(cis, trans, by = "eid", all = TRUE, sort = FALSE)

  list(
    protein_name = protein_name,
    phenotype = omic_orig,
    combo = stats::na.omit(merge(omic_pgs, omic_orig, by = "eid", sort = FALSE)),
    solo = stats::na.omit(omic_pgs)
  )
}

build_protein_specific_analysis <- function(protein_id, pxs_loader, covariate_names,
                                            genotype_ids, cfg, strict_scores = FALSE,
                                            use_scores = FALSE) {
  pxs <- as_pxs(pxs_loader)
  protein_name <- prot_clean(protein_id)
  if (!protein_name %in% names(pxs$UKBprot_df)) {
    stopf("Protein %s is not present in the loader phenotype table.", protein_name)
  }

  protein_df <- as.data.frame(pxs$UKBprot_df[, c("eid", protein_name), drop = FALSE])
  protein_nonmissing_ids <- protein_df$eid[!is.na(protein_df[[protein_name]])]

  if (use_scores) {
    omicpred_map <- read_omicpred_map(cfg$omicpred_map %||% NULL)
    gs <- build_gs_struct(
      protein_id = protein_id,
      pxs_loader = pxs,
      omicpred_map = omicpred_map,
      cfg = cfg,
      strict = strict_scores
    )
    overlap_ids <- gs$combo$eid
    combo_count <- nrow(gs$combo)
  } else {
    overlap_ids <- protein_nonmissing_ids
    combo_count <- length(protein_nonmissing_ids)
  }

  covariate_names <- unique(as.character(covariate_names))
  available_covars <- names(pxs$covars_df)
  missing_covars <- setdiff(covariate_names, available_covars)
  if (length(missing_covars) > 0L) {
    stopf(
      "Protein-specific prep is missing required covariates: %s",
      paste(missing_covars, collapse = ", ")
    )
  }

  covar_df_full <- as.data.frame(pxs$covars_df)
  covar_df <- covar_df_full[, c("eid", covariate_names), drop = FALSE]
  covar_sub <- covar_df[covar_df$eid %in% overlap_ids, , drop = FALSE]
  covar_complete <- stats::complete.cases(covar_sub[, covariate_names, drop = FALSE])
  covar_ids <- covar_sub$eid[covar_complete]
  genotype_ids <- as.integer(genotype_ids)
  analysis_ids <- overlap_ids[overlap_ids %in% intersect(covar_ids, genotype_ids)]

  list(
    protein_name = protein_name,
    protein_df = protein_df,
    covar_df = covar_df,
    analysis_ids = analysis_ids,
    counts = c(
      protein_nonmissing = length(protein_nonmissing_ids),
      pre_covars_overlap = combo_count,
      after_covars = length(covar_ids),
      after_genotype = length(analysis_ids)
    )
  )
}

read_tsv <- function(path, header = TRUE, col_classes = NA) {
  if (!file.exists(path)) {
    stopf("Required file not found: %s", path)
  }
  read.table(
    path,
    header = header,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = col_classes
  )
}

write_tsv <- function(x, path, header = TRUE) {
  write.table(
    x,
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = header,
    quote = FALSE,
    na = "NA"
  )
}

read_psam <- function(pfile_prefix) {
  psam_path <- paste0(pfile_prefix, ".psam")
  psam <- read_tsv(psam_path, header = TRUE, col_classes = "character")
  names(psam)[1:2] <- c("FID", "IID")
  psam
}

read_grm_ids <- function(prefix) {
  ids <- read_tsv(paste0(prefix, ".grm.id"), header = FALSE, col_classes = "character")
  names(ids) <- c("FID", "IID")
  ids
}

run_command <- function(command, args, log_path = NULL) {
  timestamp_msg("Running:", command, paste(args, collapse = " "))
  if (is.null(log_path)) {
    status <- system2(command, args = args)
  } else {
    ensure_dir(dirname(log_path))
    status <- system2(command, args = args, stdout = log_path, stderr = log_path)
  }
  if (!identical(status, 0L)) {
    stopf("Command failed [%s]: %s", status, command)
  }
  invisible(status)
}

read_protein_arg <- function(protein_arg, default_proteins = NULL) {
  if (is.null(protein_arg) || identical(protein_arg, "")) {
    return(default_proteins)
  }
  if (file.exists(protein_arg)) {
    vals <- trimws(readLines(protein_arg, warn = FALSE))
    vals <- vals[nzchar(vals)]
    return(vals)
  }
  trimws(unlist(strsplit(protein_arg, ",", fixed = TRUE)))
}

resolve_proteins <- function(pxs_loader, protein_arg = NULL, default_proteins = NULL) {
  proteins <- read_protein_arg(protein_arg, default_proteins = default_proteins)
  available <- if (isS4(pxs_loader)) pxs_loader@protIDs else pxs_loader$protIDs
  if (is.null(proteins) || length(proteins) == 0L) {
    proteins <- available
  }
  proteins <- unique(gsub("-", "_", proteins))
  missing <- proteins[!proteins %in% available]
  if (length(missing) > 0L) {
    stopf("Requested proteins are not present in the loader: %s", paste(missing, collapse = ", "))
  }
  proteins
}

merge_exposure_list <- function(elist, cohort_ids) {
  filtered <- lapply(elist, function(df) {
    df[df$eid %in% cohort_ids, , drop = FALSE]
  })
  Reduce(function(left, right) merge(left, right, by = "eid", all = TRUE, sort = FALSE), filtered)
}

missing_cols <- function(df, miss_rate = 0.2, ignore = character()) {
  feature_names <- setdiff(names(df), ignore)
  if (length(feature_names) == 0L || !is.finite(miss_rate)) {
    return(character())
  }
  na_rate <- colMeans(is.na(df[, feature_names, drop = FALSE]))
  names(na_rate[na_rate > miss_rate])
}

drop_missing_cols <- function(df, miss_rate = 0.2, ignore = character()) {
  drop <- missing_cols(df, miss_rate = miss_rate, ignore = ignore)
  keep <- setdiff(names(df), drop)
  df[, keep, drop = FALSE]
}

make_numeric_matrix <- function(df, ordinal_ids = character(), known_factor_vars = character(),
                                impute_method = "mean",
                                center = TRUE, scale_columns = TRUE,
                                categorical_missing = "missing_level") {
  out_parts <- list()
  feature_meta <- data.frame(
    input_feature = character(),
    encoded_feature = character(),
    kind = character(),
    stringsAsFactors = FALSE
  )

  categorical_ids <- unique(c(ordinal_ids, known_factor_vars))
  if (!impute_method %in% c("mean", "median", "none")) {
    stopf("Unsupported impute_method: %s", impute_method)
  }
  if (!categorical_missing %in% c("missing_level", "error")) {
    stopf("Unsupported categorical_missing mode: %s", categorical_missing)
  }

  for (col_name in names(df)) {
    x <- df[[col_name]]
    is_categorical <- col_name %in% categorical_ids

    if (is.logical(x) && !is_categorical) {
      x <- as.integer(x)
    }

    if (!is_categorical && (is.numeric(x) || is.integer(x))) {
      x <- as.numeric(x)
      if (all(is.na(x))) {
        next
      }
      if (anyNA(x)) {
        if (identical(impute_method, "none")) {
          stopf("Numeric feature %s still has missing values under complete-case E mode.", col_name)
        }
        fill_value <- if (impute_method == "median") {
          stats::median(x, na.rm = TRUE)
        } else {
          mean(x, na.rm = TRUE)
        }
        if (!is.finite(fill_value)) {
          fill_value <- 0
        }
        x[is.na(x)] <- fill_value
      }
      out_parts[[col_name]] <- matrix(x, ncol = 1L, dimnames = list(NULL, col_name))
      feature_meta <- rbind(
        feature_meta,
        data.frame(
          input_feature = col_name,
          encoded_feature = col_name,
          kind = "numeric",
          stringsAsFactors = FALSE
        )
      )
      next
    }

    x <- as.character(x)
    missing_mask <- is.na(x) | !nzchar(x)
    if (any(missing_mask)) {
      if (identical(categorical_missing, "error")) {
        stopf("Categorical feature %s still has missing values under complete-case E mode.", col_name)
      }
      x[missing_mask] <- "Missing"
    }
    factor_x <- if (col_name %in% ordinal_ids) {
      factor(x, levels = sort(unique(x)), ordered = FALSE)
    } else {
      factor(x)
    }
    mm <- model.matrix(~ factor_x - 1L)
    colnames(mm) <- paste0(col_name, "__", make.names(levels(factor_x)))
    out_parts[[col_name]] <- mm
    feature_meta <- rbind(
      feature_meta,
      data.frame(
        input_feature = rep(col_name, ncol(mm)),
        encoded_feature = colnames(mm),
        kind = if (col_name %in% ordinal_ids) "ordinal_factor" else "factor",
        stringsAsFactors = FALSE
      )
    )
  }

  if (length(out_parts) == 0L) {
    stopf("No features were available after encoding.")
  }

  mat <- do.call(cbind, out_parts)
  raw_for_scale <- mat
  raw_centered <- sweep(raw_for_scale, 2L, colMeans(raw_for_scale), "-")
  sd_vals <- apply(raw_centered, 2L, stats::sd)
  keep <- is.finite(sd_vals) & sd_vals > 0
  mat <- mat[, keep, drop = FALSE]
  feature_meta <- feature_meta[feature_meta$encoded_feature %in% colnames(mat), , drop = FALSE]
  sd_vals <- sd_vals[keep]

  if (center) {
    mat <- sweep(mat, 2L, colMeans(mat), "-")
  }
  if (scale_columns) {
    mat <- sweep(mat, 2L, sd_vals, "/")
  }

  list(
    matrix = unname(mat),
    colnames = colnames(mat),
    feature_meta = feature_meta,
    n_features = ncol(mat)
  )
}

write_gcta_grm_from_design <- function(design, ids, prefix, block_size = 500L,
                                       n_contributors = ncol(design)) {
  if (nrow(design) != nrow(ids)) {
    stopf("Design matrix rows (%s) do not match ID rows (%s).", nrow(design), nrow(ids))
  }
  if (ncol(design) == 0L) {
    stopf("Design matrix has zero columns for prefix %s.", prefix)
  }

  ensure_dir(dirname(prefix))
  id_path <- paste0(prefix, ".grm.id")
  grm_path <- paste0(prefix, ".grm.bin")
  n_path <- paste0(prefix, ".grm.N.bin")

  write_tsv(ids[, c("FID", "IID")], id_path, header = FALSE)

  con_grm <- file(grm_path, open = "wb")
  con_n <- file(n_path, open = "wb")
  on.exit(close(con_grm), add = TRUE)
  on.exit(close(con_n), add = TRUE)

  n <- nrow(design)
  m <- ncol(design)
  n_contrib_vec <- as.numeric(n_contributors)

  for (start in seq.int(1L, n, by = block_size)) {
    end <- min(start + block_size - 1L, n)
    block <- (design[start:end, , drop = FALSE] %*% t(design[1:end, , drop = FALSE])) / m
    for (row_idx in seq_len(nrow(block))) {
      abs_idx <- start + row_idx - 1L
      vals <- as.numeric(block[row_idx, seq_len(abs_idx)])
      writeBin(as.numeric(vals), con_grm, size = 4L)
      writeBin(as.numeric(rep.int(n_contrib_vec, length(vals))), con_n, size = 4L)
    }
  }

  invisible(prefix)
}

stream_hadamard_grm <- function(prefix_a, prefix_b, prefix_out, contributor_count = 1,
                                chunk_size = 5e6L) {
  ids_a <- read_grm_ids(prefix_a)
  ids_b <- read_grm_ids(prefix_b)
  if (!identical(ids_a, ids_b)) {
    stopf("GRM IDs are not aligned: %s vs %s", prefix_a, prefix_b)
  }

  ensure_dir(dirname(prefix_out))
  write_tsv(ids_a, paste0(prefix_out, ".grm.id"), header = FALSE)

  con_a <- file(paste0(prefix_a, ".grm.bin"), open = "rb")
  con_b <- file(paste0(prefix_b, ".grm.bin"), open = "rb")
  con_a_n <- file(paste0(prefix_a, ".grm.N.bin"), open = "rb")
  con_b_n <- file(paste0(prefix_b, ".grm.N.bin"), open = "rb")
  con_out <- file(paste0(prefix_out, ".grm.bin"), open = "wb")
  con_out_n <- file(paste0(prefix_out, ".grm.N.bin"), open = "wb")

  on.exit(close(con_a), add = TRUE)
  on.exit(close(con_b), add = TRUE)
  on.exit(close(con_a_n), add = TRUE)
  on.exit(close(con_b_n), add = TRUE)
  on.exit(close(con_out), add = TRUE)
  on.exit(close(con_out_n), add = TRUE)

  repeat {
    vals_a <- readBin(con_a, what = numeric(), n = chunk_size, size = 4L)
    if (length(vals_a) == 0L) {
      break
    }
    vals_b <- readBin(con_b, what = numeric(), n = length(vals_a), size = 4L)
    n_a <- readBin(con_a_n, what = numeric(), n = length(vals_a), size = 4L)
    n_b <- readBin(con_b_n, what = numeric(), n = length(vals_a), size = 4L)

    if (length(vals_b) != length(vals_a)) {
      stopf("Mismatched GRM chunk lengths between %s and %s.", prefix_a, prefix_b)
    }

    out_vals <- as.numeric(vals_a * vals_b)
    out_n <- if (length(n_a) == length(vals_a) && length(n_b) == length(vals_a)) {
      as.numeric(pmax(1, pmin(n_a, n_b)))
    } else {
      as.numeric(rep.int(contributor_count, length(vals_a)))
    }
    writeBin(out_vals, con_out, size = 4L)
    writeBin(out_n, con_out_n, size = 4L)
  }

  invisible(prefix_out)
}

parse_hsq <- function(path) {
  lines <- readLines(path, warn = FALSE)
  blank_idx <- which(!nzchar(trimws(lines)))
  if (length(blank_idx) > 0L) {
    lines <- lines[seq_len(blank_idx[1L] - 1L)]
  }
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2L) {
    stopf("HSQ file did not contain a parsable variance table: %s", path)
  }
  read.table(
    text = lines,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

extract_hsq_value <- function(hsq_df, source_name, value_col = "Variance") {
  idx <- which(hsq_df$Source == source_name)
  if (length(idx) == 0L) {
    return(NA_real_)
  }
  as.numeric(hsq_df[idx[1L], value_col])
}

detect_convergence <- function(log_path) {
  if (!file.exists(log_path)) {
    return(FALSE)
  }
  lines <- readLines(log_path, warn = FALSE)
  has_finish <- any(grepl("Analysis finished", lines, ignore.case = TRUE))
  has_error <- any(grepl("error:|an error occurs|failed", lines, ignore.case = TRUE))
  has_constraint_stop <- any(grepl("more than half of the variance components are constrained", lines, ignore.case = TRUE))
  isTRUE(has_finish && !has_error && !has_constraint_stop)
}

collect_log_warnings <- function(log_path) {
  if (!file.exists(log_path)) {
    return(NA_character_)
  }
  lines <- readLines(log_path, warn = FALSE)
  hits <- lines[grepl("warning|error|fail", lines, ignore.case = TRUE)]
  if (length(hits) == 0L) {
    return(NA_character_)
  }
  paste(unique(trimws(hits)), collapse = " | ")
}

calc_r2 <- function(y, pred) {
  keep <- is.finite(y) & is.finite(pred)
  if (sum(keep) < 3L) {
    return(NA_real_)
  }
  y <- y[keep]
  pred <- pred[keep]
  1 - (sum((y - pred)^2) / sum((y - mean(y))^2))
}

compute_fixed_covar_share <- function(y, discrete_df, quantitative_df) {
  dat <- data.frame(y = y, stringsAsFactors = FALSE)

  if (!is.null(discrete_df) && ncol(discrete_df) > 2L) {
    disc <- discrete_df[, -c(1L, 2L), drop = FALSE]
    for (nm in names(disc)) {
      disc[[nm]] <- factor(disc[[nm]])
    }
    dat <- cbind(dat, disc)
  }

  if (!is.null(quantitative_df) && ncol(quantitative_df) > 2L) {
    qdat <- quantitative_df[, -c(1L, 2L), drop = FALSE]
    for (nm in names(qdat)) {
      qdat[[nm]] <- as.numeric(qdat[[nm]])
    }
    dat <- cbind(dat, qdat)
  }

  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  if (nrow(dat) < 3L || ncol(dat) < 2L) {
    return(NA_real_)
  }

  fit <- stats::lm(y ~ ., data = dat)
  stats::var(stats::fitted(fit)) / stats::var(dat$y)
}

compute_predictive_gxe_delta <- function(protein_id, covar_spec, predictive_root, phenotype_df) {
  if (is.null(predictive_root) || !nzchar(predictive_root)) {
    return(list(base_r2 = NA_real_, full_r2 = NA_real_, delta_r2 = NA_real_))
  }

  rds_path <- file.path(
    predictive_root,
    covar_spec,
    paste0("OOF_components_", protein_id, "_", covar_spec, ".rds")
  )
  if (!file.exists(rds_path)) {
    return(list(base_r2 = NA_real_, full_r2 = NA_real_, delta_r2 = NA_real_))
  }

  oof <- readRDS(rds_path)
  keep_cols <- c("Covars", "Gcis", "Gtrans")
  e_cols <- grep("^E_", names(oof), value = TRUE)
  gxe_cols <- grep("^GxE", names(oof), value = TRUE)
  if (length(gxe_cols) == 0L) {
    return(list(base_r2 = NA_real_, full_r2 = NA_real_, delta_r2 = NA_real_))
  }

  oof$pred_base <- rowSums(oof[, intersect(c(keep_cols, e_cols), names(oof)), drop = FALSE], na.rm = TRUE)
  oof$pred_full <- rowSums(oof[, intersect(c(keep_cols, e_cols, gxe_cols), names(oof)), drop = FALSE], na.rm = TRUE)
  merged <- merge(phenotype_df[, c("eid", protein_id), drop = FALSE], oof[, c("eid", "pred_base", "pred_full"), drop = FALSE], by = "eid")
  y <- as.numeric(merged[[protein_id]])
  base_r2 <- calc_r2(y, merged$pred_base)
  full_r2 <- calc_r2(y, merged$pred_full)
  list(base_r2 = base_r2, full_r2 = full_r2, delta_r2 = full_r2 - base_r2)
}
