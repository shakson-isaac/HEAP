#!/usr/bin/env Rscript

# ============================================================================
# HEAP workflow status dashboard
# ============================================================================
# One command to answer: what's running, what's done, what still needs to run.
#
#   module load gcc/14.2.0 R/4.4.2
#   Rscript workflow/heap_status.R                 # CLI table + writes STATUS.md + status.html
#   Rscript workflow/heap_status.R --no-write      # CLI only
#   Rscript workflow/heap_status.R --module module2  # restrict to one module
#
# Combines THREE signals per experiment/run:
#   1. expected   — how many array tasks the run needs (from manifests / list files)
#   2. completed  — how many per-task output marker files already exist on disk
#   3. slurm      — how many matching jobs are currently R(unning)/PD(pending) in squeue
#
# Status is derived from those:
#   done       completed >= expected
#   running    matching jobs in the queue (R or PD)
#   partial    0 < completed < expected, nothing queued  (look here for problems)
#   empty      completed == 0, nothing queued            (not started yet)
#
# Reads the same sources of truth the jobs do (config/io_map.yml, the experiment
# YAMLs, the manifests, the canonical IGLOO outputs), so it never drifts from the
# real workflow. Degrades gracefully if squeue is absent or a run has not started.
# ============================================================================

suppressWarnings(suppressMessages({
  # Locate and source the path config (defines heap_project_output(), etc.)
  .self <- tryCatch(normalizePath(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
            error = function(e) NA_character_)
  .wfdir <- if (!is.na(.self)) dirname(.self) else file.path(getwd(), "workflow")
  source(file.path(.wfdir, "00_paths.R"))
  source(file.path(.wfdir, "config_helpers.R"))
}))

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ---------------------------------------------------------------------------
# CLI args
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
opt_no_write <- "--no-write" %in% args
opt_module   <- { i <- which(args == "--module"); if (length(i)) args[i + 1] else NULL }

DOCS_DIR <- heap_project_docs("status")
dir.create(DOCS_DIR, recursive = TRUE, showWarnings = FALSE)
STATUS_MD   <- file.path(DOCS_DIR, "STATUS.md")
STATUS_HTML <- file.path(DOCS_DIR, "status.html")

# ---------------------------------------------------------------------------
# Helpers: count expected tasks and completed markers
# ---------------------------------------------------------------------------

# Data-row count of a manifest TSV (excludes header).
manifest_nrows <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  n <- length(readLines(path, warn = FALSE))
  max(0L, n - 1L)
}

# Count files under `root` (recursively) whose basename matches `marker` regex.
# marker = NULL counts every file (used for modules without a known marker).
count_markers <- function(root, marker = NULL) {
  if (!dir.exists(root)) return(0L)
  ff <- list.files(root, recursive = TRUE, full.names = FALSE, all.files = FALSE)
  ff <- ff[!grepl("/$", ff)]
  if (!is.null(marker)) ff <- ff[grepl(marker, basename(ff))]
  length(ff)
}

# Number of non-blank lines in a list file (expected array size for list modules).
list_len <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  sum(nzchar(trimws(readLines(path, warn = FALSE))))
}

# ---------------------------------------------------------------------------
# squeue: one snapshot, parsed into a data.frame
# ---------------------------------------------------------------------------
# HEAP is a shared, group-writable workflow — teammates routinely launch
# modules (e.g. module6) against the same outputs — so we scope squeue to the
# whole `patel` account by default. queue_counts() then attributes jobs by name,
# so unrelated account jobs are ignored. Override with e.g.
# HEAP_STATUS_SQUEUE_ARGS="--me" (only mine) or "-u user1,user2".
read_squeue <- function() {
  empty <- data.frame(jobid = character(), name = character(),
                      state = character(), user = character(),
                      time = character(), stringsAsFactors = FALSE)
  if (nzchar(Sys.which("squeue")) == FALSE) {
    attr(empty, "available") <- FALSE
    return(empty)
  }
  extra <- Sys.getenv("HEAP_STATUS_SQUEUE_ARGS", unset = "--me")
  cmd <- paste("squeue", extra, '-h -o "%i|%j|%t|%u|%M"')
  out <- tryCatch(system(cmd, intern = TRUE, ignore.stderr = TRUE),
                  error = function(e) character(0),
                  warning = function(w) character(0))
  if (length(out) == 0) { attr(empty, "available") <- TRUE; return(empty) }
  parts <- strsplit(out, "|", fixed = TRUE)
  parts <- parts[lengths(parts) >= 5]
  if (length(parts) == 0) { attr(empty, "available") <- TRUE; return(empty) }
  df <- data.frame(
    jobid = vapply(parts, `[`, "", 1),
    name  = vapply(parts, `[`, "", 2),
    state = vapply(parts, `[`, "", 3),
    user  = vapply(parts, `[`, "", 4),
    time  = vapply(parts, `[`, "", 5),
    stringsAsFactors = FALSE
  )
  attr(df, "available") <- TRUE
  df
}

SQ <- read_squeue()
SQ_AVAILABLE <- isTRUE(attr(SQ, "available"))
SQ_SCOPE <- Sys.getenv("HEAP_STATUS_SQUEUE_ARGS", unset = "--me")

# Running/pending job counts for queue rows whose NAME matches `pattern` (regex).
queue_counts <- function(pattern) {
  if (nrow(SQ) == 0) return(c(R = 0L, PD = 0L))
  hit <- grepl(pattern, SQ$name, ignore.case = TRUE)
  c(R  = sum(hit & SQ$state == "R"),
    PD = sum(hit & SQ$state %in% c("PD", "CF")))
}

# Distinct owners of active (R/PD/CF) jobs matching `pattern`, comma-separated.
# Used to show when a teammate — not you — is running a module.
queue_who <- function(pattern) {
  if (nrow(SQ) == 0) return("")
  hit <- grepl(pattern, SQ$name, ignore.case = TRUE) &
         SQ$state %in% c("R", "PD", "CF")
  paste(sort(unique(SQ$user[hit])), collapse = ",")
}

# Owner suffix for the queue column: names anyone *other than* the current user
# who has matching jobs, so "125R/1PD (dil402)" reads at a glance. Empty when
# only you (or nobody) are running it.
ME <- Sys.getenv("USER")
own_suffix <- function(who) {
  if (!nzchar(who)) return("")
  others <- setdiff(strsplit(who, ",", fixed = TRUE)[[1]], ME)
  if (length(others) == 0) "" else paste0(" (", paste(others, collapse = ","), ")")
}

# ---------------------------------------------------------------------------
# Per-module enumeration of tracked runs
# Each row: module, run, expected, done, q_run, q_pend, status, note
# ---------------------------------------------------------------------------

rows <- list()
add_row <- function(...) rows[[length(rows) + 1]] <<- list(...)

derive_status <- function(expected, done, qR, qPD) {
  if (!is.na(expected) && expected > 0 && done >= expected) return("done")
  if (qR > 0 || qPD > 0) return("running")
  if (done > 0)          return("partial")
  "empty"
}

# ---- Manifest modules (1, 2, 3, 5): one run per manifest TSV present ---------
manifest_module <- function(module, marker, jobprefix,
                            default_subdir, output_kind = c("output", "manifest_subdir")) {
  mdir <- heap_manifest(module)
  tsvs <- if (dir.exists(mdir)) list.files(mdir, pattern = "\\.tsv$", full.names = TRUE) else character(0)
  if (length(tsvs) == 0) {
    add_row(module = module, run = "(no manifests)", expected = NA_integer_,
            done = 0L, q_run = 0L, q_pend = 0L, q_who = "", status = "empty",
            note = "generate a manifest to start")
    return(invisible())
  }
  for (tsv in sort(tsvs)) {
    exp <- sub("\\.tsv$", "", basename(tsv))
    expected <- manifest_nrows(tsv)
    # output_subdir from experiment YAML if available, else module default
    subdir <- tryCatch(load_experiment_config(module, exp)$output_subdir,
                       error = function(e) NULL) %||% default_subdir
    root <- heap_project_output(subdir, exp)
    done <- count_markers(root, marker)
    pat <- paste0(jobprefix, ".*", exp)
    qc <- queue_counts(pat)
    # also try a looser match on the experiment name alone (job-name doubling)
    if (qc["R"] + qc["PD"] == 0) { pat <- exp; qc <- queue_counts(pat) }
    add_row(module = module, run = exp, expected = expected, done = done,
            q_run = unname(qc["R"]), q_pend = unname(qc["PD"]), q_who = queue_who(pat),
            status = derive_status(expected, done, qc["R"], qc["PD"]),
            note = "")
  }
}

run_module1 <- function()
  manifest_module("module1", marker = "^predictive_r2_coarse_[0-9]+\\.txt$",
                  jobprefix = "M1", default_subdir = "module1_predictive_r2_score_partition")

run_module2 <- function()
  manifest_module("module2", marker = "^univar_assoc_[0-9]+\\.rds$",
                  jobprefix = "M2", default_subdir = "module2")

run_module3 <- function()
  manifest_module("module3", marker = "^MDres_[0-9]+\\.txt$",
                  jobprefix = "M3", default_subdir = "module3")

run_module5 <- function()
  # One `chunk_<id>.log` is written per array task (the manifest unit) under
  # logs/<edge_type>/, so it tracks task completion against the manifest row count.
  # Each chunk ALSO writes ~9 per-edge result files (E_to_P_*.tsv); the old
  # marker=NULL counted all ~110k of those -> bogus >1000%. (The log is created at
  # task start, so a handful of in-flight tasks count as done — fine for progress.)
  manifest_module("module5", marker = "^chunk_[0-9]+\\.log$",
                  jobprefix = "MR", default_subdir = "mr_edges")

# ---- Module 6 (list): prod stage, per covariate set ------------------------
run_module6 <- function() {
  n_exp <- list_len(heap_config("exposure_sets", "module6_exposures.txt")) %||% 172L
  if (is.na(n_exp)) n_exp <- 172L
  # The M6long_array job name does NOT encode the covariate type, so squeue alone
  # can't tell a `base` run from a `base_clinical` one. Attribute the queue to the
  # covar(s) with on-disk activity (output subdir created by a running task); if
  # neither has started, fall back to the primary `base`. Otherwise the same jobs
  # would wrongly show on both rows. (Durable fix: encode covar in the job name.)
  covars <- c("base", "base_clinical")
  qc  <- queue_counts("M6long")
  who <- queue_who("M6long")
  active <- covars[vapply(covars, function(cv)
    dir.exists(heap_project_output("module6_pes_longitudinal", cv)), logical(1))]
  if (length(active) == 0) active <- "base"
  for (covar in covars) {
    root <- heap_project_output("module6_pes_longitudinal", covar)
    marker <- paste0("^PESlong_", covar, "_.*_FinalModelArtifact\\.rds$")
    done <- count_markers(root, marker)
    on  <- covar %in% active
    rR  <- if (on) unname(qc["R"])  else 0L
    rPD <- if (on) unname(qc["PD"]) else 0L
    add_row(module = "module6", run = paste0("prod_", covar),
            expected = n_exp, done = done,
            q_run = rR, q_pend = rPD, q_who = if (on) who else "",
            status = derive_status(n_exp, done, rR, rPD),
            note = if (covar == "base") "primary" else "sensitivity")
  }
}

# ---- Population architecture (list): GREML per SPEC ------------------------
run_poparch <- function() {
  pf <- heap_script("population_architecture", "config", "protein_sets",
                    "all_proteins_from_loader.txt")
  n_prot <- list_len(pf) %||% NA_integer_
  for (spec in c("base", "base_clinical")) {
    root <- heap_project_output("population_architecture", spec)
    done <- count_markers(root, "_summary\\.tsv$")
    qc <- queue_counts("poparch")
    is_primary <- spec == "base"
    # Only show base_clinical row if it has data or is queued (it's optional sens).
    if (!is_primary && done == 0 && qc["R"] + qc["PD"] == 0) next
    add_row(module = "population_architecture", run = spec,
            expected = n_prot, done = done,
            q_run = unname(qc["R"]), q_pend = unname(qc["PD"]), q_who = queue_who("poparch"),
            status = derive_status(n_prot, done, qc["R"], qc["PD"]),
            note = if (is_primary) "primary" else "sensitivity")
  }
}

# ---- Exposure GWAS (list): continuous + binary ----------------------------
gwas_done_count <- function(evars_path, regenie_root) {
  exps <- if (file.exists(evars_path))
    trimws(readLines(evars_path, warn = FALSE)) else character(0)
  exps <- exps[nzchar(exps)]
  if (length(exps) == 0 || !dir.exists(regenie_root)) return(0L)
  sum(vapply(exps, function(e) {
    d <- file.path(regenie_root, e)
    # accept simplified name (<exp>.regenie) or legacy doubled name
    file.exists(file.path(d, paste0(e, ".regenie"))) ||
      file.exists(file.path(d, paste0("regenie_step2_", e, "_", e, ".regenie")))
  }, logical(1)))
}

run_gwas <- function() {
  root <- heap_gwas("regenie_step2")
  specs <- list(
    continuous = heap_path("slurm", "gwas_regenie", "evars_continuous_heap.txt"),
    binary     = heap_path("slurm", "gwas_regenie", "evars_binary_heap.txt")
  )
  for (kind in names(specs)) {
    evars <- specs[[kind]]
    expected <- list_len(evars)
    done <- gwas_done_count(evars, root)
    pat <- if (kind == "continuous") "gwas.*cont|regenie.*cont|exposures_continuous"
           else "gwas.*bin|regenie.*bin|exposures_binary"
    qc <- queue_counts(pat)
    if (qc["R"] + qc["PD"] == 0) { pat <- "gwas|regenie"; qc <- queue_counts(pat) }
    add_row(module = "gwas_regenie", run = kind,
            expected = expected, done = done,
            q_run = unname(qc["R"]), q_pend = unname(qc["PD"]), q_who = queue_who(pat),
            status = derive_status(expected, done, qc["R"], qc["PD"]),
            note = if (is.na(expected)) "run prepare_gwas_exposures.R" else "")
  }
}

# ---------------------------------------------------------------------------
# Foundation: HEAP.rds present?
# ---------------------------------------------------------------------------
run_foundation <- function() {
  rds <- heap_loader_rds
  ok <- file.exists(rds)
  fpat <- "HEAP_loader|run_HEAP"
  fq <- queue_counts(fpat)
  add_row(module = "foundation", run = "HEAP.rds",
          expected = 1L, done = as.integer(ok),
          q_run = unname(fq["R"]), q_pend = unname(fq["PD"]), q_who = queue_who(fpat),
          status = if (ok) "done" else if (sum(fq) > 0) "running" else "empty",
          note = if (ok) "everything depends on this" else "BUILD FIRST")
}

# ---------------------------------------------------------------------------
# Run the enumerators
# ---------------------------------------------------------------------------
all_runners <- list(
  foundation = run_foundation,
  module1 = run_module1, module2 = run_module2, module3 = run_module3,
  module5 = run_module5, module6 = run_module6,
  population_architecture = run_poparch, gwas_regenie = run_gwas
)
sel <- if (is.null(opt_module)) names(all_runners) else intersect(opt_module, names(all_runners))
for (m in sel) all_runners[[m]]()

DF <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
DF$pct <- ifelse(is.na(DF$expected) | DF$expected == 0, NA_real_,
                 round(100 * DF$done / DF$expected, 1))

# ---------------------------------------------------------------------------
# Render: CLI
# ---------------------------------------------------------------------------
status_glyph <- c(done = "[done]", running = "[run ]",
                  partial = "[part]", empty = "[ -- ]")
# ANSI SGR codes per status; applied only on an interactive terminal so the
# escapes never pollute pipes, log files, or the squeue-less batch case.
status_color <- c(done = "32", running = "33", partial = "1;31", empty = "90")
USE_COLOR <- isatty(stdout()) && Sys.getenv("NO_COLOR", "") == "" &&
             !identical(Sys.getenv("TERM"), "dumb")
ansi <- function(s, status)
  if (USE_COLOR) paste0("\033[", status_color[[status]], "m", s, "\033[0m") else s

# left-justify into width w, truncating with a trailing '~' if too long
pad <- function(x, w) {
  s <- as.character(x)
  s <- ifelse(nchar(s) > w, paste0(substr(s, 1, w - 1), "~"), s)
  formatC(s, width = -w, flag = " ")
}

# Fixed-width progress bar. Unicode blocks when the locale supports them,
# ASCII otherwise, so it stays readable on any terminal.
UTF8 <- grepl("UTF-8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)
bar_fill  <- if (UTF8) "█" else "#"
bar_empty <- if (UTF8) "░" else "."
progress_bar <- function(pct, width = 10) {
  if (is.na(pct)) return(strrep(" ", width))
  fill <- max(0L, min(width, as.integer(round(pct / 100 * width))))
  paste0(strrep(bar_fill, fill), strrep(bar_empty, width - fill))
}

# Summary counts — computed once, reused by the CLI header and the MD/HTML render.
n_done  <- sum(DF$status == "done");    n_run   <- sum(DF$status == "running")
n_part  <- sum(DF$status == "partial"); n_empty <- sum(DF$status == "empty")
attn    <- DF[DF$status == "partial", , drop = FALSE]   # problems, surfaced first

RULE <- strrep("=", 92)
cat("\n", RULE, "\n", sep = "")
cat("HEAP WORKFLOW STATUS   ", format(Sys.time(), "%Y-%m-%d %H:%M"),
    "   user=", Sys.getenv("USER"), "\n", sep = "")
cat("  ", ansi(paste0(n_done, " done"), "done"), " | ",
    ansi(paste0(n_run, " running"), "running"), " | ",
    ansi(paste0(n_part, " partial"), "partial"), " | ",
    ansi(paste0(n_empty, " not-started"), "empty"),
    "   (", nrow(DF), " tracked runs)\n", sep = "")
if (SQ_AVAILABLE)
  cat("  queue scope: ", SQ_SCOPE,
      "  (HEAP_STATUS_SQUEUE_ARGS to change; '--me' = only yours)\n", sep = "")
if (!SQ_AVAILABLE) cat("  (squeue unavailable — showing file-based progress only)\n")

# Needs-attention block: partial runs first, where the problems are.
if (nrow(attn) > 0) {
  cat(RULE, "\n", sep = "")
  cat(ansi("  ! NEEDS ATTENTION", "partial"),
      " — output present but nothing queued (likely failed/interrupted):\n", sep = "")
  for (i in seq_len(nrow(attn))) {
    r <- attn[i, ]
    exp_str <- if (is.na(r$expected)) "?" else as.character(r$expected)
    cat("      ", pad(paste0(r$module, "/", r$run), 40),
        r$done, "/", exp_str, "  — inspect logs\n", sep = "")
  }
}

cat(RULE, "\n", sep = "")
cat(pad("STATUS", 8), pad("MODULE", 22), pad("RUN", 26),
    pad("DONE/EXP", 11), "PROGRESS", strrep(" ", 9), "QUEUE\n", sep = "")
cat(strrep("-", 92), "\n", sep = "")

cur_mod <- ""
for (i in seq_len(nrow(DF))) {
  r <- DF[i, ]
  mod <- if (r$module == cur_mod) "" else r$module
  cur_mod <- r$module
  exp_str <- if (is.na(r$expected)) "?" else as.character(r$expected)
  pct_str <- if (is.na(r$pct)) "    -" else formatC(paste0(r$pct, "%"), width = 5)
  q_str <- if (r$q_run + r$q_pend > 0)
    paste0(r$q_run, "R/", r$q_pend, "PD", own_suffix(r$q_who)) else ""
  cat(ansi(pad(status_glyph[[r$status]], 8), r$status),
      pad(mod, 22), pad(r$run, 26),
      pad(paste0(r$done, "/", exp_str), 11),
      ansi(progress_bar(r$pct), r$status), " ", pct_str, "  ", q_str, "\n", sep = "")
}
cat(strrep("-", 92), "\n", sep = "")
cat("Dependency gates: module3 needs module1 (matching covar+family); module5 needs gwas_regenie.\n")
cat(RULE, "\n\n", sep = "")

# ---------------------------------------------------------------------------
# Render: Markdown + HTML
# ---------------------------------------------------------------------------
if (!opt_no_write) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  badge <- c(done = "✅ done", running = "🟡 running",
             partial = "🟠 partial", empty = "⬜ not started")

  md <- c(
    "# HEAP workflow status",
    "",
    paste0("_Generated ", ts, " by `workflow/heap_status.R` — user `",
           Sys.getenv("USER"), "`._"),
    if (!SQ_AVAILABLE) "\n> squeue unavailable; progress shown is file-based only.\n" else "",
    "",
    paste0("**", n_done, "** done · **", n_run, "** running · **", n_part,
           "** partial · **", n_empty, "** not started · ", nrow(DF), " tracked runs"),
    "",
    if (nrow(attn) > 0) paste0(
      "> ⚠️ **Needs attention** — output present but nothing queued ",
      "(likely failed/interrupted): ",
      paste(sprintf("`%s/%s` (%d/%s)", attn$module, attn$run, attn$done,
                    ifelse(is.na(attn$expected), "?", attn$expected)),
            collapse = ", "),
      ". Inspect logs.") else "",
    "",
    "| Status | Module | Run | Done / Expected | % | Queue | Note |",
    "|---|---|---|---|---|---|---|"
  )
  for (i in seq_len(nrow(DF))) {
    r <- DF[i, ]
    exp_str <- if (is.na(r$expected)) "?" else as.character(r$expected)
    pct_str <- if (is.na(r$pct)) "-" else paste0(r$pct, "%")
    q_str <- if (r$q_run + r$q_pend > 0)
      paste0(r$q_run, "R / ", r$q_pend, "PD", own_suffix(r$q_who)) else ""
    md <- c(md, paste0("| ", badge[[r$status]], " | ", r$module, " | `", r$run,
                       "` | ", r$done, " / ", exp_str, " | ", pct_str, " | ",
                       q_str, " | ", r$note, " |"))
  }
  md <- c(md, "",
    "**Dependency gates:** `module3` needs the matching `module1` run (same covariate set + family); `module5` needs `gwas_regenie`.",
    "",
    "Regenerate with `Rscript workflow/heap_status.R`. See `docs/REPRODUCIBILITY.md` for how to run each module.")
  writeLines(md, STATUS_MD)

  # HTML
  color <- c(done = "#1a7f37", running = "#9a6700",
             partial = "#bc4c00", empty = "#6e7781")
  bg    <- c(done = "#dafbe1", running = "#fff8c5",
             partial = "#ffec99", empty = "#f6f8fa")
  htmlrows <- character(0)
  cur <- ""
  for (i in seq_len(nrow(DF))) {
    r <- DF[i, ]
    new_group <- r$module != cur          # show module name once per group
    cur <- r$module
    mod_cell  <- if (new_group) r$module else ""
    grp_border <- if (new_group && i > 1) "border-top:2px solid #d0d7de;" else ""
    exp_str <- if (is.na(r$expected)) "?" else as.character(r$expected)
    pct_val <- if (is.na(r$pct)) 0 else r$pct
    pct_str <- if (is.na(r$pct)) "-" else paste0(r$pct, "%")
    osfx <- own_suffix(r$q_who)
    q_str <- if (r$q_run + r$q_pend > 0) {
      paste0("<b>", r$q_run, "</b>R / <b>", r$q_pend, "</b>PD",
             if (nzchar(osfx)) paste0(" <span class='sub'>", osfx, "</span>") else "")
    } else "&mdash;"
    bar <- paste0(
      "<div style='background:#eaeef2;border-radius:4px;height:14px;width:120px;",
      "display:inline-block;vertical-align:middle;overflow:hidden'>",
      "<div style='background:", color[[r$status]], ";height:14px;width:",
      round(pct_val), "%'></div></div>")
    htmlrows <- c(htmlrows, paste0(
      "<tr style='background:", bg[[r$status]], ";", grp_border,
      "border-left:4px solid ", color[[r$status]], "'>",
      "<td style='color:", color[[r$status]], ";font-weight:600'>", r$status, "</td>",
      "<td style='font-weight:600'>", mod_cell, "</td><td><code>", r$run, "</code></td>",
      "<td style='text-align:right;font-variant-numeric:tabular-nums'>", r$done, " / ", exp_str, "</td>",
      "<td>", bar, " <span class='sub'>", pct_str, "</span></td>",
      "<td style='text-align:center'>", q_str, "</td>",
      "<td style='color:#57606a'>", r$note, "</td></tr>"))
  }
  html <- c(
    "<!doctype html><html><head><meta charset='utf-8'>",
    "<meta http-equiv='refresh' content='120'>",
    "<title>HEAP status</title>",
    "<style>body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;",
    "margin:2rem;color:#1f2328}h1{margin-bottom:.2rem}",
    "table{border-collapse:collapse;width:100%;font-size:14px}",
    "th,td{padding:6px 10px;border-bottom:1px solid #d0d7de;text-align:left}",
    "th{background:#f6f8fa;position:sticky;top:0}code{background:#eff1f3;",
    "padding:1px 5px;border-radius:4px;font-size:13px}.sub{color:#57606a;font-size:13px}",
    ".pill{display:inline-block;padding:2px 8px;border-radius:10px;font-size:13px;margin-right:6px}</style>",
    "</head><body>",
    "<h1>HEAP workflow status</h1>",
    paste0("<div class='sub'>Generated ", ts, " &middot; user ", Sys.getenv("USER"),
           " &middot; auto-refresh 120s",
           if (SQ_AVAILABLE) paste0(" &middot; queue scope <code>", SQ_SCOPE, "</code>") else "",
           if (!SQ_AVAILABLE) " &middot; <b>squeue unavailable (file-based only)</b>" else "",
           "</div><p>",
           "<span class='pill' style='background:#dafbe1;color:#1a7f37'>", n_done, " done</span>",
           "<span class='pill' style='background:#fff8c5;color:#9a6700'>", n_run, " running</span>",
           "<span class='pill' style='background:#ffec99;color:#bc4c00'>", n_part, " partial</span>",
           "<span class='pill' style='background:#f6f8fa;color:#6e7781'>", n_empty, " not started</span>",
           "</p>"),
    if (nrow(attn) > 0) paste0(
      "<div style='background:#fff1e5;border:1px solid #ffb366;",
      "border-left:4px solid #bc4c00;border-radius:6px;padding:10px 14px;margin:0 0 14px'>",
      "<b style='color:#bc4c00'>&#9888; Needs attention</b> &mdash; output present but ",
      "nothing queued (likely failed/interrupted): ",
      paste(sprintf("<code>%s/%s</code> (%d/%s)", attn$module, attn$run, attn$done,
                    ifelse(is.na(attn$expected), "?", attn$expected)), collapse = ", "),
      ". Inspect logs.</div>") else "",
    "<table><thead><tr><th>Status</th><th>Module</th><th>Run</th>",
    "<th>Done/Expected</th><th>Progress</th><th>Queue</th><th>Note</th></tr></thead><tbody>",
    htmlrows,
    "</tbody></table>",
    "<p class='sub'>Dependency gates: <code>module3</code> needs the matching <code>module1</code> run (same covariate set + family); <code>module5</code> needs <code>gwas_regenie</code>. ",
    "Regenerate: <code>Rscript workflow/heap_status.R</code>. Run instructions: <code>docs/REPRODUCIBILITY.md</code>.</p>",
    "</body></html>")
  writeLines(html, STATUS_HTML)

  cat("Wrote:\n  ", STATUS_MD, "\n  ", STATUS_HTML, "\n\n", sep = "")
}
