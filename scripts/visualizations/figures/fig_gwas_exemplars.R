#!/usr/bin/env Rscript

# ============================================================================
# fig_gwas_exemplars.R  [figure_id: fig_gwas_exemplars]
# ----------------------------------------------------------------------------
# CONSOLIDATED exposure-GWAS exemplars: the two extremes that bracket the whole
# instrument range. Merges four figures that were always read as two pairs:
#   fig_gwas_manhattan       + fig_gwas_qq         (alcohol-intake frequency)
#   fig_gwas_manhattan_wheat + fig_gwas_qq_wheat   (wheat avoidance)
#
# Supplementary Note 5 already narrates them as one bracketed pair:
#   DIFFUSE  alcohol-intake frequency -- highly polygenic, signal spread over many
#            sub-genome-wide peaks with mild residual inflation => weak,
#            distributed instruments.
#   SHARP    wheat avoidance -- one strong HLA / chromosome-6 peak on a near-null
#            background => a clean, strong instrument.
# Putting them in one figure makes the contrast visible instead of asking the
# reader to flip between four supplementary pages.
#
# Each row is one exposure: Manhattan (left) + its QQ / calibration panel (right).
#
# Input : gwas/regenie_step2/<exposure>/<exposure>.regenie via load_exposure_gwas()
# Output: figures/supplement/gwas/fig_gwas_exemplars.{pdf,png} + data tsv
#
# Authored at 6.5in = the supplement's \textwidth (scale 1.0).
# ============================================================================

local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  common <- cand[dir.exists(cand)][1]
  if (is.na(common)) stop("cannot locate common/ helpers")
  for (f in c("figure_paths", "load_heap_results", "plot_theme",
              "label_helpers", "export_helpers"))
    source(file.path(common, paste0(f, ".R")))
})
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_gwas_exemplars")
BS <- 7.5

GW_THR  <- -log10(5e-8)
SUG_THR <- -log10(1e-5)

EXEMPLARS <- list(
  list(key = "alcohol",
       exposure = "alcohol_intake_frequency_f1558_0_0",
       blurb    = "diffuse, highly polygenic"),
  list(key = "wheat",
       exposure = "never_eat_eggs_dairy_wheat_sugar_f6144_0_0.multi_Wheat_products",
       blurb    = "sharply localized (HLA)"))

completed <- list_exposure_gwas(completed_only = TRUE)
for (e in EXEMPLARS)
  if (!e$exposure %in% completed)
    stop("Exemplar exposure has no completed GWAS: ", e$exposure)

# --------------------------------------------------------------- builders ----
manhattan_panel <- function(g, ttl) {
  chr_info <- g[, .(chr_len = max(pos)), by = chr][order(chr)]
  chr_info[, offset := c(0, head(cumsum(as.numeric(chr_len)), -1))]
  chr_info[, center := offset + chr_len / 2]
  off <- setNames(chr_info$offset, chr_info$chr)

  leads <- heap_gwas_lead_variants(g, log10p_thresh = GW_THR, window = 5e5)
  d <- heap_gwas_thin(g, keep_full = 4, pos_bin = 2e5, y_bin = 0.1)
  d[, cum := pos + off[as.character(chr)]]
  d[, band := factor(chr %% 2)]
  d[, sig := LOG10P >= GW_THR]
  ymax <- max(d$LOG10P, na.rm = TRUE)

  # label the strongest loci only; a handful is enough at this size
  lab <- if (nrow(leads)) {
    leads[, cum := pos + off[as.character(chr)]]
    leads[order(-LOG10P)][seq_len(min(6L, .N))]
  } else leads[0]

  p <- ggplot(d, aes(cum, LOG10P)) +
    geom_hline(yintercept = SUG_THR, linetype = "dashed", colour = "grey75", linewidth = .25) +
    geom_hline(yintercept = GW_THR,  linetype = "dashed", colour = "#C0392B", linewidth = .3) +
    geom_point(data = d[sig == FALSE], aes(colour = band), size = .22, alpha = .6,
               show.legend = FALSE) +
    scale_colour_manual(values = c(`0` = "#9FB6CD", `1` = "#34557A")) +
    geom_point(data = d[sig == TRUE], colour = "#E07B39", size = .38, alpha = .9) +
    scale_x_continuous(breaks = chr_info$center, labels = chr_info$chr,
                       expand = expansion(mult = .01)) +
    scale_y_continuous(limits = c(0, ymax * 1.12), expand = expansion(mult = c(0, .02))) +
    labs(x = "Chromosome", y = expression(-log[10](italic(p)))) +
    theme_heap(base_size = BS) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(size = BS - 3),
          plot.margin = margin(10, 4, 2, 2))
  if (nrow(lab))
    p <- p + geom_text_repel(data = lab, aes(label = ID), size = 1.7,
                             max.overlaps = Inf, min.segment.length = 0,
                             box.padding = .3, segment.size = .18,
                             segment.colour = "grey65", colour = "grey20", seed = 7)
  # millions of variants: rasterize the point cloud, keep axes/labels/lines vector
  heap_rasterize(p)
}

qq_panel <- function(g, lambda) {
  obs <- sort(g$LOG10P, decreasing = TRUE)
  n   <- length(obs)
  ex  <- -log10((seq_len(n) - 0.5) / n)

  keep   <- obs >= 2
  idx_lo <- which(!keep)
  if (length(idx_lo) > 20000L) {
    step   <- ceiling(length(idx_lo) / 20000L)
    idx_lo <- idx_lo[seq(1L, length(idx_lo), by = step)]
  }
  sel <- sort(c(which(keep), idx_lo))
  qq <- data.table(expected = ex[sel], observed = obs[sel],
                   ci_lo = -log10(qbeta(0.975, sel, n - sel + 1)),
                   ci_hi = -log10(qbeta(0.025, sel, n - sel + 1)))
  qq <- qq[qq[, .I[1L], by = .(a = round(expected, 2), b = round(observed, 1))]$V1]
  setorder(qq, expected)

  ggplot(qq, aes(expected, observed)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "grey87", alpha = .7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "#C0392B", linewidth = .3) +
    geom_point(size = .35, alpha = .75, colour = "#34557A") +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("lambda[GC] == %.3f", lambda), parse = TRUE,
             hjust = -0.15, vjust = 1.6, size = 2.1, colour = "grey25") +
    labs(x = expression(Expected~-log[10](italic(p))),
         y = expression(Observed~-log[10](italic(p)))) +
    theme_heap(base_size = BS) +
    theme(panel.grid = element_blank(), plot.margin = margin(10, 4, 2, 2))
}

# ------------------------------------------------------------------ build ----
panels <- list(); out <- list()
for (e in EXEMPLARS) {
  g <- load_exposure_gwas(e$exposure)
  g <- g[is.finite(LOG10P) & chr %in% 1:22]
  setorder(g, chr, pos)
  lambda <- heap_gwas_lambda(g)
  n_gw   <- sum(g$LOG10P >= GW_THR)
  leads  <- heap_gwas_lead_variants(g, log10p_thresh = GW_THR, window = 5e5)
  message(sprintf("  %-8s %-28s variants=%d lambda=%.3f gw-sig=%d lead loci=%d",
                  e$key, heap_exposure_label(e$exposure), nrow(g), lambda, n_gw, nrow(leads)))

  panels[[paste0(e$key, "_man")]] <- manhattan_panel(g, e$key)
  panels[[paste0(e$key, "_qq")]]  <- qq_panel(g, lambda)
  out[[e$key]] <- data.table(exemplar = e$key,
                             exposure = e$exposure,
                             label = as.character(heap_exposure_label(e$exposure)),
                             blurb = e$blurb,
                             n_variants = nrow(g), lambda_gc = lambda,
                             n_gwsig = n_gw, n_lead = nrow(leads))
}

p <- (panels$alcohol_man | panels$alcohol_qq) /
     (panels$wheat_man   | panels$wheat_qq) +
  plot_layout(widths = c(2.05, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

heap_emit_figure(p, figure_id, data = rbindlist(out), category = "supplement",
                 subdir = "gwas", formats = c("pdf", "png"),
                 width = 6.5, height = 5.4, website = TRUE)

message("fig_gwas_exemplars: done.")
