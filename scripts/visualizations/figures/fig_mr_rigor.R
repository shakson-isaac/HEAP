#!/usr/bin/env Rscript

# ============================================================================
# fig_mr_rigor.R  [figure_id: fig_mr_rigor]
# ----------------------------------------------------------------------------
# REPLACES fig_mr_supplement (2026-07-11), which had become half-redundant. That
# figure was a six-panel raster composite; three of its panels have since been
# promoted or merged elsewhere and were being PRINTED TWICE in the supplement:
#   a coloc gate            -> now fig_mr_coloc panel a
#   d MR-refines-M3 funnel  -> now fig_mr_refines_mediation panel a
#   f causal-core chains    -> now fig_mr_refines_mediation panel b
#
# What is left is the part nothing else carries -- does the MR agree with itself,
# with the second pQTL platform, and with the observational data?
#
#   a  replication   cis P->D effects, UKB Olink vs deCODE SomaScan
#   b  sensitivity   share of Tier-1 hits passing Steiger direction, MR-Egger
#                    (no directional pleiotropy) and weighted-median (outlier-robust)
#   c  concordance   MR exposure->protein effects vs the Module-2 observational
#                    estimates -- the genetic and observational reads agree in sign
#
# Panel b is also what fig_mr_sensitivity should have been: that figure was a
# 10k-point IVW-vs-Egger blob whose "sensitivity-flagged" cloud sat OFF the
# diagonal, visually contradicting the sentence citing it ("IVW and MR-Egger
# slopes in close agreement"). It is retired; this states the pass rates instead.
#
# Input : mr_edges/summary/{MRmotifs.tsv, DECODE/MRmotifs.tsv, mr_sensitivity_long.tsv}
#         module2 test-split statE (observational E->P) via load_module2_results()
# Output: figures/supplement/module5/fig_mr_rigor.{pdf,png} + data tsv
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
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel); library(scales)
})

figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_mr_rigor")
BS <- 7.5
sd <- heap_resolve_output(file.path("mr_edges", "summary"), must_exist = TRUE)

# ---- a: two-arm replication (cis P->D) --------------------------------------
keep <- c("Protein", "Disease", "beta_PDcis", "se_PDcis", "padj_PDcis")
ukb <- unique(fread(file.path(sd, "MRmotifs.tsv"), select = keep), by = c("Protein", "Disease"))
dec <- unique(fread(file.path(sd, "DECODE", "MRmotifs.tsv"), select = keep),
              by = c("Protein", "Disease"))
rp <- merge(ukb, dec, by = c("Protein", "Disease"), suffixes = c("_ukb", "_dec"))
rp <- rp[is.finite(beta_PDcis_ukb) & is.finite(beta_PDcis_dec) & padj_PDcis_ukb < 0.05]
rp[, concord := sign(beta_PDcis_ukb) == sign(beta_PDcis_dec)]
r_rep   <- rp[, cor(beta_PDcis_ukb, beta_PDcis_dec)]
rho_rep <- rp[, cor(beta_PDcis_ukb, beta_PDcis_dec, method = "spearman")]
p_rep <- round(100 * mean(rp$concord))
lab_a <- rp[abs(beta_PDcis_ukb) > 0.15][order(-abs(beta_PDcis_ukb))][seq_len(min(7, .N))]
message(sprintf("replication: n=%d | r=%.2f | %d%% sign-concordant", nrow(rp), r_rep, p_rep))

pa <- ggplot(rp, aes(beta_PDcis_ukb, beta_PDcis_dec)) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = .25) +
  geom_vline(xintercept = 0, colour = "grey88", linewidth = .25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = .3) +
  geom_point(aes(colour = concord), size = 1.2, alpha = .85) +
  geom_text_repel(data = lab_a, aes(label = Protein), size = 1.7, colour = "grey20",
                  segment.size = .18, segment.colour = "grey65", min.segment.length = 0,
                  box.padding = .35, max.overlaps = Inf, seed = 1, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.6, size = 2.0, colour = "grey25",
           label = sprintf("Spearman rho = %.2f", rho_rep)) +
  scale_colour_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "grey70"),
                      labels = c(`TRUE` = "concordant sign", `FALSE` = "discordant"),
                      name = NULL) +
  labs(x = expression("UKB Olink  "*beta*"  (cis P"%->%"D)"),
       y = expression("deCODE SomaScan  "*beta)) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2))

# ---- b: does a robust estimator agree with the primary one? ------------------
# REPLACED 2026-08-29 (author review of Supp Fig 29). The panel used to give the
# share of Tier-1 hits "passing" Steiger, the MR-Egger intercept and the weighted
# median. Two of those three are 100% BY CONSTRUCTION -- Tier 1 requires
# steiger_ok and robust_pass (build_mr_tables.R:272-281) -- so the bars restated
# the tier definition rather than testing it. The third sat at 94%, and that 6%
# is not a failure rate either: of the 219 Tier-1 hits with a flagged Egger
# intercept, ALL 219 were carried through by MR-PRESSO (164) or the weighted
# median (173). The panel was measuring how many needed rescuing while presenting
# it as a pass rate.
#
# What it shows instead is the question the check was meant to answer: does an
# estimator that is robust to a few bad instruments land in the same place as the
# primary IVW estimate? That is the same "do two independent estimates agree?"
# form as panels (a) and (c), so the figure now reads as one argument.
#
# wm_b is computed by build_mr_tables.R but not written to mr_sensitivity_long,
# so it is read straight from the per-edge robust-estimator files and cached.
# ~3.5k small reads, about 20 s cold, instant thereafter.
te   <- fread(file.path(sd, "mr_tiered_edges.tsv"))[dataset == "UKB" &
              mr_tier_final %in% c("Tier1", "Tier1plus")]
CACHE <- file.path(sd, ".tier1_robust_estimators.tsv")
if (file.exists(CACHE)) {
  rb <- fread(CACHE)
} else {
  PE <- heap_resolve_output(file.path("mr_edges", "MR_UKB_primary"), must_exist = TRUE)
  fp <- file.path(PE, te$edge_dir, te$src_id, te$tgt_id,
                  paste0(te$edge_dir, "_mr_methods.tsv"))
  rb <- rbindlist(lapply(seq_along(fp), function(i) {
    if (!file.exists(fp[i])) return(NULL)
    d <- tryCatch(fread(fp[i], showProgress = FALSE), error = function(e) NULL)
    if (is.null(d) || !"method" %in% names(d)) return(NULL)
    wm <- d[method == "Weighted median"]
    if (!nrow(wm)) return(NULL)
    data.table(edge_dir = te$edge_dir[i], src_id = te$src_id[i], tgt_id = te$tgt_id[i],
               wm_b = as.numeric(wm$b[1]), wm_se = as.numeric(wm$se[1]))
  }), fill = TRUE)
  fwrite(rb, CACHE, sep = "\t")
}

sens <- fread(file.path(sd, "mr_sensitivity_long.tsv"))[dataset == "UKB"]
ag <- merge(rb, sens[, .(edge_dir, src_id, tgt_id, b, clean, rescued_presso, rescued_median)],
            by = c("edge_dir", "src_id", "tgt_id"))
ag <- ag[is.finite(b) & is.finite(wm_b)]
ag[, status := fifelse(clean, "clean", "rescued")]
r_wm    <- cor(ag$b, ag$wm_b, use = "complete.obs")
rho_wm  <- cor(ag$b, ag$wm_b, method = "spearman", use = "complete.obs")
pct_sgn <- 100 * mean(sign(ag$b) == sign(ag$wm_b), na.rm = TRUE)
message(sprintf("estimator agreement: %s Tier-1 hits | r = %.2f | %.0f%% sign-concordant | %s rescued",
                comma(nrow(ag)), r_wm, pct_sgn, comma(ag[status == "rescued", .N])))

LIM <- range(c(ag$b, ag$wm_b), na.rm = TRUE)
pb <- ggplot(ag, aes(b, wm_b)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey55", linewidth = .3) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = .25) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = .25) +
  geom_point(aes(colour = status), size = .55, alpha = .55) +
  scale_colour_manual(values = c(clean = "#9EC9C9", rescued = "#E8A33D"), name = NULL,
                      labels = c(clean = "clean", rescued = "rescued by MR-PRESSO\nor weighted median")) +
  coord_cartesian(xlim = LIM, ylim = LIM) +
  annotate("text", x = LIM[1], y = LIM[2], hjust = 0, vjust = 1, size = 1.9, colour = "grey25",
           label = sprintf("Spearman rho = %.2f", rho_wm)) +
  labs(x = expression("Primary IVW"~beta), y = expression("Weighted-median"~beta)) +
  theme_heap(base_size = BS) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2.5, lineheight = .9),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2))

# ---- c: MR vs observational (E->P) ------------------------------------------
mm <- fread(file.path(sd, "MRmotifs.tsv"))
ep <- unique(mm[is.finite(beta_EP), .(Exposure, Protein, beta_EP, padj_EP)])[padj_EP < 0.05]
obs <- load_module2_results(covarType = "base", split = "test")$statE[
         , .(Exposure = ID, Protein = omicID, obs_beta = Estimate)]
d <- merge(ep, obs, by = c("Exposure", "Protein"))[is.finite(beta_EP) & is.finite(obs_beta)]
d[, concord := sign(beta_EP) == sign(obs_beta)]
p_con <- round(100 * mean(d$concord))
rho   <- cor(d$beta_EP, d$obs_beta, method = "spearman")
message(sprintf("E->P concordance: n=%s | %d%% sign-concordant | Spearman rho=%.2f",
                comma(nrow(d)), p_con, rho))

pc <- ggplot(d, aes(obs_beta, beta_EP)) +
  geom_hline(yintercept = 0, colour = "grey88", linewidth = .25) +
  geom_vline(xintercept = 0, colour = "grey88", linewidth = .25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = .3) +
  geom_point(aes(colour = concord), size = .7, alpha = .55, stroke = 0) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.6, size = 2.0, colour = "grey25",
           label = sprintf("Spearman rho = %.2f", rho)) +
  scale_colour_manual(values = c(`TRUE` = "#1A6B30", `FALSE` = "#D95F0E"),
                      labels = c(`TRUE` = "same direction", `FALSE` = "opposite"),
                      name = NULL) +
  labs(x = expression("Observational  "*beta*"  (HEAP exposure"%->%"protein)"),
       y = expression("MR  "*beta)) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2),
        legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(10, 4, 2, 2))

p <- ((pa | pb) / pc) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  rp[, .(panel = "a", key = paste(Protein, Disease, sep = " | "), value = beta_PDcis_dec)],
  ag[, .(panel = "b", key = paste(edge_dir, src_id, tgt_id, status, sep = " | "), value = wm_b)],
  d[,  .(panel = "c", key = paste(Exposure, Protein, sep = " | "), value = beta_EP)]),
  use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module5",
                 formats = c("pdf", "png"), width = 6.5, height = 5.6, website = TRUE)

message("fig_mr_rigor: done.")
