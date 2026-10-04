#!/usr/bin/env Rscript

# ============================================================================
# fig_mediation_flows.R  [figure_id: fig_mediation_flows]
# ----------------------------------------------------------------------------
# CONSOLIDATED (2026-07-11). Merges fig_disease_mediation + fig_mediation_alluvial,
# which were cited by ONE sentence: "LEP, ADM, and FABP4 are a few examples of
# these hubs, linking exposure categories like physical activity to multiple
# metabolic, renal, and circulatory disease outcomes."
#
# Neither predecessor showed that. fig_mediation_alluvial showed flows but not
# breadth; fig_disease_mediation was a partitioned forest up to 28in tall that
# mixed genetic and exposomic drivers and was not about hubs at all. This makes
# the sentence's claim directly readable:
#
#   a  breadth  the most pleiotropic mediator proteins, by the NUMBER of diseases
#               each mediates, stacked by the EXPOSURE CATEGORY driving them -- so
#               "LEP links physical activity to many outcomes" is read off the bar
#   b  flows    exposure -> protein -> disease alluvial for the strongest links,
#               showing what those hubs actually connect
#
# NB distinct from fig_mediation_hubs, which breaks the same pleiotropic proteins
# down by ORGAN SYSTEM (the reporter-spectrum view feeding the main figure). This
# one is by EXPOSURE CATEGORY, which is what the citing sentence is about.
# Per-disease effect sizes for these proteins are in fig_mediation_prioritization
# panel b and are deliberately not repeated here.
#
# Input : module3 partitioned_categories via load_module3_results()
# Output: figures/supplement/module3/fig_mediation_flows.{pdf,png} + data tsv
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
  library(data.table); library(ggplot2); library(patchwork); library(ggalluvial)
})

a <- commandArgs(trailingOnly = TRUE); a <- a[!startsWith(a, "--")]
a <- a[!a %in% c("fig_mediation_flows", "all_main", "all_supplement", "all", "website")]
covarType <- if (length(a) >= 1) a[1] else "base"
family    <- if (length(a) >= 2) a[2] else "lasso"
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_mediation_flows")
BS <- 7.5
ALPHA <- 0.05
N_HUB <- 12L
N_PER_DRIVER <- 2L

sel <- c("protID", "DZ_ID", "predictor", "predictor_class", "effect_type",
         "effect_logHR", "effect_HR", "delta_p")
md <- load_module3_results(covarType = covarType, family = family,
                           mode = "partitioned_categories", select = sel)
nie <- heap_md_fdr(md[effect_type == "NIE"], "delta_p", "delta_q")
nie <- nie[is.finite(effect_HR) & is.finite(delta_q) & delta_q < ALPHA]
nie[, category := heap_md_category(predictor)]
nie <- nie[!is.na(category)]                    # modifiable exposure drivers only
nie[, catf := heap_category_pretty(category)]
nie[, absLog := abs(effect_logHR)]

pal_cat <- HEAP_ECAT_COLORS
names(pal_cat) <- heap_category_pretty(names(pal_cat))

# ---- a: hub breadth ----------------------------------------------------------
# NB bar length = the number of DISTINCT diseases the protein mediates. Stacking by
# category would double-count (one disease can be driven by several categories), so
# the bar is a single length and the fill is the protein's DOMINANT driving category.
u   <- unique(nie[, .(protID, DZ_ID, catf)])
brd <- u[, .(n_dz = uniqueN(DZ_ID)), by = protID][order(-n_dz)]
brd[, rank := .I]
dom <- nie[, .N, by = .(protID, catf)][order(protID, -N)][, .SD[1], by = protID][, .(protID, dom = catf)]
brd <- merge(brd, dom, by = "protID")
setorder(brd, rank)          # merge() re-sorts by the join key; restore breadth order

# the sentence names LEP, ADM and FABP4 as hub examples -- they are NOT the top hubs
# (LEP is 227th of 1,236 by breadth), so show the leaders AND the named three, marked,
# rather than a top-N that silently omits the proteins the text points at.
NAMED <- c("LEP", "ADM", "FABP4")
show <- unique(c(head(brd$protID, 9L), NAMED))
sb <- brd[protID %in% show][order(-n_dz)]
sb[, protID := factor(protID, levels = rev(protID))]
sb[, named := as.character(protID) %in% NAMED]

message(sprintf("hubs: %d mediator proteins | top: %s | named-in-text ranks: %s",
                nrow(brd), paste(sprintf("%s(%d)", brd$protID[1:3], brd$n_dz[1:3]), collapse = " "),
                paste(sprintf("%s=#%d", NAMED, brd[match(NAMED, protID), rank]), collapse = " ")))

pa <- ggplot(sb, aes(n_dz, protID, fill = dom)) +
  geom_col(width = .68) +
  geom_text(aes(label = sprintf("%d%s", n_dz, fifelse(named, "  *", ""))),
            hjust = -0.2, size = 1.9, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = pal_cat, name = "Dominant driver", drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, .16))) +
  labs(x = "Distinct diseases the protein mediates", y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = BS - 1, face = "italic"),
        legend.position = "right", legend.key.size = unit(6, "pt"),
        legend.text = element_text(size = BS - 2),
        legend.title = element_text(size = BS - 1),
        plot.margin = margin(10, 4, 2, 2)) +
  guides(fill = guide_legend(ncol = 1))

# ---- b: the flows -------------------------------------------------------------
fl <- copy(nie)
fl[, driver_lab := heap_md_predictor_label(predictor)]
fl[, protein := as.character(protID)]
fl[, disease := stringr::str_trunc(heap_pretty_disease(DZ_ID), 24)]
setorder(fl, -delta_p)
fl  <- unique(fl, by = c("driver_lab", "protein", "disease"), fromLast = TRUE)
setorder(fl, driver_lab, -absLog)
pdt <- fl[, head(.SD, N_PER_DRIVER), by = driver_lab]
setorder(pdt, -absLog)
pdt[, weight := 1]
pdt[, catf := factor(catf, levels = names(pal_cat))]
message(sprintf("  alluvial: %d links | %d exposures x %d proteins x %d diseases",
                nrow(pdt), uniqueN(pdt$driver_lab), uniqueN(pdt$protein),
                uniqueN(pdt$disease)))

# ggalluvial 0.12.5 geoms call a ggplot2 internal removed in ggplot2 4.0 -- their
# draw_key errors whenever a legend is built. Draw flows/strata with
# show.legend = FALSE; panel a already carries the category legend.
pb <- ggplot(pdt, aes(axis1 = driver_lab, axis2 = protein, axis3 = disease, y = weight)) +
  geom_alluvium(aes(fill = catf), alpha = .85, width = 1/9,
                curve_type = "sigmoid", show.legend = FALSE) +
  geom_stratum(width = 1/5, fill = "white", colour = "grey45", linewidth = .25,
               show.legend = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 1.45, lineheight = .9) +
  scale_x_discrete(limits = c("Exposure", "Protein", "Disease"),
                   expand = expansion(mult = c(.10, .16))) +
  scale_fill_manual(values = pal_cat, na.value = "grey80", guide = "none") +
  labs(y = "Mediation links") +
  theme_heap(base_size = BS) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = BS - 1, face = "bold"),
        plot.margin = margin(10, 4, 2, 2))

p <- (pa / pb) + plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9),
        plot.tag.position = c(0, 1))

out <- rbindlist(list(
  brd[, .(panel = "a", key = protID, value = as.numeric(n_dz))],
  pdt[, .(panel = "b", key = paste(driver_lab, protein, disease, sep = " | "),
          value = effect_HR)]), use.names = TRUE)

heap_emit_figure(p, figure_id, data = out, category = "supplement", subdir = "module3",
                 formats = c("pdf", "png"), width = 6.5, height = 7.2, website = TRUE)

message("fig_mediation_flows: done.")
