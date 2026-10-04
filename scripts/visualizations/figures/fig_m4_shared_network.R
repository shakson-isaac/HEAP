#!/usr/bin/env Rscript
# ============================================================================
# fig_m4_shared_network.R  [figure_id: fig_m4_shared_network]  -- Fig5 panel d
# ----------------------------------------------------------------------------
# The plasma proteome as a SHARED LANGUAGE, split into causal intermediates and
# disease reporters. Tripartite network read left->right:
#   (1) lifestyle exposures (observational) + trials (HERITAGE, GLP1-RA)
#   (2) the shared plasma proteins they all move, in two blocks:
#         CAUSAL INTERMEDIATES  -- genetically causal for disease (forward MR)
#         DISEASE REPORTERS     -- downstream markers of disease (reverse MR)
#   (3) cardiometabolic disease (central hub: causes feed in, markers branch out)
# Uniform signed encoding: edge COLOUR = direction (red raises / blue lowers the
# node it points to); genetic LINETYPE = evidence type (solid colocalized,
# dashed MR, dotted reporter-reverse). Perturbation edges: thin=observational,
# bold=trial. This directly draws the paper's minority-causal / majority-reporter
# thesis instead of hiding it in a node colour.
#
# Pure drawer: curated cast + all statistics come from the support tables written
# by scripts/support/build_shared_language_network.R.
#
# Renders two ways (per MULTIPANEL_FIGURE_GUIDE):
#   standalone (default) -- title + reading subtitle, for solo review
#   CELL (HEAP_CELL=1)   -- panel size for the Fig5 composite (no big title)
#
# Input : support/intervention_compare/shared_language_network_{nodes,edges}.tsv
# Output: standalone -> figures/exploratory/module2/fig_m4_shared_network.{pdf,png}
#         CELL       -> figures/exploratory/module4/fig_m4_shared_network_cell.{png,pdf}
# ============================================================================
local({
  cand <- c(file.path(getwd(), "scripts", "visualizations", "common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cand[dir.exists(cand)][1]
  for (f in c("figure_paths.R", "load_heap_results.R", "plot_theme.R",
              "label_helpers.R", "export_helpers.R")) source(file.path(cm, f))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

CELL      <- nzchar(Sys.getenv("HEAP_CELL"))
figure_id <- Sys.getenv("HEAP_FIGURE_ID", unset = "fig_m4_shared_network")
FIGDIR    <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module4")
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
indir     <- heap_project_output("support", "intervention_compare")

nodes <- fread(file.path(indir, "shared_language_network_nodes.tsv"))
edges <- fread(file.path(indir, "shared_language_network_edges.tsv"))

## ---- layout -----------------------------------------------------------------
srcord <- c("Strenuous sports", "Processed meat", "Usual walking pace",
            "Current smoking", "HERITAGE", "GLP1-RA")
src <- nodes[kind %in% c("exp_obs", "exp_rct")]
src[, ord := match(id, srcord)]; setorder(src, ord)
src[, sy := seq(2.80, -0.90, length.out = .N)][, sx := -3.9]
src[, fl := fifelse(kind == "exp_obs", "exposure", "RCT")]

dord <- c("Hypertension", "Lipid disorder", "Type-2 diabetes", "Obesity")
dz <- nodes[kind == "disease"]; dz[, ord := match(label, dord)]; setorder(dz, ord)
dz[, dy := seq(1.90, -0.30, length.out = .N)][, dx := 3.2]

efw <- edges[etype == "gen_fwd"]                       # from = protein, to = disease
erv <- edges[etype == "gen_rev"]                       # from = disease, to = protein
pr  <- nodes[kind == "protein"]

# causal block (top): order by forward disease then breadth
prim_fwd <- efw[, .(dzl = to[which.max(weight)]), by = .(protein = from)]
prc <- merge(pr[class == "causal"], prim_fwd, by.x = "id", by.y = "protein", all.x = TRUE)
prc <- merge(prc, dz[, .(dzl = label, ddy = dy)], by = "dzl", all.x = TRUE)
setorder(prc, -ddy, -breadth); prc[, py := seq(2.90, 1.81, length.out = .N)][, px := 0]
# reporter block (bottom): order by reverse disease then breadth
prim_rev <- erv[, .(dzl = from[1]), by = .(protein = to)]
prr <- merge(pr[class == "reporter"], prim_rev, by.x = "id", by.y = "protein", all.x = TRUE)
prr <- merge(prr, dz[, .(dzl = label, rdy = dy)], by = "dzl", all.x = TRUE)
setorder(prr, -rdy, -breadth); prr[, py := seq(1.25, -1.15, length.out = .N)][, px := 0]
prall <- rbind(prc[, .(id, label, breadth, px, py, R2_E)], prr[, .(id, label, breadth, px, py, R2_E)])

## ---- edges with coordinates -------------------------------------------------
eo <- merge(edges[etype == "obs"],    src[, .(id, sx, sy)], by.x = "from", by.y = "id")
eo <- merge(eo, prall[, .(id, px, py)], by.x = "to", by.y = "id")
ei <- merge(edges[etype == "interv"], src[, .(id, sx, sy)], by.x = "from", by.y = "id")
ei <- merge(ei, prall[, .(id, px, py)], by.x = "to", by.y = "id")
efw <- merge(efw, prc[, .(id, px, py)], by.x = "from", by.y = "id")
efw <- merge(efw, dz[, .(id = label, dx, dy)], by.x = "to", by.y = "id")
erv <- merge(erv, dz[, .(id = label, dx, dy)], by.x = "from", by.y = "id")
erv <- merge(erv, prr[, .(id, px, py)], by.x = "to", by.y = "id")
# reverse edges ON the causal proteins (they too are broad obesity/T2D markers) -- drawn faint
erc <- merge(edges[etype == "gen_rev_causal"], dz[, .(id = label, dx, dy)], by.x = "from", by.y = "id")
erc <- merge(erc, prc[, .(id, px, py)], by.x = "to", by.y = "id")

## ---- palettes / encodings ---------------------------------------------------
rk  <- c(`-1` = "#2166AC", `1` = "#B2182B")
ltv <- c(colocalized = "solid", `cis (coloc pending)` = "44", `reporter (reverse)` = "12")
ar  <- arrow(length = unit(0.06, "cm"), type = "closed")

p <- ggplot() +
  ## flow header
  annotate("text", x = -3.9, y = 3.72, label = "1.  lifestyle & trials", size = 3, fontface = "bold", colour = "grey30") +
  annotate("text", x = 0,    y = 3.72, label = "2.  shared proteins",    size = 3, fontface = "bold", colour = "grey30") +
  annotate("text", x = 3.2,  y = 3.72, label = "3.  disease",            size = 3, fontface = "bold", colour = "grey30") +
  annotate("segment", x = -2.75, xend = -1.15, y = 3.72, yend = 3.72, arrow = arrow(length = unit(0.11, "cm"), type = "closed"), colour = "grey55") +
  annotate("segment", x = 1.20, xend = 2.45, y = 3.72, yend = 3.72, arrow = arrow(length = unit(0.11, "cm"), type = "closed"), colour = "grey55") +
  annotate("segment", x = -0.75, xend = 1.35, y = -0.02, yend = -0.02, linewidth = 0.3, colour = "grey82", linetype = "22") +
  ## perturbation -> protein
  geom_segment(data = eo, aes(sx + 0.25, sy, xend = px - 0.14, yend = py, colour = factor(sign)), linewidth = 0.25, alpha = 0.4) +
  geom_segment(data = ei, aes(sx + 0.25, sy, xend = px - 0.14, yend = py, colour = factor(sign)), linewidth = 0.55, alpha = 0.85) +
  ## faint reverse edges ON causal proteins (they are ALSO broad obesity/T2D markers) -- background layer
  geom_segment(data = erc, aes(dx - 0.5, dy, xend = px + 0.14, yend = py, colour = factor(sign), linetype = tier), linewidth = 0.22, alpha = 0.25, arrow = arrow(length = unit(0.04, "cm"), type = "closed")) +
  ## forward causal: protein -> disease (solid, on top)
  geom_segment(data = efw, aes(px + 0.14, py, xend = dx - 0.5, yend = dy, colour = factor(sign), linetype = tier), linewidth = 0.5, alpha = 0.95, arrow = ar) +
  ## reverse reporter: disease -> protein (thin dotted, subordinate but coloured by direction)
  geom_segment(data = erv, aes(dx - 0.5, dy, xend = px + 0.14, yend = py, colour = factor(sign), linetype = tier), linewidth = 0.38, alpha = 0.8, arrow = ar) +
  ## nodes: fill = exposome R2 (how exposure-responsive the protein is); size = breadth
  geom_point(data = prall, aes(px, py, size = breadth, fill = R2_E), shape = 21, colour = "grey35", stroke = 0.5) +
  geom_text(data = prall, aes(px + 0.2, py, label = label), hjust = 0, size = 2.15, colour = "grey15") +
  ## source boxes drawn with CONSTANT fills (kept off the fill scale, which is used by node R2_E)
  geom_label(data = src[fl == "exposure"], aes(sx, sy, label = label), fill = "#607D8B", size = 2.5, fontface = "bold", label.size = 0, colour = "white", hjust = 1) +
  geom_label(data = src[fl == "RCT"],      aes(sx, sy, label = label), fill = "#CC7722", size = 2.5, fontface = "bold", label.size = 0, colour = "white", hjust = 1) +
  geom_label(data = dz, aes(dx, dy, label = label), size = 2.6, fontface = "bold", fill = "#37474F", colour = "white", hjust = 0) +
  ## group cues
  annotate("text", x = -3.9, y = 3.32, label = "OBSERVATIONAL  (lifestyle)", size = 2.5, colour = "#607D8B", fontface = "bold", hjust = 1) +
  annotate("text", x = -3.9, y = 0.30, label = "INTERVENTIONAL  (trials)", size = 2.5, colour = "#CC7722", fontface = "bold", hjust = 1) +
  annotate("label", x = 0.05, y = 3.30, label = "CAUSAL INTERMEDIATES   (protein → disease)", size = 2.6, colour = "#B2182B", fontface = "bold", hjust = 0.5, fill = "white", label.size = 0, label.padding = unit(0.5, "mm")) +
  annotate("label", x = 0.05, y = 1.56, label = "DISEASE REPORTERS   (disease → protein, reverse)", size = 2.6, colour = "#546E7A", fontface = "bold", hjust = 0.5, fill = "white", label.size = 0, label.padding = unit(0.5, "mm")) +
  scale_colour_manual(values = rk, na.value = "grey80", name = "edge: pushes the next node", labels = c("down (lowers)", "up (raises)"),
                      guide = guide_legend(order = 1, override.aes = list(linewidth = 1.3))) +
  scale_linetype_manual(values = ltv, name = "genetic link", breaks = c("colocalized", "cis (coloc pending)", "reporter (reverse)"),
                        guide = guide_legend(order = 2, override.aes = list(colour = "grey30"))) +
  scale_size_area("exposures reading it", max_size = 5.0, breaks = c(10, 30, 60), guide = guide_legend(order = 3)) +
  scale_fill_gradient(low = "#E5F5E0", high = "#1B7837", name = "exposome R²", na.value = "grey85",
                      guide = guide_colourbar(order = 4, barheight = grid::unit(0.9, "cm"), barwidth = grid::unit(0.35, "cm"))) +
  coord_cartesian(xlim = c(-6.0, 5.0), ylim = c(-1.5, 4.0), clip = "off") +
  theme_void(base_size = 10) +
  theme(plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        legend.position = "right",
        legend.title = element_text(size = 7.3, face = "bold"), legend.text = element_text(size = 7),
        plot.margin = margin(6, 6, 6, 10))

if (!CELL) {
  p <- p + labs(title = "Shared proteins split into causal intermediates and disease reporters",
    subtitle = paste0("Lifestyle and trials move a shared protein set (red = raises it, blue = lowers).\n",
      "Top block is genetically causal for disease by Tier-1 colocalized cis-pQTL MR (solid); bottom block is a downstream marker (dotted, disease → protein).")) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 7.2, hjust = 0.5, colour = "grey35", margin = margin(b = 4)))
}

if (CELL) {
  ggsave(file.path(FIGDIR, "fig_m4_shared_network_cell.png"), p, width = 9.7, height = 4.6, dpi = 400, bg = "white")
  ggsave(file.path(FIGDIR, "fig_m4_shared_network_cell.pdf"), p, width = 9.7, height = 4.6, bg = "white")
  message("fig_m4_shared_network CELL done (", nrow(prc), " causal / ", nrow(prr), " reporter)")
} else {
  heap_emit_figure(p, figure_id, data = edges, category = "exploratory",
                   formats = c("pdf", "png"), width = 11, height = 6.8, website = FALSE)
  message("fig_m4_shared_network standalone done (", nrow(prc), " causal / ", nrow(prr), " reporter; ",
          nrow(efw), " forward / ", nrow(erv), " reverse edges)")
}
