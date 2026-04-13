# plot_HEAP_GLP1_MR_onepanel <- function(
    #     EXPOSURE_TO_PLOT,
#     GLP_ARM = c("GLP1_1", "GLP1_2"),
#     disease_for_arm = NULL,
#     mr_alpha = 0.05,
#     
#     # text controls
#     title = "HEAP vs GLP1 protein shifts (MR edge encoded)",
#     subtitle = NULL,
#     xlab = "HEAP beta (Exposure → Protein)",
#     ylab = NULL,
#     show_title = TRUE,
#     show_subtitle = TRUE,
#     show_axis_titles = TRUE,
#     
#     # correlation annotation controls
#     corr_model = "Model6",
#     corr_source = c("GLP1", "HERITAGE"),   # which precomputed correlation to annotate
#     corr_loc = c("topleft","topright","bottomleft","bottomright"),
#     corr_digits = 2,
#     
#     # label controls
#     label_n = 12,
#     label_only_mr = FALSE,
#     
#     # legend / compactness
#     legend_position = c("right","bottom","none"),
#     base_size = 12
# ) {
#   
#   # -----------------------------
#   # Helpers (self-contained)
#   # -----------------------------
#   fmt_p <- function(p) {
#     if (is.na(p)) return("NA")
#     if (p < 1e-300) return("<1e-300")
#     format.pval(p, digits = 2, eps = 1e-300)
#   }
#   
#   # canonicalize exposure IDs for relaxed matching:
#   # - only removes a trailing ordinal level like "_04", "_1", etc.
#   # - does NOT touch ".multi_*" suffixes
#   canon_exposure <- function(x) {
#     sub("_[0-9]+$", "", x)
#   }
#   
#   cap_n <- function(x, n) x[seq_len(min(length(x), n))]
#   
#   # -----------------------------
#   # Arg parsing
#   # -----------------------------
#   GLP_ARM <- match.arg(GLP_ARM)
#   corr_source <- match.arg(corr_source)
#   corr_loc <- match.arg(corr_loc)
#   legend_position <- match.arg(legend_position)
#   
#   # if (is.null(disease_for_arm)) {
#   #   disease_for_arm <- if (GLP_ARM == "GLP1_1") "finngen_R12_E4_OBESITY" else "finngen_R12_T2D"
#   # }
#   # 
#   # glp_beta_col <- if (GLP_ARM == "GLP1_1"){
#   #                             "GLP1_effect1"
#   #                   } else if (GLP_ARM == "GLP1_2"){
#   #                             "GLP1_effect2"
#   #                   } #else "HERITAGE_effect"
#   # if (is.null(ylab)) ylab <- paste0("GLP1 effect (", GLP_ARM, ")")
#   
#   GLP_ARM <- match.arg(GLP_ARM, choices = c("GLP1_1", "GLP1_2", "HERITAGE"))
#   
#   y_col <- if (GLP_ARM == "GLP1_1") {
#     "GLP1_effect1"
#   } else if (GLP_ARM == "GLP1_2") {
#     "GLP1_effect2"
#   } else {
#     "HERITAGE_effect"
#   }
#   
#   if (is.null(ylab)) ylab <- paste0(y_col)
#   
#   exposure_key_plot <- canon_exposure(EXPOSURE_TO_PLOT)
#   
#   # -----------------------------
#   # 1) Build HEAP subset (RELAXED exposure matching)
#   # -----------------------------
#   heap_dt <- as.data.table(HEAPint@sList[[corr_model]])
#   if (!all(c("ID","EntrezGeneSymbol") %in% names(heap_dt))) {
#     stop("HEAPint@sList[[modelType]] must contain ID and EntrezGeneSymbol.")
#   }
#   setnames(heap_dt,
#            old = c("ID","EntrezGeneSymbol","Estimate","Std. Error"),
#            new = c("Exposure","Protein","beta_HEAP","se_HEAP"),
#            skip_absent = TRUE)
#   
#   heap_dt[, Exposure_key := canon_exposure(Exposure)]
#   
#   # prefer exact exposure if present; else fall back to key match
#   if (any(heap_dt$Exposure == EXPOSURE_TO_PLOT, na.rm = TRUE)) {
#     heap_pick <- heap_dt[Exposure == EXPOSURE_TO_PLOT]
#   } else {
#     heap_pick <- heap_dt[Exposure_key == exposure_key_plot]
#   }
#   
#   # collapse to 1 row per protein if multiple exposure strings map to the same key
#   heap_sub <- heap_pick[
#     , .(
#       Exposure  = EXPOSURE_TO_PLOT,
#       beta_HEAP = mean(beta_HEAP, na.rm = TRUE),
#       beta_GLP1 = mean(get(y_col), na.rm = TRUE)
#     ),
#     by = .(Protein)
#   ][!is.na(beta_GLP1)]
#   
#   if (nrow(heap_sub) == 0) stop("No HEAP rows for exposure with ", GLP_ARM, " available after relaxed matching.")
#   
#   # -----------------------------
#   # 2) Merge MR (restricted to ONE disease) using relaxed exposure matching
#   # -----------------------------
#   MR <- as.data.table(MRres)
#   MR[, Exposure_key := canon_exposure(Exposure)]
#   
#   mr_sub <- MR[Exposure_key == exposure_key_plot & Disease == disease_for_arm]
#   
#   # Merge by Protein only (exposure already fixed by heap_sub)
#   dt <- merge(heap_sub, mr_sub, by = "Protein", all.x = TRUE)
#   dt[, Exposure := EXPOSURE_TO_PLOT]
#   
#   # -----------------------------
#   # 3) Merge Soma–Olink correlation + size metric
#   # -----------------------------
#   prot_rel_dt <- as.data.table(prot_rel)
#   if ("EntrezGeneSymbol" %in% names(prot_rel_dt)) {
#     setnames(prot_rel_dt, "EntrezGeneSymbol", "Protein", skip_absent = TRUE)
#   }
#   if (!"r_crossv2" %in% names(prot_rel_dt)) {
#     if ("r_cross" %in% names(prot_rel_dt)) prot_rel_dt[, r_crossv2 := r_cross]
#   }
#   
#   dt <- merge(dt, prot_rel_dt[, .(Protein, r_crossv2)], by = "Protein", all.x = TRUE)
#   
#   # r_crossv2 is SomaScan–Olink correlation (can be negative)
#   dt[, olink_soma_r := r_crossv2]
#   
#   # size: use r directly, but only positive values contribute to size (negative/0 -> 0)
#   dt[, size_r := pmax(0, olink_soma_r)]
#   dt[is.na(size_r), size_r := 0]
#   
#   # alpha: r<=0 (or NA) extremely faint
#   dt[, alpha_r := fifelse(!is.na(olink_soma_r) & olink_soma_r > 0, 0.9, 0.01)]
#   
#   # -----------------------------
#   # 4) MR edge significance + priority assignment
#   # PDcis > PDtrans > DP > None
#   # -----------------------------
#   dt[, sig_PDcis   := !is.na(padj_PDcis)   & padj_PDcis   < mr_alpha]
#   dt[, sig_PDtrans := !is.na(padj_PDtrans) & padj_PDtrans < mr_alpha]
#   dt[, sig_DP      := !is.na(padj_DP)      & padj_DP      < mr_alpha]
#   
#   dt[, mr_edge_sig := fifelse(
#     sig_PDcis, "PDcis",
#     fifelse(sig_PDtrans, "PDtrans",
#             fifelse(sig_DP, "DP", "None"))
#   )]
#   dt[, mr_edge_sig := factor(mr_edge_sig, levels = c("None","PDcis","PDtrans","DP"))]
#   
#   dt[, n_sig_edges := as.integer(sig_PDcis) + as.integer(sig_PDtrans) + as.integer(sig_DP)]
#   
#   # -----------------------------
#   # 5) Labels: prioritize MR edge category (PDcis > PDtrans > DP) then fill extremes
#   # -----------------------------
#   pd_cis_prots   <- dt[mr_edge_sig=="PDcis"][order(-abs(beta_GLP1))]$Protein |> unique()
#   pd_trans_prots <- dt[mr_edge_sig=="PDtrans"][order(-abs(beta_GLP1))]$Protein |> unique()
#   dp_prots       <- dt[mr_edge_sig=="DP"][order(-abs(beta_GLP1))]$Protein |> unique()
#   
#   n_cis   <- ceiling(label_n * 0.5)
#   n_trans <- ceiling(label_n * 0.3)
#   n_dp    <- ceiling(label_n * 0.2)
#   
#   lab_mr <- unique(c(
#     cap_n(pd_cis_prots, n_cis),
#     cap_n(pd_trans_prots, n_trans),
#     cap_n(dp_prots, n_dp)
#   ))
#   
#   n_left <- max(0, label_n - length(lab_mr))
#   
#   fill_prots <- unique(c(
#     dt[order(-abs(beta_HEAP))]$Protein,
#     dt[order(-abs(beta_GLP1))]$Protein
#   ))
#   
#   lab_fill <- setdiff(fill_prots, lab_mr)
#   lab_prots <- unique(c(lab_mr, cap_n(lab_fill, n_left)))
#   
#   # optional: only label MR hits (if requested) but still cap to label_n
#   if (isTRUE(label_only_mr)) {
#     lab_prots <- cap_n(lab_mr, label_n)
#   }
#   
#   dt[, label_me := Protein %in% lab_prots]
#   dt[, label_once := label_me & !duplicated(Protein)]
#   
#   # -----------------------------
#   # 6) Pull precomputed correlation + p-value for this exposure
#   # NOTE: here we also relax matching on the rowname exposure key
#   # -----------------------------
#   ctab <- as.data.table(HEAPint@cList[[corr_model]], keep.rownames = "Exposure")
#   ptab <- as.data.table(HEAPint@pList[[corr_model]], keep.rownames = "Exposure")
#   
#   ctab[, Exposure_key := canon_exposure(Exposure)]
#   ptab[, Exposure_key := canon_exposure(Exposure)]
#   
#   corr_col <- if (corr_source == "HERITAGE") "HERITAGE_effect" else y_col
#   
#   # prefer exact exposure if present; else key match
#   if (any(ctab$Exposure == EXPOSURE_TO_PLOT, na.rm = TRUE)) {
#     r_val <- ctab[Exposure == EXPOSURE_TO_PLOT, get(corr_col)]
#     p_val <- ptab[Exposure == EXPOSURE_TO_PLOT, get(corr_col)]
#   } else {
#     r_val <- ctab[Exposure_key == exposure_key_plot, get(corr_col)]
#     p_val <- ptab[Exposure_key == exposure_key_plot, get(corr_col)]
#   }
#   
#   r_txt <- ifelse(is.na(r_val), "NA", formatC(r_val, digits = corr_digits, format = "f"))
#   p_txt <- fmt_p(p_val)
#   
#   corr_label <- paste0(
#     corr_source, " corr: r = ", r_txt, "\n",
#     "p = ", p_txt
#   )
#   
#   if (is.null(subtitle)) {
#     subtitle <- paste0("Exposure: ", EXPOSURE_TO_PLOT, " | Disease: ", disease_for_arm)
#   }
#   
#   # -----------------------------
#   # 7) Compact plot theme + correlation placement
#   # -----------------------------
#   xr <- range(dt$beta_HEAP, na.rm = TRUE)
#   yr <- range(dt$beta_GLP1, na.rm = TRUE)
#   xpad <- diff(xr) * 0.02
#   ypad <- diff(yr) * 0.04
#   
#   ann_x <- switch(corr_loc,
#                   "topleft"     = xr[1] + xpad,
#                   "bottomleft"  = xr[1] + xpad,
#                   "topright"    = xr[2] - xpad,
#                   "bottomright" = xr[2] - xpad)
#   ann_y <- switch(corr_loc,
#                   "topleft"     = yr[2] - ypad,
#                   "topright"    = yr[2] - ypad,
#                   "bottomleft"  = yr[1] + ypad,
#                   "bottomright" = yr[1] + ypad)
#   ann_hjust <- if (grepl("right", corr_loc)) 1 else 0
#   ann_vjust <- if (grepl("bottom", corr_loc)) 0 else 1
#   
#   # -----------------------------
#   # 8) Draw plot
#   # -----------------------------
#   p <- ggplot(dt, aes(x = beta_HEAP, y = beta_GLP1)) +
#     geom_hline(yintercept = 0, linewidth = 0.25, color = "grey75") +
#     geom_vline(xintercept = 0, linewidth = 0.25, color = "grey75") +
#     
#     geom_point(aes(color = mr_edge_sig, size = size_r, alpha = alpha_r)) +
#     
#     ggrepel::geom_text_repel(
#       data = dt[label_once == TRUE],
#       aes(label = Protein),
#       size = 3,
#       box.padding = 0.2,
#       point.padding = 0.12,
#       min.segment.length = 0,
#       max.overlaps = 60
#     ) +
#     
#     annotate("label",
#              x = ann_x, y = ann_y,
#              label = corr_label,
#              hjust = ann_hjust, vjust = ann_vjust,
#              size = 3.2,
#              label.size = 0.25,
#              fill = "white") +
#     
#     scale_alpha_identity(guide = "none") +
#     scale_size_continuous(
#       range = c(1.6, 5.2),
#       breaks = c(0.25, 0.5, 0.75),
#       limits = c(0, 1),
#       name = "SomaScan–Olink\nr (pos only)"
#     ) +
#     
#     scale_color_manual(
#       values = c(
#         "None"    = "grey75",
#         "PDcis"   = "#1b9e77",
#         "PDtrans" = "#7570b3",
#         "DP"      = "#d95f02"
#       ),
#       name = paste0("MR significant edge\n(adj.p<", mr_alpha, ")")
#     ) +
#     
#     labs(
#       title = if (show_title) title else NULL,
#       subtitle = if (show_subtitle) subtitle else NULL,
#       x = if (show_axis_titles) xlab else NULL,
#       y = if (show_axis_titles) ylab else NULL
#     ) +
#     
#     theme_classic(base_size = base_size) +
#     theme(
#       legend.position = if (legend_position == "none") "none" else legend_position,
#       legend.title = element_text(face = "bold"),
#       legend.key.height = unit(0.8, "lines"),
#       legend.spacing.y = unit(0.2, "lines"),
#       plot.title.position = "plot",
#       plot.margin = margin(6, 6, 6, 6),
#       plot.subtitle = element_text(size = base_size - 1)
#     ) +
#     coord_cartesian(clip = "off")
#   
#   return(p)
# }