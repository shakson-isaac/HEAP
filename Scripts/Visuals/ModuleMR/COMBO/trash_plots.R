#, bg = "white")

# comp_long <- motif_label_counts %>%
#   select(motif_label, n_shared_same, n_unique_UKB, n_unique_DECODE) %>%
#   tidyr::pivot_longer(-motif_label, names_to = "component", values_to = "n") %>%
#   mutate(component = recode(component,
#                             n_shared_same = "Shared",
#                             n_unique_UKB = "UKB pQTLs",
#                             n_unique_DECODE = "deCODE pQTLs"))
# 
# p_components <- ggplot(comp_long,
#                        aes(x = reorder(motif_label, n, FUN = sum), y = n + 1, fill = component)) +
#   geom_col(position = position_dodge(width = 0.8)) +
#   coord_flip() +
#   scale_y_log10(labels = comma,
#                 expand = expansion(mult = c(0, 0.10))) +
#   scale_x_discrete(limits = positions) +
#   theme_classic() +
#   labs(x = NULL, y = "Triad count (log10)", title = "Shared vs Unique Triad Hits (UKB vs deCODE)") +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(size = 10),
#         axis.text.y = element_text(size = 12))
# 
# p_components
# 
# ggsave(file.path(OUTDIR, "MotifLabel_Separate_DECODEUKB.png"),
#        p_components, width = 6, height = 4, units = "in", dpi = 400)

# #Make 0's NA:
# motif_label_counts[motif_label_counts == 0] <- NA
# motif_label_counts[motif_label_counts == 1] <- 2
# 
# 
# plot_long <- motif_label_counts %>%
#   select(motif_label, n_shared_same, n_unique_UKB, n_unique_DECODE) %>%
#   pivot_longer(-motif_label, names_to="component", values_to="n") %>%
#   mutate(component = recode(component,
#                             n_shared_same="Shared (same label)",
#                             n_unique_UKB="UKB only",
#                             n_unique_DECODE="deCODe only"
#   ))
# 
# ggplot(plot_long, aes(x = reorder(motif_label, n, FUN=sum), y = n, fill = component)) +
#   geom_col(position = position_dodge(width = 0.8)) +
#   coord_flip() +
#   scale_y_log10() +
#   theme_bw() +
#   labs(x=NULL, y="Triplets (log10)")
# 
# 
# 
# 
# # Plot: overlap/unique stacked by motif_label (LOG SCALE + smart labels)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(scales)
# 
# plot_lab_overlap <- motif_label_counts %>%
#   select(motif_label, n_shared_same, n_unique_UKB, n_unique_DECODE) %>%
#   pivot_longer(-motif_label, names_to = "category", values_to = "n") %>%
#   mutate(
#     category = recode(category,
#                       n_shared_same   = "Shared (same label)",
#                       n_unique_UKB    = "UKB only (label differs)",
#                       n_unique_DECODE = "deCODe only (label differs)"),
#     n = as.numeric(n)
#   )
# 
# # totals for % labels (raw scale)
# totals_by_label <- plot_lab_overlap %>%
#   group_by(motif_label) %>%
#   summarise(total = sum(n, na.rm = TRUE), .groups = "drop")
# 
# stack_order <- c("deCODe only (label differs)",
#                  "UKB only (label differs)",
#                  "Shared (same label)")  # left->right after flip
# 
# plot_lab_overlap_lab <- plot_lab_overlap %>%
#   left_join(totals_by_label, by = "motif_label") %>%
#   mutate(
#     category = factor(category, levels = stack_order),
#     pct = ifelse(total > 0, 100 * n / total, NA_real_),
#     label = case_when(
#       n <= 0 ~ "",
#       n < 100 ~ as.character(n),
#       TRUE ~ paste0(round(pct), "%")
#     ),
#     n_plot = n + 1  # <-- pseudo-count ONLY for log scale
#   ) %>%
#   arrange(motif_label, category) %>%
#   group_by(motif_label) %>%
#   mutate(
#     xmin = cumsum(lag(n_plot, default = 0)),
#     xmax = cumsum(n_plot),
#     xmid = (xmin + xmax) / 2
#   ) %>%
#   ungroup() %>%
#   mutate(
#     # optional: hide labels for very thin segments (raw n)
#     label = ifelse(n < 5, "", label)
#   )
# 
# p1 <- ggplot(plot_lab_overlap_lab,
#              aes(y = reorder(motif_label, total), x = n_plot, fill = category)) +
#   geom_col() +
#   geom_text(aes(x = xmid, label = label), size = 3) +
#   scale_x_continuous(
#     trans = scales::log10_trans(),
#     breaks = c(1, 10, 100, 1e3, 1e4, 1e5),
#     labels = function(x) scales::comma(x - 1)  # show ticks in original n units
#   ) +
#   theme_bw() +
#   labs(y = NULL, x = "# triplets (log10 scale)",
#        title = "Motif-label overlap vs unique (UKB vs deCODe)") +
#   theme(legend.position = "bottom")
# 
# p1


#ggsave(file.path(OUTDIR, "MotifLabel_overlap_unique_stacked_log.png"),
#       p1, width = 8, height = 5, units = "in", dpi = 400)
# plot_lab_overlap <- motif_label_counts %>%
#   select(motif_label, n_shared_same, n_unique_UKB, n_unique_DECODE) %>%
#   pivot_longer(-motif_label, names_to = "category", values_to = "n") %>%
#   mutate(category = recode(category,
#                            n_shared_same = "Shared (same label)",
#                            n_unique_UKB = "UKB only (label differs)",
#                            n_unique_DECODE = "deCODe only (label differs)"))
# 
# p1 <- ggplot(plot_lab_overlap, aes(x = reorder(motif_label, n), y = n, fill = category)) +
#   geom_col() +
#   coord_flip() +
#   theme_bw() +
#   labs(x = NULL, y = "# triplets", title = "Motif-label overlap vs unique (UKB vs deCODe)")
# 
# p1
# 
# ggsave(file.path(OUTDIR, "MotifLabel_overlap_unique_stacked.png"),
#        p1, width = 8, height = 5, units = "in", dpi = 400)
