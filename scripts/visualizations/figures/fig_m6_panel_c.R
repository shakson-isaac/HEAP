#!/usr/bin/env Rscript
# ============================================================================
# fig_m6_panel_c.R  -> Module 6 panel c (Q2): does the PES track WITHIN-PERSON change?
# ----------------------------------------------------------------------------
# Left (breadth): within-person Delta-correlation (change in proteome score vs
#   change in actual exposure) per exposure, with a 95% bootstrap CI (resampling
#   held-out people; support/module6_within_ci.R). Colored = proteome, grey =
#   covariate benchmark. Dashed line = no tracking (0). One exemplar labeled per
#   category. Pack-years smoking is labeled at ~0 on purpose: it reads almost
#   perfectly cross-sectionally (panel b) yet its CI straddles zero here -- it is
#   cumulative, so it does not change within a person.
# Right, two concrete exemplars of "tracking change":
#   top    -- alcohol (continuous dose): when reported drinking changes between
#             visits, the proteome alcohol-score moves with it (r with 95% CI).
#   bottom -- smoking (discrete reversibility): the current-smoker score drops
#             ~1.8 SD when people quit and rises when they start; stayers are flat.
# ============================================================================
local({
  cand <- c(file.path(getwd(),"scripts","visualizations","common"),
            "/n/groups/patel/shakson_ukb/HEAP/scripts/visualizations/common")
  cm <- cand[dir.exists(cand)][1]
  for (f in c("figure_paths.R","load_heap_results.R","plot_theme.R","label_helpers.R","export_helpers.R")) source(file.path(cm,f))
})
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(patchwork) })

GREEN <- "#1B7837"; GREY <- "grey60"
# --- render profile: standalone (default) vs small composite CELL ----------
CELL  <- nzchar(Sys.getenv("HEAP_CELL"))
BS    <- if (CELL) 7   else 10     # base_size
PSZ_S <- if (CELL) 0.9 else 1.4    # cloud point
PSZ_L <- if (CELL) 2.3 else 2.9    # exemplar point (bolder than the cloud)
LBL   <- if (CELL) 1.8 else 2.5    # exemplar label (mm)
TTL   <- if (CELL) 8   else 13     # title (pt)
TTLS  <- if (CELL) 7.5 else 9.6    # exemplar sub-panel title (pt)
ANN   <- if (CELL) 2.4 else 3.4    # alcohol r annotation
ANN2  <- if (CELL) 2.1 else 2.9    # smoking slab annotation
C_W   <- 4.42; C_H <- 2.85         # cell size (in); narrower than d; taller (grow bottom row for vertical spacing -- keep C_H==D_H so fonts stay 1:1)
od <- heap_project_output("module6_pes_longitudinal","base")
CAST  <- c("alcohol_intake_frequency_f1558_0_0","pack_years_of_smoking_f20161_0_0","oily_fish_intake_f1329_0_0",
           "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0","processed_meat_intake_f1349_0_0")
# one labeled exemplar per category; for Smoking we force BOTH the high current-status
# tracker and pack-years (~0) so the cumulative-vs-current contrast is explicit.
PREFER <- c(Diet_Weekly = "processed_meat_intake_f1349_0_0",
            Exercise_Freq = "number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0",
            Smoking = "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days",
            Sleep = "sleep_duration_f1160_0_0",
            Alcohol = "alcohol_intake_frequency_f1558_0_0")
EX_FORCE <- c("pack_years_of_smoking_f20161_0_0")   # always labeled (the cautionary case)
EX_LABEL <- c(
  alcohol_intake_frequency_f1558_0_0 = "Alcohol frequency",
  pack_years_of_smoking_f20161_0_0 = "Pack-years smoking",
  current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days = "Current smoking",
  processed_meat_intake_f1349_0_0 = "Processed meat",
  number_of_days_week_of_vigorous_physical_activity_10_plus_minutes_f904_0_0 = "Vigorous activity",
  oily_fish_intake_f1329_0_0 = "Oily fish", coffee_intake_f1498_0_0 = "Coffee",
  sleep_duration_f1160_0_0 = "Sleep duration",
  `vitamin_and_mineral_supplements_f6155_0_0.multi_Multivitamins_..._minerals` = "Multivitamins",
  frequency_of_solarium_sunlamp_use_f2277_0_0 = "Solarium use",
  plays_computer_games_f2237_0_0 = "Computer games",
  average_total_household_income_before_tax_f738_0_0 = "Household income")
PRETTY <- c(Deprivation_Indices="Deprivation / income", Exercise_Freq="Exercise (frequency)",
            Sun_Exposure="Sun exposure", Diet_Weekly="Diet", Internet_Usage="Internet use")
pcat <- function(x) ifelse(x %in% names(PRETTY), PRETTY[x], x)
shorten <- function(x){ x <- gsub("_"," ",x); ifelse(nchar(x)>24, paste0(substr(x,1,22),"…"), x) }
pretty_ex <- function(id) ifelse(id %in% names(EX_LABEL), EX_LABEL[id], shorten(heap_exposure_label(id)))

# ---------------------------------------------------------------------------
# LEFT: breadth -- within-person Delta-correlation per exposure (+ bootstrap CI)
# ---------------------------------------------------------------------------
wc <- fread(file.path(od, "PESlong_base_WithinDeltaCorCI.tsv"))
wc <- wc[!category %in% c("Sexual_Factors")]            # immutable within person
wc[, is_cast := exposure_id %in% CAST]
ord <- wc[, .(m = median(dcor_prot, na.rm = TRUE)), by = category][order(m)]$category
wc[, category := factor(category, levels = ord)]
wc[, ycat := as.integer(category)]
set.seed(1); wc[, yj := ycat + runif(.N, -0.16, 0.16)]
wc[, prefid := ifelse(as.character(category) %in% names(PREFER), PREFER[as.character(category)], NA_character_)]
wc[, pref := !is.na(prefid) & exposure_id == prefid]
ex <- wc[order(category, -pref, -is_cast, -dcor_prot)][, .SD[1], by = category]
ex <- rbind(ex, wc[exposure_id %in% EX_FORCE & !exposure_id %in% ex$exposure_id], fill = TRUE)  # force cautionary labels
ex[, lab := pretty_ex(exposure_id)]
ex[, yj := ycat]    # centre the highlighted exemplar on its row
gp <- 0.045 * (0.98 - (-0.12))                    # fixed gap: label sits this far RIGHT of the dot, on-row
ex[, `:=`(lx = dcor_prot + gp, ly = as.numeric(ycat), hj = 0)]   # ly MUST be double -- integer ycat truncates fractional offsets (ycat+0.62 -> 2), which silently dropped labels onto the wrong row
# two Smoking-row exemplars at very different dcor (0.10 vs 0.77): label them on
# DIFFERENT bands so neither crowds the other or the Internet-use row below.
# pack-years lifts ABOVE-left (clears the dense low-dcor cluster); current smoking sits
# ON-row to the LEFT of its far-right dot (the empty 0.3-0.7 gap), leader pointing right
# (no room to its right at dcor ~0.77; below would hit "Computer games"; on the dot would
# bury the text under the large point).
ex[exposure_id == "pack_years_of_smoking_f20161_0_0", `:=`(lx = dcor_prot + gp, ly = ycat + 0.62, hj = 0)]
ex[exposure_id == "current_tobacco_smoking_f1239_0_0_Yes._on_most_or_all_days", `:=`(lx = dcor_prot - 0.05, ly = ycat, hj = 1)]
ex[exposure_id == "pack_years_of_smoking_f20161_0_0", lab := "Pack-years"]   # drop "smoking" (axis row already says Smoking) so the above-left label clears "Current smoking" on the row

pL <- ggplot(wc, aes(y = yj)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4, linetype = "22") +
  geom_errorbarh(aes(xmin = prot_lo, xmax = prot_hi, colour = category), height = 0, linewidth = 0.32, alpha = 0.30) +
  geom_point(aes(x = dcor_prot, colour = category), size = PSZ_S, alpha = if (CELL) 0.6 else 0.9) +
  geom_errorbarh(data = ex, aes(xmin = prot_lo, xmax = prot_hi, colour = category), height = 0, linewidth = 0.6) +
  geom_point(data = ex, aes(x = dcor_prot, colour = category), size = PSZ_L) +
  geom_point(data = ex, aes(x = dcor_prot, y = yj), size = PSZ_L, shape = 1, colour = "grey20",
             stroke = if (CELL) 0.4 else 0.5, inherit.aes = FALSE) +  # dark outline ring -> exemplar reads bolder
  {if (CELL) list(
      geom_segment(data = ex, aes(x = dcor_prot, xend = lx, y = ycat, yend = ly), colour = "grey55", linewidth = 0.25, inherit.aes = FALSE),
      geom_text(data = ex, aes(x = lx, y = ly, label = lab), hjust = ex$hj, colour = "#222222", size = LBL, inherit.aes = FALSE))
   else ggrepel::geom_text_repel(data = ex, aes(x = dcor_prot, label = lab), colour = "#222222", size = LBL,
      box.padding = 0.5, point.padding = 0.2, min.segment.length = 0, segment.colour = "grey65",
      max.overlaps = 40, force = 1.4, seed = 1)} +
  scale_colour_exposure(drop = TRUE, guide = "none") +
  scale_y_continuous(breaks = seq_along(ord), labels = pcat(ord), expand = expansion(add = 0.6)) +
  scale_x_continuous(limits = c(-0.12, 0.98), expand = expansion(mult = c(0.01, 0.04))) +
  labs(title = if (CELL) "Within-person tracking" else "The proteome score tracks within-person change in the exposure",
       x = if (CELL) expression(Delta*"-correlation  ("*Delta*"score vs "*Delta*"exposure)")
           else expression("within-person "*Delta*"-correlation  ("*Delta*"score vs "*Delta*"exposure)"), y = NULL) +
  theme_heap(base_size = BS) +
  theme(panel.grid.major.y = if (CELL) element_blank() else element_line(colour = "grey93"),
        panel.grid.major.x = if (CELL) element_blank() else element_line(colour = "grey92"),
        axis.text.y = element_text(size = rel(0.85)),
        axis.title = element_text(face = if (CELL) "plain" else "bold"),
        plot.title = element_text(size = if (CELL) TTL else rel(1.05), hjust = if (CELL) 0.5 else 0),
        plot.title.position = "panel", plot.margin = margin(4,3,2,5))

# ---------------------------------------------------------------------------
# RIGHT-TOP: alcohol -- continuous dose tracking (person-level), r with 95% CI
# ---------------------------------------------------------------------------
hsA <- fread(file.path(od, "PESlong_base_alcohol_intake_frequency_f1558_0_0_HoldoutScores.tsv"))
bA <- hsA[instance == 0, .(eid, y0 = y_raw, p0 = pes_prot_z)]
fA <- hsA[instance %in% c(2,3), .(eid, y1 = y_raw, p1 = pes_prot_z)]
chA <- merge(fA, bA, by = "eid"); chA[, `:=`(dY = y1 - y0, dP = p1 - p0)]
chA <- chA[is.finite(dY) & is.finite(dP)]
ctA <- cor.test(chA$dY, chA$dP)
rlab <- sprintf("r = %.2f  [%.2f, %.2f]", ctA$estimate, ctA$conf.int[1], ctA$conf.int[2])
set.seed(1)
pTR <- ggplot(chA, aes(dY, dP)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.4) +
  geom_jitter(width = 0.13, height = 0, alpha = 0.12, colour = GREEN, size = if (CELL) 0.5 else 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = GREEN, fill = GREEN, alpha = 0.18, linewidth = if (CELL) 0.7 else 0.9) +
  annotate("text", x = -4.7, y = 2.6, hjust = 0, label = rlab, colour = GREEN, fontface = "bold", size = ANN) +
  scale_x_continuous(breaks = seq(-4, 4, 2)) + coord_cartesian(ylim = c(-3, 3)) +
  labs(title = if (CELL) "Alcohol dose" else "Alcohol — continuous dose: score follows reported drinking",
       x = if (CELL) expression(Delta*" drinking") else expression(Delta*" drinking frequency  (+ = drinks more)"),
       y = if (CELL) expression(Delta*" score (z)") else expression(Delta*" alcohol-score (z)")) +
  theme_heap(base_size = BS) + theme(plot.title = element_text(size = TTLS, hjust = if (CELL) 0.5 else 0),
        plot.title.position = "panel", axis.title = element_text(face = if (CELL) "plain" else "bold"),
        plot.margin = margin(7,6,2,6), panel.grid.major = element_blank())   # .major (not .grid) -- theme_heap sets the child, parent won't override

# ---------------------------------------------------------------------------
# RIGHT-BOTTOM: smoking -- discrete reversibility (current-smoker score)
# ---------------------------------------------------------------------------
hsS <- fread(file.path(od, "PESlong_base_smoking_status_f20116_0_0_Current_HoldoutScores.tsv"))
bS <- hsS[instance == 0, .(eid, y0 = y_raw, p0 = pes_prot_z)]
fS <- hsS[instance %in% c(2,3), .(eid, y1 = y_raw, p1 = pes_prot_z)]
chS <- merge(fS, bS, by = "eid")
chS[, grp := fcase(y0 == 1 & y1 == 0, "Quit", y0 == 0 & y1 == 1, "Started",
                   y0 == 1 & y1 == 1, "Stayed smoker", y0 == 0 & y1 == 0, "Stayed non-smoker")]
# within-person Delta-correlation (status change vs score change) + 95% CI for the annotation
ctS <- cor.test(chS$y1 - chS$y0, chS$p1 - chS$p0)
slab <- sprintf("r = %.2f  [%.2f, %.2f]", ctS$estimate, ctS$conf.int[1], ctS$conf.int[2])  # plain r (Delta glyph won't render); shown INSIDE the panel
# before -> after LEVELS (not the change): smokers stay high, non-smokers low, quit falls, start rises.
gL  <- chS[!is.na(grp), .(n=.N, b_m=mean(p0), b_se=sd(p0)/sqrt(.N), f_m=mean(p1), f_se=sd(p1)/sqrt(.N)), by=grp]
gLL <- rbindlist(list(gL[, .(grp, n, visit="Baseline",  m=b_m, se=b_se)],
                      gL[, .(grp, n, visit="Follow-up", m=f_m, se=f_se)]))
gLL[, visit := factor(visit, levels=c("Baseline","Follow-up"))]
gLL[, `:=`(lo = m-1.96*se, hi = m+1.96*se)]
gLL[, changed := grp %in% c("Quit","Started")]
# short labels for CELL; Delta-corr moves to the subtitle so the plot interior is clean.
# stayers (flat grey lines) read as plain "Smoker"/"Non-smoker"; changers keep Quit/Started.
gLL[, glab := if (CELL) fcase(grp=="Stayed smoker","Smoker", grp=="Stayed non-smoker","Non-smoker", default=grp) else grp]
gLL[, flab := if (CELL) fifelse(changed, sprintf("%s (n=%d)", glab, n), glab) else sprintf("%s (n=%d)", grp, n)]  # cell: n on the TRANSITION groups (Started/Quit); stayers stay short
# DETERMINISTIC labels: spread the 4 groups evenly down the clear right margin,
# ordered by their follow-up score, each tied back to its endpoint by a leader line.
fu <- copy(gLL[visit == "Follow-up"])[order(-m)]
.yr <- range(c(gLL$lo, gLL$hi))
fu[, `:=`(ytar = seq(.yr[2], .yr[1], length.out = .N), xtar = 2.24)]
pBR <- ggplot(gLL, aes(visit, m, group=grp, colour=changed)) +
  geom_segment(data=fu, aes(x=2.05, xend=xtar-0.04, y=m, yend=ytar, colour=changed),
               linewidth=0.25, alpha=0.55, inherit.aes=FALSE) +
  geom_line(linewidth=if (CELL) 0.8 else 1.0, alpha=0.9) +
  geom_pointrange(aes(ymin=lo, ymax=hi), size=if (CELL) 0.28 else 0.45) +
  geom_text(data=fu, aes(x=xtar, y=ytar, label=flab, colour=changed), hjust=0, size=LBL, inherit.aes=FALSE) +
  annotate("text", x=0.5, y=Inf, hjust=0, vjust=1.3, label=slab, colour=GREEN, fontface="bold", size=ANN2) +  # stat INSIDE the panel, far top-left so it clears the "Smoker" label
  scale_colour_manual(values=c(`TRUE`=GREEN, `FALSE`=GREY), guide="none") +
  scale_x_discrete(expand=expansion(add=c(0.55, if (CELL) 1.95 else 1.8))) +   # right room for the "Started (n=15)" labels; keep the Baseline<->Follow-up gap
  scale_y_continuous(expand=expansion(mult=c(0.12, if (CELL) 0.20 else 0.10))) +
  labs(title=if (CELL) "Smoking reversibility" else "Smoking — discrete reversibility: the score follows status change",
       x=NULL, y=if (CELL) "score (z)" else "current-smoker score (z)") +
  theme_heap(base_size=BS) + theme(plot.title=element_text(size=TTLS, hjust=if (CELL) 0.5 else 0),
        plot.title.position="panel", axis.title=element_text(face=if (CELL) "plain" else "bold"),
        axis.text.x=element_text(size=if (CELL) rel(0.82) else rel(1)),
        plot.margin=margin(7,6,2,6), panel.grid.major=element_blank())

# ---------------------------------------------------------------------------
p <- pL | (pTR / pBR)
p <- p + plot_layout(widths = c(1.3, 1))   # narrower within-person tracking (pL); frees width for the alcohol/smoking stack + panel d
FIGDIR <- file.path(Sys.getenv("HEAP_PROJECT_ROOT", "/n/groups/patel/IGLOO/UKB/HEAP"), "figures", "exploratory/module6")
if (CELL) {
  # NB: overall title removed -> lives in the separate figure legend file.
  ggsave(file.path(FIGDIR,"fig_m6_panel_c_cell.png"), p, width = C_W, height = C_H, dpi = 400, bg = "white")
  ggsave(file.path(FIGDIR,"fig_m6_panel_c_cell.pdf"), p, width = C_W, height = C_H, bg = "white")
  message("panel c CELL done")
} else {
  p <- p + plot_annotation(title = "The exposure scores track change within the same person over time",
          subtitle = "left: Δ-correlation per exposure (proteome colored + 95% CI, covariate benchmark grey); right: two ways the score follows change",
          theme = theme(plot.title = element_text(face = "bold", size = 13),
                        plot.subtitle = element_text(size = 9.5, colour = "grey35")))
  ggsave(file.path(FIGDIR,"fig_m6_panel_c.png"), p, width = 12.5, height = 6.0, dpi = 175, bg = "white")
  message("panel c done; alcohol ", rlab, "; quit baseline->follow ",
          round(gL[grp=="Quit"]$b_m,2), "->", round(gL[grp=="Quit"]$f_m,2))
}
