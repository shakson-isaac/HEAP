#!/usr/bin/env bash
# ============================================================================
# make_figures.sh — rebuild the main (Figs 1-6) and supplementary figures of the
# HEAP manuscript from completed pipeline outputs.
#
# Usage (from anywhere):
#   bash scripts/visualizations/make_figures.sh            # main + supp
#   bash scripts/visualizations/make_figures.sh main       # Figs 1-6 only
#   bash scripts/visualizations/make_figures.sh supp       # supplementary only
#   bash scripts/visualizations/make_figures.sh fig3       # one main figure
#
# Inputs : $HEAP_PROJECT_ROOT/output/...   (module outputs + figure summaries;
#          see scripts/visualizations/README.md for the summary scripts that
#          must run first)
# Outputs: $HEAP_PROJECT_ROOT/figures/{exploratory,supplement}/<module>/...
#          plus a copy of each finished main figure at
#          $HEAP_PROJECT_ROOT/figures/manuscript/Fig<N>.pdf
#
# Needs  : R 4.4.2 (packages in config/r_packages.tsv), tectonic (TikZ
#          schematics and the vector composites of Figs 3-4; set $TECTONIC or put
#          it on PATH), and pdftocairo (poppler).
# ============================================================================
set -euo pipefail

export HEAP_ROOT="${HEAP_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export HEAP_PATHS_FILE="${HEAP_PATHS_FILE:-$HEAP_ROOT/workflow/00_paths.R}"
: "${HEAP_PROJECT_ROOT:?set HEAP_PROJECT_ROOT to the directory holding output/ (figures/ is written beside it)}"
export HEAP_PROJECT_ROOT

V="$HEAP_ROOT/scripts/visualizations"
FIG="$HEAP_PROJECT_ROOT/figures"
EXP="$FIG/exploratory"
FINAL="$FIG/manuscript"
cd "$HEAP_ROOT"   # several panel scripts source common/ relative to the repo root
# Not every panel script creates its output directory; on a fresh tree they must exist.
mkdir -p "$EXP"/module{1,3,4,5,6} "$FIG/main/module2" "$FINAL"

TECTONIC="${TECTONIC:-$(command -v tectonic || true)}"
[ -x "$HEAP_ROOT/.toolchain/pandocenv/bin/tectonic" ] && [ -z "$TECTONIC" ] && TECTONIC="$HEAP_ROOT/.toolchain/pandocenv/bin/tectonic"
need() { command -v "$1" >/dev/null 2>&1 || [ -x "$1" ] || { echo "ERROR: $1 not found" >&2; exit 1; }; }

say()  { printf '\n== %s ==\n' "$*"; }
r()    { echo "  Rscript $*"; Rscript "$@"; }
cell() { echo "  HEAP_CELL=1 Rscript $*"; HEAP_CELL=1 Rscript "$@"; }

# Compile a standalone TikZ file to <outdir>/<stem>.pdf and a PNG at <dpi>.
# The composites measure and place the PNG, so its resolution is part of the figure.
schematic() {  # schematic <tex> <outdir> <dpi>
  local tex="$1" out="$2" dpi="$3" stem; stem="$(basename "$tex" .tex)"
  need "$TECTONIC"; need pdftocairo; mkdir -p "$out"
  echo "  tectonic $tex -> $out/$stem.{pdf,png} @${dpi}dpi"
  "$TECTONIC" -X compile --outdir "$out" "$tex" >/dev/null
  pdftocairo -png -r "$dpi" -singlefile "$out/$stem.pdf" "$out/$stem"
}

# Compile a TeX file that a builder wrote, in its own directory.
compile_here() {  # compile_here <dir> <stem>
  need "$TECTONIC"
  echo "  tectonic $1/$2.tex"
  ( cd "$1" && "$TECTONIC" -X compile "$2.tex" >/dev/null )
}

publish() {  # publish <N> <source pdf>
  mkdir -p "$FINAL"; cp -f "$2" "$FINAL/Fig$1.pdf"; echo "  -> $FINAL/Fig$1.pdf"
}

fig1() {
  say "Fig 1 — genetic vs. exposomic spectrum"
  schematic "$V/figures/fig1_exp_S1.tex"          "$EXP/module1" 400
  schematic "$V/figures/fig1_partition_block.tex" "$EXP/module1" 380
  r    "$V/build_panelA_S1.R"
  HEAP_VERT=1 cell "$V/figures/cmp_greml_heap_spectrum.R"
  cell "$V/figures/fig_expo_signatures.R"
  r    "$V/build_module1_fig1_horizontal.R"
  publish 1 "$EXP/module1/fig_module1_fig1_horizontal.pdf"
}

fig2() {
  say "Fig 2 — exposure-protein associations"
  r "$V/build_module2_fig3_composite.R"
  publish 2 "$FIG/main/module2/fig_exposomic_composite.pdf"
}

fig3() {
  say "Fig 3 — mediation"
  schematic "$V/figures/fig_module3_schematic_compact.tex" "$EXP/module3" 300
  cell "$V/figures/fig_mediation_scale_main.R"
  cell "$V/figures/fig_mediation_pleiotropy.R"
  cell "$V/figures/fig_mediation_forest.R"
  r    "$V/build_module3_composite_vector.R"
  compile_here "$EXP/module3" Module3_composite_vector
  publish 3 "$EXP/module3/Module3_composite_vector.pdf"
}

fig4() {
  say "Fig 4 — Mendelian randomization"
  # panel a-right: the MR evidence DAG fragment, wrapped with its shared style
  local schem="$HEAP_ROOT/mr_schematics" tmp
  tmp="$(mktemp -d)"
  cp "$schem/style.tex" "$schem/mr_evidence_dag.tex" "$tmp/"
  printf '%s\n' '\documentclass[border=8pt]{standalone}' '\input{style.tex}' \
    '\begin{document}' '\input{mr_evidence_dag.tex}' '\end{document}' > "$tmp/fig_mr_evidence_dag.tex"
  schematic "$tmp/fig_mr_evidence_dag.tex" "$schem" 300
  rm -rf "$tmp"
  cell "$V/build_mr_tier_ladder.R"
  r    "$V/build_mr_panelb_folded.R"
  r    "$V/build_mr_panelc.R"
  r    "$V/build_mr_paneld.R"
  cell "$V/figures/fig_mr_main_mediators_dag.R"
  r    "$V/build_mr_composite_vector.R"
  compile_here "$EXP/module5" Module5_composite_vector
  publish 4 "$EXP/module5/Module5_composite_vector.pdf"
}

fig5() {
  say "Fig 5 — intervention trials"
  schematic "$V/figures/fig_m4_schematic_compact.tex" "$EXP/module4" 300
  cell "$V/figures/fig_m4_panel_b.R"
  cell "$V/figures/fig_m4_panel_d.R"
  cell "$V/figures/fig_m4_shared_network.R"
  r    "$V/build_module4_composite_v4.R"
  publish 5 "$EXP/module4/Fig5_intervention_composite_v4.pdf"
}

fig6() {
  say "Fig 6 — proteome-based exposure scores"
  schematic "$V/figures/fig_module6_schematic_compact.tex" "$EXP/module6" 200
  cell "$V/figures/fig_m6_panel_b.R"
  cell "$V/figures/fig_m6_panel_c.R"
  cell "$V/figures/fig_m6_panel_d.R"
  r    "$V/build_module6_composite_landscape.R"
  publish 6 "$EXP/module6/Module6_composite_landscape.pdf"
}

# Supplementary figures, in the order they appear in the supplement. Each script
# writes $FIG/supplement/<module>/<id>.{pdf,png} and its plotted data to
# $FIG/data/<module>/<id>.tsv.
SUPP=(
  fig_greml_vs_r2 fig_tissue_enrichment fig_pathway_enrichment
  fig_pes_incremental_value fig_pes_covariate_robustness fig_instrument_diagnostics
  fig_gwas_exemplars fig_ldsc_rg fig_category_biology
  fig_traintest_stability_categories fig_exwas_miami fig_gxe_architecture
  fig_pathway_themes fig_tissue_themes fig_mediation_flows
  fig_mediation_prioritization fig_mediation_proportion fig_mr_coloc
  fig_mr_motif_overview fig_mr_refines_mediation fig_mr_attrition
  fig_pes_panel_size fig_pes_specificity fig_pes_tracking_design
  fig_exposure_cormatrix fig_dose_response fig_mr_rigor
  fig_module2_spec_sensitivity fig_module3_spec_sensitivity fig_pes_imaging_tracking
  fig_module2_gxe_spec_sensitivity fig_variance_architecture fig_pathways_by_component
  fig_gxe_noise_floor fig_module1_robustness
)
supp() {
  say "Supplementary figures (${#SUPP[@]})"
  local id failed=()
  for id in "${SUPP[@]}"; do r "$V/figures/$id.R" || failed+=("$id"); done
  [ ${#failed[@]} -eq 0 ] || { echo "FAILED: ${failed[*]}" >&2; return 1; }
}

main() { fig1; fig2; fig3; fig4; fig5; fig6; }

[ $# -eq 0 ] && set -- main supp
for target in "$@"; do
  case "$target" in
    main|supp|fig1|fig2|fig3|fig4|fig5|fig6) "$target" ;;
    *) echo "unknown target: $target (main | supp | fig1..fig6)" >&2; exit 2 ;;
  esac
done
say "done"
