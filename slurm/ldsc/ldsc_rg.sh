#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-08:00
#SBATCH --mem=24G
#SBATCH -p short
###############################################################################
# ldsc_rg.sh — pairwise genetic correlation (rg) between exposure GWAS via LDSC.
#
# Consumes the munged sumstats produced by the h2 stage
# (${IGLOO}/UKB/HEAP/output/gwas/ldsc/munged/<exposure>.sumstats.gz), so run the
# h2 array first. For each exposure i (in list order) it runs one ldsc.py --rg
# call of  exposure_i  vs every LATER exposure (upper triangle — each unordered
# pair computed once), writing one log per root exposure:
#     ${OUT_ROOT}/rg/<exposure_i>.rg.log
# The rg estimate + s.e. + p live in the "Summary of Genetic Correlation Results"
# table at the bottom of each log.
#
# Which exposures to correlate:
#   slurm/ldsc/ldsc_rg_exposures.txt   (one exposure id per line) if present,
#   otherwise ALL munged exposures found in the munged dir.
# rg is noisy for low-heritability traits, so curating that list is recommended.
#
# Submit:   sbatch slurm/ldsc/ldsc_rg.sh
# Override resources/paths with the LDSC_* env vars (see run_ldsc_h2.sh).
###############################################################################
set -euo pipefail
umask 0002

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
IGLOO="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}"
LDSC_ENV="${LDSC_ENV:-${IGLOO}/envs/ldsc_env}"
[[ -x "${LDSC_ENV}/bin/python" ]] || LDSC_ENV="${CONDA_ENVS:-$HOME/.conda/envs}/ldsc"
LDSC_HOME="${LDSC_HOME:-${IGLOO}/LDSC/ldsc}"
LDSC_LD_DIR="${LDSC_LD_DIR:-${IGLOO}/LDSC/eur_w_ld_chr}"
OUT_ROOT="${LDSC_OUT_ROOT:-${IGLOO}/UKB/HEAP/output/gwas/ldsc}"
PY="${LDSC_ENV}/bin/python"

MUNGED_DIR="${OUT_ROOT}/munged"
RG_DIR="${OUT_ROOT}/rg"
mkdir -p "${RG_DIR}" "/n/groups/patel/IGLOO/UKB/HEAP/logs/ldsc"

[[ -x "${PY}" ]] || { echo "ERROR: python not found in LDSC_ENV=${LDSC_ENV}" >&2; exit 1; }
[[ -d "${MUNGED_DIR}" ]] || { echo "ERROR: munged dir not found: ${MUNGED_DIR} (run h2 array first)" >&2; exit 1; }

# --- assemble the exposure list ----------------------------------------------
LIST="${LDSC_RG_LIST:-${HEAP_ROOT}/slurm/ldsc/ldsc_rg_exposures.txt}"
declare -a EXPS=()
if [[ -f "${LIST}" ]]; then
  while IFS= read -r line; do
    e=$(echo "${line}" | tr -d '[:space:]'); [[ -z "${e}" ]] && continue
    if [[ -f "${MUNGED_DIR}/${e}.sumstats.gz" ]]; then EXPS+=("${e}"); else
      echo "WARN: no munged sumstats for '${e}', skipping" >&2
    fi
  done < "${LIST}"
else
  echo "No ${LIST}; using ALL munged exposures."
  for f in "${MUNGED_DIR}"/*.sumstats.gz; do
    [[ -e "${f}" ]] || continue
    EXPS+=("$(basename "${f}" .sumstats.gz)")
  done
fi

N=${#EXPS[@]}
if [[ "${N}" -lt 2 ]]; then
  echo "ERROR: need >=2 munged exposures for rg; found ${N}." >&2
  exit 1
fi
echo "[ldsc rg] ${N} exposures; computing upper-triangle genetic correlations."

# --- upper-triangle rg: exposure_i vs exposures_{i+1..N} in one call each -----
for (( i=0; i<N-1; i++ )); do
  root="${EXPS[$i]}"
  rest=()
  for (( j=i+1; j<N; j++ )); do rest+=("${MUNGED_DIR}/${EXPS[$j]}.sumstats.gz"); done
  others=$(IFS=,; echo "${rest[*]}")
  echo ">>> rg root [${i}/${N}]: ${root}  vs $((N-i-1)) others"
  "${PY}" "${LDSC_HOME}/ldsc.py" \
    --rg "${MUNGED_DIR}/${root}.sumstats.gz,${others}" \
    --ref-ld-chr "${LDSC_LD_DIR}/" \
    --w-ld-chr "${LDSC_LD_DIR}/" \
    --out "${RG_DIR}/${root}.rg" || echo "WARN: rg failed for root ${root}" >&2
  chmod g+w "${RG_DIR}/${root}.rg.log" 2>/dev/null || true
done

echo "[ldsc rg] done -> ${RG_DIR}/"
