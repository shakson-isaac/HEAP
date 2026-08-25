#!/bin/bash
###############################################################################
# run_ldsc_h2.sh  —  LD Score Regression heritability + intercept for ONE exposure
#
# Runs the three-step LDSC pipeline on a HEAP exposure GWAS (REGENIE step 2):
#   1. preprocess : regenie summary stats -> LDSC-ready columns (SNP A1 A2 N BETA SE P)
#   2. munge      : ldsc/munge_sumstats.py (HapMap3 SNP filter + allele merge)
#   3. h2         : ldsc/ldsc.py --h2  -> SNP heritability AND the LDSC intercept
#                   (the intercept is the model-based genomic-inflation estimate:
#                    it separates true polygenicity from confounding/structure,
#                    unlike lambda_GC which conflates the two).
#
# Usage:
#   run_ldsc_h2.sh <exposure_id>
#   e.g. run_ldsc_h2.sh alcohol_intake_frequency_f1558_0_0
#
# All locations default to the shared IGLOO copies so a teammate can run it as-is;
# every path is overridable via the environment variables below.
#
#   LDSC_ENV          conda env prefix with python2.7 + LDSC deps
#                       default: /n/groups/patel/IGLOO/envs/ldsc_env
#                       (falls back to ${CONDA_ENVS:-$HOME/.conda/envs}/ldsc)
#   LDSC_HOME         the LDSC package (ldsc.py, munge_sumstats.py)
#                       default: /n/groups/patel/IGLOO/LDSC/ldsc
#   LDSC_LD_DIR       --ref-ld-chr / --w-ld-chr reference (eur_w_ld_chr)
#                       default: /n/groups/patel/IGLOO/LDSC/eur_w_ld_chr
#   LDSC_SNPLIST      --merge-alleles HapMap3 snplist
#                       default: ${LDSC_LD_DIR}/w_hm3.snplist
#                       (falls back to decompressing LDSCORE-w_hm3.snplist.bz2)
#   HEAP_REGENIE_STEP2_DIR  exposure GWAS root
#                       default: /n/groups/patel/IGLOO/UKB/HEAP/output/gwas/regenie_step2
#   LDSC_OUT_ROOT     LDSC output root
#                       default: /n/groups/patel/IGLOO/UKB/HEAP/output/gwas/ldsc
#   LDSC_FORCE=1      recompute even if the h2 log already exists (default: skip)
#   LDSC_KEEP_TMP=1   keep the node-local working dir (default: clean on exit)
###############################################################################
set -euo pipefail
umask 0002   # group-writable outputs for hpc_patel team runs

PHENO="${1:-}"
if [[ -z "${PHENO}" ]]; then
  echo "Usage: $0 <exposure_id>" >&2
  exit 2
fi

# --- Resolve resources (shared IGLOO defaults, all overridable) --------------
IGLOO="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}"
LDSC_ENV="${LDSC_ENV:-${IGLOO}/envs/ldsc_env}"
[[ -x "${LDSC_ENV}/bin/python" ]] || LDSC_ENV="${CONDA_ENVS:-$HOME/.conda/envs}/ldsc"
LDSC_HOME="${LDSC_HOME:-${IGLOO}/LDSC/ldsc}"
LDSC_LD_DIR="${LDSC_LD_DIR:-${IGLOO}/LDSC/eur_w_ld_chr}"
LDSC_SNPLIST="${LDSC_SNPLIST:-${LDSC_LD_DIR}/w_hm3.snplist}"
STEP2_DIR="${HEAP_REGENIE_STEP2_DIR:-${IGLOO}/UKB/HEAP/output/gwas/regenie_step2}"
OUT_ROOT="${LDSC_OUT_ROOT:-${IGLOO}/UKB/HEAP/output/gwas/ldsc}"

PY="${LDSC_ENV}/bin/python"   # call the env python directly (no conda activation needed)

# --- Sanity checks -----------------------------------------------------------
[[ -x "${PY}" ]]                         || { echo "ERROR: python not found in LDSC_ENV=${LDSC_ENV}" >&2; exit 1; }
[[ -f "${LDSC_HOME}/ldsc.py" ]]          || { echo "ERROR: ldsc.py not found in LDSC_HOME=${LDSC_HOME}" >&2; exit 1; }
[[ -f "${LDSC_HOME}/munge_sumstats.py" ]]|| { echo "ERROR: munge_sumstats.py not found in LDSC_HOME=${LDSC_HOME}" >&2; exit 1; }
[[ -d "${LDSC_LD_DIR}" ]]                || { echo "ERROR: LD reference dir not found: ${LDSC_LD_DIR}" >&2; exit 1; }

# Source regenie file: canonical <exp>.regenie, then legacy doubled name.
SRC="${STEP2_DIR}/${PHENO}/${PHENO}.regenie"
if [[ ! -f "${SRC}" ]]; then
  LEGACY="${STEP2_DIR}/${PHENO}/regenie_step2_${PHENO}_${PHENO}.regenie"
  if [[ -f "${LEGACY}" ]]; then SRC="${LEGACY}"; else
    echo "ERROR: exposure GWAS not found for '${PHENO}':" >&2
    echo "       ${SRC}" >&2
    echo "       (run the gwas_regenie stage for this exposure first)" >&2
    exit 1
  fi
fi

# HapMap3 snplist: use the plain file if present, else decompress the bz2 once.
if [[ ! -f "${LDSC_SNPLIST}" ]]; then
  BZ2="${IGLOO}/LDSC/LDSCORE-w_hm3.snplist.bz2"
  if [[ -f "${BZ2}" ]]; then
    LDSC_SNPLIST="${IGLOO}/LDSC/LDSCORE-w_hm3.snplist"
    [[ -f "${LDSC_SNPLIST}" ]] || bunzip2 -kc "${BZ2}" > "${LDSC_SNPLIST}"
  else
    echo "ERROR: HapMap3 snplist not found (set LDSC_SNPLIST): ${LDSC_SNPLIST}" >&2
    exit 1
  fi
fi

# --- Output layout -----------------------------------------------------------
MUNGED_DIR="${OUT_ROOT}/munged"
H2_DIR="${OUT_ROOT}/h2"
mkdir -p "${MUNGED_DIR}" "${H2_DIR}"

MUNGED="${MUNGED_DIR}/${PHENO}.sumstats.gz"   # reused by ldsc_rg.sh (genetic correlation)
H2_PREFIX="${H2_DIR}/${PHENO}"                # ldsc.py writes ${H2_PREFIX}.log

# Idempotent: skip a completed exposure unless LDSC_FORCE=1.
if [[ -f "${H2_PREFIX}.log" && "${LDSC_FORCE:-0}" != "1" ]]; then
  echo "[${PHENO}] h2 log already exists (${H2_PREFIX}.log); skipping (LDSC_FORCE=1 to redo)."
  exit 0
fi

# Node-local working dir for the preprocessed intermediate (cleaned on exit).
WORKDIR="${SLURM_TMPDIR:-/tmp}/ldsc_${PHENO}_${SLURM_JOB_ID:-$$}"
mkdir -p "${WORKDIR}"
cleanup() { [[ "${LDSC_KEEP_TMP:-0}" == "1" ]] || rm -rf "${WORKDIR}"; }
trap cleanup EXIT
PRE="${WORKDIR}/${PHENO}.pre.sumstats.gz"

echo "============================================================"
echo " LDSC h2  |  exposure: ${PHENO}"
echo "   src     : ${SRC}"
echo "   env     : ${LDSC_ENV}"
echo "   ldsc    : ${LDSC_HOME}"
echo "   ld ref  : ${LDSC_LD_DIR}"
echo "   snplist : ${LDSC_SNPLIST}"
echo "   out     : ${OUT_ROOT}"
echo "============================================================"

# --- Step 1: preprocess regenie -> LDSC columns ------------------------------
# regenie step2 cols: 1 CHROM 2 GENPOS 3 ID 4 ALLELE0 5 ALLELE1 6 A1FREQ 7 INFO
#                     8 N 9 TEST 10 BETA 11 SE 12 CHISQ 13 LOG10P 14 EXTRA
# LDSC A1 = effect (tested) allele = regenie ALLELE1 ($5); A2 = ALLELE0 ($4).
# P reconstructed from LOG10P: P = 10^(-LOG10P) = exp(-LOG10P*ln10).
echo ">>> [${PHENO}] Step 1/3: preprocessing regenie -> ${PRE}"
READER="cat"; gzip -t "${SRC}" >/dev/null 2>&1 && READER="zcat"
${READER} "${SRC}" | awk 'BEGIN{OFS="\t"}
  NR==1 { print "SNP","A1","A2","N","BETA","SE","P"; next }
  $13!="NA" && $13!="" {
    p = exp(-$13*log(10));
    print $3, $5, $4, $8, $10, $11, p
  }' | gzip -c > "${PRE}"
echo "[${PHENO}] preprocessed ($(zcat "${PRE}" | wc -l) lines incl header)"

# --- Step 2: munge -----------------------------------------------------------
echo ">>> [${PHENO}] Step 2/3: munge_sumstats.py"
"${PY}" "${LDSC_HOME}/munge_sumstats.py" \
  --sumstats "${PRE}" \
  --snp SNP --a1 A1 --a2 A2 --p P --N-col N \
  --signed-sumstats BETA,0 \
  --merge-alleles "${LDSC_SNPLIST}" \
  --chunksize 500000 \
  --out "${MUNGED_DIR}/${PHENO}"
chmod g+w "${MUNGED}" "${MUNGED_DIR}/${PHENO}.log" 2>/dev/null || true
echo "[${PHENO}] munged -> ${MUNGED}"

# --- Step 3: LDSC heritability + intercept -----------------------------------
echo ">>> [${PHENO}] Step 3/3: ldsc.py --h2"
"${PY}" "${LDSC_HOME}/ldsc.py" \
  --h2 "${MUNGED}" \
  --ref-ld-chr "${LDSC_LD_DIR}/" \
  --w-ld-chr "${LDSC_LD_DIR}/" \
  --out "${H2_PREFIX}"
chmod g+w "${H2_PREFIX}.log" 2>/dev/null || true

echo "[${PHENO}] DONE -> ${H2_PREFIX}.log"
echo "----- LDSC h2 summary -----"
grep -E "Total Observed scale h2|Lambda GC|Mean Chi\^2|Intercept|Ratio" "${H2_PREFIX}.log" || true
