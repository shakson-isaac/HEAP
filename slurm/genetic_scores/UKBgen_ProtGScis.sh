#!/bin/bash
#SBATCH -c 1                              # 1 cores
#SBATCH -t 0-4:00                         # Runtime of 8 hours, in D-HH:MM format
#SBATCH --mem=10G                          # Memory total in 10 GB
#SBATCH -p short                          # Run in a partition: short
#SBATCH --array=1-2704

HEAP_ROOT="${HEAP_ROOT:-/n/groups/patel/shakson_ukb/HEAP}"
umask 0002   # group-writable outputs for hpc_patel team runs
LEGACY_UKB_ROOT="${HEAP_LEGACY_UKB_ROOT:-/n/groups/patel/shakson_ukb/UK_Biobank}"
SCRATCH_ROOT="${HEAP_SCRATCH_ROOT:-/n/scratch/users/${USER:0:1}/${USER}}"
IGLOO_ROOT="${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}"
NO_BACKUP_UKB_ROOT="${HEAP_NO_BACKUP_UKB_ROOT:-/n/no_backup2/patel/uk_biobank}"
NO_BACKUP_IGLOO_ROOT="${HEAP_NO_BACKUP_IGLOO_ROOT:-/n/no_backup2/patel/IGLOO}"
export HEAP_ROOT LEGACY_UKB_ROOT SCRATCH_ROOT IGLOO_ROOT NO_BACKUP_UKB_ROOT NO_BACKUP_IGLOO_ROOT
mkdir -p "${HEAP_IGLOO_ROOT:-/n/groups/patel/IGLOO}/UKB/HEAP/logs/genetic_scores"

#NOTES
#2704 proteins for UKB multi-ancestry:
#1G RAM is NOT enough caused a memory error w/ PLINK2
#2G is NOT enough memory for 4-5k variant PGS calculation
#12 hour time was to much limit this to 3-4 hours.

#Load modules
module load gcc/9.2.0
module load plink2/2.0.20220814
module load sqlite3/3.43.1

#give enough temp file space
ulimit -n 10000

#Obtain filename for each job array:
fn=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${HEAP_ROOT}/config/genetic_scores/UKBprotGScis.txt" | tr -d '"')
protID=$(echo "${fn}" | cut -d'_' -f1)

# Specify file directories:
# Genotype bgen: canonical cold-storage copy (IGLOO), fall back to scratch.
BGEN_FILE="${NO_BACKUP_IGLOO_ROOT}/UKB/Genetics/UKBallchr.bgen"
[ -f "${BGEN_FILE}" ] || BGEN_FILE="${SCRATCH_ROOT}/UKB_intermediate/Genetics/UKBallchr.bgen"
BGEN_SAMP="${NO_BACKUP_UKB_ROOT}/ukb_genetics/52887/project_52887_genetics/ukb52887_imp_chr1_v3_s487296.sample"
# OmicsPred cis/trans SNP scores: prefer shared IGLOO copy, fall back to legacy.
OMICSPRED_ROOT="${IGLOO_ROOT}/UKB/OMICSPRED"
[ -d "${OMICSPRED_ROOT}" ] || OMICSPRED_ROOT="${LEGACY_UKB_ROOT}/Data/OMICSPRED"
SCORE_FILE="${OMICSPRED_ROOT}/cis_trans_snps/cis/${fn}"
OUTPUT_FILE="${IGLOO_ROOT}/UKB/ProtGScis/${protID}" 


# Temporary files for RSID, CHRPOS, BGEN_FILT, NEW_BGEN, and TEMP_ModelFile:
RSID_FILE=$(mktemp)
CHRPOS_FILE=$(mktemp)
BGEN_FILT=$(mktemp)
NEW_BGEN=$(mktemp)
TEMP_MODEL_FILE=$(mktemp)


# Run bgenix command with variables
# bgenix tool: prefer shared IGLOO bin, fall back to legacy UK_Biobank bin.
BGENIX_DIR="${IGLOO_ROOT}/bin"
[ -x "${BGENIX_DIR}/bgenix" ] || BGENIX_DIR="${LEGACY_UKB_ROOT}/bin"
cd "${BGENIX_DIR}"


##### CODE to FILTER OUT VARIANTS FROM SCORE FILE: ####
# Extract rsID:
awk -F'\t' 'NR>1 && $2 ~ /^[0-9]+$/ { print $1 }' "$SCORE_FILE" > "$RSID_FILE"


# Extract CHR POS:
awk -F'\t' 'NR>1 && $2 ~ /^[0-9]+$/ { print sprintf("%02d", $2)":"$3"-"$3 }' "$SCORE_FILE" > "$CHRPOS_FILE"


# Cleaned Model File: Remove comment lines and save to temp file
sed '/^#/d' "$SCORE_FILE" > "$TEMP_MODEL_FILE"


#Method: Filter SNPs from UKB allchr file based of score file
#Remember to put a space in cmd output
./bgenix \
    -g "$BGEN_FILE" \
    -incl-rsids "$RSID_FILE" \
    -incl-range "$CHRPOS_FILE" \
    > "$BGEN_FILT"


# Write index file for bgen (bgen.bgi)
./bgenix \
    -g "$BGEN_FILT"\
    -index -clobber


##### CODE to Align and Account for Multiallelic SNPs ####
# Join Variant IDs to develop score properly:
# REMEMBER: UKB bgen files: have 'allele1' column being the reference allele

# Rename variables for clarity:
BGEN_BGI="${BGEN_FILT}.bgi"

# Import the betas into the sqlite database as a table called Betas
sqlite3 "$BGEN_BGI" "DROP TABLE IF EXISTS Betas;"
sqlite3 -separator $'\t' "$BGEN_BGI" ".import $TEMP_MODEL_FILE Betas"


sqlite3 "$BGEN_BGI" "DROP TABLE IF EXISTS Joined;"
# And inner join it to the index table (Variants), making a new table (Joined)
# By joining on alleles as well as chromosome and position 
# we can ensure only the relevant alleles from any multi-allelic SNPs are retained
sqlite3 -header "$BGEN_BGI" \
"CREATE TABLE Joined AS 
  SELECT Variant.*, Betas.chr_name, Betas.effect_weight FROM Variant INNER JOIN Betas 
    ON Variant.chromosome = printf('%02d', Betas.chr_name) 
    AND Variant.position = Betas.chr_position 
    AND Variant.allele1 = Betas.other_allele 
    AND Variant.allele2 = Betas.effect_allele 
  UNION 
  SELECT Variant.*, Betas.chr_name, -Betas.effect_weight FROM Variant INNER JOIN Betas 
    ON Variant.chromosome = printf('%02d', Betas.chr_name) 
    AND Variant.position = Betas.chr_position 
    AND Variant.allele1 = Betas.effect_allele AND 
    Variant.allele2 = Betas.other_allele;"

echo "SQLITE database aligned variants properly: ${protID}"


# Filter the .bgen file to include only the alleles specified in the Betas for each SNP 
./bgenix \
    -g "$BGEN_FILT" \
    -table Joined  \
    > "$NEW_BGEN"

# And produce an index file for the new .bgen
./bgenix \
    -g "$NEW_BGEN" \
    -index -clobber

echo "new bgen file produces and ready to use: ${protID}"


#### CODE to Output Polygenic Scores ####
# Run plink on filtered bgen, original sample file, and protein or metabolite score.
plink2 \
  --bgen "$NEW_BGEN" 'ref-first'\
  --sample "$BGEN_SAMP" \
  --score  "$TEMP_MODEL_FILE" 1 4 6 header list-variants cols=scoresums \
  --out "$OUTPUT_FILE"

# Clean up the temporary files
rm "$TEMP_MODEL_FILE"
rm "$RSID_FILE"
rm "$CHRPOS_FILE"
rm "$BGEN_FILT"
rm "$NEW_BGEN"

echo "Done with example polygenic score: ${protID}"
