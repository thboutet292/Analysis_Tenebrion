#!/bin/bash
#SBATCH --job-name=KRAKEN_DUAL
#SBATCH --output=log/TAXONOMY_DUAL_%A_%a.out
#SBATCH --error=log/TAXONOMY_DUAL_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G                    
#SBATCH --partition=normal,lmge
#SBATCH --time=24:00:00

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   shotgun_kraken2_bracken_dual.slurm
#   TOOLS   :   rclone ; Apptainer (Kraken2 / Bracken)
#
#   GOAL    :   Profilage taxonomique sur DEUX bases de données (PlusPFP & Addon).
#               0. Importation des 2 bases Kraken2 depuis le S3 vers le Scratch.
#               1. Importation des FASTQ "clean" depuis le S3 vers le Scratch.
#               2. Assignation Kraken2 (Double passage).
#               3. Ré-estimation Bracken (Double passage).
#               4. Exportation vers le dossier local results/.
#######################################################################################

set -uo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] DÉMARRAGE DU PROFILAGE TAXONOMIQUE DOUBLE"
echo "Job Master ID : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}} | Task : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "====================================================================="

module purge
module load rclone/1.55.1

# =====================================================================
# CONFIGURATION
# =====================================================================
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_IN_DIR="shotgun_cleaned/fastq"

# Chemins des deux bases sur le S3
S3_DB_PFP="ref/kraken2_pluspfp"
S3_DB_ADDON="ref/kraken2_insect_addon"

PROJECT_DIR="${PWD}"
APPTAINER_BIN="/usr/bin/apptainer"
SIF_IMG="/home/thboutet/projets/16_tenebrion/containers/kraken2_bracken.sif"

# Création des 4 dossiers de résultats
LOCAL_OUT_DIR="${PROJECT_DIR}/results"
mkdir -p "${LOCAL_OUT_DIR}/kraken2_pluspfp" "${LOCAL_OUT_DIR}/bracken_pluspfp"
mkdir -p "${LOCAL_OUT_DIR}/kraken2_addon" "${LOCAL_OUT_DIR}/bracken_addon"

THREADS="${SLURM_CPUS_PER_TASK}"
READ_LEN=150

# =====================================================================
# GESTION DE L'ARRAY SLURM DYNAMIQUE
# =====================================================================
mapfile -t SAMPLES < <(rclone ls "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/" | awk '{print $2}' | grep "_clean_R1.fastq.gz" | sed 's/_clean_R1.fastq.gz//' | sort)

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] ${N_SAMPLES} échantillons détectés."
    SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | awk '{print $1}')
    sbatch --array=0-$((N_SAMPLES - 1)) "${SCRIPT_PATH}"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# =====================================================================
# ÉTAPE 0 : IMPORTATION DES DEUX BASES DEPUIS LE S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Préparation des bases sur le Scratch..."

# 0.A - BASE PLUSPFP
KRAKEN_DB_PFP="/storage/scratch/${USER}/kraken2_pfp_db_shared"
mkdir -p "${KRAKEN_DB_PFP}"
FLAG_PFP="${KRAKEN_DB_PFP}/.sync_done"

exec 9> "/storage/scratch/${USER}/kraken_pfp_sync.lock"
flock -x 9
if [ ! -f "${FLAG_PFP}" ]; then
    echo " -> Importation de la base PlusPFP..."
    rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_DB_PFP}/" "${KRAKEN_DB_PFP}/" -P --transfers 8
    touch "${FLAG_PFP}"
fi
flock -u 9

# 0.B - BASE ADDON
KRAKEN_DB_ADDON="/storage/scratch/${USER}/kraken2_addon_db_shared"
mkdir -p "${KRAKEN_DB_ADDON}"
FLAG_ADDON="${KRAKEN_DB_ADDON}/.sync_done"

exec 8> "/storage/scratch/${USER}/kraken_addon_sync.lock"
flock -x 8
if [ ! -f "${FLAG_ADDON}" ]; then
    echo " -> Importation de la base Insect Addon..."
    rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_DB_ADDON}/" "${KRAKEN_DB_ADDON}/" -P --transfers 8
    touch "${FLAG_ADDON}"
fi
flock -u 8

# =====================================================================
# GESTION DU SCRATCH ÉCHANTILLON
# =====================================================================
JOB_SCRATCH="/storage/scratch/${USER}/kraken_sample_${SAMPLE}_${SLURM_JOB_ID}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

cleanup() {
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Nettoyage de l'espace de calcul..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
    echo "====================================================================="
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 1 : IMPORTATION DES READS
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Importation de ${SAMPLE} depuis S3..."
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_clean_R1.fastq.gz" ./
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_clean_R2.fastq.gz" ./

# =====================================================================
# ÉTAPE 2 : ASSIGNATION KRAKEN2 (DOUBLE PASSAGE)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2A] Kraken2 sur PlusPFP..."
${APPTAINER_BIN} exec --bind "${KRAKEN_DB_PFP}:${KRAKEN_DB_PFP}" --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2 \
    --db "${KRAKEN_DB_PFP}" --threads "${THREADS}" --paired --gzip-compressed --use-names \
    --report "${SAMPLE}_pluspfp.kreport" --output /dev/null \
    "${SAMPLE}_clean_R1.fastq.gz" "${SAMPLE}_clean_R2.fastq.gz"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2B] Kraken2 sur Insect Addon..."
${APPTAINER_BIN} exec --bind "${KRAKEN_DB_ADDON}:${KRAKEN_DB_ADDON}" --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2 \
    --db "${KRAKEN_DB_ADDON}" --threads "${THREADS}" --paired --gzip-compressed --use-names \
    --report "${SAMPLE}_addon.kreport" --output /dev/null \
    "${SAMPLE}_clean_R1.fastq.gz" "${SAMPLE}_clean_R2.fastq.gz"

# =====================================================================
# ÉTAPE 3 : RÉ-ESTIMATION BRACKEN (DOUBLE PASSAGE)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3A] Bracken sur PlusPFP..."
${APPTAINER_BIN} exec --bind "${KRAKEN_DB_PFP}:${KRAKEN_DB_PFP}" --bind "${PWD}:${PWD}" "${SIF_IMG}" bracken \
    -d "${KRAKEN_DB_PFP}" -i "${SAMPLE}_pluspfp.kreport" -o "${SAMPLE}_pluspfp.bracken" \
    -w "${SAMPLE}_pluspfp_bracken_species.kreport" -r "${READ_LEN}" -l S

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3B] Bracken sur Insect Addon..."
# /!\ Cette commande peut échouer si bracken-build n'a pas été lancé sur la base Addon.
# L'option "|| true" permet au script de continuer même en cas d'erreur ici.
${APPTAINER_BIN} exec --bind "${KRAKEN_DB_ADDON}:${KRAKEN_DB_ADDON}" --bind "${PWD}:${PWD}" "${SIF_IMG}" bracken \
    -d "${KRAKEN_DB_ADDON}" -i "${SAMPLE}_addon.kreport" -o "${SAMPLE}_addon.bracken" \
    -w "${SAMPLE}_addon_bracken_species.kreport" -r "${READ_LEN}" -l S || true

# =====================================================================
# ÉTAPE 4 : EXPORTATION VERS LE /HOME
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 4] Sauvegarde des résultats..."

# Sauvegarde PlusPFP
cp "${SAMPLE}_pluspfp.kreport" "${LOCAL_OUT_DIR}/kraken2_pluspfp/"
cp "${SAMPLE}_pluspfp.bracken" "${LOCAL_OUT_DIR}/bracken_pluspfp/" 2>/dev/null || true
cp "${SAMPLE}_pluspfp_bracken_species.kreport" "${LOCAL_OUT_DIR}/bracken_pluspfp/" 2>/dev/null || true

# Sauvegarde Addon
cp "${SAMPLE}_addon.kreport" "${LOCAL_OUT_DIR}/kraken2_addon/"
cp "${SAMPLE}_addon.bracken" "${LOCAL_OUT_DIR}/bracken_addon/" 2>/dev/null || true
cp "${SAMPLE}_addon_bracken_species.kreport" "${LOCAL_OUT_DIR}/bracken_addon/" 2>/dev/null || true

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] Analyse taxonomique terminée pour ${SAMPLE}."
