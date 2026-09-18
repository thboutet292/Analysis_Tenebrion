#!/bin/bash
#SBATCH --job-name=build_matam_db
#SBATCH --output=log/BUILD_MATAM_DB_%j.out
#SBATCH --error=log/BUILD_MATAM_DB_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=long,lmge
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/home/thboutet/projets/16_tenebrion

#######################################################################################
#   GOAL : Préparation de la base SILVA 138.2 avec l'outil OFFICIEL de MATAM
#######################################################################################

set -euo pipefail

# =====================================================================
# 1. CONFIGURATION
# =====================================================================
module purge
module load rclone/1.55.1

DB_URL="https://www.arb-silva.de/fileadmin/silva_databases/release_138.2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"
DB_FASTA_GZ="SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
CONTAINERS_DIR="${PROJECT_DIR}/containers"
# Utilisation du conteneur MATAM pour avoir accès à matam_db_preprocessing.py
MATAM_SIF="${CONTAINERS_DIR}/matam_1.6.1--hdcf5f25_1.sif"
SINGULARITY="/usr/bin/singularity"

S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
# Nouveau dossier pour différencier de la base manuelle
S3_DEST="${S3_REMOTE}:${S3_BUCKET}/resources/matam_db_138_2_official"

mkdir -p "${CONTAINERS_DIR}" "log"

# =====================================================================
# 2. GESTION DU SCRATCH
# =====================================================================
SCRATCH_BASE="/storage/scratch/${USER}"
SCRATCH_JOB="${SCRATCH_BASE}/build_matam_db_${SLURM_JOB_ID}"
WORKDIR="${SCRATCH_JOB}/work"
TMP_WORK="${SCRATCH_JOB}/tmp"
MATAM_DB_OUT="${WORKDIR}/matam_db_138_2"

mkdir -p "${WORKDIR}" "${TMP_WORK}" "${MATAM_DB_OUT}"
export TMPDIR="${TMP_WORK}"
export SINGULARITY_TMPDIR="${TMP_WORK}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Nettoyage du scratch local..."
    rm -rf "${SCRATCH_JOB}" 2>/dev/null || true
}
trap cleanup EXIT

# =====================================================================
# 3. TÉLÉCHARGEMENT DE SILVA 138.2 NR99
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ÉTAPE 1] Téléchargement de SILVA..."
cd "${WORKDIR}"
wget --tries=3 --waitretry=5 -O "${DB_FASTA_GZ}" "${DB_URL}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Décompression..."
gunzip "${DB_FASTA_GZ}"
RAW_FASTA="SILVA_138.2_SSURef_NR99_tax_silva.fasta"

# =====================================================================
# 4. PRÉPARATION OFFICIELLE MATAM
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ÉTAPE 2] Préparation avec matam_db_preprocessing.py..."

# Le script officiel gère le clustering et l'indexation SortMeRNA automatiquement
${SINGULARITY} exec --cleanenv \
    --env "LC_ALL=C" --env "LANG=C" --env "TMPDIR=/work/tmp" \
    -B "${WORKDIR}:/work" -B "${TMP_WORK}:/work/tmp" \
    "${MATAM_SIF}" bash -lc "
        matam_db_preprocessing.py \
            -i '/work/${RAW_FASTA}' \
            -d '/work/matam_db_138_2' \
            --cpu ${SLURM_CPUS_PER_TASK} \
            --max_memory 60000 \
            -v
    "

echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Base de données prête. Contenu :"
ls -lh "${MATAM_DB_OUT}"

# =====================================================================
# 5. TRANSFERT S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ÉTAPE 3] Transfert vers S3 (${S3_DEST})..."
rclone sync "${MATAM_DB_OUT}" "${S3_DEST}" -P --transfers=4 --checkers=8

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] OPÉRATION TERMINÉE AVEC SUCCÈS"
echo "====================================================================="