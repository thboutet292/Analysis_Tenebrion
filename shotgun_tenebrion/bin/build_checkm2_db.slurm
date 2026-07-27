#!/bin/bash
#SBATCH --job-name=BUILD_CHECKM2_DB
#SBATCH --output=log/BUILD_CHECKM2_DB_%j.out
#SBATCH --error=log/BUILD_CHECKM2_DB_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal,lmge
#SBATCH --time=02:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   build_checkm2_db.slurm
#
#   GOAL    :   Téléchargement de la base de données CheckM2 (uniref100.KO.dmnd)
#               via le binaire interne à l'image Apptainer, puis upload sur le S3.
#               À ne lancer qu'une seule fois.
#
#   USAGE   :   sbatch build_checkm2_db_s3.slurm
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] TÉLÉCHARGEMENT BASE CHECKM2"
echo "Job ID : ${SLURM_JOB_ID}"
echo "====================================================================="

module purge
module load rclone/1.55.1
module load apptainer 2>/dev/null || true

# =====================================================================
# CONFIGURATION
# =====================================================================
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
IMG_CHECKM2="${PROJECT_DIR}/containers/checkm2_v1.0.2.sif"

S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_DEST="${S3_REMOTE}:${S3_BUCKET}/ref/checkm2_db"

# Espace de travail scratch (la base fait ~3.5 Go)
JOB_SCRATCH="/storage/scratch/${USER}/checkm2_db_${SLURM_JOB_ID}"
DB_LOCAL="${JOB_SCRATCH}/CheckM2_database"
mkdir -p "${DB_LOCAL}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression du scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 1 : TÉLÉCHARGEMENT VIA CHECKM2
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Téléchargement de la base (~3.5 Go)..."

apptainer exec \
    --writable-tmpfs \
    --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
    "${IMG_CHECKM2}" \
    checkm2 database --download \
        --path "${DB_LOCAL}" || true

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Téléchargement terminé."

# Vérification
DB_FILE=$(find "${DB_LOCAL}" -name "*.dmnd" | head -1)
if [ -z "${DB_FILE}" ]; then
    echo "[ERREUR] Fichier .dmnd introuvable après téléchargement." >&2
    exit 1
fi
echo " -> Base trouvée : ${DB_FILE} ($(du -sh ${DB_FILE} | cut -f1))"

# =====================================================================
# ÉTAPE 2 : UPLOAD VERS LE S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Upload vers le S3..."

rclone copy "${DB_LOCAL}" "${S3_DEST}" --transfers 4 -P

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Upload terminé."
echo " -> Chemin S3 : ${S3_DEST}"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] TERMINÉ. Base CheckM2 disponible sur le S3."
echo "====================================================================="