#!/bin/bash
#SBATCH --job-name=CHECKM2
#SBATCH --output=log/PHASE4_CHECKM2_%j.out
#SBATCH --error=log/PHASE4_CHECKM2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal,lmge
#SBATCH --time=12:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   shotgun_checkm2.slurm
#
#   GOAL    :   PHASE 4 - CONTRÔLE QUALITÉ DES MAGs (CheckM2).
#               Évalue la complétude et la contamination des génomes générés par
#               les TROIS binners (MetaBAT2, SemiBin2, CONCOCT) via des modèles
#               de Machine Learning, en une seule soumission SLURM.
#
#   BDD :   La base de données (~3.5 Go) est rapatriée depuis le S3
#               vers le scratch avant le lancement, puis supprimée après.
#               Elle est téléchargée UNE SEULE FOIS et réutilisée pour les 3 runs.
#               Lancer build_checkm2_db.slurm une première fois avant.
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] QC DES MAGs (CheckM2) -- 3 BINNERS"
echo "Job ID : ${SLURM_JOB_ID}"
echo "====================================================================="

module purge
module load rclone/1.55.1
module load apptainer 2>/dev/null || true

# =====================================================================
# CONFIGURATION DES CHEMINS
# =====================================================================
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

IMG_CHECKM2="${PROJECT_DIR}/containers/checkm2_v1.0.2.sif"

S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_DB="s3_uca:lmge-tenebrion/ref/checkm2_db"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p log/

# Association nom du binner -> dossier contenant les bins .fa
declare -A BIN_DIRS
BIN_DIRS["metabat2"]="${PROJECT_DIR}/results/metabat2_bins"
BIN_DIRS["semibin2"]="${PROJECT_DIR}/results/semibin2_bins/bins_fasta"
BIN_DIRS["concoct"]="${PROJECT_DIR}/results/concoct_bins/bins_fasta"

DIR_CHECKM2_ROOT="${PROJECT_DIR}/results/checkm2_qc"
mkdir -p "${DIR_CHECKM2_ROOT}"

# =====================================================================
# ÉTAPE 1 : RAPATRIEMENT DE LA BASE DEPUIS LE S3 (UNE SEULE FOIS)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Rapatriement de la base CheckM2 depuis le S3..."

JOB_SCRATCH="/storage/scratch/${USER}/checkm2_${SLURM_JOB_ID}"
DB_LOCAL="${JOB_SCRATCH}/CheckM2_database"
mkdir -p "${DB_LOCAL}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression du scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

rclone copy "${S3_DB}" "${DB_LOCAL}" --transfers 4 -P

# Détection automatique du fichier .dmnd (nom peut varier selon la version)
DB_FILE=$(find "${DB_LOCAL}" -name "*.dmnd" | head -1)
if [ -z "${DB_FILE}" ]; then
    echo "[ERREUR] Fichier .dmnd introuvable après rapatriement depuis le S3." >&2
    echo "         Vérifier que build_checkm2_db.slurm a bien été lancé." >&2
    exit 1
fi
echo " -> Base disponible : ${DB_FILE}"

export TMPDIR="${JOB_SCRATCH}/tmp"
mkdir -p "${TMPDIR}"
export APPTAINER_TMPDIR="${TMPDIR}"

# =====================================================================
# ÉTAPE 2 : BOUCLE CHECKM2 SUR LES TROIS JEUX DE BINS
# =====================================================================
NB_TOOLS_RUN=0

for TOOL in "${!BIN_DIRS[@]}"; do
    DIR_BINS="${BIN_DIRS[${TOOL}]}"
    DIR_OUT="${DIR_CHECKM2_ROOT}/${TOOL}"

    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${TOOL}] Vérification des bins dans ${DIR_BINS}..."

    NB_BINS=$(ls -1 "${DIR_BINS}/"*.fa 2>/dev/null | wc -l || echo 0)
    if [ "${NB_BINS}" -eq 0 ]; then
        echo "[AVERTISSEMENT] Aucun bin trouvé pour ${TOOL} dans ${DIR_BINS}. Ce binner est ignoré." >&2
        continue
    fi
    echo " -> ${NB_BINS} MAGs détectés pour ${TOOL}."

    mkdir -p "${DIR_OUT}"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${TOOL}] Exécution CheckM2 (Diamond BlastP)..."

    apptainer exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
        --bind "${TMPDIR}:/tmp" \
        "${IMG_CHECKM2}" \
        checkm2 predict \
            --threads "${THREADS}" \
            --input "${DIR_BINS}" \
            --extension fa \
            --output-directory "${DIR_OUT}" \
            --database_path "${DB_FILE}"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${TOOL}] Rapport généré : ${DIR_OUT}/quality_report.tsv"
    NB_TOOLS_RUN=$((NB_TOOLS_RUN + 1))
done

if [ "${NB_TOOLS_RUN}" -eq 0 ]; then
    echo "[ERREUR] Aucun des trois binners n'a produit de bins exploitables." >&2
    exit 1
fi

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FIN DU CONTRÔLE QUALITÉ (${NB_TOOLS_RUN}/3 binners traités)."
echo " -> Rapports : ${DIR_CHECKM2_ROOT}/<binner>/quality_report.tsv"
echo "====================================================================="