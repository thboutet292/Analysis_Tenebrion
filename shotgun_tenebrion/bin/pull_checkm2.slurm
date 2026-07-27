#!/bin/bash
#SBATCH --job-name=pull_checkm2
#SBATCH --output=log/PULL_CHECKM2_%j.out
#SBATCH --error=log/PULL_CHECKM2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal,lmge
#SBATCH --time=02:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   pull_checkm2.slurm
#
#   GOAL    :   1. Téléchargement de l'image Singularity/Apptainer de CheckM2.
#               2. Téléchargement de la base de données Machine Learning (Diamond).
#
#               CheckM2 évalue la qualité (Complétude et Contamination) des MAGs
#               en utilisant des modèles prédictifs de Machine Learning.
#######################################################################################

set -euo pipefail

module purge
module load apptainer 2>/dev/null || true

cd "${SLURM_SUBMIT_DIR}"
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
CONTAINERS_DIR="${PROJECT_DIR}/containers"
DB_DIR="${PROJECT_DIR}/databases/CheckM2_database"

mkdir -p "${CONTAINERS_DIR}" "${DB_DIR}" log

# Isolation du cache Apptainer sur le scratch
SCRATCH_BASE="/storage/scratch/${USER}"
APPTAINER_CACHEDIR="${SCRATCH_BASE}/apptainer_cache_checkm2"
APPTAINER_TMPDIR="${SCRATCH_BASE}/apptainer_tmp_checkm2"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"
export APPTAINER_CACHEDIR
export APPTAINER_TMPDIR
export TMPDIR="${APPTAINER_TMPDIR}"

IMG_NAME="checkm2_v1.0.2.sif"
IMG_PATH="${CONTAINERS_DIR}/${IMG_NAME}"
DOCKER_URL="docker://quay.io/biocontainers/checkm2:1.0.2--pyh7cba7a3_0"

echo "========================================================"
echo "START -- $(date)"
echo "HOSTNAME -- $(hostname)"
echo "========================================================"

# ---------------------------------------------------------
# 1. TÉLÉCHARGEMENT DE L'IMAGE
# ---------------------------------------------------------
if [ -f "${IMG_PATH}" ]; then
    echo "[OK] Image déjà présente : ${IMG_PATH}"
else
    echo "-> Téléchargement de l'image CheckM2 : ${DOCKER_URL} ..."
    if command -v apptainer &> /dev/null; then
        apptainer pull "${IMG_PATH}" "${DOCKER_URL}"
    else
        singularity pull "${IMG_PATH}" "${DOCKER_URL}"
    fi
fi

# ---------------------------------------------------------
# 2. TÉLÉCHARGEMENT DE LA BASE DE DONNÉES
# ---------------------------------------------------------
echo "--------------------------------------------------------"
if [ -f "${DB_DIR}/uniref100.KO.1.dmnd" ]; then
    echo "[OK] Base de données CheckM2 déjà présente dans : ${DB_DIR}"
else
    echo "-> Téléchargement de la base de données CheckM2 (~3 Go)..."
    
    # On utilise l'outil intégré dans le conteneur pour télécharger la DB officielle
    APPTAINER_BIN=$(command -v apptainer || command -v singularity)
    "${APPTAINER_BIN}" exec "${IMG_PATH}" checkm2 database --download --path "${DB_DIR}"
    
    echo "[OK] Base de données téléchargée avec succès !"
fi

# ---------------------------------------------------------
# 3. VÉRIFICATION
# ---------------------------------------------------------
echo "--------------------------------------------------------"
echo "TEST -- Vérification de l'installation :"
APPTAINER_BIN=$(command -v apptainer || command -v singularity)
"${APPTAINER_BIN}" exec "${IMG_PATH}" checkm2 --version || true
echo "--------------------------------------------------------"

echo "Nettoyage du cache Apptainer/Singularity..."
rm -rf "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

echo "END -- $(date)"