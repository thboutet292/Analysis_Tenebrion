#!/bin/bash
#SBATCH --job-name=pull_coverm_sif
#SBATCH --output=log/PULL_COVERM_%j.out
#SBATCH --error=log/PULL_COVERM_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal,lmge
#SBATCH --time=01:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   pull_coverm.slurm
#
#   GOAL    :   Téléchargement de l'image Singularity/Apptainer de CoverM.
#               Dernière version stable : 0.7.0 (Docker Hub : wwood/coverm).
#
#               CoverM calcule la profondeur de couverture des reads sur les contigs
#               assemblés. Il gère en interne l'indexation et l'alignement (minimap2/BWA).
#               Sortie : tableau contig × profondeur (mean, trimmed_mean, rpkm...).
#######################################################################################

set -euo pipefail

module purge

cd "${SLURM_SUBMIT_DIR}"
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
CONTAINERS_DIR="${PROJECT_DIR}/containers"
mkdir -p "${CONTAINERS_DIR}" log

# Isolation du cache Apptainer sur le scratch
SCRATCH_BASE="/storage/scratch/${USER}"
APPTAINER_CACHEDIR="${SCRATCH_BASE}/apptainer_cache"
APPTAINER_TMPDIR="${SCRATCH_BASE}/apptainer_tmp"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"
export APPTAINER_CACHEDIR
export APPTAINER_TMPDIR
export TMPDIR="${APPTAINER_TMPDIR}"

IMG_NAME="coverm_v0.7.0.sif"
IMG_PATH="${CONTAINERS_DIR}/${IMG_NAME}"
DOCKER_URL="docker://quay.io/biocontainers/coverm:0.7.0--h07ea13f_0"

echo "START -- $(date)"
echo "HOSTNAME -- $(hostname)"
echo "IMG_PATH -- ${IMG_PATH}"

if [ -f "${IMG_PATH}" ]; then
    echo "Image déjà présente : ${IMG_PATH}"
else
    echo "Téléchargement de l'image CoverM : ${DOCKER_URL} ..."
    if command -v apptainer &> /dev/null; then
        apptainer pull "${IMG_PATH}" "${DOCKER_URL}"
    else
        singularity pull "${IMG_PATH}" "${DOCKER_URL}"
    fi
fi

echo "--------------------------------------------------------"
echo "TEST -- Vérification de l'installation de CoverM :"

APPTAINER_BIN=$(command -v apptainer || command -v singularity)
"${APPTAINER_BIN}" exec "${IMG_PATH}" coverm --version || true

echo "--------------------------------------------------------"

echo "Nettoyage du cache Apptainer/Singularity..."
rm -rf "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

echo "END -- $(date)"
