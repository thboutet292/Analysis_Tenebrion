#!/bin/bash
#SBATCH --job-name=pull_metabat2_sif
#SBATCH --output=log/PULL_METABAT2_%j.out
#SBATCH --error=log/PULL_METABAT2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal,lmge
#SBATCH --time=01:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   pull_metabat2.slurm
#
#   GOAL    :   Téléchargement de l'image Singularity/Apptainer de MetaBAT2.
#               Dernière version stable : 2.15 (Biocontainers).
#
#               MetaBAT2 reconstruit des génomes (MAGs) à partir de contigs métagénomiques.
#               Il utilise les fréquences de tétramères (composition) et la 
#               co-abondance différentielle issue des fichiers BAM (via jgi_summarize).
#######################################################################################

set -euo pipefail

module purge
module load apptainer 2>/dev/null || true

cd "${SLURM_SUBMIT_DIR}"
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
CONTAINERS_DIR="${PROJECT_DIR}/containers"
mkdir -p "${CONTAINERS_DIR}" log

# Isolation du cache Apptainer sur le scratch
SCRATCH_BASE="/storage/scratch/${USER}"
APPTAINER_CACHEDIR="${SCRATCH_BASE}/apptainer_cache_metabat2"
APPTAINER_TMPDIR="${SCRATCH_BASE}/apptainer_tmp_metabat2"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"
export APPTAINER_CACHEDIR
export APPTAINER_TMPDIR
export TMPDIR="${APPTAINER_TMPDIR}"

IMG_NAME="metabat2_v2.15.sif"
IMG_PATH="${CONTAINERS_DIR}/${IMG_NAME}"
# L'image officielle la plus stable et légère sur Biocontainers
DOCKER_URL="docker://quay.io/biocontainers/metabat2:2.15--h986a166_1"

echo "START -- $(date)"
echo "HOSTNAME -- $(hostname)"
echo "IMG_PATH -- ${IMG_PATH}"

if [ -f "${IMG_PATH}" ]; then
    echo "Image déjà présente : ${IMG_PATH}"
else
    echo "Téléchargement de l'image MetaBAT2 : ${DOCKER_URL} ..."
    if command -v apptainer &> /dev/null; then
        apptainer pull "${IMG_PATH}" "${DOCKER_URL}"
    else
        singularity pull "${IMG_PATH}" "${DOCKER_URL}"
    fi
fi

echo "--------------------------------------------------------"
echo "TEST -- Vérification de l'installation de MetaBAT2 :"

APPTAINER_BIN=$(command -v apptainer || command -v singularity)
# On utilise -h (help) car metabat2 n'a pas de flag --version standard
"${APPTAINER_BIN}" exec "${IMG_PATH}" metabat2 -h || true

echo "--------------------------------------------------------"

echo "Nettoyage du cache Apptainer/Singularity..."
rm -rf "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

echo "END -- $(date)"