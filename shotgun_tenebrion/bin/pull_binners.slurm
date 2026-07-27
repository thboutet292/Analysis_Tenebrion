#!/bin/bash
#SBATCH --job-name=pull_binners
#SBATCH --output=log/PULL_BINNERS_%j.out
#SBATCH --error=log/PULL_BINNERS_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal,lmge
#SBATCH --time=02:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   pull_binners.slurm
#
#   GOAL    :   Téléchargement des images Apptainer pour les trois outils de binning
#               utilisés en production : MetaBAT2, CONCOCT et SemiBin2.
#######################################################################################

set -euo pipefail

module purge
module load apptainer 2>/dev/null || true

cd "${SLURM_SUBMIT_DIR}"
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
CONTAINERS_DIR="${PROJECT_DIR}/containers"
mkdir -p "${CONTAINERS_DIR}" log

# Isolation du cache
SCRATCH_BASE="/storage/scratch/${USER}"
export APPTAINER_CACHEDIR="${SCRATCH_BASE}/apptainer_cache_binners"
export APPTAINER_TMPDIR="${SCRATCH_BASE}/apptainer_tmp_binners"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"
export TMPDIR="${APPTAINER_TMPDIR}"

echo "========================================================"
echo "START -- $(date)"
echo "========================================================"

# --- IMAGE 1 : METABAT2 ---
IMG_METABAT="${CONTAINERS_DIR}/metabat2_v2.15.sif"
URL_METABAT="docker://quay.io/biocontainers/metabat2:2.15--h986a166_1"

if [ ! -f "${IMG_METABAT}" ]; then
    echo "-> Téléchargement de MetaBAT2..."
    apptainer pull "${IMG_METABAT}" "${URL_METABAT}"
fi

# --- IMAGE 2 : CONCOCT ---
IMG_CONCOCT="${CONTAINERS_DIR}/concoct_1.1.0.sif"
URL_CONCOCT="docker://quay.io/biocontainers/concoct:1.1.0--py38h7be5676_2"

if [ ! -f "${IMG_CONCOCT}" ]; then
    echo "-> Téléchargement de CONCOCT..."
    apptainer pull "${IMG_CONCOCT}" "${URL_CONCOCT}"
fi

# --- IMAGE 3 : SEMIBIN2 ---
IMG_SEMIBIN="${CONTAINERS_DIR}/semibin_2.1.0.sif"
URL_SEMIBIN="docker://quay.io/biocontainers/semibin:2.1.0--pyhdfd78af_0"

if [ ! -f "${IMG_SEMIBIN}" ]; then
    echo "-> Téléchargement de SemiBin2..."
    apptainer pull "${IMG_SEMIBIN}" "${URL_SEMIBIN}"
fi

echo "========================================================"
echo "TEST D'INSTALLATION"
echo "========================================================"
# metabat2 n'a pas de flag --version standard, on utilise -h (help)
apptainer exec "${IMG_METABAT}" metabat2 -h || true
apptainer exec "${IMG_CONCOCT}" concoct --version || true
apptainer exec "${IMG_SEMIBIN}" SemiBin2 --version || true

rm -rf "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"
echo "END -- $(date)"