#!/bin/bash
#SBATCH --job-name=install_picrust2
#SBATCH --output=log/INSTALL_PICRUST2_%j.out
#SBATCH --error=log/INSTALL_PICRUST2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=normal
#SBATCH --time=01:00:00
#SBATCH --chdir=/home/thboutet/projets/16_tenebrion

set -euo pipefail
module purge

# Chargement des pré-requis de base
module load gcc/8.1.0 2>/dev/null || true

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
LOG_DIR="${PROJECT_DIR}/log"
CONTAINERS_DIR="${PROJECT_DIR}/containers"
mkdir -p "${LOG_DIR}" "${CONTAINERS_DIR}"

SINGULARITY="/usr/bin/singularity"
if [ ! -x "${SINGULARITY}" ]; then
    SINGULARITY="$(command -v singularity || true)"
fi
if [ -z "${SINGULARITY}" ] || [ ! -x "${SINGULARITY}" ]; then
    echo "ERROR: singularity/apptainer introuvable" >&2
    exit 1
fi

################################################################################
#   SCRATCH (TMP + CACHE) + CLEANUP robuste
################################################################################

SCRATCH_BASE="/storage/scratch/${USER}"
SCRATCH_JOB="${SCRATCH_BASE}/install_picrust2_${SLURM_JOB_ID}"
SCRATCH_TMP="${SCRATCH_JOB}/tmp"
mkdir -p "${SCRATCH_JOB}" "${SCRATCH_TMP}"

export SINGULARITY_TMPDIR="${SCRATCH_TMP}"
export SINGULARITY_CACHEDIR="${SCRATCH_BASE}/singularity_cache"
mkdir -p "${SINGULARITY_CACHEDIR}" 2>/dev/null || true

cleanup() {
    # Nettoyage scratch du job uniquement
    rm -rf "${SCRATCH_JOB}" 2>/dev/null || true
}
trap cleanup EXIT

################################################################################
#   LOG HEADER
################################################################################

echo "START -- $(date)"
echo "HOSTNAME -- $(hostname)"
echo "PROJECT_DIR -- ${PROJECT_DIR}"
echo "SINGULARITY -- ${SINGULARITY}"
echo "SINGULARITY_VERSION -- $(${SINGULARITY} --version 2>/dev/null || true)"
echo "SINGULARITY_CACHEDIR -- ${SINGULARITY_CACHEDIR}"
echo "SINGULARITY_TMPDIR -- ${SINGULARITY_TMPDIR}"
echo "SCRATCH_JOB -- ${SCRATCH_JOB}"
df -h "${SCRATCH_BASE}" || true
echo

################################################################################
#   1) BUILD PICRUST2 IMAGE (depuis biocontainers)
################################################################################

PICRUST2_VERSION="2.5.2"
# Tag officiel de Biocontainers pour la version stable 2.5.2
PICRUST2_DOCKER_IMAGE="docker://quay.io/biocontainers/picrust2:2.5.2--pyhdfd78af_0"

PICRUST2_SIF="${CONTAINERS_DIR}/picrust2.sif"
PICRUST2_TMP="${PICRUST2_SIF}.tmp"

echo "################################################################################"
echo "# 1) BUILD PICRUST2 IMAGE"
echo "################################################################################"
echo "PICRUST2_VERSION -- ${PICRUST2_VERSION}"
echo "PICRUST2_DOCKER_IMAGE -- ${PICRUST2_DOCKER_IMAGE}"
echo "PICRUST2_SIF -- ${PICRUST2_SIF}"

if [ -s "${PICRUST2_SIF}" ]; then
    echo "PICRUST2_SIF -- déjà présent -> skip build"
    ls -lh "${PICRUST2_SIF}" | sed 's/^/   /'
else
    rm -f "${PICRUST2_TMP}"
    echo "BUILD -- ${SINGULARITY} build ${PICRUST2_TMP} ${PICRUST2_DOCKER_IMAGE}"
    ${SINGULARITY} build "${PICRUST2_TMP}" "${PICRUST2_DOCKER_IMAGE}"
    mv -f "${PICRUST2_TMP}" "${PICRUST2_SIF}"
    echo "BUILD OK -- ${PICRUST2_SIF}"
    ls -lh "${PICRUST2_SIF}" | sed 's/^/   /'
fi
echo

################################################################################
#   2) CHECK PICRUST2 BINARIES
################################################################################

echo "CHECK PICRUST2 BINARIES"
${SINGULARITY} exec --cleanenv "${PICRUST2_SIF}" bash -lc '
    set -e
    echo "PATH=$PATH"
    echo "--- Version Python ---"
    python --version 2>/dev/null || true
    echo "--- Validation des binaires PICRUSt2 ---"
    for x in picrust2_pipeline.py place_seqs.py hsp.py pathway_pipeline.py; do
        echo "--- which $x"
        command -v "$x" || true
    done
    echo "--- Test execution picrust2_pipeline.py ---"
    picrust2_pipeline.py -v || true
'
echo

echo "INSTALL DONE -- $(date)"
