#!/bin/bash
#SBATCH --job-name=REBUILD_VP
#SBATCH --output=log/REBUILD_VP_%j.out
#SBATCH --error=log/REBUILD_VP_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=normal,lmge
#SBATCH --time=01:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE : REBUILD VIRAL + PEPTIDASE"
echo "====================================================================="

module purge
module load apptainer 2>/dev/null || true
module load rclone/1.55.1 2>/dev/null || true

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
S3_DEST="s3_uca:lmge-tenebrion/ref/dram_db"
IMG_DRAM="${PROJECT_DIR}/containers/dram_v1.4.6.sif"
APPTAINER_BIN=$(command -v apptainer || command -v singularity)
TARGET_CONFIG="/usr/local/lib/python3.10/site-packages/mag_annotator/CONFIG"

mkdir -p log/

JOB_SCRATCH="/storage/scratch/${USER}/dram_vp_${SLURM_JOB_ID}"
DB_DIR="${JOB_SCRATCH}/dram_db"
RAW_DIR="${JOB_SCRATCH}/raw_files"
FAKE_CONFIG="${JOB_SCRATCH}/dram_config.json"

mkdir -p "${DB_DIR}" "${RAW_DIR}" "${JOB_SCRATCH}/tmp"
export APPTAINER_TMPDIR="${JOB_SCRATCH}/tmp"
export TMPDIR="${JOB_SCRATCH}/tmp"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Upload S3..."
    rclone copy "${DB_DIR}/" "${S3_DEST}/" --transfers 8 -P || true
    echo " -> Upload terminé."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Rapatriement depuis S3..."
rclone copy "${S3_DEST}/database_files/viral.merged.protein.faa.gz" "${RAW_DIR}/" -P
rclone copy "${S3_DEST}/database_files/merops_peptidases_nr.faa" "${RAW_DIR}/" -P

echo "{}" > "${FAKE_CONFIG}"
chmod 777 "${FAKE_CONFIG}"

# TA fonction pour contourner le bug des arguments multiples !
run_dram_select() {
    local SELECT="$1"
    local EXTRA_ARG="$2"
    local EXTRA_VAL="$3"

    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Compilation : ${SELECT}..."
    
    rm -rf "${DB_DIR}/database_files" 2>/dev/null || true

    "${APPTAINER_BIN}" exec \
        --writable-tmpfs \
        --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
        --bind "${FAKE_CONFIG}:${TARGET_CONFIG}" \
        "${IMG_DRAM}" \
        DRAM-setup.py prepare_databases \
            --output_dir "${DB_DIR}" \
            "${EXTRA_ARG}" "${EXTRA_VAL}" \
            --select_db "${SELECT}" \
            --threads "${SLURM_CPUS_PER_TASK}" \
            --skip_uniref \
            --verbose

    echo " -> ${SELECT} terminé."
}

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Compilation séquentielle..."
run_dram_select "viral" "--viral_loc" "${RAW_DIR}/viral.merged.protein.faa.gz"
run_dram_select "peptidase" "--peptidase_loc" "${RAW_DIR}/merops_peptidases_nr.faa"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] COMPILATION TERMINÉE."
echo "====================================================================="