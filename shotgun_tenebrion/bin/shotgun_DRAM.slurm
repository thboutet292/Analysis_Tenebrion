#!/bin/bash
#SBATCH --job-name=DRAM
#SBATCH --output=log/DRAM_%j.out
#SBATCH --error=log/DRAM_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --partition=lmge,long
#SBATCH --time=76:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   GOAL    :   DRAM annotate + distill sur un bin
#   USAGE   :   sbatch bin/shotgun_DRAM.slurm 
#######################################################################################

set -euo pipefail

BIN_NAME="${1}"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DRAM — ${BIN_NAME} — Job ${SLURM_JOB_ID}"
echo "====================================================================="

module purge
module load apptainer 2>/dev/null || true
module load rclone/1.55.1 2>/dev/null || true

APPTAINER_BIN=$(command -v apptainer || command -v singularity)
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
IMG_DRAM="${PROJECT_DIR}/containers/dram_v1.4.6.sif"
S3_DB="s3_uca:lmge-tenebrion/ref/dram_db"
TARGET_CONFIG="/usr/local/lib/python3.10/site-packages/mag_annotator/CONFIG"
TARGET_ANNOTATE="/usr/local/lib/python3.10/site-packages/mag_annotator/annotate_bins.py"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

BIN_FA="${PROJECT_DIR}/results/filtered_bins/${BIN_NAME}.fa"
DIR_OUT_RAW="${PROJECT_DIR}/results/dram_annotation_raw/${BIN_NAME}"
DIR_OUT_DISTILL="${PROJECT_DIR}/results/dram_annotation_distill/${BIN_NAME}"

mkdir -p log/ "${DIR_OUT_RAW}" "${DIR_OUT_DISTILL}"

JOB_SCRATCH="/storage/scratch/${USER}/dram_${SLURM_JOB_ID}"
DB_LOCAL="${JOB_SCRATCH}/dram_db"
PATCH_DIR="${JOB_SCRATCH}/patch"
FAKE_TMP="${JOB_SCRATCH}/fake_tmp"
CONFIG_FILE="${JOB_SCRATCH}/dram_config.json"

mkdir -p "${DB_LOCAL}" "${PATCH_DIR}" "${FAKE_TMP}" "${JOB_SCRATCH}/apptainer_tmp"
export APPTAINER_TMPDIR="${JOB_SCRATCH}/apptainer_tmp"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 1 : TÉLÉCHARGEMENT DB
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Téléchargement base DRAM..."

rclone copy "${S3_DB}/" "${DB_LOCAL}/" \
    --exclude "database_files/**" \
    --exclude "kofam_profiles/profiles/**" \
    --exclude "tmp/**" \
    --transfers 16 -P

echo " -> DB disponible : $(du -sh ${DB_LOCAL} | cut -f1)"

# =====================================================================
# ÉTAPE 2 : CONFIG DRAM
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Config DRAM..."

echo "{}" > "${CONFIG_FILE}"
chmod 777 "${CONFIG_FILE}"

PEPTIDASE_MMSEQS=$(find "${DB_LOCAL}" -name "peptidases.*" ! -name "*.h3?" 2>/dev/null | grep -v database_files | head -1 || echo "")
VIRAL_MMSEQS=$(find "${DB_LOCAL}" -name "refseq_viral*.mmsdb" ! -name "*.dbtype" ! -name "*.idx" ! -name "*.index" ! -name "*_h*" ! -name "*.lookup" ! -name "*.source" 2>/dev/null | grep -v database_files | head -1 || echo "")
VOG_HMM=$(find "${DB_LOCAL}" -name "vog_latest_hmms.txt" 2>/dev/null | head -1 || echo "")
VOG_ANNOT=$(find "${DB_LOCAL}" -name "vog_annotations_latest.tsv.gz" 2>/dev/null | head -1 || echo "")
PFAM_MMSPRO=$(find "${DB_LOCAL}" -name "pfam.mmspro" 2>/dev/null | head -1 || echo "")
DESCRIPTION_DB=$(find "${DB_LOCAL}" -name "description_db.sqlite" 2>/dev/null | head -1 || echo "")

echo " -> Peptidase : ${PEPTIDASE_MMSEQS:-NON TROUVÉ}"
echo " -> Viral     : ${VIRAL_MMSEQS:-NON TROUVÉ}"
echo " -> VOG HMM   : ${VOG_HMM:-NON TROUVÉ}"
echo " -> Pfam      : ${PFAM_MMSPRO:-NON TROUVÉ}"

PEPTIDASE_ARG=""; [ -n "${PEPTIDASE_MMSEQS}" ] && PEPTIDASE_ARG="--peptidase_loc ${PEPTIDASE_MMSEQS}"
VIRAL_ARG="";     [ -n "${VIRAL_MMSEQS}" ]     && VIRAL_ARG="--viral_loc ${VIRAL_MMSEQS}"
VOG_ARG="";       [ -n "${VOG_HMM}" ]           && VOG_ARG="--vogdb_loc ${VOG_HMM}"
VOG_ANNOT_ARG=""; [ -n "${VOG_ANNOT}" ]         && VOG_ANNOT_ARG="--vog_annotations ${VOG_ANNOT}"
PFAM_ARG="";      [ -n "${PFAM_MMSPRO}" ]       && PFAM_ARG="--pfam_loc ${PFAM_MMSPRO}"
DESC_ARG="";      [ -n "${DESCRIPTION_DB}" ]    && DESC_ARG="--description_db_loc ${DESCRIPTION_DB}"

"${APPTAINER_BIN}" exec \
    --writable-tmpfs \
    --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
    --bind "${FAKE_TMP}:/tmp" \
    --bind "${CONFIG_FILE}:${TARGET_CONFIG}" \
    "${IMG_DRAM}" \
    DRAM-setup.py set_database_locations \
        --kofam_hmm_loc             "${DB_LOCAL}/kofam_profiles.hmm" \
        --kofam_ko_list_loc         "${DB_LOCAL}/kofam_ko_list.tsv" \
        --dbcan_loc                 "${DB_LOCAL}/dbCAN-HMMdb-V11.txt" \
        --dbcan_fam_activities_loc  "${DB_LOCAL}/CAZyDB.08062022.fam-activities.txt" \
        --dbcan_subfam_ec_loc       "${DB_LOCAL}/CAZyDB.08062022.fam.subfam.ec.txt" \
        --genome_summary_form_loc   "${DB_LOCAL}/genome_summary_form.20260521.tsv" \
        --module_step_form_loc      "${DB_LOCAL}/module_step_form.20260521.tsv" \
        --etc_module_database_loc   "${DB_LOCAL}/etc_mdoule_database.20260521.tsv" \
        --function_heatmap_form_loc "${DB_LOCAL}/function_heatmap_form.20260521.tsv" \
        --amg_database_loc          "${DB_LOCAL}/amg_database.20260521.tsv" \
        ${PEPTIDASE_ARG} \
        ${VIRAL_ARG} \
        ${VOG_ARG} \
        ${VOG_ANNOT_ARG} \
        ${PFAM_ARG} \
        ${DESC_ARG}

# =====================================================================
# ÉTAPE 3 : PATCH annotate_bins.py
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] Patch annotate_bins.py..."

"${APPTAINER_BIN}" exec "${IMG_DRAM}" \
    cat "${TARGET_ANNOTATE}" > "${PATCH_DIR}/annotate_bins.py"

sed -i \
    's/annotation_str += "; %s (db=%s)" % (annotation\.pfam_hits, "pfam")/if hasattr(annotation, "pfam_hits"): annotation_str += "; %s (db=%s)" % (annotation.pfam_hits, "pfam")/' \
    "${PATCH_DIR}/annotate_bins.py"

echo " -> Patch OK."

# =====================================================================
# ÉTAPE 4 : DRAM.py annotate
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 4] Annotate — ${BIN_NAME}..."

rm -rf "${DIR_OUT_RAW}"

"${APPTAINER_BIN}" exec \
    --writable-tmpfs \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
    --bind "${FAKE_TMP}:/tmp" \
    --bind "${CONFIG_FILE}:${TARGET_CONFIG}" \
    --bind "${PATCH_DIR}/annotate_bins.py:${TARGET_ANNOTATE}" \
    --env TMPDIR="/tmp" --env TEMP="/tmp" --env TMP="/tmp" \
    "${IMG_DRAM}" \
    DRAM.py annotate \
        -i "${BIN_FA}" \
        -o "${DIR_OUT_RAW}" \
        --threads "${THREADS}" \
        --min_contig_size 1000 \
        --verbose

# =====================================================================
# ÉTAPE 5 : DRAM.py distill
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 5] Distill — ${BIN_NAME}..."

rm -rf "${DIR_OUT_DISTILL}"

RRNA_ARG=""; [ -f "${DIR_OUT_RAW}/rrnas.tsv" ] && RRNA_ARG="--rrna_path ${DIR_OUT_RAW}/rrnas.tsv"
TRNA_ARG=""; [ -f "${DIR_OUT_RAW}/trnas.tsv" ] && TRNA_ARG="--trna_path ${DIR_OUT_RAW}/trnas.tsv"

"${APPTAINER_BIN}" exec \
    --writable-tmpfs \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
    --bind "${FAKE_TMP}:/tmp" \
    --bind "${CONFIG_FILE}:${TARGET_CONFIG}" \
    --env TMPDIR="/tmp" --env TEMP="/tmp" --env TMP="/tmp" \
    "${IMG_DRAM}" \
    DRAM.py distill \
        -i "${DIR_OUT_RAW}/annotations.tsv" \
        -o "${DIR_OUT_DISTILL}" \
        ${RRNA_ARG} \
        ${TRNA_ARG}

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] TERMINÉ — ${BIN_NAME}"
echo "====================================================================="