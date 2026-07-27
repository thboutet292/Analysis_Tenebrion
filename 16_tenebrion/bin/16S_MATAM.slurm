#!/bin/bash
#SBATCH --job-name=16S_MATAM
#SBATCH --output=log/PROD_HYBRID_%A_%a.out
#SBATCH --error=log/PROD_HYBRID_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --partition=long,lmge
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/home/thboutet/projets/16_tenebrion

#######################################################################################

#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   DATE    :   03 - 03 - 2026
#   SCRIPT  :   16S_MATAM.slurm
#   TOOLS   :   USEARCH 9.2.64 ; MATAM 1.6.1 ; Salmon 1.10.2 (via Conda)
#
#   GOAL    :   Massive production run on 44 samples (full Hybrid Pipeline).
#               - Python cleaning (Quality filtering & trimming, streamed reading)
#               - USEARCH dereplication (Algorithmic acceleration)
#               - MATAM assembly (de novo 16S reconstruction, LCA bypass)
#               - Salmon quantification (Probabilistic EM algorithm, offline mode)

#######################################################################################

# sbatch bin/16S_MATAM.slurm

######################################################################################
set -uo pipefail

# =====================================================================
# INITIALIZATION AND TOOL LOADING
# =====================================================================
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] START UP OF THE PIPELINE"
echo "Job Master ID  : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}"
echo "Task Array ID  : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Nœud de calcul : ${SLURMD_NODENAME}"
echo "====================================================================="

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] Purging and loading modules..."
module purge
module load usearch/9.2.64 2>/dev/null || true

# Activate the Conda environment (library stability)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] Enabling the Conda Environment (Salmon) ..."
module load conda/23.3.1
eval "$(conda shell.bash hook)"
conda activate /home/thboutet/projets/16_tenebrion/resources/env_salmon

# =====================================================================
# PATH AND VARIABLE CONFIGURATION
# =====================================================================

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
RAW_DIR="${PROJECT_DIR}/data/raw"
OUT_BASE="${PROJECT_DIR}/results/PRODUCTION_HYBRID"
RES_DIR="${PROJECT_DIR}/resources"
CONTAINERS_DIR="${PROJECT_DIR}/containers"

MATAM_SIF="${CONTAINERS_DIR}/matam_1.6.1--hdcf5f25_1.sif"
MATAM_DB_PREFIX="${RES_DIR}/matam_db_138/SILVA_138_1"
SINGULARITY="/usr/bin/singularity"
THREADS="${SLURM_CPUS_PER_TASK}"

# =====================================================================
# SLURM ARRAY MANAGEMENT 
# =====================================================================
mapfile -t SAMPLES < <(ls "${RAW_DIR}"/*_R1.fastq.gz | xargs -n 1 basename | sed 's/_R1.fastq.gz//' | sort)

# Automatic Array launch if the script is run without --array
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] Detection of ${N_SAMPLES} samples."
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] Submission of the Array job (10 tasks max in parallel)..."
    sbatch --array=0-$((N_SAMPLES - 1))%11 "$0"
    echo "---------------------------------------------------------------------"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
OUT_DIR="${OUT_BASE}/${SAMPLE}"
mkdir -p "${OUT_DIR}" "log"

# =====================================================================
# LOCAL SCRATCH MANAGEMENT (WORKING SPACE)
# =====================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] Configure Local Scratch Space..."
JOB_SCRATCH="/storage/scratch/${USER}/prod_hyb_${SAMPLE}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
WORKDIR="${JOB_SCRATCH}/work"
TMP_WORK="${JOB_SCRATCH}/tmp"

mkdir -p "${WORKDIR}" "${TMP_WORK}"
export TMPDIR="${TMP_WORK}"
export SINGULARITY_TMPDIR="${TMP_WORK}"

cleanup() {
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Cleaning the scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [END] Task completed for ${SAMPLE}."
    echo "---------------------------------------------------------------------"
}
trap cleanup EXIT

# Wrapper function for Singularity (MATAM / Python)
matam_exec() {
    "${SINGULARITY}" exec --cleanenv \
        --env "LC_ALL=C" --env "LANG=C" --env "TMPDIR=/work/tmp" \
        -B "${WORKDIR}:/work" -B "${TMP_WORK}:/work/tmp" -B "${RES_DIR}:${RES_DIR}" \
        "${MATAM_SIF}" bash -lc "$*"
}

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉBUT DU TRAITEMENT : ${SAMPLE}"
echo "====================================================================="

# =====================================================================
# STEP 1: EXACT CLEANING 
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 1] Python Cleanup (Q>30, L>50, Clip=15)..."

CLIPPED_FASTQ_CONT="/work/${SAMPLE}_clipped.fastq" 
PY_CLIP="${WORKDIR}/clip_fastq_strict.py"

cat > "${PY_CLIP}" <<'PY'
import sys, gzip
inp = sys.argv[1]
outp = sys.argv[2]
CLIP_LEN = 15
CAP = 74
MIN_Q = 34
with gzip.open(inp, "rt", encoding="utf-8", errors="replace") as fi, open(outp, "w", encoding="utf-8") as fo:
    while True:
        h = fi.readline(); s = fi.readline(); p = fi.readline(); q = fi.readline()
        if not (s and p and q): break
        s_cut = s.rstrip("\r\n")[CLIP_LEN:]
        q_cut = q.rstrip("\r\n")[CLIP_LEN:]
        if len(s_cut) >= 50:
            q2 = "".join(chr(max(MIN_Q, min(ord(ch), CAP))) for ch in q_cut)
            fo.write(f"{h.rstrip()}\n{s_cut}\n{p.rstrip()}\n{q2}\n")
PY

matam_exec "python /work/clip_fastq_strict.py '${RAW_DIR}/${SAMPLE}_R1.fastq.gz' '${CLIPPED_FASTQ_CONT}'"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 1] Cleaning completed."

# =====================================================================
# STEP 2: DEREPLICATION
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 2] USEARCH Dereplication..."

DEREP_FASTQ_WORK="${WORKDIR}/${SAMPLE}_derep.fastq"
DEREP_FASTQ_CONT="/work/${SAMPLE}_derep.fastq"

usearch -fastx_uniques "${WORKDIR}/${SAMPLE}_clipped.fastq" -fastqout "${DEREP_FASTQ_WORK}" -sizeout -quiet
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 2] Dereplication completed."

# =====================================================================
# STEP 3: MATAM ASSEMBLY
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 3] MATAM Assembly..."

MATAM_CONT_OUT="/work/02_matam"
MATAM_SCRATCH_OUT="${WORKDIR}/02_matam"
mkdir -p "${MATAM_SCRATCH_OUT}"

matam_exec "\
matam_assembly.py \
    --verbose \
    --cpu '${THREADS}' \
    --max_memory 8000 \
    --assembler SGA \
    --ref_db '${MATAM_DB_PREFIX}' \
    --out_dir '${MATAM_CONT_OUT}' \
    --input_fastq '${DEREP_FASTQ_CONT}' \
"
SCAFFOLDS_FA="${MATAM_SCRATCH_OUT}/final_assembly.fa"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 3] Assembly completed."

# =====================================================================
# STEP 4: SALMON QUANTIFICATION
# =====================================================================

echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 4] Salmon Quantification Salmon..."

SALMON_INDEX="${WORKDIR}/salmon_index"
SALMON_OUT="${WORKDIR}/salmon_quant"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 4.1] Building the index of scaffolds..."
salmon index -t "${SCAFFOLDS_FA}" -i "${SALMON_INDEX}" -k 31 --no-version-check

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 4.2] Quantification (Mapping-back)..."
salmon quant -i "${SALMON_INDEX}" -l U -r "${WORKDIR}/${SAMPLE}_clipped.fastq" -p "${THREADS}" --validateMappings --no-version-check -o "${SALMON_OUT}" --quiet
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [STEP 4] Quantification completed."

# =====================================================================
# STEP 5: DATA RETRIEVAL
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 5] Data retrieval..."

cp -L "${SCAFFOLDS_FA}" "${OUT_DIR}/${SAMPLE}_scaffolds.fasta"
cp -L "${SALMON_OUT}/quant.sf" "${OUT_DIR}/${SAMPLE}_abundance_salmon.tsv"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 5] Data successfully copied in: ${OUT_DIR}"
