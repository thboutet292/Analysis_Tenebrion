#!/bin/bash
#SBATCH --job-name=16S_fastqc
#SBATCH --output=log/FASTQC_%A_%a.out
#SBATCH --error=log/FASTQC_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=normal,lmge
#SBATCH --time=02:00:00
#SBATCH --chdir=/home/thboutet/projets/16_tenebrion

#######################################################################################

#   AUTHORS :   Thomas BOUTET (IHP,LMGE)
#   DATE    :   03 - 02 - 2026
#   SCRIPT  :   16S_fastqc.slurm
#   TOOLS   :   fastqc (v0.11.7)

#######################################################################################

# sbatch bin/16S_fastqc.slurm

#######################################################################################

set -euo pipefail

# Load modules
module purge
module load java
module load fastqc/0.11.7

# Working directory path
cd "${SLURM_SUBMIT_DIR}"
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

# Definition of directory tree paths
RAW_DIR="${PROJECT_DIR}/data/raw"
FASTQC_OUT="${PROJECT_DIR}/results/qc/fastqc"
LOG_DIR="${PROJECT_DIR}/log"

# Use cluster scratch space for writing files
WORKDIR="/storage/scratch/${USER}"

# Create directories if they don't already exist
mkdir -p "${FASTQC_OUT}" "${LOG_DIR}" "${WORKDIR}"

if [ ! -d "${RAW_DIR}" ]; then
    echo "ERROR: folder not found: ${RAW_DIR}" >&2
    exit 1
fi

# Count the number of unique samples by listing the R1 files
# sed extracts the exact sample prefix for counting
N=$(ls -1 "${RAW_DIR}"/*_R1.fastq.gz 2>/dev/null \
    | sed 's|.*/||' \
    | sed 's/_R1\.fastq\.gz$//' \
    | sort -u \
    | wc -l)

if [ "${N}" -eq 0 ]; then
    echo "ERROR: no *_R1.fastq.gz file found in ${RAW_DIR}" >&2
    exit 1
fi

# Auto-submission
if [ -z "${SLURM_ARRAY_TASK_COUNT:-}" ]; then
    echo "Auto-submission: 0-$((${N}-1))"
    sbatch --array=0-$((${N}-1)) "$0" "$@"
    exit 0
fi

# Store the sample prefixes in an indexed bash array
mapfile -t SAMPLES < <(ls -1 "${RAW_DIR}"/*_R1.fastq.gz \
    | sed 's|.*/||' \
    | sed 's/_R1\.fastq\.gz$//' \
    | sort -u)

# Assign the sample to process for this specific sub-job via its array ID
SAMPLE="${SAMPLES[${SLURM_ARRAY_TASK_ID}]}"
R1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"

echo "START job -- $(date)"
echo "HOSTNAME -- $(hostname)"
echo "SAMPLE -- ${SAMPLE}"
echo "R1 -- ${R1}"
echo "R2 -- ${R2}"
echo "SLURM_JOB_ID -- ${SLURM_JOB_ID}"
echo "SLURM_ARRAY_TASK_ID -- ${SLURM_ARRAY_TASK_ID}"

# Check the integrity of the read1/read2 pair before running the computation
if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "ERROR: incomplete R1/R2 pair for the sample ${SAMPLE}" >&2
    exit 1
fi

# Create a unique temporary working directory for this sub-job in the scratch space.
# Prevents write conflicts between the different tasks of the array.
TASKDIR="${WORKDIR}/16S_fastqc_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TASKDIR}"

# Run FastQC on the pair.
fastqc -t "${SLURM_CPUS_PER_TASK}" -o "${TASKDIR}" "${R1}" "${R2}"

# Move the generated reports to the project's persistent storage
mv "${TASKDIR}"/* "${FASTQC_OUT}/"

# Clean up the temporary space on the compute node
rm -rf "${TASKDIR}"

echo "END job -- $(date)"
