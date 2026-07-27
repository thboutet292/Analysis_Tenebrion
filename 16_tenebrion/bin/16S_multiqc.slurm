#!/bin/bash
#SBATCH --job-name=16S_multiqc
#SBATCH --output=log/MULTIQC_%j.out
#SBATCH --error=log/MULTIQC_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --partition=normal,lmge
#SBATCH --time=00:30:00
#SBATCH --chdir=/home/thboutet/projets/16_tenebrion

#######################################################################################

#   AUTHORS :   Thomas BOUTET (IHP,LMGE)
#   DATE    :   03 - 02 - 2026
#   SCRIPT  :   16S_multiqc.slurm
#   TOOLS   :   MultiQC (v1.7) + gcc/8.1.0

#######################################################################################

# sbatch bin/16S_multiqc.slurm

######################################################################################

set -euo pipefail

# Load modules
module purge
module load gcc/8.1.0
module load python/3.7.1
module load MultiQC/1.7

# Working directory pathcd "${SLURM_SUBMIT_DIR}"
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

# Definition of input and output paths 
FASTQC_DIR="${PROJECT_DIR}/results/qc/fastqc"
MULTIQC_OUT="${PROJECT_DIR}/results/qc/multiqc"
LOG_DIR="${PROJECT_DIR}/log"

# Silently create output directories if they don't already exist
mkdir -p "${MULTIQC_OUT}" "${LOG_DIR}"

# Preliminary check that the folder containing the raw results exists
if [ ! -d "${FASTQC_DIR}" ]; then
    echo "ERROR: FastQC folder not found: ${FASTQC_DIR}" >&2
    exit 1
fi

# Count FastQC HTML reports to make sure the previous step ran correctly
NHTML=$(ls -1 "${FASTQC_DIR}"/*_fastqc.html 2>/dev/null | wc -l)
if [ "${NHTML}" -eq 0 ]; then
    echo "ERROR: no *_fastqc.html file found in ${FASTQC_DIR}" >&2
    exit 1
fi


echo "START job -- $(date)"
echo "HOSTNAME -- $(hostname)"
echo "FASTQC_DIR -- ${FASTQC_DIR}"
echo "MULTIQC_OUT -- ${MULTIQC_OUT}"
echo "NB_FASTQC_HTML -- ${NHTML}"
echo "SLURM_JOB_ID -- ${SLURM_JOB_ID}"

# Run MultiQC on the target folder
# The --force argument overwrites the previous report if it exists
multiqc "${FASTQC_DIR}" \
    --outdir "${MULTIQC_OUT}" \
    --force

echo "END job -- $(date)"
