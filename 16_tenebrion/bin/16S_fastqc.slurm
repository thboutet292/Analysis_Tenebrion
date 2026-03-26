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

# sbatch /bin/16S_fastqc.slurm

#######################################################################################

set -euo pipefail

# Chargement des modules
module purge
module load java
module load fastqc/0.11.7

# Chemin repertoire de travail
cd "${SLURM_SUBMIT_DIR}"
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

# Définition des chemins d'arborescence
RAW_DIR="${PROJECT_DIR}/data/raw"
FASTQC_OUT="${PROJECT_DIR}/results/qc/fastqc"
LOG_DIR="${PROJECT_DIR}/log"

# Utilisation du scratch du cluster pour écriture des fichiers
WORKDIR="/storage/scratch/${USER}"

# Création des dossiers si jamais besoin de les créer 
mkdir -p "${FASTQC_OUT}" "${LOG_DIR}" "${WORKDIR}"

if [ ! -d "${RAW_DIR}" ]; then
    echo "ERROR: dossier introuvable: ${RAW_DIR}" >&2
    exit 1
fi

# Compte le nombre d'échantillons uniques en listant les fichiers R1
# Le sed extrait le préfixe exact de l'échantillon pour le comptage
N=$(ls -1 "${RAW_DIR}"/*_R1.fastq.gz 2>/dev/null \
    | sed 's|.*/||' \
    | sed 's/_R1\.fastq\.gz$//' \
    | sort -u \
    | wc -l)

if [ "${N}" -eq 0 ]; then
    echo "ERROR: aucun fichier *_R1.fastq.gz trouvé dans ${RAW_DIR}" >&2
    exit 1
fi

# Auto-soumission 
if [ -z "${SLURM_ARRAY_TASK_COUNT:-}" ]; then
    echo "Auto-submission en array: 0-$((${N}-1))"
    sbatch --array=0-$((${N}-1)) "$0" "$@"
    exit 0
fi

# Stocke les préfixes des échantillons dans un tableau bash indexé
mapfile -t SAMPLES < <(ls -1 "${RAW_DIR}"/*_R1.fastq.gz \
    | sed 's|.*/||' \
    | sed 's/_R1\.fastq\.gz$//' \
    | sort -u)

# Assigne l'échantillon à traiter pour ce sous-job spécifique via son ID d'array
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

# Vérification de l'intégrité de la paire read1/read2 avant de lancer le calcul
if [ ! -f "${R1}" ] || [ ! -f "${R2}" ]; then
    echo "ERROR: paire R1/R2 incomplète pour l’échantillon ${SAMPLE}" >&2
    exit 1
fi

# Création d'un dossier de travail temporaire et unique pour ce sous-job dans le scratch.
# Prévient les conflits d'écriture entre les différentes tâches de l'array.
TASKDIR="${WORKDIR}/16S_fastqc_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TASKDIR}"

# Exécution de FastQC sur la paire.
fastqc -t "${SLURM_CPUS_PER_TASK}" -o "${TASKDIR}" "${R1}" "${R2}"

# Rapatriement des rapports générés vers le stockage persistant du projet
mv "${TASKDIR}"/* "${FASTQC_OUT}/"

# Nettoyage de l'espace temporaire sur le noeud de calcul
rm -rf "${TASKDIR}"

echo "END job -- $(date)"