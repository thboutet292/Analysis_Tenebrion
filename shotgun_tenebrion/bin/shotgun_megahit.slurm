#!/bin/bash
#SBATCH --job-name=MEGAHIT_FULL
#SBATCH --output=log/MEGAHIT_FULL_%A_%a.out
#SBATCH --error=log/MEGAHIT_FULL_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G                    
#SBATCH --partition=long,lmge
#SBATCH --time=48:00:00              

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   shotgun_megahit_full.slurm
#   TOOLS   :   rclone ; MEGAHIT
#
#   GOAL    :   Assemblage de novo complet du microbiome via MEGAHIT.
#
#   STEPS   :   1. Importation dynamique S3 -> Scratch des fichiers entiers.
#               2. Assemblage par graphes de de Bruijn via MEGAHIT (--presets meta-sensitive).
#               3. Nettoyage et rapatriement des métriques et contigs finaux.
#######################################################################################

# Mode strict : arrêt sur erreur ou variable non définie
set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] DÉMARRAGE DE L'ASSEMBLAGE MEGAHIT (FULL)"
echo "Job Master ID : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}} | Task : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "====================================================================="

# Chargement exclusif des modules fonctionnels validés sur l'environnement
module purge
module load rclone/1.55.1
module load gcc/8.1.0
module load megahit/1.2.9

# Configuration S3 & Chemins locaux
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_IN_DIR="shotgun_cleaned/fastq"

PROJECT_DIR="${PWD}"
LOCAL_OUT_DIR="${PROJECT_DIR}/results/assembly_megahit"
mkdir -p "${LOCAL_OUT_DIR}"

# Paramètres d'exécution
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# =====================================================================
# GESTION DE L'ARRAY SLURM DYNAMIQUE
# =====================================================================
mapfile -t SAMPLES < <(rclone ls "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/" | awk '{print $2}' | grep "_clean_R1.fastq.gz" | sed 's/_clean_R1.fastq.gz//' | sort)

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    if [ ${N_SAMPLES} -eq 0 ]; then
        echo "[ERREUR] Aucun échantillon _clean_R1 détecté sur le S3." >&2
        exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] ${N_SAMPLES} échantillons détectés. Lancement."
    SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | awk '{print $1}')
    sbatch --array=0-$((N_SAMPLES - 1)) "${SCRIPT_PATH}"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# =====================================================================
# GESTION DU SCRATCH ISOLÉ
# =====================================================================
JOB_SCRATCH="/storage/scratch/${USER}/megahit_full_${SAMPLE}_${SLURM_JOB_ID:-local}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression du répertoire temporaire..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 1 : IMPORTATION DES DONNÉES COMPLÈTES
# =====================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Rapatriement des FASTQ complets..."

# Plus de manipulation de flux : copie directe des fichiers entiers vers le scratch
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_clean_R1.fastq.gz" ./
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_clean_R2.fastq.gz" ./

# =====================================================================
# ÉTAPE 2 : ASSEMBLAGE MEGAHIT
# =====================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Lancement de l'assemblage (meta-sensitive)..."

# Exécution sur les fichiers entiers
megahit \
    -1 "${SAMPLE}_clean_R1.fastq.gz" \
    -2 "${SAMPLE}_clean_R2.fastq.gz" \
    --num-cpu-threads "${THREADS}" \
    --presets meta-sensitive \
    --out-dir megahit_out \
    --out-prefix "${SAMPLE}_megahit"

# =====================================================================
# ÉTAPE 3 : SAUVEGARDE DES RÉSULTATS
# =====================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] Exportation des contigs vers le répertoire projet..."

cp "megahit_out/${SAMPLE}_megahit.contigs.fa" "${LOCAL_OUT_DIR}/"
cp "megahit_out/${SAMPLE}_megahit.log" "${LOCAL_OUT_DIR}/"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] Succès du pipeline pour l'échantillon ${SAMPLE}."