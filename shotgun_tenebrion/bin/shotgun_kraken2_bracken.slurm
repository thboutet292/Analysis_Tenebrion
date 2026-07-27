#!/bin/bash
#SBATCH --job-name=KRAKEN_PFP
#SBATCH --output=log/TAXONOMY_PFP_%A_%a.out
#SBATCH --error=log/TAXONOMY_PFP_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G                    
#SBATCH --partition=normal,lmge
#SBATCH --time=24:00:00

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   DATE    :   18 - 03 - 2026
#   SCRIPT  :   shotgun_kraken2_bracken.slurm
#   TOOLS   :   rclone ; Apptainer (Kraken2 / Bracken)
#
#   GOAL    :   Profilage taxonomique du microbiome sur données décontaminées.
#               0. Importation de la base Kraken2 depuis le S3 vers le Scratch.
#               1. Importation des FASTQ "clean" depuis le S3 vers le Scratch.
#               2. Assignation par k-mers via Kraken2 (avec lecture des .gz).
#               3. Ré-estimation bayésienne de l'abondance (Espèce) via Bracken.
#               4. Rapatriement des rapports vers le dossier local results/.
#######################################################################################

set -uo pipefail

# =====================================================================
# INITIALISATION ET VARIABLES
# =====================================================================
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] DÉMARRAGE DU PROFILAGE TAXONOMIQUE"
echo "Job Master ID : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}} | Task : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "====================================================================="

module purge
module load rclone/1.55.1

# Configuration S3
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_DB_DIR="ref/kraken2_pluspfp"
S3_IN_DIR="shotgun_cleaned/fastq"

# Configuration Locale
PROJECT_DIR="${PWD}" # S'adapte dynamiquement à ton dossier courant
APPTAINER_BIN="/usr/bin/apptainer"
SIF_IMG="/home/thboutet/projets/16_tenebrion/containers/kraken2_bracken.sif"

# Dossier de sortie pointant vers "results/"
LOCAL_OUT_DIR="${PROJECT_DIR}/results"
mkdir -p "${LOCAL_OUT_DIR}/kraken2_pluspfp" "${LOCAL_OUT_DIR}/bracken_pluspfp"

THREADS="${SLURM_CPUS_PER_TASK}"
READ_LEN=150 # Longueur attendue des reads pour Bracken

# =====================================================================
# GESTION DE L'ARRAY SLURM DYNAMIQUE (Depuis le S3)
# =====================================================================
mapfile -t SAMPLES < <(rclone ls "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/" | awk '{print $2}' | grep "_clean_R1.fastq.gz" | sed 's/_clean_R1.fastq.gz//' | sort)

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] ${N_SAMPLES} échantillons détectés sur le S3."
    
    # On récupère le chemin exact du script soumis pour que le master relance l'array proprement
    SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | awk '{print $1}')
    sbatch --array=0-$((N_SAMPLES - 1)) "${SCRIPT_PATH}"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# =====================================================================
# ÉTAPE 0 : VÉRIFICATION ET IMPORTATION DE LA BASE DEPUIS LE S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Préparation de la base Kraken2 sur le Scratch..."

# MODIFICATION: Dossier exclusif pour la base PlusPFP
KRAKEN_DB_SCRATCH="/storage/scratch/${USER}/kraken2_pfp_db_shared"
mkdir -p "${KRAKEN_DB_SCRATCH}"
FLAG_FILE="${KRAKEN_DB_SCRATCH}/.sync_done"

# Le verrou englobe tout le bloc (Nouveau verrou pour la base PFP).
exec 9> "/storage/scratch/${USER}/kraken_pfp_s3_sync.lock"
flock -x 9

# Si le fichier témoin (.sync_done) n'existe pas, c'est que la base n'est pas encore copiée à 100%
if [ ! -f "${FLAG_FILE}" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Base absente du Scratch. Importation depuis le S3..."
    rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_DB_DIR}/" "${KRAKEN_DB_SCRATCH}/" -P --transfers 8
    
    # On pose le fichier drapeau UNIQUEMENT quand rclone a fini son travail sans erreur
    touch "${FLAG_FILE}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Importation de la base réussie."
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Base de données trouvée et complète sur le Scratch."
fi

flock -u 9

# =====================================================================
# GESTION DU SCRATCH (Espace de travail ISOLÉ par échantillon)
# =====================================================================
JOB_SCRATCH="/storage/scratch/${USER}/kraken_sample_${SAMPLE}_${SLURM_JOB_ID}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

cleanup() {
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Nettoyage des FASTQ temporaires sur le Scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
    echo "====================================================================="
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 1 : IMPORTATION DES DONNÉES DEPUIS S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Importation de ${SAMPLE} depuis le S3..."

rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_clean_R1.fastq.gz" ./
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_clean_R2.fastq.gz" ./

# =====================================================================
# ÉTAPE 2 : ASSIGNATION TAXONOMIQUE (KRAKEN2)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Exécution de Kraken2..."

${APPTAINER_BIN} exec --bind "${KRAKEN_DB_SCRATCH}:${KRAKEN_DB_SCRATCH}" --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2 \
    --db "${KRAKEN_DB_SCRATCH}" \
    --threads "${THREADS}" \
    --paired \
    --gzip-compressed \
    --use-names \
    --report "${SAMPLE}_pluspfp.kreport" \
    --output /dev/null \
    "${SAMPLE}_clean_R1.fastq.gz" "${SAMPLE}_clean_R2.fastq.gz"

# =====================================================================
# ÉTAPE 3 : RÉ-ESTIMATION DES ABONDANCES (BRACKEN)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] Exécution de Bracken (Niveau Espèce)..."

${APPTAINER_BIN} exec --bind "${KRAKEN_DB_SCRATCH}:${KRAKEN_DB_SCRATCH}" --bind "${PWD}:${PWD}" "${SIF_IMG}" bracken \
    -d "${KRAKEN_DB_SCRATCH}" \
    -i "${SAMPLE}_pluspfp.kreport" \
    -o "${SAMPLE}_pluspfp.bracken" \
    -w "${SAMPLE}_pluspfp_bracken_species.kreport" \
    -r "${READ_LEN}" \
    -l S

# =====================================================================
# ÉTAPE 4 : EXPORTATION VERS LE /HOME
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 4] Sauvegarde des résultats dans le dossier local results/..."

cp "${SAMPLE}_pluspfp.kreport" "${LOCAL_OUT_DIR}/kraken2_pluspfp/"
cp "${SAMPLE}_pluspfp.bracken" "${LOCAL_OUT_DIR}/bracken_pluspfp/"
# MODIFICATION: Correction du nom du fichier pour correspondre à l'Étape 3
cp "${SAMPLE}_pluspfp_bracken_species.kreport" "${LOCAL_OUT_DIR}/bracken_pluspfp/"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] Analyse taxonomique terminée pour ${SAMPLE}."
