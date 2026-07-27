#!/bin/bash
#SBATCH --job-name=KRAKEN_RAW
#SBATCH --output=log/TAXONOMY_RAW_%A_%a.out
#SBATCH --error=log/TAXONOMY_RAW_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G                    
#SBATCH --partition=normal,lmge
#SBATCH --time=06:00:00

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE) / Assistant
#   DATE    :   26 - 03 - 2026
#   SCRIPT  :   shotgun_kraken2_raw.slurm
#   TOOLS   :   rclone ; Apptainer (Kraken2 / Bracken)
#
#   GOAL    :   Profilage taxonomique du microbiome sur données BRUTES.
#               0. Importation de la base Kraken2 depuis le S3 vers le Scratch.
#               1. Importation des FASTQ bruts depuis le S3 vers le Scratch.
#               2. Assignation par k-mers via Kraken2 (avec lecture des .gz).
#               3. Ré-estimation bayésienne de l'abondance (Espèce) via Bracken.
#               4. Rapatriement des rapports vers le dossier local results/.
#######################################################################################

set -uo pipefail

# =====================================================================
# INITIALISATION ET VARIABLES
# =====================================================================
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] DÉMARRAGE DU PROFILAGE TAXONOMIQUE (BRUT)"
echo "Job Master ID : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}} | Task : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "====================================================================="

module purge
module load rclone/1.55.1

# --- Configuration S3 (Pour les inputs) ---
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_DB_DIR="ref/kraken2_standard"
S3_IN_DIR="shotgun_bruts" 

# --- Configuration Locale ---
PROJECT_DIR="${PWD}"
APPTAINER_BIN="/usr/bin/apptainer"
SIF_IMG="/home/thboutet/projets/16_tenebrion/containers/kraken2_bracken.sif"

# --- Dossiers de sortie ---

LOCAL_OUT_DIR="${PROJECT_DIR}/results"
mkdir -p "${LOCAL_OUT_DIR}/kraken2_raw" "${LOCAL_OUT_DIR}/bracken_raw"

THREADS="${SLURM_CPUS_PER_TASK}"
READ_LEN=150 # Longueur attendue des reads pour Bracken

# =====================================================================
# GESTION DE L'ARRAY SLURM DYNAMIQUE (Depuis le S3)
# =====================================================================
# MODIFICATION de la logique de parsing :
# 1. 'grep "_R1.fastq.gz"' filtre uniquement les reads R1 bruts (évite de compter en double avec R2).
# 2. 'sed 's/_R1.fastq.gz//'' supprime la fin de la chaîne pour ne garder que l'ID (ex: AR_E40).
mapfile -t SAMPLES < <(rclone ls "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/" | awk '{print $2}' | grep "_R1.fastq.gz" | sed 's/_R1.fastq.gz//' | sort)

# Bloc de soumission du Master Job (si pas déjà dans un array)
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] ${N_SAMPLES} échantillons bruts détectés sur le S3."
    
    # Récupère le chemin exact du script pour relancer l'array
    SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | awk '{print $1}')
    sbatch --array=0-$((N_SAMPLES - 1)) "${SCRIPT_PATH}"
    exit 0
fi

# Récupération du nom de l'échantillon pour la tâche courante
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# =====================================================================
# ÉTAPE 0 : VÉRIFICATION ET IMPORTATION DE LA BASE DEPUIS LE S3
# =====================================================================
# Cette étape utilise un système de verrouillage (flock). C'est une excellente
# pratique d'optimisation (complexité spatiale/réseau) : cela évite que 20 jobs 
# téléchargent la même énorme base Kraken2 en même temps sur le scratch.
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Préparation de la base Kraken2 sur le Scratch..."

KRAKEN_DB_SCRATCH="/storage/scratch/${USER}/kraken2_db_shared"
mkdir -p "${KRAKEN_DB_SCRATCH}"
FLAG_FILE="${KRAKEN_DB_SCRATCH}/.sync_done"

exec 9> "/storage/scratch/${USER}/kraken_s3_sync.lock"
flock -x 9

if [ ! -f "${FLAG_FILE}" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Base absente du Scratch. Importation depuis le S3..."
    rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_DB_DIR}/" "${KRAKEN_DB_SCRATCH}/" -P --transfers 8
    touch "${FLAG_FILE}"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Importation de la base réussie."
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Base de données trouvée et complète sur le Scratch."
fi

flock -u 9 # Déverrouillage pour les autres jobs de l'array

# =====================================================================
# GESTION DU SCRATCH (Espace de travail ISOLÉ par échantillon)
# =====================================================================
# Création d'un dossier temporaire unique par job pour éviter les collisions
JOB_SCRATCH="/storage/scratch/${USER}/kraken_sample_${SAMPLE}_${SLURM_JOB_ID}_RAW"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

cleanup() {
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Nettoyage des FASTQ bruts temporaires sur le Scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
    echo "====================================================================="
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 1 : IMPORTATION DES DONNÉES DEPUIS S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Importation des bruts de ${SAMPLE} depuis le S3..."

# MODIFICATION : Retrait de "_clean" dans les noms de fichiers recherchés
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_R1.fastq.gz" ./
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_R2.fastq.gz" ./

# =====================================================================
# ÉTAPE 2 : ASSIGNATION TAXONOMIQUE (KRAKEN2)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Exécution de Kraken2..."

# Rappel pédagogique : Kraken2 utilise des k-mers exacts. S'il y a ambiguïté 
# (un k-mer appartient à plusieurs espèces), il remonte à l'ancêtre commun le plus récent (LCA).
${APPTAINER_BIN} exec --bind "${KRAKEN_DB_SCRATCH}:${KRAKEN_DB_SCRATCH}" --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2 \
    --db "${KRAKEN_DB_SCRATCH}" \
    --threads "${THREADS}" \
    --paired \
    --gzip-compressed \
    --use-names \
    --report "${SAMPLE}_raw.kreport" \
    --output /dev/null \
    "${SAMPLE}_R1.fastq.gz" "${SAMPLE}_R2.fastq.gz" # MODIFICATION : inputs sans "_clean"

# =====================================================================
# ÉTAPE 3 : RÉ-ESTIMATION DES ABONDANCES (BRACKEN)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] Exécution de Bracken (Niveau Espèce)..."

# Rappel pédagogique : Bracken prend le rapport de Kraken2 et redistribue 
# mathématiquement (approche bayésienne) les reads assignés à des rangs supérieurs 
# (Genre, Famille) vers le rang Espèce, en se basant sur la probabilité de trouver 
# chaque espèce dans la base de données.
${APPTAINER_BIN} exec --bind "${KRAKEN_DB_SCRATCH}:${KRAKEN_DB_SCRATCH}" --bind "${PWD}:${PWD}" "${SIF_IMG}" bracken \
    -d "${KRAKEN_DB_SCRATCH}" \
    -i "${SAMPLE}_raw.kreport" \
    -o "${SAMPLE}_raw.bracken" \
    -w "${SAMPLE}_raw_species.kreport" \
    -r "${READ_LEN}" \
    -l S

# =====================================================================
# ÉTAPE 4 : EXPORTATION VERS LE /HOME
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 4] Sauvegarde des résultats dans le dossier local..."

# MODIFICATION : Envoi vers les sous-dossiers spécifiques "_raw" créés en haut de script
cp "${SAMPLE}_raw.kreport" "${LOCAL_OUT_DIR}/kraken2_raw/"
cp "${SAMPLE}_raw.bracken" "${LOCAL_OUT_DIR}/bracken_raw/"
cp "${SAMPLE}_raw_species.kreport" "${LOCAL_OUT_DIR}/bracken_raw/"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] Analyse taxonomique terminée pour ${SAMPLE} (Données brutes)."
