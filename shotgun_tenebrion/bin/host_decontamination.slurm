#!/bin/bash
#SBATCH --job-name=SHOTGUN_DECONTAM
#SBATCH --output=log/DECONTAM_%A_%a.out
#SBATCH --error=log/DECONTAM_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal,lmge
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   DATE    :   17 - 03 - 2026
#   SCRIPT  :   host_decontamination_v2.slurm
#   TOOLS   :   rclone 1.55.1 ; minimap2 2.24 ; samtools 1.15
#
#   GOAL    :   Retrait des lectures de l'hôte (Tenebrio molitor) pour le shotgun.
#               - Gestion dynamique du génome de référence via S3.
#               - Mapping short-reads via Minimap2 (algorithme optimisé).
#               - Extraction stricte du microbiome (lectures non-mappées).
#               - Rapport statistique automatisé (Comptage et Ratio).
#######################################################################################

set -euo pipefail

# =====================================================================
# INITIALISATION ET CHARGEMENT DES OUTILS
# =====================================================================
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] DÉMARRAGE DE LA DÉCONTAMINATION"
echo "Job Master ID  : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}"
echo "Task Array ID  : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Nœud de calcul : ${SLURMD_NODENAME}"
echo "====================================================================="

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] Chargement des modules HPC..."
module purge
module load rclone/1.55.1
module load minimap2/2.24
module load gcc/8.1.0
module load samtools/1.16.1

# =====================================================================
# CONFIGURATION DES CHEMINS ET VARIABLES (S3 & LOCAL)
# =====================================================================
# Paramètres S3
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_RAW_DIR="shotgun_bruts"
S3_OUT_DIR="shotgun_cleaned"
S3_REF_FNA="ref/genome/GCA_963966145.1_icTenMoli1.1_genomic.fna"

# Paramètres Projets
PROJECT_DIR="/home/thboutet/projets/16_tenebrion"
RES_DIR="${PROJECT_DIR}/resources"
REF_INDEX="${RES_DIR}/tenebrio_molitor_v1.1.mmi"
THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "${RES_DIR}" "log"

# =====================================================================
# GESTION DE L'ARRAY SLURM (DÉTECTION DYNAMIQUE SUR S3)
# =====================================================================
# On liste le S3 pour identifier les paires R1
mapfile -t SAMPLES < <(rclone ls "${S3_REMOTE}:${S3_BUCKET}/${S3_RAW_DIR}/" | awk '{print $2}' | grep "_R1.fastq.gz" | sed 's/_R1.fastq.gz//' | sort)

# Auto-soumission de l'Array si lancé en mode script simple
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] ${N_SAMPLES} échantillons détectés."
    sbatch --array=0-$((N_SAMPLES - 1)) "$0"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# =====================================================================
# GESTION DU SCRATCH LOCAL (ACCÉLÉRATION I/O)
# =====================================================================
JOB_SCRATCH="/storage/scratch/${USER}/decontam_${SAMPLE}_${SLURM_JOB_ID}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

# Fonction de nettoyage automatique
cleanup() {
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression du scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] Tâche terminée."
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 0 : PRÉPARATION DU GÉNOME HÔTE (S3 -> INDEX MMI)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Vérification de l'index hôte..."

# Cette section évite que 44 jobs n'indexent le même fichier en même temps
if [ ! -f "${REF_INDEX}" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Index non trouvé. Importation et création..."
    rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_REF_FNA}" "./"
    GENOME_LOCAL=$(basename "${S3_REF_FNA}")
    
    # Création de l'index Minimap2 (.mmi) pour optimiser la RAM et le temps CPU
    minimap2 -t "${THREADS}" -d "${REF_INDEX}" "${GENOME_LOCAL}"
    rm "${GENOME_LOCAL}"
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Index déjà disponible : ${REF_INDEX}"
fi

# =====================================================================
# ÉTAPE 1 : IMPORTATION DES DONNÉES DEPUIS LE S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Importation des FASTQ de l'échantillon ${SAMPLE}..."

rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_RAW_DIR}/${SAMPLE}_R1.fastq.gz" ./
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_RAW_DIR}/${SAMPLE}_R2.fastq.gz" ./

# =====================================================================
# ÉTAPE 2 : MAPPING CONTRE L'HÔTE (Minimap2)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Mapping contre le génome Ténébrion (Minimap2 -sr)..."

# -ax sr : Mode optimisé pour short-reads (mapping global)
# samtools view -bS : Conversion directe en flux vers le format BAM compressé
minimap2 -ax sr -t "${THREADS}" "${REF_INDEX}" "${SAMPLE}_R1.fastq.gz" "${SAMPLE}_R2.fastq.gz" \
    | samtools view -@ "${THREADS}" -bS -o host_mapped.bam

# =====================================================================
# ÉTAPE 3 : FILTRAGE ET EXTRACTION DU MICROBIOME
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] Extraction des reads non-mappés (Microbiome)..."

# Utilisation du flag SAM -f 12 :
# - 4 : Read unmapped
# - 8 : Mate unmapped
# On ne garde que les paires où AUCUN des deux reads n'est lié à l'hôte.
samtools view -@ "${THREADS}" -b -f 12 host_mapped.bam > microbiome_clean.bam

# Conversion du BAM filtré en fichiers FASTQ compressés
samtools fastq -@ "${THREADS}" \
    -1 "${SAMPLE}_clean_R1.fastq.gz" \
    -2 "${SAMPLE}_clean_R2.fastq.gz" \
    -0 /dev/null -s /dev/null -n microbiome_clean.bam

# =====================================================================
# ÉTAPE 4 : CALCUL DES STATISTIQUES DE DÉCONTAMINATION
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 4] Calcul des métriques de filtrage..."

# Comptage rapide des reads via Samtools view -c
TOTAL_READS=$(samtools view -c host_mapped.bam)
MICRO_READS=$(samtools view -c microbiome_clean.bam)
HOST_READS=$((TOTAL_READS - MICRO_READS))

# Calcul des pourcentages via AWK (gestion des flottants)
PERC_HOST=$(awk "BEGIN {printf \"%.2f\", ($HOST_READS/$TOTAL_READS)*100}")
PERC_MICRO=$(awk "BEGIN {printf \"%.2f\", ($MICRO_READS/$TOTAL_READS)*100}")

# Génération du rapport de décontamination
REPORT="${SAMPLE}_decontam_report.txt"
echo -e "Sample\tTotal_Reads\tHost_Reads\tHost_%\tMicrobiome_Reads\tMicrobiome_%" > "${REPORT}"
echo -e "${SAMPLE}\t${TOTAL_READS}\t${HOST_READS}\t${PERC_HOST}\t${MICRO_READS}\t${PERC_MICRO}" >> "${REPORT}"

cat "${REPORT}"

# =====================================================================
# ÉTAPE 5 : EXPORTATION DES LIVRABLES VERS LE S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 5] Exportation des données propres vers le S3..."

rclone copy "${SAMPLE}_clean_R1.fastq.gz" "${S3_REMOTE}:${S3_BUCKET}/${S3_OUT_DIR}/fastq/"
rclone copy "${SAMPLE}_clean_R2.fastq.gz" "${S3_REMOTE}:${S3_BUCKET}/${S3_OUT_DIR}/fastq/"
rclone copy "${REPORT}" "${S3_REMOTE}:${S3_BUCKET}/${S3_OUT_DIR}/stats/"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] PIPELINE TERMINÉ AVEC SUCCÈS"
echo "====================================================================="