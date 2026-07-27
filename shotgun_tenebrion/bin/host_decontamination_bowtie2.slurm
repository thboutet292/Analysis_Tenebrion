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
#   DATE    :   Mai 2026
#   SCRIPT  :   host_decontamination_bowtie2.slurm
#   TOOLS   :   rclone 1.55.1 ; bowtie2 2.3.4.3 ; samtools 1.16.1
#
#   GOAL    :   Retrait des lectures de l'hôte (Tenebrio molitor) pour le shotgun.
#               - Construction "Safe" de l'index Bowtie2 (avec Mutex File Lock).
#               - Alignement short-reads via Bowtie2 (--very-sensitive).
#               - Extraction stricte du microbiome (lectures non-mappées, flag -f 12).
#               - Exportation vers S3 (shotgun_cleaned/fastq).
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] DÉMARRAGE DE LA DÉCONTAMINATION v5 (Auto-Index)"
echo "Job Master ID  : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}"
echo "Task Array ID  : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Nœud de calcul : ${SLURMD_NODENAME}"
echo "====================================================================="

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] Chargement des modules HPC..."
module purge
module load rclone/1.55.1
module load bowtie2/2.3.4.3
module load gcc/8.1.0
module load samtools/1.16.1

# =====================================================================
# CONFIGURATION DES CHEMINS ET VARIABLES
# =====================================================================
THREADS="${SLURM_CPUS_PER_TASK:-16}"

S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_RAW_DIR="shotgun_bruts"
S3_OUT_DIR="shotgun_cleaned"
S3_REF_FNA="ref/genome/GCA_963966145.1_icTenMoli1.1_genomic.fna"

PROJECT_DIR="/home/thboutet/projets/shotgun_tenebrion"
RES_DIR="${PROJECT_DIR}/resources"
REF_INDEX="${RES_DIR}/tenebrio_molitor_v1.1"

# Dossier servant de verrou atomique pour l'indexation
LOCK_DIR="${RES_DIR}/bowtie2_index.lock"

mkdir -p "${RES_DIR}"

# =====================================================================
# GESTION DE L'ARRAY SLURM (DÉTECTION DYNAMIQUE)
# =====================================================================
mapfile -t SAMPLES < <(rclone ls "${S3_REMOTE}:${S3_BUCKET}/${S3_RAW_DIR}/" \
    | awk '{print $2}' \
    | grep "_R1.fastq.gz" \
    | sed 's/_R1.fastq.gz//' \
    | sort)

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    if [ ${N_SAMPLES} -eq 0 ]; then
        echo "[ERREUR] Aucun échantillon détecté sur le S3." >&2
        exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] ${N_SAMPLES} échantillons détectés."
    SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | awk '{print $1}')
    sbatch --array=0-$((N_SAMPLES - 1)) "${SCRIPT_PATH}"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# =====================================================================
# GESTION DU SCRATCH LOCAL
# =====================================================================
JOB_SCRATCH="/storage/scratch/${USER}/decontam_${SAMPLE}_${SLURM_JOB_ID}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

cleanup() {
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression du scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] Tâche terminée."
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 0 : VÉRIFICATION ET CONSTRUCTION SÉCURISÉE DE L'INDEX
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Vérification de l'index Bowtie2..."

if [ ! -f "${REF_INDEX}.1.bt2" ]; then
    # Tentative de création du verrou. mkdir est atomique : un seul processus réussira.
    if mkdir "${LOCK_DIR}" 2>/dev/null; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Index absent. Ce processus prend le lead pour l'indexation..."
        
        echo " -> Téléchargement du génome hôte depuis S3..."
        rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_REF_FNA}" "${RES_DIR}/"
        GENOME_LOCAL="${RES_DIR}/$(basename ${S3_REF_FNA})"
        
        echo " -> Construction de l'index Bowtie2 (patientez)..."
        bowtie2-build --threads "${THREADS}" "${GENOME_LOCAL}" "${REF_INDEX}"
        
        echo " -> Nettoyage du FASTA et libération du verrou..."
        rm "${GENOME_LOCAL}"
        rmdir "${LOCK_DIR}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Indexation terminée avec succès."
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Un autre processus construit actuellement l'index."
        echo " -> Mise en attente de ce job..."
        
        # Boucle d'attente : on vérifie toutes les 30 secondes si le verrou a disparu
        while [ -d "${LOCK_DIR}" ] || [ ! -f "${REF_INDEX}.1.bt2" ]; do
            sleep 30
        done
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] L'index est maintenant disponible ! Poursuite de l'analyse."
    fi
else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 0] Index déjà disponible."
fi

# =====================================================================
# ÉTAPE 1 : IMPORTATION DES DONNÉES DEPUIS LE S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Importation des FASTQ de l'échantillon ${SAMPLE}..."

rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_RAW_DIR}/${SAMPLE}_R1.fastq.gz" ./
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_RAW_DIR}/${SAMPLE}_R2.fastq.gz" ./

# =====================================================================
# ÉTAPE 2 : MAPPING CONTRE L'HÔTE (Bowtie2)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Mapping contre le génome Ténébrion..."

bowtie2 \
    -x "${REF_INDEX}" \
    -1 "${SAMPLE}_R1.fastq.gz" \
    -2 "${SAMPLE}_R2.fastq.gz" \
    --very-sensitive \
    --threads "${THREADS}" \
    | samtools view -@ "${THREADS}" -b -o host_mapped.bam

# =====================================================================
# ÉTAPE 3 : FILTRAGE ET EXTRACTION DU MICROBIOME
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] Extraction des reads non-mappés (Microbiome)..."

samtools view -@ "${THREADS}" -b -f 12 host_mapped.bam > microbiome_clean.bam

samtools fastq -@ "${THREADS}" \
    -1 "${SAMPLE}_clean_R1.fastq.gz" \
    -2 "${SAMPLE}_clean_R2.fastq.gz" \
    -0 /dev/null -s /dev/null -n \
    microbiome_clean.bam

# =====================================================================
# ÉTAPE 4 : CALCUL DES STATISTIQUES
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 4] Calcul des métriques de filtrage..."

TOTAL_READS=$(samtools view -c host_mapped.bam)
MICRO_READS=$(samtools view -c microbiome_clean.bam)
HOST_READS=$((TOTAL_READS - MICRO_READS))

PERC_HOST=$(awk "BEGIN {printf \"%.2f\", ($HOST_READS/$TOTAL_READS)*100}")
PERC_MICRO=$(awk "BEGIN {printf \"%.2f\", ($MICRO_READS/$TOTAL_READS)*100}")

REPORT="${SAMPLE}_decontam_report.txt"
echo -e "Sample\tTotal_Reads\tHost_Reads\tHost_%\tMicrobiome_Reads\tMicrobiome_%" > "${REPORT}"
echo -e "${SAMPLE}\t${TOTAL_READS}\t${HOST_READS}\t${PERC_HOST}\t${MICRO_READS}\t${PERC_MICRO}" >> "${REPORT}"

cat "${REPORT}"

# =====================================================================
# ÉTAPE 5 : EXPORTATION VERS LE S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 5] Exportation des données propres vers le S3..."

rclone copy "${SAMPLE}_clean_R1.fastq.gz" "${S3_REMOTE}:${S3_BUCKET}/${S3_OUT_DIR}/fastq/"
rclone copy "${SAMPLE}_clean_R2.fastq.gz" "${S3_REMOTE}:${S3_BUCKET}/${S3_OUT_DIR}/fastq/"
rclone copy "${REPORT}" "${S3_REMOTE}:${S3_BUCKET}/${S3_OUT_DIR}/stats/"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] PIPELINE TERMINÉ AVEC SUCCÈS — ${SAMPLE}"
echo "====================================================================="