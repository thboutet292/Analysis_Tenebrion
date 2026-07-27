#!/bin/bash
#SBATCH --job-name=COVERM_PROD
#SBATCH --output=log/COVERM_PROD_%A_%a.out
#SBATCH --error=log/COVERM_PROD_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal,lmge
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   shotgun_coverm.slurm
#
#   GOAL    :   PHASE 2 EN PRODUCTION (JOB ARRAY).
#               1. [MUTEX] Création du catalogue global non-redondant (MMseqs2).
#                  (Fusion de tous les contigs MEGAHIT >= 500 bp).
#               2. [ARRAY] Alignement des reads complets sur le catalogue.
#               3. Calcul des abondances (trimmed_mean, count) et sauvegarde des BAMs.
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE PHASE 2 (PRODUCTION)"
echo "Job Master ID : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}} | Task : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "====================================================================="

module purge
module load MMseqs2/13-45111
module load rclone/1.55.1
module load apptainer 2>/dev/null || true

# =====================================================================
# CONFIGURATION DES CHEMINS
# =====================================================================
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_CLEAN_DIR="shotgun_cleaned/fastq"

DIR_MEGAHIT="${PROJECT_DIR}/results/assembly_megahit"
DIR_CATALOG="${PROJECT_DIR}/results/global_catalog"
DIR_COVERM="${PROJECT_DIR}/results/coverm"
IMG_COVERM="${PROJECT_DIR}/containers/coverm_v0.7.0.sif"

mkdir -p "${DIR_CATALOG}" "${DIR_COVERM}/bams" log/

THREADS="${SLURM_CPUS_PER_TASK:-16}"

MERGED_FASTA="${DIR_CATALOG}/all_megahit_contigs.fasta"
CATALOG_FASTA="${DIR_CATALOG}/non_redundant_catalog.fasta"
LOCK_FILE="${DIR_CATALOG}/.catalog_done.lock"

# =====================================================================
# DÉTECTION S3 ET SOUMISSION DE L'ARRAY
# =====================================================================
mapfile -t SAMPLES < <(rclone ls "${S3_REMOTE}:${S3_BUCKET}/${S3_CLEAN_DIR}/" \
    | awk '{print $2}' \
    | grep "_clean_R1.fastq.gz" \
    | sed 's/_clean_R1.fastq.gz//' \
    | sort)

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    if [ ${N_SAMPLES} -eq 0 ]; then
        echo "[ERREUR] Aucun échantillon trouvé sur le S3." >&2
        exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] Soumission de la quantification pour ${N_SAMPLES} échantillons."
    SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | awk '{print $1}')
    sbatch --array=0-$((N_SAMPLES - 1)) "${SCRIPT_PATH}"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# =====================================================================
# ÉTAPE 1 : CRÉATION DU CATALOGUE GLOBAL (Tâche Leader / Mutex)
# =====================================================================
if [ ! -f "${LOCK_FILE}" ]; then
    if mkdir "${DIR_CATALOG}/.lockdir" 2>/dev/null; then
        echo "---------------------------------------------------------------------"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Construction du Catalogue..."
        
        # 1.A : Concaténation avec RENOMMAGE pour éviter les doublons
        echo " -> Concaténation, Renommage et filtrage (>= 1500 bp)..."
        
        # On vide le fichier au cas où
        > "${MERGED_FASTA}"
        
        for f in "${DIR_MEGAHIT}"/*_megahit.contigs.fa; do
            sname=$(basename "$f" _megahit.contigs.fa)
            # On ajoute le nom de l'échantillon au nom du contig : >k141_1 devient >AR_E40_k141_1
            sed "s/>/>${sname}_/" "$f" >> "${MERGED_FASTA}"
        done

        # Filtrage par taille sur le fichier fusionné
        # (On utilise un fichier temporaire pour ne pas saturer la RAM)
        mv "${MERGED_FASTA}" "${MERGED_FASTA}.tmp"
        awk 'BEGIN{RS=">";FS="\n"} NR>1{seq=""; for(i=2;i<=NF;i++) seq=seq$i; if(length(seq)>=1500) printf ">%s\n%s\n", $1, seq}' "${MERGED_FASTA}.tmp" > "${MERGED_FASTA}"
        rm "${MERGED_FASTA}.tmp"
        
        # 1.B : Clustering MMseqs2
        echo " -> Déréplication avec MMseqs2..."
        MMSEQS_TMP="/storage/scratch/${USER}/mmseqs_tmp_${SLURM_JOB_ID}"
        mkdir -p "${MMSEQS_TMP}"

        mmseqs easy-linclust "${MERGED_FASTA}" "${DIR_CATALOG}/nr_catalog" "${MMSEQS_TMP}" \
            --min-seq-id 0.95 -c 0.90 --cov-mode 1 --threads "${THREADS}"
            
        mv "${DIR_CATALOG}/nr_catalog_rep_seq.fasta" "${CATALOG_FASTA}"
        
        # Nettoyage
        rm -rf "${MMSEQS_TMP}" "${DIR_CATALOG}/nr_catalog_all_seqs.fasta" "${DIR_CATALOG}/nr_catalog_cluster.tsv" "${MERGED_FASTA}"
        rmdir "${DIR_CATALOG}/.lockdir"
        touch "${LOCK_FILE}"
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Catalogue unique créé."
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Attente du catalogue..."
        while [ ! -f "${LOCK_FILE}" ]; do sleep 30; done
    fi
fi

# =====================================================================
# ÉTAPE 2 : TÉLÉCHARGEMENT COMPLET DES READS ET QUANTIFICATION
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Traitement complet de l'échantillon : ${SAMPLE}"

# Scratch local massif (les reads font plusieurs Go)
JOB_SCRATCH="/storage/scratch/${USER}/coverm_${SAMPLE}_${SLURM_JOB_ID}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Nettoyage du scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

echo " -> Téléchargement des reads complets de ${SAMPLE} depuis le S3..."
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_CLEAN_DIR}/${SAMPLE}_clean_R1.fastq.gz" ./
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_CLEAN_DIR}/${SAMPLE}_clean_R2.fastq.gz" ./

R1="${SAMPLE}_clean_R1.fastq.gz"
R2="${SAMPLE}_clean_R2.fastq.gz"

echo " -> Lancement de CoverM (avec redirection du /tmp vers le scratch)..."

# On force les variables d'environnement temporaires
export TMPDIR="${JOB_SCRATCH}"
export APPTAINER_TMPDIR="${JOB_SCRATCH}"

apptainer exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
    --bind "${JOB_SCRATCH}:/tmp" \
    "${IMG_COVERM}" \
    coverm contig \
        --reference "${CATALOG_FASTA}" \
        -1 "${R1}" -2 "${R2}" \
        --mapper minimap2-sr \
        --methods trimmed_mean count \
        --min-read-percent-identity 0.95 \
        --threads "${THREADS}" \
        --output-file "${DIR_COVERM}/${SAMPLE}_coverage.tsv" \
        --bam-file-cache-directory "${DIR_COVERM}/bams"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FIN DU TRAITEMENT POUR : ${SAMPLE}"
echo "====================================================================="