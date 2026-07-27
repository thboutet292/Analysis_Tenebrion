#!/bin/bash
#SBATCH --job-name=COVERM_QUANT
#SBATCH --output=log/COVERM_QUANT_%A_%a.out
#SBATCH --error=log/COVERM_QUANT_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal,lmge
#SBATCH --time=08:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   coverm_quantification.slurm
#   TOOLS   :   rclone ; Apptainer ; CoverM v0.7.0
#
#   GOAL    :   Quantification des reads sur les contigs assemblés (MEGAHIT et
#               metaSPAdes) pour chaque échantillon du métagénome de Tenebrio molitor.
#
#               CoverM remappes les reads bruts (clean) sur les contigs et calcule
#               la profondeur de couverture par contig — proxy direct de l'abondance
#               relative de chaque séquence dans l'échantillon.
#
#               Deux seuils de filtrage des contigs :
#               - ≥ 500 pb  → pour CoverM (maximise le signal d'abondance)
#               - ≥ 1000 pb → copie filtrée pour Kraken2 (meilleure classification)
#
#               Métriques CoverM calculées :
#               - mean            → profondeur moyenne par contig
#               - trimmed_mean    → profondeur robuste (exclut les 5% extrêmes)
#               - covered_fraction → fraction du contig effectivement couverte
#               - rpkm            → abondance normalisée (Reads Per Kilobase per Million)
#
#   STEPS   :   1. Détection dynamique des reads clean sur S3.
#               2. Filtrage des contigs (≥ 500 pb et ≥ 1000 pb) via awk.
#               3. CoverM sur les contigs MEGAHIT filtrés (≥ 500 pb).
#               4. CoverM sur les contigs metaSPAdes filtrés (≥ 500 pb).
#               5. Export des tableaux de couverture + contigs filtrés ≥ 1000 pb.
#
#   NOTE    :   Lancer pull_coverm.slurm avant ce script si l'image n'est pas
#               encore présente dans containers/coverm_v0.7.0.sif
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] DÉMARRAGE COVERM QUANTIFICATION"
echo "Job Master ID : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}} | Task : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "====================================================================="

module purge
module load rclone/1.55.1
module load apptainer 2>/dev/null || true

# =====================================================================
# CONFIGURATION ET CHEMINS
# =====================================================================
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_IN_DIR="shotgun_cleaned/fastq"

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
CONTAINER="${PROJECT_DIR}/containers/coverm_v0.7.0.sif"
MEGAHIT_DIR="${PROJECT_DIR}/results/assembly_megahit"
SPADES_DIR="${PROJECT_DIR}/results/assembly_metaspades"
OUT_DIR="${PROJECT_DIR}/results/coverm"
FILTERED_DIR="${PROJECT_DIR}/results/contigs_filtered"
mkdir -p "${OUT_DIR}/megahit" "${OUT_DIR}/metaspades" \
         "${FILTERED_DIR}/megahit" "${FILTERED_DIR}/metaspades" log/

THREADS="${SLURM_CPUS_PER_TASK:-16}"
MIN_COVERM=500    # Seuil CoverM : maximise le signal d'abondance
MIN_KRAKEN=1000   # Seuil Kraken2 : meilleure classification taxonomique

APPTAINER_BIN=$(command -v apptainer || command -v singularity)

# Vérification de l'image
if [ ! -f "${CONTAINER}" ]; then
    echo "[ERREUR] Image CoverM introuvable : ${CONTAINER}" >&2
    echo "         Lancer pull_coverm.slurm d'abord." >&2
    exit 1
fi

# =====================================================================
# DÉTECTION DYNAMIQUE DES ÉCHANTILLONS SUR S3
# =====================================================================
mapfile -t SAMPLES < <(rclone ls "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/" \
    | awk '{print $2}' \
    | grep "_clean_R1.fastq.gz" \
    | sed 's/_clean_R1.fastq.gz//' \
    | sort)

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    if [ ${N_SAMPLES} -eq 0 ]; then
        echo "[ERREUR] Aucun échantillon détecté sur le S3." >&2
        exit 1
    fi
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] ${N_SAMPLES} échantillons détectés. Soumission en lot."
    SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | awk '{print $1}')
    sbatch --array=0-$((N_SAMPLES - 1)) "${SCRIPT_PATH}"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SAMPLE] Traitement de : ${SAMPLE}"

# =====================================================================
# GESTION DU SCRATCH ISOLÉ
# =====================================================================
JOB_SCRATCH="/storage/scratch/${USER}/coverm_${SAMPLE}_${SLURM_JOB_ID:-local}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression du répertoire temporaire..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 1 : IMPORTATION DES READS CLEAN DEPUIS S3
# =====================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Rapatriement des reads clean de ${SAMPLE}..."

rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_clean_R1.fastq.gz" ./
rclone copy "${S3_REMOTE}:${S3_BUCKET}/${S3_IN_DIR}/${SAMPLE}_clean_R2.fastq.gz" ./

R1="${JOB_SCRATCH}/${SAMPLE}_clean_R1.fastq.gz"
R2="${JOB_SCRATCH}/${SAMPLE}_clean_R2.fastq.gz"

# =====================================================================
# ÉTAPE 2 : FILTRAGE DES CONTIGS PAR LONGUEUR (awk)
# =====================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Filtrage des contigs par longueur..."

filter_contigs() {
    local INPUT="$1"
    local OUTPUT="$2"
    local MINLEN="$3"
    awk -v minlen="${MINLEN}" \
        '/^>/{header=$0; seq=""; next}
         {seq=seq$0}
         END{print header"\n"seq}  # flush last
         /^>/{if(seq!="" && length(seq)>=minlen) print prev_header"\n"seq; prev_header=header; seq=""}
         END{if(seq!="" && length(seq)>=minlen) print prev_header"\n"seq}' \
        "${INPUT}" > "${OUTPUT}" 2>/dev/null || \
    # Méthode alternative plus robuste avec awk multi-ligne
    awk -v minlen="${MINLEN}" 'BEGIN{RS=">"; FS="\n"}
        NR>1{
            seq=""; for(i=2;i<=NF;i++) seq=seq$i;
            if(length(seq)>=minlen) printf ">%s\n%s\n", $1, seq
        }' "${INPUT}" > "${OUTPUT}"
}

# MEGAHIT
MEGAHIT_FA="${MEGAHIT_DIR}/${SAMPLE}_megahit.contigs.fa"
MEGAHIT_500="${JOB_SCRATCH}/${SAMPLE}_megahit_min500.fa"
MEGAHIT_1000="${FILTERED_DIR}/megahit/${SAMPLE}_megahit_min1000.fa"

if [ -f "${MEGAHIT_FA}" ]; then
    echo "    [MEGAHIT] Filtrage ≥ ${MIN_COVERM} pb..."
    filter_contigs "${MEGAHIT_FA}" "${MEGAHIT_500}" "${MIN_COVERM}"
    echo "    [MEGAHIT] Filtrage ≥ ${MIN_KRAKEN} pb (pour Kraken2)..."
    filter_contigs "${MEGAHIT_FA}" "${MEGAHIT_1000}" "${MIN_KRAKEN}"
    echo "    [MEGAHIT] Contigs ≥${MIN_COVERM}pb : $(grep -c '^>' ${MEGAHIT_500}) | ≥${MIN_KRAKEN}pb : $(grep -c '^>' ${MEGAHIT_1000})"
else
    echo "    [WARN] Contigs MEGAHIT introuvables pour ${SAMPLE} : ${MEGAHIT_FA}"
fi

# metaSPAdes
SPADES_FA="${SPADES_DIR}/${SAMPLE}_metaspades_contigs.fa"
SPADES_500="${JOB_SCRATCH}/${SAMPLE}_metaspades_min500.fa"
SPADES_1000="${FILTERED_DIR}/metaspades/${SAMPLE}_metaspades_min1000.fa"

if [ -f "${SPADES_FA}" ]; then
    echo "    [metaSPAdes] Filtrage ≥ ${MIN_COVERM} pb..."
    filter_contigs "${SPADES_FA}" "${SPADES_500}" "${MIN_COVERM}"
    echo "    [metaSPAdes] Filtrage ≥ ${MIN_KRAKEN} pb (pour Kraken2)..."
    filter_contigs "${SPADES_FA}" "${SPADES_1000}" "${MIN_KRAKEN}"
    echo "    [metaSPAdes] Contigs ≥${MIN_COVERM}pb : $(grep -c '^>' ${SPADES_500}) | ≥${MIN_KRAKEN}pb : $(grep -c '^>' ${SPADES_1000})"
else
    echo "    [WARN] Contigs metaSPAdes introuvables pour ${SAMPLE} : ${SPADES_FA}"
fi

# =====================================================================
# ÉTAPE 3 : COVERM SUR CONTIGS MEGAHIT
# =====================================================================
if [ -f "${MEGAHIT_500}" ] && [ -s "${MEGAHIT_500}" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] CoverM sur contigs MEGAHIT ≥${MIN_COVERM}pb..."

    "${APPTAINER_BIN}" exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
        "${CONTAINER}" \
        coverm contig \
            -1 "${R1}" \
            -2 "${R2}" \
            --reference "${MEGAHIT_500}" \
            --methods mean trimmed_mean covered_fraction rpkm \
            --min-read-percent-identity 0.95 \
            --threads "${THREADS}" \
            --output-file "${OUT_DIR}/megahit/${SAMPLE}_megahit_coverage.tsv"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] CoverM MEGAHIT terminé."
    echo "    Sortie : ${OUT_DIR}/megahit/${SAMPLE}_megahit_coverage.tsv"
else
    echo "[WARN] Fichier MEGAHIT filtré vide ou absent, CoverM ignoré pour ${SAMPLE}."
fi

# =====================================================================
# ÉTAPE 4 : COVERM SUR CONTIGS METASPADES
# =====================================================================
if [ -f "${SPADES_500}" ] && [ -s "${SPADES_500}" ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 4] CoverM sur contigs metaSPAdes ≥${MIN_COVERM}pb..."

    "${APPTAINER_BIN}" exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
        "${CONTAINER}" \
        coverm contig \
            -1 "${R1}" \
            -2 "${R2}" \
            --reference "${SPADES_500}" \
            --methods mean trimmed_mean covered_fraction rpkm \
            --min-read-percent-identity 0.95 \
            --threads "${THREADS}" \
            --output-file "${OUT_DIR}/metaspades/${SAMPLE}_metaspades_coverage.tsv"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 4] CoverM metaSPAdes terminé."
    echo "    Sortie : ${OUT_DIR}/metaspades/${SAMPLE}_metaspades_coverage.tsv"
else
    echo "[WARN] Fichier metaSPAdes filtré vide ou absent, CoverM ignoré pour ${SAMPLE}."
fi

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] CoverM terminé pour ${SAMPLE}."
echo "  Tableaux de couverture : ${OUT_DIR}/"
echo "  Contigs ≥1000pb (Kraken2) : ${FILTERED_DIR}/"
echo "====================================================================="