#!/bin/bash
#SBATCH --job-name=METAQUAST_CMP
#SBATCH --output=log/METAQUAST_CMP_%A_%a.out
#SBATCH --error=log/METAQUAST_CMP_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=normal,lmge
#SBATCH --time=04:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   metaquast_comparison.slurm
#   TOOLS   :   Apptainer ; MetaQUAST v5.2.0
#
#   GOAL    :   Comparaison de la qualité des assemblages MEGAHIT vs metaSPAdes
#               pour chaque échantillon du métagénome de Tenebrio molitor.
#
#               Métriques clés évaluées :
#               - N50 / L50                    → longueur et nb de contigs "utiles"
#               - Nombre et taille des contigs → quantité vs qualité
#               - Misassemblies                → contigs mal assemblés (jonctions incorrectes)
#               - Contigs ≥ 500 / ≥ 1000 pb   → fraction biologiquement exploitable
#               - GC content                   → cohérence avec le microbiote attendu
#
#   STEPS   :   1. Détection dynamique des échantillons communs aux deux assembleurs.
#               2. Lancement MetaQUAST par échantillon (MEGAHIT vs metaSPAdes).
#               3. Lancement MetaQUAST global sur tous les échantillons combinés.
#               4. Export des rapports dans results/metaquast/
#
#   NOTE    :   Lancer pull_metaquast.slurm avant ce script si l'image n'est pas
#               encore présente dans containers/metaquast_v5.2.0.sif
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] DÉMARRAGE METAQUAST COMPARISON"
echo "Job Master ID : ${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}} | Task : ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "====================================================================="

module purge
module load apptainer 2>/dev/null || true

# =====================================================================
# CONFIGURATION ET CHEMINS
# =====================================================================
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
CONTAINER="${PROJECT_DIR}/containers/metaquast_v5.2.0.sif"
MEGAHIT_DIR="${PROJECT_DIR}/results/assembly_megahit"
SPADES_DIR="${PROJECT_DIR}/results/assembly_metaspades"
OUT_DIR="${PROJECT_DIR}/results/metaquast"
mkdir -p "${OUT_DIR}" log/

THREADS="${SLURM_CPUS_PER_TASK:-16}"
MIN_CONTIG=500  

APPTAINER_BIN=$(command -v apptainer || command -v singularity)

# Vérification de l'image
if [ ! -f "${CONTAINER}" ]; then
    echo "[ERREUR] Image MetaQUAST introuvable : ${CONTAINER}" >&2
    echo "         Lancer pull_metaquast.slurm d'abord." >&2
    exit 1
fi

# =====================================================================
# DÉTECTION DYNAMIQUE DES ÉCHANTILLONS COMMUNS AUX DEUX ASSEMBLEURS
# =====================================================================
mapfile -t MEGAHIT_SAMPLES < <(ls "${MEGAHIT_DIR}"/*_megahit.contigs.fa 2>/dev/null | \
    xargs -n1 basename | sed 's/_megahit\.contigs\.fa//' | sort)

mapfile -t SPADES_SAMPLES < <(ls "${SPADES_DIR}"/*_metaspades_contigs.fa 2>/dev/null | \
    xargs -n1 basename | sed 's/_metaspades_contigs\.fa//' | sort)

# Intersection : échantillons présents dans les deux assembleurs uniquement
mapfile -t SAMPLES < <(comm -12 \
    <(printf '%s\n' "${MEGAHIT_SAMPLES[@]}") \
    <(printf '%s\n' "${SPADES_SAMPLES[@]}"))

if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "[ERREUR] Aucun échantillon commun détecté entre MEGAHIT et metaSPAdes." >&2
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] ${#SAMPLES[@]} échantillons communs détectés :"
printf '    - %s\n' "${SAMPLES[@]}"

# =====================================================================
# GESTION DE L'ARRAY SLURM DYNAMIQUE
# =====================================================================
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    N_SAMPLES=${#SAMPLES[@]}
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ARRAY] Soumission en lot pour ${N_SAMPLES} échantillons."
    SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | awk -F= '/Command=/{print $2}' | awk '{print $1}')
    sbatch --array=0-$((N_SAMPLES - 1)) "${SCRIPT_PATH}"
    exit 0
fi

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

MEGAHIT_FA="${MEGAHIT_DIR}/${SAMPLE}_megahit.contigs.fa"
SPADES_FA="${SPADES_DIR}/${SAMPLE}_metaspades_contigs.fa"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [SAMPLE] Traitement de : ${SAMPLE}"
echo "    MEGAHIT    : ${MEGAHIT_FA}"
echo "    metaSPAdes : ${SPADES_FA}"

# =====================================================================
# GESTION DU SCRATCH ISOLÉ
# =====================================================================
JOB_SCRATCH="/storage/scratch/${USER}/metaquast_${SAMPLE}_${SLURM_JOB_ID:-local}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression du répertoire temporaire..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 1 : METAQUAST PAR ÉCHANTILLON — MEGAHIT vs metaSPAdes
# =====================================================================
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Lancement MetaQUAST pour ${SAMPLE}..."

SAMPLE_OUT="${OUT_DIR}/per_sample/${SAMPLE}"
mkdir -p "${SAMPLE_OUT}"

"${APPTAINER_BIN}" exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
    "${CONTAINER}" \
    metaquast.py \
        "${MEGAHIT_FA}" \
        "${SPADES_FA}" \
        --labels "MEGAHIT,metaSPAdes" \
        --min-contig "${MIN_CONTIG}" \
        --threads "${THREADS}" \
        --output-dir "${SAMPLE_OUT}" \
        --no-icarus \
        --max-ref-number 0   

echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Terminé pour ${SAMPLE}."
echo "    Rapport : ${SAMPLE_OUT}/report.html"

# =====================================================================
# ÉTAPE 2 : RAPPORT GLOBAL (tous échantillons combinés)
# Lancé uniquement par le task 0
# =====================================================================
if [ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Génération du rapport global multi-échantillons..."

    MEGAHIT_ALL=()
    SPADES_ALL=()
    LABELS_ALL=()

    for S in "${SAMPLES[@]}"; do
        MEGAHIT_ALL+=("${MEGAHIT_DIR}/${S}_megahit.contigs.fa")
        SPADES_ALL+=("${SPADES_DIR}/${S}_metaspades_contigs.fa")
        LABELS_ALL+=("${S}_MEGAHIT" "${S}_metaSPAdes")
    done

    GLOBAL_OUT="${OUT_DIR}/global_comparison"
    mkdir -p "${GLOBAL_OUT}"

    LABEL_STR=$(IFS=','; echo "${LABELS_ALL[*]}")
    ALL_FA=("${MEGAHIT_ALL[@]}" "${SPADES_ALL[@]}")

    "${APPTAINER_BIN}" exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        "${CONTAINER}" \
        metaquast.py \
            "${ALL_FA[@]}" \
            --labels "${LABEL_STR}" \
            --min-contig "${MIN_CONTIG}" \
            --threads "${THREADS}" \
            --output-dir "${GLOBAL_OUT}" \
            --no-icarus \

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Rapport global : ${GLOBAL_OUT}/report.html"
fi

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] MetaQUAST terminé pour ${SAMPLE}."
echo "====================================================================="
