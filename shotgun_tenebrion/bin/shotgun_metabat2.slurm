#!/bin/bash
#SBATCH --job-name=METABAT2
#SBATCH --output=log/PHASE3_METABAT2_%j.out
#SBATCH --error=log/PHASE3_METABAT2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal,lmge
#SBATCH --time=12:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   shotgun_metabat2.slurm
#
#   GOAL    :   PHASE 3 - BINNING (Reconstruction des Génomes/MAGs).
#               1. Création de la matrice de profondeur de couverture globale (JGI).
#               2. Regroupement des contigs en MAGs avec MetaBAT2.
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE PHASE 3 : BINNING (MetaBAT2)"
echo "Job ID : ${SLURM_JOB_ID}"
echo "====================================================================="

module purge
module load apptainer 2>/dev/null || true

# =====================================================================
# CONFIGURATION DES CHEMINS
# =====================================================================
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

DIR_CATALOG="${PROJECT_DIR}/results/global_catalog"
CATALOG_FASTA="${DIR_CATALOG}/non_redundant_catalog.fasta"

DIR_COVERM="${PROJECT_DIR}/results/coverm"
BAM_DIR="${DIR_COVERM}/bams"

DIR_BINS="${PROJECT_DIR}/results/metabat2_bins"
IMG_METABAT="${PROJECT_DIR}/containers/metabat2_v2.15.sif"

mkdir -p "${DIR_BINS}" log/

THREADS="${SLURM_CPUS_PER_TASK:-16}"
DEPTH_FILE="${DIR_BINS}/global_depth.txt"

# --- Vérification des fichiers d'entrée ---
if [ ! -f "${CATALOG_FASTA}" ]; then
    echo "[ERREUR] Le catalogue global est introuvable à : ${CATALOG_FASTA}" >&2
    exit 1
fi

NB_BAMS=$(ls -1 "${BAM_DIR}/"*.bam 2>/dev/null | wc -l)
if [ "${NB_BAMS}" -eq 0 ]; then
    echo "[ERREUR] Aucun fichier BAM trouvé dans ${BAM_DIR}." >&2
    exit 1
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${NB_BAMS} fichiers BAM détectés. Analyse de co-abondance prête."

# =====================================================================
# ÉTAPE 1 : CALCUL DE LA PROFONDEUR DE COUVERTURE (jgi_summarize)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Résumé des profondeurs de contigs..."

# On trie la liste des BAMs pour garantir que l'ordre des colonnes est constant
BAM_LIST=$(ls -1 "${BAM_DIR}/"*.bam | sort)

apptainer exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    "${IMG_METABAT}" \
    jgi_summarize_bam_contig_depths \
        --outputDepth "${DEPTH_FILE}" \
        ${BAM_LIST}

echo "[$(date '+%Y-%m-%d %H:%M:%S')] -> Matrice de profondeur générée : ${DEPTH_FILE}"

# =====================================================================
# ÉTAPE 2 : BINNING (CRÉATION DES MAGs)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Lancement de MetaBAT2..."

# 2.A : Création d'une copie du catalogue filtrée à >= 1500 pb
# (Obligatoire car MetaBAT2 refuse de lire un fichier contenant des petits contigs)
CATALOG_1500="${DIR_BINS}/catalog_min1500.fasta"

echo " -> Filtrage du catalogue pour MetaBAT2 (>= 1500 bp)..."
awk 'BEGIN{RS=">";FS="\n"} NR>1{seq=""; for(i=2;i<=NF;i++) seq=seq$i; if(length(seq)>=1500) printf ">%s\n%s\n", $1, seq}' "${CATALOG_FASTA}" > "${CATALOG_1500}"

# 2.B : Lancement de l'algorithme
echo " -> Exécution de l'algorithme de clustering..."
apptainer exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    "${IMG_METABAT}" \
    metabat2 \
        -i "${CATALOG_1500}" \
        -a "${DEPTH_FILE}" \
        -o "${DIR_BINS}/bin" \
        -m 1500 \
        -t "${THREADS}" \
        -v

# Nettoyage du fichier temporaire
rm "${CATALOG_1500}"

# Statistiques finales
NB_BINS=$(ls -1 "${DIR_BINS}/"bin*.fa 2>/dev/null | wc -l || echo 0)

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FIN DU BINNING."
echo "MetaBAT2 a généré ${NB_BINS} bins (MAGs) potentiels dans ${DIR_BINS}/"
echo "====================================================================="