#!/bin/bash
#SBATCH --job-name=PROD_MAXBIN2
#SBATCH --output=log/PROD_MAXBIN2_%j.out
#SBATCH --error=log/PROD_MAXBIN2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --partition=normal,lmge
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   shotgun_maxbin2.slurm
#   GOAL    :   PRODUCTION : Binning complet du catalogue avec MaxBin2.
#               Utilisation d'un algorithme Expectation-Maximization basé sur la
#               composition en tétranucléotides et l'abondance (profondeur des BAM).
#######################################################################################

set -euo pipefail

module purge
module load apptainer 2>/dev/null || true

# =====================================================================
# CONFIGURATION DES CHEMINS
# =====================================================================
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

CATALOG_FASTA="${PROJECT_DIR}/results/global_catalog/non_redundant_catalog.fasta"
BAM_DIR="${PROJECT_DIR}/results/coverm/bams"
IMG_MAXBIN="${PROJECT_DIR}/containers/maxbin2_2.2.7.sif"
IMG_METABAT="${PROJECT_DIR}/containers/metabat2_v2.15.sif"

DIR_BINS="${PROJECT_DIR}/results/maxbin2_bins"
mkdir -p "${DIR_BINS}" log/

export APPTAINER_CACHEDIR="/storage/scratch/${USER}/apptainer_cache_maxbin2"
export APPTAINER_TMPDIR="/storage/scratch/${USER}/apptainer_tmp_maxbin2"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
BIND_OPT="--bind ${PROJECT_DIR}:${PROJECT_DIR},${APPTAINER_TMPDIR}:/tmp"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE PRODUCTION MAXBIN2"
echo "Job ID : ${SLURM_JOB_ID}"
echo "====================================================================="

# ---------------------------------------------------------------------
# CONTRÔLE D'INTÉGRITÉ DES FICHIERS D'ENTRÉE
# ---------------------------------------------------------------------
echo "[ETAPE 0] Vérification des prérequis d'alignement..."

NB_BAMS=$(ls -1 "${BAM_DIR}/"*.bam 2>/dev/null | wc -l || echo 0)
NB_BAIS=$(ls -1 "${BAM_DIR}/"*.bai 2>/dev/null | wc -l || echo 0)

echo " -> Détection de ${NB_BAMS} fichiers BAM et ${NB_BAIS} fichiers d'index (.bai)."
if [ "${NB_BAMS}" -gt "${NB_BAIS}" ]; then
    echo "[ERREUR FATALE] Il manque des fichiers d'index (.bai) dans ${BAM_DIR}."
    echo "La génération de la matrice de profondeur nécessite un index pour chaque fichier BAM. Veuillez les générer (samtools index)."
    exit 1
elif [ "${NB_BAMS}" -eq 0 ]; then
    echo "[ERREUR FATALE] Aucun fichier BAM détecté dans ${BAM_DIR}."
    exit 1
fi

BAM_LIST=$(ls -1 "${BAM_DIR}/"*.bam | tr '\n' ' ')

# =====================================================================
# ÉTAPE 1 : PRÉPARATION DU CATALOGUE
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Préparation du catalogue..."

CATALOG_1500="${DIR_BINS}/catalog_min1500.fasta"

if [ ! -f "${CATALOG_1500}" ]; then
    echo " -> Filtrage du catalogue pour retenir les contigs >= 1500 bp..."
    awk 'BEGIN{RS=">";FS="\n"} NR>1{seq=""; for(i=2;i<=NF;i++) seq=seq$i; if(length(seq)>=1500) printf ">%s\n%s\n", $1, seq}' "${CATALOG_FASTA}" > "${CATALOG_1500}"
else
    echo " -> Catalogue filtré existant détecté : ${CATALOG_1500}"
fi

# =====================================================================
# ÉTAPE 2 : GÉNÉRATION DE LA MATRICE DE PROFONDEUR (jgi_summarize_bam_contig_depths)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Calcul de la profondeur par BAM (outil MetaBAT2)..."
echo " -> MaxBin2 n'a pas d'outil de calcul de profondeur intégré, on réutilise jgi_summarize_bam_contig_depths."

DEPTH_MATRIX="${DIR_BINS}/depth_matrix.txt"

if [ ! -f "${DEPTH_MATRIX}" ]; then
    apptainer exec ${BIND_OPT} "${IMG_METABAT}" jgi_summarize_bam_contig_depths \
        --outputDepth "${DEPTH_MATRIX}" \
        ${BAM_LIST}
else
    echo " -> Matrice de profondeur existante détectée : ${DEPTH_MATRIX}"
fi

# =====================================================================
# ÉTAPE 3 : CONVERSION EN FICHIERS D'ABONDANCE (FORMAT MAXBIN2)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] Génération des fichiers d'abondance par échantillon..."
echo " -> MaxBin2 attend un fichier à 2 colonnes (contig, profondeur) par échantillon, listés dans -abund_list."

DIR_ABUND="${DIR_BINS}/abundance_files"
mkdir -p "${DIR_ABUND}"
ABUND_LIST="${DIR_BINS}/abund_list.txt"

NB_SAMPLES=$(( ($(head -1 "${DEPTH_MATRIX}" | awk '{print NF}') - 3) / 2 ))
echo " -> ${NB_SAMPLES} échantillon(s) détecté(s) dans la matrice de profondeur."

> "${ABUND_LIST}"
for i in $(seq 1 "${NB_SAMPLES}"); do
    COL=$(( 2 * i + 2 ))
    ABUND_FILE="${DIR_ABUND}/sample_${i}.abund"
    awk -v col="${COL}" 'BEGIN{FS=OFS="\t"} NR>1{print $1, $col}' "${DEPTH_MATRIX}" > "${ABUND_FILE}"
    echo "${ABUND_FILE}" >> "${ABUND_LIST}"
done

# =====================================================================
# ÉTAPE 4 : EXÉCUTION DU BINNING (MAXBIN2)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 4] Lancement du binning MaxBin2..."

MAXBIN_PREFIX="${DIR_BINS}/maxbin2_out"

apptainer exec ${BIND_OPT} "${IMG_MAXBIN}" run_MaxBin.pl \
    -contig "${CATALOG_1500}" \
    -abund_list "${ABUND_LIST}" \
    -out "${MAXBIN_PREFIX}" \
    -thread "${THREADS}"

# =====================================================================
# ÉTAPE 5 : FORMATAGE DES RÉSULTATS
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 5] Validation et formatage des sorties..."

DIR_BINS_FASTA="${DIR_BINS}/bins_fasta"
mkdir -p "${DIR_BINS_FASTA}"

# MaxBin2 nomme ses bins <prefix>.001.fasta, <prefix>.002.fasta, etc. On les renomme
# en .fa pour rester cohérent avec les sorties des autres binners (MetaBAT2, CONCOCT, SemiBin2).
for BIN_FASTA in "${MAXBIN_PREFIX}".*.fasta; do
    [ -e "${BIN_FASTA}" ] || continue
    BASENAME=$(basename "${BIN_FASTA}" .fasta)
    cp "${BIN_FASTA}" "${DIR_BINS_FASTA}/${BASENAME}.fa"
done

NB_BINS_MAXBIN=$(ls -1 "${DIR_BINS_FASTA}/"*.fa 2>/dev/null | wc -l || echo 0)

# =====================================================================
# NETTOYAGE
# =====================================================================
rm -rf "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"
echo "---------------------------------------------------------------------"
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FIN DE LA PRODUCTION MAXBIN2."
echo "Nombre de Bins (MAGs) générés : ${NB_BINS_MAXBIN}"
echo "Les génomes sont disponibles dans : ${DIR_BINS_FASTA}"
echo "====================================================================="