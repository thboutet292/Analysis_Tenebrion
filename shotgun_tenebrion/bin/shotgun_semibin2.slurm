#!/bin/bash
#SBATCH --job-name=PROD_SEMIBIN2
#SBATCH --output=log/PROD_SEMIBIN2_%j.out
#SBATCH --error=log/PROD_SEMIBIN2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --partition=normal,lmge
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   run_semibin2.slurm
#   GOAL    :   PRODUCTION : Binning complet du catalogue avec SemiBin2.
#               Utilisation de l'apprentissage contrastif auto-supervisé (Deep Learning)
#               via le mode single_easy_bin.
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
IMG_SEMIBIN="${PROJECT_DIR}/containers/semibin_2.1.0.sif"

DIR_BINS="${PROJECT_DIR}/results/semibin2_bins"
mkdir -p "${DIR_BINS}" log/

export APPTAINER_CACHEDIR="/storage/scratch/${USER}/apptainer_cache_semibin2"
export APPTAINER_TMPDIR="/storage/scratch/${USER}/apptainer_tmp_semibin2"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
BIND_OPT="--bind ${PROJECT_DIR}:${PROJECT_DIR},${APPTAINER_TMPDIR}:/tmp"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE PRODUCTION SEMIBIN2"
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
    echo "SemiBin2 nécessite un index pour chaque fichier BAM. Veuillez les générer (samtools index)."
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
# ÉTAPE 2 : EXÉCUTION DU RÉSEAU DE NEURONES (SEMIBIN2)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Lancement du binning SemiBin2..."
echo " -> Cette étape comprend la génération des caractéristiques (k-mers, profondeur) et l'entraînement du modèle."

apptainer exec ${BIND_OPT} "${IMG_SEMIBIN}" SemiBin2 single_easy_bin \
    --input-fasta "${CATALOG_1500}" \
    --input-bam ${BAM_LIST} \
    --output "${DIR_BINS}" \
    --threads "${THREADS}" \
    --compression none

# =====================================================================
# ÉTAPE 3 : FORMATAGE DES RÉSULTATS
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] Validation et formatage des sorties..."

DIR_BINS_FASTA="${DIR_BINS}/bins_fasta"
mkdir -p "${DIR_BINS_FASTA}"

# SemiBin2 place ses résultats dans le sous-dossier output_bins. On les déplace pour standardisation.
if [ -d "${DIR_BINS}/output_bins" ]; then
    cp "${DIR_BINS}/output_bins/"*.fa "${DIR_BINS_FASTA}/" 2>/dev/null || true
fi

NB_BINS_SEMIBIN=$(ls -1 "${DIR_BINS_FASTA}/"*.fa 2>/dev/null | wc -l || echo 0)

# =====================================================================
# NETTOYAGE
# =====================================================================
rm -rf "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"
echo "---------------------------------------------------------------------"
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FIN DE LA PRODUCTION SEMIBIN2."
echo "Nombre de Bins (MAGs) générés : ${NB_BINS_SEMIBIN}"
echo "Les génomes sont disponibles dans : ${DIR_BINS_FASTA}"
echo "====================================================================="