#!/bin/bash
#SBATCH --job-name=TEST_BINNERS
#SBATCH --output=log/TEST_BINNERS_%j.out
#SBATCH --error=log/TEST_BINNERS_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=normal,lmge
#SBATCH --time=04:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   test_binners.slurm
#   GOAL    :   Validation technique (Dry-Run) de CONCOCT et SemiBin2
#               sur un jeu de données minimal (avec indexation BAM).
#######################################################################################

set -euo pipefail

module purge
module load apptainer 2>/dev/null || true

# =====================================================================
# CONFIGURATION DES CHEMINS
# =====================================================================
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

CATALOG_FASTA="${PROJECT_DIR}/results/global_catalog/non_redundant_catalog.fasta"
BAM_DIR_SRC="${PROJECT_DIR}/results/coverm/bams"

IMG_CONCOCT="${PROJECT_DIR}/containers/concoct_1.1.0.sif"
IMG_SEMIBIN="${PROJECT_DIR}/containers/semibin_2.1.0.sif"

DIR_TEST="${PROJECT_DIR}/results/test_binning"
DIR_TEST_BAMS="${DIR_TEST}/bams"
mkdir -p "${DIR_TEST_BAMS}" log/

export APPTAINER_CACHEDIR="/storage/scratch/${USER}/apptainer_cache_test_bin"
export APPTAINER_TMPDIR="/storage/scratch/${USER}/apptainer_tmp_test_bin"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE DU DRY-RUN (CONCOCT & SEMIBIN2)"
echo "====================================================================="

BIND_OPT="--bind ${PROJECT_DIR}:${PROJECT_DIR},${APPTAINER_TMPDIR}:/tmp"

# =====================================================================
# ÉTAPE 0 : PRÉPARATION DES DONNÉES RÉDUITES
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 0] Préparation des données de test..."

CATALOG_TEST="${DIR_TEST}/test_catalog.fasta"
echo " -> Extraction du catalogue réduit..."
head -n 20000 "${CATALOG_FASTA}" > "${CATALOG_TEST}"

CATALOG_1500="${DIR_TEST}/catalog_min1500.fasta"
awk 'BEGIN{RS=">";FS="\n"} NR>1{seq=""; for(i=2;i<=NF;i++) seq=seq$i; if(length(seq)>=1500) printf ">%s\n%s\n", $1, seq}' "${CATALOG_TEST}" > "${CATALOG_1500}"

BAM_ARRAY=("${BAM_DIR_SRC}"/*.bam)

if [ ${#BAM_ARRAY[@]} -lt 2 ]; then
    echo "[ERREUR] Moins de 2 fichiers BAM trouvés dans ${BAM_DIR_SRC}."
    exit 1
fi

echo " -> Copie des fichiers BAM vers le dossier de test..."
cp "${BAM_ARRAY[0]}" "${DIR_TEST_BAMS}/"
cp "${BAM_ARRAY[1]}" "${DIR_TEST_BAMS}/"

# CORRECTION : Création des index (.bai) indispensables pour la lecture
echo " -> Indexation des BAMs (samtools index)..."
BAM_TEST_1="${DIR_TEST_BAMS}/$(basename "${BAM_ARRAY[0]}")"
BAM_TEST_2="${DIR_TEST_BAMS}/$(basename "${BAM_ARRAY[1]}")"

apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" samtools index "${BAM_TEST_1}"
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" samtools index "${BAM_TEST_2}"

BAM_LIST=$(ls -1 "${DIR_TEST_BAMS}/"*.bam | tr '\n' ' ')

# =====================================================================
# ÉTAPE 1 : TEST CONCOCT
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [TEST CONCOCT] Démarrage..."
DIR_CONCOCT="${DIR_TEST}/concoct_output"
mkdir -p "${DIR_CONCOCT}"

echo " -> [1/5] cut_up_fasta..."
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" cut_up_fasta.py "${CATALOG_1500}" -c 10000 -o 0 --merge_last -b "${DIR_CONCOCT}/contigs_10K.bed" > "${DIR_CONCOCT}/contigs_10K.fa"

echo " -> [2/5] concoct_coverage_table..."
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" concoct_coverage_table.py "${DIR_CONCOCT}/contigs_10K.bed" ${BAM_LIST} > "${DIR_CONCOCT}/coverage_table.tsv"

echo " -> [3/5] concoct clustering..."
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" concoct --composition_file "${DIR_CONCOCT}/contigs_10K.fa" --coverage_file "${DIR_CONCOCT}/coverage_table.tsv" -b "${DIR_CONCOCT}/" -t 8

echo " -> [4/5] merge_cutup..."
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" merge_cutup_clustering.py "${DIR_CONCOCT}/clustering_gt1000.csv" > "${DIR_CONCOCT}/clustering_merged.csv"

echo " -> [5/5] extract_fasta_bins..."
mkdir -p "${DIR_CONCOCT}/bins"
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" extract_fasta_bins.py "${CATALOG_1500}" "${DIR_CONCOCT}/clustering_merged.csv" --output_path "${DIR_CONCOCT}/bins"

echo "[SUCCÈS] CONCOCT s'est exécuté complètement."

# =====================================================================
# ÉTAPE 2 : TEST SEMIBIN2 (Corrigé pour les catalogues globaux)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [TEST SEMIBIN2] Démarrage..."
DIR_SEMIBIN="${DIR_TEST}/semibin2_output"
mkdir -p "${DIR_SEMIBIN}"

# On utilise single_easy_bin, qui est beaucoup plus tolérant sur les noms de contigs.
# Il acceptera le catalogue global tel quel, tout en utilisant la liste des BAMs
# pour le profilage de co-abondance (multi-mapping).
echo " -> SemiBin2 single_easy_bin..."
apptainer exec ${BIND_OPT} "${IMG_SEMIBIN}" SemiBin2 single_easy_bin \
    --input-fasta "${CATALOG_1500}" \
    --input-bam ${BAM_LIST} \
    --output "${DIR_SEMIBIN}" \
    --threads 8 \
    --compression none

echo "[SUCCÈS] SemiBin2 s'est exécuté complètement."
# =====================================================================
# NETTOYAGE
# =====================================================================
rm -rf "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FIN DU DRY-RUN."
echo "====================================================================="