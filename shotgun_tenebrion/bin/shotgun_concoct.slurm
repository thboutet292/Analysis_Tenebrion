#!/bin/bash
#SBATCH --job-name=PROD_CONCOCT
#SBATCH --output=log/PROD_CONCOCT_%j.out
#SBATCH --error=log/PROD_CONCOCT_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=normal, lmge
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   run_concoct_prod.slurm
#   GOAL    :   PRODUCTION : Binning complet du catalogue avec CONCOCT.
#               Basé sur l'abondance (BAMs) et la composition k-mer.
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

IMG_CONCOCT="${PROJECT_DIR}/containers/concoct_1.1.0.sif"

DIR_BINS="${PROJECT_DIR}/results/concoct_bins"
mkdir -p "${DIR_BINS}" log/

export APPTAINER_CACHEDIR="/storage/scratch/${USER}/apptainer_cache_concoct"
export APPTAINER_TMPDIR="/storage/scratch/${USER}/apptainer_tmp_concoct"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
BIND_OPT="--bind ${PROJECT_DIR}:${PROJECT_DIR},${APPTAINER_TMPDIR}:/tmp"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE PRODUCTION CONCOCT"
echo "Job ID : ${SLURM_JOB_ID}"
echo "====================================================================="

# Vérification des index BAM (CRITIQUE POUR CONCOCT)
NB_BAMS=$(find "${BAM_DIR}" -maxdepth 1 -name "*.bam" 2>/dev/null | wc -l)
NB_BAIS=$(find "${BAM_DIR}" -maxdepth 1 -name "*.bai" 2>/dev/null | wc -l)

echo " -> Détection de ${NB_BAMS} fichiers BAM et ${NB_BAIS} fichiers d'index (.bai)."
if [ "${NB_BAMS}" -gt "${NB_BAIS}" ]; then
    echo "Veuillez indexer vos BAMs (samtools index) avant de lancer ce script."
    exit 1
fi

BAM_LIST=$(ls -1 "${BAM_DIR}/"*.bam | tr '\n' ' ')

# =====================================================================
# ÉTAPE 1 : PRÉPARATION DU CATALOGUE
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Préparation du catalogue complet..."

CATALOG_1500="${DIR_BINS}/catalog_min1500.fasta"

# Filtrage à >= 1500 bp (si le fichier n'existe pas déjà pour éviter de le refaire)
if [ ! -f "${CATALOG_1500}" ]; then
    echo " -> Filtrage du catalogue (>= 1500 bp)..."
    awk 'BEGIN{RS=">";FS="\n"} NR>1{seq=""; for(i=2;i<=NF;i++) seq=seq$i; if(length(seq)>=1500) printf ">%s\n%s\n", $1, seq}' "${CATALOG_FASTA}" > "${CATALOG_1500}"
else
    echo " -> Catalogue filtré déjà présent."
fi

# =====================================================================
# ÉTAPE 2 : EXÉCUTION DE CONCOCT
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Lancement du workflow CONCOCT..."

# 2.A Découpage des contigs longs
echo " -> [1/5] Découpage des contigs (cut_up_fasta)..."
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" cut_up_fasta.py "${CATALOG_1500}" -c 10000 -o 0 --merge_last -b "${DIR_BINS}/contigs_10K.bed" > "${DIR_BINS}/contigs_10K.fa"

# 2.B Création de la matrice d'abondance
# (C'est l'étape qui peut prendre le plus de temps)
echo " -> [2/5] Calcul de la table de couverture globale (patientez...)..."
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" concoct_coverage_table.py "${DIR_BINS}/contigs_10K.bed" ${BAM_LIST} > "${DIR_BINS}/coverage_table.tsv"

# 2.C L'algorithme de Clustering
echo " -> [3/5] Apprentissage et Clustering (Modèles Gaussiens)..."
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" concoct \
    --composition_file "${DIR_BINS}/contigs_10K.fa" \
    --coverage_file "${DIR_BINS}/coverage_table.tsv" \
    -b "${DIR_BINS}/" \
    -t "${THREADS}"

# 2.D Re-fusion des fragments virtuels
echo " -> [4/5] Re-fusion des fragments..."
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" merge_cutup_clustering.py "${DIR_BINS}/clustering_gt1000.csv" > "${DIR_BINS}/clustering_merged.csv"

# 2.E Génération finale des fichiers FASTA (les Bins/MAGs)
echo " -> [5/5] Extraction des Bins en fichiers FASTA..."
DIR_BINS_FASTA="${DIR_BINS}/bins_fasta"
mkdir -p "${DIR_BINS_FASTA}"
apptainer exec ${BIND_OPT} "${IMG_CONCOCT}" extract_fasta_bins.py "${CATALOG_1500}" "${DIR_BINS}/clustering_merged.csv" --output_path "${DIR_BINS_FASTA}"

# Renommage des .fasta en .fa pour être homogène avec MetaBAT2 et CheckM2
for f in "${DIR_BINS_FASTA}/"*.fasta; do
    [ -f "$f" ] && mv -- "$f" "${f%.fasta}.fa"
done

# Statistiques de sortie
NB_BINS_CONCOCT=$(ls -1 "${DIR_BINS_FASTA}/"*.fa 2>/dev/null | wc -l || echo 0)

# =====================================================================
# NETTOYAGE
# =====================================================================
rm -rf "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"
echo "---------------------------------------------------------------------"
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FIN DE LA PRODUCTION CONCOCT."
echo "Nombre de Bins (MAGs) générés : ${NB_BINS_CONCOCT}"
echo "Les génomes sont disponibles dans : ${DIR_BINS_FASTA}"
echo "====================================================================="