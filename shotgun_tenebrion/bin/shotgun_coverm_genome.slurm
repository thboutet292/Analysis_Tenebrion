#!/bin/bash
#SBATCH --job-name=PHASE6_ABUND
#SBATCH --output=log/PHASE6_ABUND_%j.out
#SBATCH --error=log/PHASE6_ABUND_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=normal,lmge
#SBATCH --time=04:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   phase6_coverm_genome.slurm
#
#   GOAL    :   PHASE 6 - CALCUL D'ABONDANCE DES MAGs.
#               - Calcul des abondances avec CoverM.
#               - Remplacement automatique des "bin.X" par la taxonomie GTDB-Tk.
#               - Nettoyage des noms de colonnes d'échantillons.
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE PHASE 6 : ABONDANCE DES MAGS"
echo "Job ID : ${SLURM_JOB_ID}"
echo "====================================================================="

module purge
module load apptainer 2>/dev/null || true

# =====================================================================
# CONFIGURATION DES CHEMINS
# =====================================================================
PROJECT_DIR="${SLURM_SUBMIT_DIR}"

DIR_BINS="${PROJECT_DIR}/results/filtered_bins"
DIR_BAMS="${PROJECT_DIR}/results/coverm/bams"
DIR_OUT="${PROJECT_DIR}/results/abundance_mags"
GTDB_SUMMARY="${PROJECT_DIR}/results/gtdbtk_taxo/classify/gtdbtk.bac120.summary.tsv"

mkdir -p "${DIR_OUT}" log/

IMG_COVERM="${PROJECT_DIR}/containers/coverm_v0.7.0.sif"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

JOB_SCRATCH="/storage/scratch/${USER}/coverm_genome_${SLURM_JOB_ID}"
mkdir -p "${JOB_SCRATCH}"
export TMPDIR="${JOB_SCRATCH}"

# =====================================================================
# ÉTAPE 1 : LANCEMENT DE COVERM
# =====================================================================
echo "---------------------------------------------------------------------"
echo "-> Calcul des abondances avec CoverM (mode genome)..."

RAW_TSV="${DIR_OUT}/MAGs_raw_summary.tsv"
FINAL_TSV="${DIR_OUT}/MAGs_abundance_taxo.tsv"

apptainer exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    --bind "${JOB_SCRATCH}:${JOB_SCRATCH}" \
    --bind "${JOB_SCRATCH}:/tmp" \
    "${IMG_COVERM}" \
    coverm genome \
        --genome-fasta-directory "${DIR_BINS}" \
        -x fa \
        --bam-files "${DIR_BAMS}/"*.bam \
        --methods relative_abundance count \
        --min-covered-fraction 0 \
        --output-file "${RAW_TSV}" \
        --threads "${THREADS}"

# =====================================================================
# ÉTAPE 2 : FUSION AVEC LA TAXONOMIE GTDB-TK ET NETTOYAGE DES COLONNES
# =====================================================================
echo "---------------------------------------------------------------------"
echo "-> Ajout de la taxonomie et nettoyage de l'en-tête..."

if [ -f "${GTDB_SUMMARY}" ]; then
    awk -F'\t' -v OFS='\t' '
        # 1er passage : on lit le fichier GTDB-Tk et on mémorise la taxonomie
        NR==FNR {
            if (NR > 1) { tax[$1] = $2 }
            next
        }
        # 2ème passage : on lit le tableau CoverM
        {
            if (FNR == 1) {
                # Nettoyage de la toute première cellule
                $1 = "Taxonomy"
                
                # Nettoyage du nom des échantillons dans toutes les colonnes
                gsub(/non_redundant_catalog\.fasta\./, "", $0)
                gsub(/_clean_R1\.fastq\.gz/, "", $0)
                gsub(/\.bam/, "", $0) # Au cas où CoverM ajouterait .bam
                
                print $0
            } else if ($1 == "unmapped") {
                print $0
            } else {
                # Remplacement du nom du bin par sa taxonomie
                if ($1 in tax) { $1 = tax[$1] }
                print $0
            }
        }
    ' "${GTDB_SUMMARY}" "${RAW_TSV}" > "${FINAL_TSV}"
    
    echo "Taxonomie ajoutée et colonnes nettoyées avec succès !"
    rm "${RAW_TSV}"
else
    echo "[ATTENTION] Fichier GTDB introuvable (${GTDB_SUMMARY})."
    mv "${RAW_TSV}" "${FINAL_TSV}"
fi

# =====================================================================
# NETTOYAGE
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Nettoyage du scratch..."
rm -rf "${JOB_SCRATCH}"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FIN DU CALCUL D'ABONDANCE."
echo "Tableau final généré : ${FINAL_TSV}"
echo "====================================================================="