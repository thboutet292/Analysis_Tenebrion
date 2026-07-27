#!/bin/bash
#SBATCH --job-name=PHASE5_GTDBTK
#SBATCH --output=log/PHASE5_GTDBTK_%j.out
#SBATCH --error=log/PHASE5_GTDBTK_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --tmp=100G
#SBATCH --partition=normal,lmge
#SBATCH --time=24:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   shotgun_gtdbtk.slurm
#
#   GOAL    :   PHASE 5 - FILTRAGE & ASSIGNATION TAXONOMIQUE DES MAGs, pour les
#               TROIS binners (MetaBAT2, SemiBin2, CONCOCT), en une seule soumission.
#
#               Pour chaque binneur :
#                 1. Filtrage des bins selon le standard MIMAG Medium-Quality :
#                    Completeness >= 100% ET Contamination < 10% (CheckM2).
#                 2. Assignation taxonomique avec GTDB-Tk (classify_wf).
#
#               La base GTDB-Tk (~85 Go) n'est rapatriée qu'UNE SEULE FOIS depuis
#               le S3 et réutilisée pour les trois runs (au lieu d'un téléchargement
#               par binneur), pour économiser temps et bande passante.
#
#               Les trois runs GTDB-Tk sont lancés SÉQUENTIELLEMENT : l'outil
#               nécessite ~150 Go de RAM et monopolise les CPUs, donc les lancer
#               en parallèle dépasserait l'enveloppe mémoire allouée (250 Go).
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE PHASE 5 : FILTRAGE & GTDB-Tk (3 binners)"
echo "Job ID : ${SLURM_JOB_ID}"
echo "====================================================================="

module purge
module load apptainer 2>/dev/null || true
module load rclone/1.55.1 2>/dev/null || true

# =====================================================================
# BLOCK 0 : CONFIGURATION GLOBALE
# =====================================================================

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

IMG_GTDBTK="${PROJECT_DIR}/containers/gtdbtk_v2.3.2.sif"
S3_DB="s3_uca:lmge-tenebrion/ref/gtdbtk_r214"

# --- Seuil de qualité (standard MIMAG Medium-Quality) ---
COMPLETENESS_MIN=100
CONTAMINATION_MAX=10

# --- MetaBAT2 ---
DIR_METABAT2_BINS="${PROJECT_DIR}/results/metabat2_bins"
CHECKM2_METABAT2="${PROJECT_DIR}/results/checkm2_qc/metabat2/quality_report.tsv"
DIR_METABAT2_FILTERED="${PROJECT_DIR}/results/filtered_bins_metabat2"
DIR_METABAT2_OUT="${PROJECT_DIR}/results/gtdbtk_taxo_metabat2"

# --- SemiBin2 ---
DIR_SEMIBIN2_BINS="${PROJECT_DIR}/results/semibin2_bins/bins_fasta"
CHECKM2_SEMIBIN2="${PROJECT_DIR}/results/checkm2_qc/semibin2/quality_report.tsv"
DIR_SEMIBIN2_FILTERED="${PROJECT_DIR}/results/filtered_bins_semibin2"
DIR_SEMIBIN2_OUT="${PROJECT_DIR}/results/gtdbtk_taxo_semibin2"

# --- CONCOCT ---
DIR_CONCOCT_BINS="${PROJECT_DIR}/results/concoct_bins/bins_fasta"
CHECKM2_CONCOCT="${PROJECT_DIR}/results/checkm2_qc/concoct/quality_report.tsv"
DIR_CONCOCT_FILTERED="${PROJECT_DIR}/results/filtered_bins_concoct"
DIR_CONCOCT_OUT="${PROJECT_DIR}/results/gtdbtk_taxo_concoct"

mkdir -p "${DIR_METABAT2_FILTERED}" "${DIR_SEMIBIN2_FILTERED}" "${DIR_CONCOCT_FILTERED}" log/

# =====================================================================
# ÉTAPE 1 : RAPATRIEMENT DE LA BASE GTDB DEPUIS LE S3 (UNE SEULE FOIS)
# =====================================================================
LOCAL_DB_DIR="/storage/scratch/${USER}/gtdbtk_db_${SLURM_JOB_ID}"
mkdir -p "${LOCAL_DB_DIR}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression de la base GTDB sur le scratch..."
    rm -rf "${LOCAL_DB_DIR}" 2>/dev/null || true
}
trap cleanup EXIT

echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Rapatriement de la base GTDB depuis le S3..."
rclone copy "${S3_DB}" "${LOCAL_DB_DIR}" --transfers 16

# =====================================================================
# FONCTION : run_gtdbtk_for_binner
#
# Arguments :
#   $1 = nom du binneur (pour les logs)         ex: "MetaBAT2"
#   $2 = rapport CheckM2 (quality_report.tsv)
#   $3 = répertoire contenant les bins bruts (.fa)
#   $4 = répertoire de destination des bins filtrés
#   $5 = répertoire de sortie GTDB-Tk
# =====================================================================
run_gtdbtk_for_binner() {

    local BINNER_NAME="$1"
    local CHECKM2_REPORT="$2"
    local DIR_ALL_BINS="$3"
    local DIR_FILTERED="$4"
    local DIR_OUT="$5"

    echo ""
    echo "###################################################################"
    echo "  BINNER : ${BINNER_NAME}"
    echo "###################################################################"

    # ------------------------------------------------------------------
    # ETAPE A : Vérification des fichiers d'entrée
    # ------------------------------------------------------------------
    if [ ! -f "${CHECKM2_REPORT}" ]; then
        echo "[AVERT.] Rapport CheckM2 introuvable pour ${BINNER_NAME} : ${CHECKM2_REPORT}. Binner ignoré." >&2
        return
    fi

    if [ ! -d "${DIR_ALL_BINS}" ]; then
        echo "[AVERT.] Répertoire de bins introuvable pour ${BINNER_NAME} : ${DIR_ALL_BINS}. Binner ignoré." >&2
        return
    fi

    # ------------------------------------------------------------------
    # ETAPE B : Filtrage des MAGs (Completeness >= 100% ET Contamination < 10%)
    # ------------------------------------------------------------------
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${BINNER_NAME}] Filtrage (Completeness >= ${COMPLETENESS_MIN}% ET Contamination < ${CONTAMINATION_MAX}%)..."

    # On repart d'un répertoire vide pour éviter de mélanger des bins
    # d'un run précédent avec ceux du run en cours.
    rm -f "${DIR_FILTERED}/"*.fa

    # La colonne 1 du rapport CheckM2 contient le nom du bin (sans extension .fa).
    # La colonne 2 contient la Completeness, la colonne 3 la Contamination.
    local NB_PASSED=0
    while IFS=$'\t' read -r bin_name completeness contamination; do
        local passes
        passes=$(awk -v comp="${completeness}" -v cont="${contamination}" \
                     -v min="${COMPLETENESS_MIN}" -v max="${CONTAMINATION_MAX}" \
                     'BEGIN { print (comp >= min && cont < max) ? "yes" : "no" }')
        if [ "${passes}" = "yes" ]; then
            local BIN_FILE="${DIR_ALL_BINS}/${bin_name}.fa"
            if [ -f "${BIN_FILE}" ]; then
                cp "${BIN_FILE}" "${DIR_FILTERED}/"
                NB_PASSED=$((NB_PASSED + 1))
            else
                echo "[AVERT.] Bin attendu mais introuvable sur disque : ${BIN_FILE}" >&2
            fi
        fi
    done < <(awk -F'\t' 'NR>1 { print $1 "\t" $2 "\t" $3 }' "${CHECKM2_REPORT}")

    echo " -> ${NB_PASSED} MAGs ont passé le filtre pour ${BINNER_NAME}."

    if [ "${NB_PASSED}" -eq 0 ]; then
        echo "[AVERT.] Aucun MAG ne respecte le seuil pour ${BINNER_NAME}. Binner ignoré." >&2
        return
    fi

    # ------------------------------------------------------------------
    # ETAPE C : Lancement de GTDB-Tk (base déjà rapatriée à l'étape 1)
    # ------------------------------------------------------------------
    echo "---------------------------------------------------------------------"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${BINNER_NAME}] Lancement de GTDB-Tk classify_wf..."

    # On repart d'un répertoire de sortie vide pour éviter les conflits.
    rm -rf "${DIR_OUT}"
    mkdir -p "${DIR_OUT}"

    apptainer exec \
        --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
        --bind "${LOCAL_DB_DIR}:${LOCAL_DB_DIR}" \
        --env GTDBTK_DATA_PATH="${LOCAL_DB_DIR}" \
        "${IMG_GTDBTK}" \
        gtdbtk classify_wf \
            --genome_dir "${DIR_FILTERED}" \
            --out_dir    "${DIR_OUT}" \
            --extension  fa \
            --skip_ani_screen \
            --cpus       "${THREADS}"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${BINNER_NAME}] Terminé. Résultats dans : ${DIR_OUT}"
}

# =====================================================================
# LANCEMENT SÉQUENTIEL DES TROIS BINNERS
# =====================================================================

run_gtdbtk_for_binner \
    "MetaBAT2" \
    "${CHECKM2_METABAT2}" \
    "${DIR_METABAT2_BINS}" \
    "${DIR_METABAT2_FILTERED}" \
    "${DIR_METABAT2_OUT}"

run_gtdbtk_for_binner \
    "SemiBin2" \
    "${CHECKM2_SEMIBIN2}" \
    "${DIR_SEMIBIN2_BINS}" \
    "${DIR_SEMIBIN2_FILTERED}" \
    "${DIR_SEMIBIN2_OUT}"

run_gtdbtk_for_binner \
    "CONCOCT" \
    "${CHECKM2_CONCOCT}" \
    "${DIR_CONCOCT_BINS}" \
    "${DIR_CONCOCT_FILTERED}" \
    "${DIR_CONCOCT_OUT}"

# =====================================================================
# RÉSUMÉ FINAL
# =====================================================================
echo ""
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FIN DE L'ASSIGNATION TAXONOMIQUE (3 binners)."
echo ""
echo "  MetaBAT2  -> bins filtrés : ${DIR_METABAT2_FILTERED}/"
echo "              résultats GTDB : ${DIR_METABAT2_OUT}/"
echo ""
echo "  SemiBin2  -> bins filtrés : ${DIR_SEMIBIN2_FILTERED}/"
echo "              résultats GTDB : ${DIR_SEMIBIN2_OUT}/"
echo ""
echo "  CONCOCT   -> bins filtrés : ${DIR_CONCOCT_FILTERED}/"
echo "              résultats GTDB : ${DIR_CONCOCT_OUT}/"
echo "====================================================================="