#!/bin/bash
#SBATCH --job-name=smart_build_db
#SBATCH --output=log/SMART_BUILD_DB_%j.out
#SBATCH --error=log/SMART_BUILD_DB_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --partition=lmge,normal
#SBATCH --time=24:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   GOAL : Téléchargement conditionnel des 4 bases nécessaires à HUMAnN 3.
#          Le script vérifie d'abord sur le S3. Si la base existe, il passe à la suite.
#######################################################################################

set -uo pipefail

echo "[$(date '+%Y-%m-%d %H:%M:%S')] --- DÉMARRAGE DU BUILD INTELLIGENT ---"

module purge
module load rclone/1.55.1

# --- CONFIGURATION ---
S3_DEST="s3_uca:lmge-tenebrion/ref/humann3_db"
IMG_PATH="/home/thboutet/projets/shotgun_tenebrion/containers/humann_v3.sif"
JOB_SCRATCH="/storage/scratch/${USER}/smart_build_${SLURM_JOB_ID}"

mkdir -p "${JOB_SCRATCH}"
APP_EXEC="apptainer exec -B ${JOB_SCRATCH}:${JOB_SCRATCH} ${IMG_PATH}"

# =====================================================================
# FONCTION DE VÉRIFICATION S3
# =====================================================================
# Cette fonction compte s'il y a des fichiers dans le dossier S3 spécifié.
# Renvoie "vrai" (0) s'il y a des fichiers, "faux" (1) sinon.
check_s3_exists() {
    local target_path="$1"
    local count=$(rclone ls "${target_path}" | head -n 5 | wc -l)
    if [ "$count" -gt 0 ]; then
        return 0
    else
        return 1
    fi
}

# =====================================================================
# ÉTAPE 1 : CHOCOPHLAN (Nucléotides)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [1/4] Vérification de ChocoPhlAn..."
if check_s3_exists "${S3_DEST}/chocophlan/"; then
    echo " -> Déjà présent sur le S3. Étape ignorée."
else
    echo " -> Absent. Téléchargement en cours..."
    ${APP_EXEC} humann_databases --download chocophlan full "${JOB_SCRATCH}/chocophlan" || echo "Config error ignored"
    echo " -> Envoi vers S3..."
    rclone copy "${JOB_SCRATCH}/chocophlan/" "${S3_DEST}/chocophlan/" -P --transfers 8
fi

# =====================================================================
# ÉTAPE 2 : UNIREF90 (Protéines)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [2/4] Vérification de UniRef90..."
if check_s3_exists "${S3_DEST}/uniref/"; then
    echo " -> Déjà présent sur le S3. Étape ignorée."
else
    echo " -> Absent. Téléchargement en cours..."
    ${APP_EXEC} humann_databases --download uniref uniref90_diamond "${JOB_SCRATCH}/uniref" || echo "Config error ignored"
    echo " -> Envoi vers S3..."
    rclone copy "${JOB_SCRATCH}/uniref/" "${S3_DEST}/uniref/" -P --transfers 8
fi

# =====================================================================
# ÉTAPE 3 : UTILITY MAPPING (Annotations)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [3/4] Vérification de Utility Mapping..."
if check_s3_exists "${S3_DEST}/utility_mapping/"; then
    echo " -> Déjà présent sur le S3. Étape ignorée."
else
    echo " -> Absent. Téléchargement en cours..."
    ${APP_EXEC} humann_databases --download utility_mapping full "${JOB_SCRATCH}/utility_mapping" || echo "Config error ignored"
    echo " -> Envoi vers S3..."
    rclone copy "${JOB_SCRATCH}/utility_mapping/" "${S3_DEST}/utility_mapping/" -P --transfers 8
fi

# =====================================================================
# ÉTAPE 4 : METAPHLAN DB (La nouvelle base manquante)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [4/4] Vérification de MetaPhlAn DB..."
if check_s3_exists "${S3_DEST}/metaphlan_db/"; then
    echo " -> Déjà présent sur le S3. Étape ignorée."
else
    echo " -> Absent. Téléchargement et indexation en cours..."
    mkdir -p "${JOB_SCRATCH}/metaphlan_db"
    ${APP_EXEC} metaphlan --install --bowtie2db "${JOB_SCRATCH}/metaphlan_db"
    echo " -> Envoi vers S3..."
    rclone copy "${JOB_SCRATCH}/metaphlan_db/" "${S3_DEST}/metaphlan_db/" -P --transfers 8
fi

# =====================================================================
# NETTOYAGE
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Nettoyage du Scratch..."
rm -rf "${JOB_SCRATCH}"

echo "--- TOUTES LES BASES SONT PRÊTES SUR LE S3 ---"