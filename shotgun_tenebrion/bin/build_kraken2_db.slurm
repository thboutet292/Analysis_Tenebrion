#!/bin/bash
#SBATCH --job-name=BUILD_K2_PFP
#SBATCH --output=log/BUILD_DB_PFP_%j.out
#SBATCH --error=log/BUILD_DB_PFP_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=normal,lmge
#SBATCH --time=12:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   DATE    :   27 - 03 - 2026
#   SCRIPT  :   build_kraken_db_pluspfp.slurm
#
#   GOAL    :   Téléchargement propre de la base Kraken2 PlusPFP.
#               (Standard + Protozoa + Fungi + Plant)
#               - Utilisation du /scratch pour éviter la saturation du /home.
#               - Décompression de l'archive.
#               - Transfert définitif vers le stockage objet (S3).
#######################################################################################

set -uo pipefail

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE DE LA CONSTRUCTION DE LA BASE KRAKEN2 PLUSPFP"

module purge
module load rclone/1.55.1

# =====================================================================
# CONFIGURATION
# =====================================================================
# URL officielle de la base PlusPFP (Ben Langmead)
DB_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20240112.tar.gz"

# Destination sur le S3 (Nouveau dossier pour ne pas écraser l'ancienne)
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_DEST="${S3_REMOTE}:${S3_BUCKET}/ref/kraken2_pluspfp"

# Espace de travail Scratch
JOB_SCRATCH="/storage/scratch/${USER}/kraken_db_build_${SLURM_JOB_ID}_PFP"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

# =====================================================================
# ÉTAPE 1 : TÉLÉCHARGEMENT SUR LE SCRATCH
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Téléchargement de l'archive (approx. 75 Go)..."
wget -c "${DB_URL}" -O pluspfp.tar.gz

# =====================================================================
# ÉTAPE 2 : DÉCOMPRESSION
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Extraction des fichiers .k2d..."

mkdir -p db_files
tar -zxvf pluspfp.tar.gz -C db_files/

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Nettoyage de l'archive compressée locale..."
rm pluspfp.tar.gz

# =====================================================================
# ÉTAPE 3 : UPLOAD VERS LE S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 3] Transfert vers le S3..."
rclone copy "db_files/" "${S3_DEST}/" -P --transfers 8

# =====================================================================
# NETTOYAGE FINAL
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Nettoyage du Scratch..."
cd /storage/scratch/${USER}
rm -rf "${JOB_SCRATCH}"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] TERMINÉ. LA BASE PLUSPFP EST DISPONIBLE SUR LE S3."
echo "Chemin : ${S3_DEST}"
echo "====================================================================="