#!/bin/bash
#SBATCH --job-name=BUILD_GTDB_S3
#SBATCH --output=log/BUILD_GTDB_S3_%j.out
#SBATCH --error=log/BUILD_GTDB_S3_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=normal,lmge
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   build_gtdbtk_db_s3.slurm
#
#   GOAL    :   Téléchargement propre de la base GTDB-Tk R214.
#######################################################################################

set -euo pipefail

echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE DE LA CONSTRUCTION DE LA BASE GTDB-TK R214"

module purge
module load rclone/1.55.1

# =====================================================================
# CONFIGURATION
# =====================================================================
# URL correcte avec le miroir européen et le sous-dossier 'auxillary_files'
DB_URL="https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxillary_files/gtdbtk_r214_data.tar.gz"

# Destination sur le S3
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_DEST="${S3_REMOTE}:${S3_BUCKET}/ref/gtdbtk_r214"

# Espace de travail Scratch
JOB_SCRATCH="/storage/scratch/${USER}/gtdb_db_build_${SLURM_JOB_ID}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

# =====================================================================
# ÉTAPE 1 : TÉLÉCHARGEMENT SUR LE SCRATCH
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 1] Téléchargement de l'archive (85 Go)..."
wget -c "${DB_URL}" -O gtdb_data.tar.gz

# =====================================================================
# ÉTAPE 2 : DÉCOMPRESSION
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ETAPE 2] Extraction des fichiers..."

mkdir -p db_files
tar -zxvf gtdb_data.tar.gz -C db_files/ --strip-components=1

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Nettoyage de l'archive compressée locale..."
rm gtdb_data.tar.gz

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
echo "[$(date '+%Y-%m-%d %H:%M:%S')] TERMINÉ. LA BASE GTDB R214 EST SUR LE S3."
echo "====================================================================="