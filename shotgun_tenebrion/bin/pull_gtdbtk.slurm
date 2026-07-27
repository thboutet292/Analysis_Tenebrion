#!/bin/bash
#SBATCH --job-name=PULL_GTDB_IMG
#SBATCH --output=log/PULL_GTDB_IMG_%j.out
#SBATCH --error=log/PULL_GTDB_IMG_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal,lmge
#SBATCH --time=01:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET
#   SCRIPT  :   pull_gtdbtk.slurm
#   GOAL    :   Télécharger uniquement le conteneur du logiciel GTDB-Tk (v2.3.2).
#######################################################################################

set -euo pipefail

echo "========================================================"
echo "START -- $(date)"
echo "========================================================"

module purge
module load apptainer 2>/dev/null || true

# On redirige le cache vers le scratch pour éviter de bloquer ton /home
export APPTAINER_CACHEDIR="/storage/scratch/${USER}/apptainer_cache_gtdb"
export APPTAINER_TMPDIR="/storage/scratch/${USER}/apptainer_tmp_gtdb"
mkdir -p "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

mkdir -p containers

echo "-> Téléchargement de l'image GTDB-Tk (Docker)..."
apptainer pull containers/gtdbtk_v2.3.2.sif docker://ecogenomic/gtdbtk:2.3.2

echo "-> Nettoyage du cache..."
rm -rf "${APPTAINER_CACHEDIR}" "${APPTAINER_TMPDIR}"

echo "========================================================"
echo "END -- $(date)"
echo "Fichier créé : containers/gtdbtk_v2.3.2.sif"
echo "========================================================"