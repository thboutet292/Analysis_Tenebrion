#!/bin/bash
#SBATCH --job-name=MERGE_DRAM
#SBATCH --output=log/MERGE_DRAM_%j.out
#SBATCH --error=log/MERGE_DRAM_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=normal,lmge
#SBATCH --time=00:30:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   GOAL : Fusionner les résultats de taxonomie GTDB-Tk avec les annotations
#          et les statistiques métaboliques générées par DRAM pour chaque MAG.
#   OUT  : results/final_bins_summary.tsv (Tableau Tidy Data formaté pour R)
#######################################################################################

set -euo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE : FUSION DRAM + GTDB-Tk"
echo "Job ID : ${SLURM_JOB_ID}"
echo "====================================================================="

module purge
module load apptainer 2>/dev/null || true

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
IMG_DRAM="${PROJECT_DIR}/containers/dram_v1.4.6.sif"
APPTAINER_BIN=$(command -v apptainer || command -v singularity)

# Vérification de l'image
if [ ! -f "${IMG_DRAM}" ]; then
    echo "[ERREUR] Image DRAM introuvable : ${IMG_DRAM}" >&2
    exit 1
fi

mkdir -p log/

# =====================================================================
# SCRIPT PYTHON EMBARQUÉ (Here-Document)
# =====================================================================
# Plutôt que d'avoir un fichier .py séparé, on l'écrit dynamiquement
# dans un fichier temporaire pour que tout soit encapsulé dans ce .slurm

PYTHON_SCRIPT="merge_temp_${SLURM_JOB_ID}.py"

cat << 'EOF' > "${PYTHON_SCRIPT}"
import pandas as pd
import glob
import os
import sys

print("Démarrage de la fusion des données par le script Python...")

# 1. Chargement de la taxonomie GTDB-Tk
print(" -> Chargement de GTDB-Tk...")
gtdb_files = glob.glob("results/gtdbtk_taxo_concoct/*.summary.tsv")
if not gtdb_files:
    print("[ERREUR] Aucun fichier GTDB-Tk trouvé dans results/gtdbtk_taxo_concoct/")
    sys.exit(1)

gtdb_list = [pd.read_csv(f, sep='\t') for f in gtdb_files]
gtdb_df = pd.concat(gtdb_list, ignore_index=True)

# Nettoyage GTDB : on garde l'ID et la taxonomie, on renomme pour matcher DRAM
gtdb_df = gtdb_df[['user_genome', 'classification']]
gtdb_df = gtdb_df.rename(columns={'user_genome': 'genome'})

# 2. Chargement des statistiques génomiques (DRAM genome_stats.tsv)
print(" -> Chargement des statistiques DRAM...")
stat_files = glob.glob("results/dram_concoct_distill/*/genome_stats.tsv")
if not stat_files:
    print("[ERREUR] Aucun fichier genome_stats.tsv trouvé.")
    sys.exit(1)

stat_list = [pd.read_csv(f, sep='\t') for f in stat_files]
stat_df = pd.concat(stat_list, ignore_index=True)

# 3. Chargement des annotations métaboliques (DRAM product.tsv)
print(" -> Chargement des fonctions DRAM...")
prod_files = glob.glob("results/dram_concoct_distill/*/product.tsv")
if not prod_files:
    print("[ERREUR] Aucun fichier product.tsv trouvé.")
    sys.exit(1)

prod_list = [pd.read_csv(f, sep='\t') for f in prod_files]
prod_df = pd.concat(prod_list, ignore_index=True)

# 4. Fusions (Merge)
print(" -> Fusion des tableaux (Outer Merge sur genome)...")
# Fusion stat + prod (si certains bins n'ont pas de fonctions, on garde quand même leurs stats)
dram_merged = pd.merge(stat_df, prod_df, on='genome', how='outer')

# Fusion finale avec GTDB-Tk (on garde tous les génomes identifiés par GTDB)
final_df = pd.merge(gtdb_df, dram_merged, on='genome', how='right')

# Remplacement des NaN par 0 pour les colonnes numériques (scores DRAM)
# Cela facilite énormément l'analyse quantitative sous R.
# On identifie les colonnes numériques avant de remplir
cols_to_fill = final_df.columns.difference(['genome', 'classification'])
final_df[cols_to_fill] = final_df[cols_to_fill].fillna(0)

# 5. Sauvegarde
output_file = "results/final_bins_summary_concoct.tsv"
final_df.to_csv(output_file, sep='\t', index=False)

print(f"Terminé ! Fichier combiné généré : {output_file}")
print(f"Nombre total de MAGs dans le tableau : {len(final_df)}")
EOF

# =====================================================================
# EXÉCUTION DU SCRIPT PYTHON DANS LE CONTENEUR
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exécution de la fusion via Pandas (Apptainer)..."

"${APPTAINER_BIN}" exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    "${IMG_DRAM}" \
    python "${PYTHON_SCRIPT}"

# Nettoyage du script python temporaire
rm -f "${PYTHON_SCRIPT}"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] FUSION TERMINÉE"
echo "Fichier de résultat : results/final_bins_summary.tsv"
echo "====================================================================="