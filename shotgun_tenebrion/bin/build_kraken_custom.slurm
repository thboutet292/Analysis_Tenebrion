#!/bin/bash
#SBATCH --job-name=BUILD_ADDON
#SBATCH --output=log/BUILD_ADDON_%j.out
#SBATCH --error=log/BUILD_ADDON_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G                    
#SBATCH --partition=normal,lmge
#SBATCH --time=06:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   build_kraken_custom.slurm
#   GOAL    :   Base Kraken2 "Add-on" ciblée (Fungi, Virus, Parasites)
#               Utilisation du BYPASS taxonomique pour accélération fulgurante.
#######################################################################################

set -uo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE DE LA CONSTRUCTION DE LA BASE ADDON"
echo "====================================================================="

module purge
module load rclone/1.55.1

# --- Configuration ---
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_DEST="${S3_REMOTE}:${S3_BUCKET}/ref/kraken2_insect_addon"

APPTAINER_BIN="/usr/bin/apptainer"
SIF_IMG="/home/thboutet/projets/16_tenebrion/containers/kraken2_bracken.sif"
THREADS="${SLURM_CPUS_PER_TASK}"

# --- Scratch ---
JOB_SCRATCH="/storage/scratch/${USER}/kraken_addon_build_${SLURM_JOB_ID}"
mkdir -p "${JOB_SCRATCH}"
cd "${JOB_SCRATCH}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Nettoyage du scratch..."
    cd /storage/scratch/${USER}
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

# =====================================================================
# ÉTAPE 1 : INSTALLATION DE NCBI DATASETS
# =====================================================================
echo "[ETAPE 1] Téléchargement de NCBI Datasets CLI..."
curl -s -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod +x datasets

# =====================================================================
# ÉTAPE 2 : TÉLÉCHARGEMENT ET INJECTION DES TAXID (BYPASS)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 2] Téléchargement des génomes d'intérêt..."
mkdir -p genomes

TAXA=(
    "Beauveria bassiana" "Beauveria" "Metarhizium anisopliae" "Metarhizium robertsii"
    "Metarhizium brunneum" "Purpureocillium lilacinum" "Cordyceps" "Ophiocordyceps"
    "Lecanicillium" "Meyerozyma guilliermondii" "Wickerhamomyces" "Starmerella bombicola"
    "Iteradensovirus" "Scindoavirus" "Iridovirus" "Gregarina" "Nosema"
    "Vairimorpha necatrix" "Crithidia fasciculata" "Crithidia bombi" "Leptomonas pyrrhocoris"
    "Methanobrevibacter cuticularis" "Methanobrevibacter curvatus" "Methanobrevibacter filiformis"
    "Baculoviridae"
    "Alphacypovirus"
    "Alphaentomopoxvirus"
    "Entomophthorales"
    "Tolypocladium"
)

SUCCESS_LIST=()
FAILED_LIST=()

for TAXON in "${TAXA[@]}"; do
    echo "-> Cible : ${TAXON}"
    SAFE_NAME=$(echo "${TAXON}" | tr ' ' '_')
    
    ./datasets download genome taxon "${TAXON}" --include genome --filename "${SAFE_NAME}.zip" --no-progressbar 2>/dev/null
    
    if [ -f "${SAFE_NAME}.zip" ]; then
        # Extraction du TaxID depuis les métadonnées
        TAXID=$(unzip -p "${SAFE_NAME}.zip" "ncbi_dataset/data/assembly_data_report.jsonl" 2>/dev/null | grep -o '"taxId": *[0-9]*' | head -1 | tr -dc '0-9' || true)
        
        if [ -z "$TAXID" ]; then
            FAILED_LIST+=("${TAXON} (Pas de TaxID trouvé)")
            continue
        fi
        
        mkdir -p tmp_extract
        unzip -q -j "${SAFE_NAME}.zip" "*.fna" -d tmp_extract/ 2>/dev/null
        
        for f in tmp_extract/*.fna; do
            if [ -f "$f" ]; then
                # Injection du TaxID dans le FASTA
                sed -i "s/^>/>kraken:taxid|${TAXID}|/g" "$f"
                mv "$f" "genomes/${SAFE_NAME}_$(basename "$f")"
            fi
        done
        rm -rf tmp_extract
        
        SUCCESS_LIST+=("${TAXON} (TaxID: ${TAXID})")
        echo "   [OK] Téléchargé et TaxID injecté."
    else
        FAILED_LIST+=("${TAXON}")
        echo "   [ECHEC] Introuvable sur NCBI."
    fi
done

# =====================================================================
# ÉTAPE 3 : TAXONOMIE LÉGÈRE (BYPASS)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 3] Configuration de la Taxonomie (Mode Bypass)..."

DB_DIR="K2_INSECT_ADDON"
mkdir -p "${DB_DIR}/taxonomy"

echo "   [3.1] Téléchargement de l'arbre taxonomique léger..."
wget -q -O taxdump.tar.gz "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
tar -xzf taxdump.tar.gz -C "${DB_DIR}/taxonomy" nodes.dmp names.dmp
rm taxdump.tar.gz

echo "   [3.2] Création de faux dictionnaires de traduction..."
touch "${DB_DIR}/taxonomy/nucl_gb.accession2taxid"
touch "${DB_DIR}/taxonomy/nucl_wgs.accession2taxid"
echo "   [OK] Taxonomie prête."

# =====================================================================
# ÉTAPE 4 : AJOUT ET CONSTRUCTION KRAKEN2
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 4] Ajout des FASTA et compilation (Build)..."

count=$(ls -1 genomes/*.fna 2>/dev/null | wc -l)
if [ "$count" -eq 0 ]; then
    echo "[ERREUR] Aucun fichier FASTA n'a pu être téléchargé. Arrêt."
    exit 1
fi

# 4.1 Ajout à la bibliothèque
for f in genomes/*.fna; do
    ${APPTAINER_BIN} exec --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2-build --add-to-library "$f" --db "${DB_DIR}" > /dev/null
done

# 4.2 Compilation de la base
${APPTAINER_BIN} exec --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2-build --build --threads "${THREADS}" --db "${DB_DIR}"

# 4.3 Nettoyage interne pour gagner de la place
${APPTAINER_BIN} exec --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2-build --clean --db "${DB_DIR}" > /dev/null

# =====================================================================
# ÉTAPE 5 : EXPORT VERS LE S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 5] Transfert de la nouvelle base vers le S3..."
rclone copy "${DB_DIR}/" "${S3_DEST}/" -P --transfers 8

echo "====================================================================="
echo " BILAN DES TÉLÉCHARGEMENTS (NCBI DATASETS) "
echo "====================================================================="
echo "✅ TAXONS TÉLÉCHARGÉS ET INTÉGRÉS :"
for t in "${SUCCESS_LIST[@]:-}"; do echo "  - $t"; done

echo " "
echo "❌ TAXONS INTROUVABLES OU ÉCHOUÉS :"
for t in "${FAILED_LIST[@]:-}"; do echo "  - $t"; done
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] Script terminé."