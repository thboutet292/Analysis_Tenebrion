#!/bin/bash
#SBATCH --job-name=BUILD_ADDON
#SBATCH --output=log/BUILD_ADDON_%j.out
#SBATCH --error=log/BUILD_ADDON_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G                    
#SBATCH --partition=normal,lmge
#SBATCH --time=04:00:00
#SBATCH --chdir=/home/thboutet/projets/shotgun_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   SCRIPT  :   build_kraken_addon_final.slurm
#   GOAL    :   Créer la base Kraken2 Add-on pour les taxons d'intérêts qui ne sont pas dans la BDD Kraken2  
#               - Utilise les noms de genres pour contourner l'ICTV/NCBI.
#               - Injection directe des TaxID pour garantir la compilation.
#######################################################################################

set -uo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] DÉMARRAGE DE LA CONSTRUCTION DE LA BASE ADDON"
echo "====================================================================="

module purge
module load rclone/1.55.1

S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_DEST="${S3_REMOTE}:${S3_BUCKET}/ref/kraken2_insect_addon"

APPTAINER_BIN="/usr/bin/apptainer"
SIF_IMG="/home/thboutet/projets/16_tenebrion/containers/kraken2_bracken.sif"
THREADS="${SLURM_CPUS_PER_TASK}"

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
# ÉTAPE 1 : INSTALLATION DES OUTILS
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 1] Téléchargement de NCBI Datasets CLI..."
curl -s -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
chmod +x datasets

# =====================================================================
# ÉTAPE 2 : TÉLÉCHARGEMENT & INJECTION DES TAXIDS
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 2] Téléchargement et formatage des en-têtes FASTA..."
mkdir -p genomes

# Liste des absents (mise à jour avec les genres pour ceux qui plantaient)
TAXA=(
    "Beauveria brongniartii"
    "Metarhizium anisopliae"
    "Metarhizium robertsii"
    "Metarhizium brunneum"
    "Purpureocillium lilacinum"
    "Cordyceps"            
    "Ophiocordyceps"       
    "Lecanicillium"        
    "Meyerozyma guilliermondii"
    "Wickerhamomyces"      
    "Starmerella bombicola"
    "Iteradensovirus"      
    "Scindoavirus"         
    "Iridovirus"           
    "Gregarina"            
    "Nosema"               
    "Vairimorpha necatrix"
    "Crithidia fasciculata"
    "Crithidia bombi"
    "Leptomonas pyrrhocoris"
    "Methanobrevibacter cuticularis"
    "Methanobrevibacter curvatus"
    "Methanobrevibacter filiformis"
)

SUCCESS_LIST=()
FAILED_LIST=()

for TAXON in "${TAXA[@]}"; do
    echo "-> Traitement de : ${TAXON}"
    SAFE_NAME=$(echo "${TAXON}" | tr ' ' '_')
    
    ./datasets download genome taxon "${TAXON}" --include genome --filename "${SAFE_NAME}.zip" --no-progressbar 2>/dev/null
    
    if [ -f "${SAFE_NAME}.zip" ]; then
        mkdir -p tmp_extract
        unzip -q -j "${SAFE_NAME}.zip" "*.fna" "ncbi_dataset/data/assembly_data_report.jsonl" -d tmp_extract/ 2>/dev/null
        
        # Récupération du TaxID
        TAXID=$(grep -o '"taxId":[0-9]*' tmp_extract/assembly_data_report.jsonl | head -n 1 | cut -d':' -f2)
        
        if [ -n "$TAXID" ]; then
            for f in tmp_extract/*.fna; do
                if [ -f "$f" ]; then
                    NEW_NAME="genomes/${SAFE_NAME}_$(basename "$f")"
                    mv "$f" "$NEW_NAME"
                    # Injection du TaxID dans l'en-tête FASTA
                    sed -i "s/^>/>kraken:taxid|${TAXID}| /" "$NEW_NAME"
                fi
            done
            SUCCESS_LIST+=("${TAXON} (TaxID: ${TAXID})")
            echo "   [OK] Injecté avec le TaxID ${TAXID}."
        else
            FAILED_LIST+=("${TAXON} (TaxID introuvable)")
            echo "   [ERREUR] Impossible de trouver le TaxID."
        fi
        rm -rf tmp_extract "${SAFE_NAME}.zip"
    else
        FAILED_LIST+=("${TAXON} (Echec téléchargement)")
        echo "   [ECHEC] Introuvable sur NCBI."
    fi
done

# =====================================================================
# ÉTAPE 3 : TAXONOMIE MANUELLE
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 3] Téléchargement manuel de la taxonomie (Squelette)..."
DB_DIR="K2_INSECT_ADDON"
mkdir -p "${DB_DIR}/taxonomy"
wget -q -O taxdump.tar.gz "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
tar -xzf taxdump.tar.gz -C "${DB_DIR}/taxonomy" nodes.dmp names.dmp
rm taxdump.tar.gz

# =====================================================================
# ÉTAPE 4 : BUILD KRAKEN2
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 4] Ajout des séquences marquées et compilation..."

for f in genomes/*.fna; do
    ${APPTAINER_BIN} exec --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2-build \
        --add-to-library "$f" \
        --db "${DB_DIR}"
done

${APPTAINER_BIN} exec --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2-build \
    --build \
    --threads "${THREADS}" \
    --db "${DB_DIR}"

${APPTAINER_BIN} exec --bind "${PWD}:${PWD}" "${SIF_IMG}" kraken2-build \
    --clean \
    --db "${DB_DIR}"

# =====================================================================
# ÉTAPE 5 : EXPORT S3
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[ETAPE 5] Transfert vers le S3..."
rclone copy "${DB_DIR}/" "${S3_DEST}/" -P

echo "====================================================================="
echo "TAXONS VALIDÉS DANS LA BASE :"
for t in "${SUCCESS_LIST[@]:-}"; do echo "  - $t"; done
echo " "
echo "TAXONS INTROUVABLES OU ÉCHOUÉS :"
for t in "${FAILED_LIST[@]:-}"; do echo "  - $t"; done
echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FIN] Script terminé."