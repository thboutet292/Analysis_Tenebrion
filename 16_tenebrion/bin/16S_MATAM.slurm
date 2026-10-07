#!/bin/bash
#SBATCH --job-name=MATAM
#SBATCH --output=log/MATAM_%A_%a.out
#SBATCH --error=log/MATAM_%A_%a.err
#SBATCH --array=0-43     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=140G       
#SBATCH --partition=long
#SBATCH --time=10-00:00:00
#SBATCH --chdir=/home/thboutet/projets/16_tenebrion

#######################################################################################

#   AUTHORS :   Thomas BOUTET (IHP,LMGE)
#   DATE    :   2026
#   SCRIPT  :   16S_MATAM.slurm
#   TOOLS   :   MATAM 1.6.1 ; Salmon 1.10.2 

#######################################################################################

# sbatch bin/16S_MATAM.slurm

######################################################################################

set -uo pipefail

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INIT] DÉMARRAGE DU RUN PRODUCTION MATAM PE"
echo "====================================================================="

module purge
module load rclone/1.55.1
module load conda/23.3.1
eval "$(conda shell.bash hook)"
conda activate /home/thboutet/projets/16_tenebrion/resources/env_salmon

PROJECT_DIR="/home/thboutet/projets/16_tenebrion"
OUT_BASE="${PROJECT_DIR}/results/MATAM_PE_FINAL"
CONTAINERS_DIR="${PROJECT_DIR}/containers"
MATAM_SIF="${CONTAINERS_DIR}/matam_1.6.1--hdcf5f25_1.sif"
SINGULARITY="/usr/bin/singularity"
THREADS="${SLURM_CPUS_PER_TASK}"

# Paramètres S3
S3_REMOTE="s3_uca"
S3_BUCKET="lmge-tenebrion"
S3_DB_PATH="${S3_REMOTE}:${S3_BUCKET}/resources/matam_db_138_2_official"
DB_PREFIX_NAME="SILVA_138.2_SSURef_NR99_tax_silva_NR95"

# =====================================================================
# 1. GESTION DU SCRATCH ET DÉCOUVERTE DYNAMIQUE SUR S3
# =====================================================================
JOB_SCRATCH="/storage/scratch/${USER}/matam_pe_${SLURM_JOB_ID}"
WORKDIR="${JOB_SCRATCH}/work"
TMP_WORK="${JOB_SCRATCH}/tmp"
LOCAL_DB_DIR="${WORKDIR}/matam_db_138_2"

mkdir -p "${WORKDIR}" "${TMP_WORK}"
export TMPDIR="${TMP_WORK}"
export SINGULARITY_TMPDIR="${TMP_WORK}"

cleanup() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [CLEANUP] Suppression du scratch..."
    rm -rf "${JOB_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

# Récupération de la liste des fichiers R1 depuis le bucket S3
echo "Recherche des fichiers sur S3..."
FILES=($(rclone lsf "${S3_REMOTE}:${S3_BUCKET}/16s_bruts/" | grep '_R1.fastq.gz' | sort))
CURRENT_R1_NAME=${FILES[$SLURM_ARRAY_TASK_ID]}
SAMPLE=$(basename "${CURRENT_R1_NAME}" "_R1.fastq.gz")
CURRENT_R2_NAME="${SAMPLE}_R2.fastq.gz"

echo "Échantillon sélectionné pour le run PE : ${SAMPLE}"

OUT_DIR="${OUT_BASE}/${SAMPLE}"
mkdir -p "${OUT_DIR}" "log"

# Rapatriement des données (Fichiers FastQ bruts + Base MATAM)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Rapatriement des fastq depuis S3..."
rclone copy "${S3_REMOTE}:${S3_BUCKET}/16s_bruts/${CURRENT_R1_NAME}" "${WORKDIR}/"
rclone copy "${S3_REMOTE}:${S3_BUCKET}/16s_bruts/${CURRENT_R2_NAME}" "${WORKDIR}/"

CURRENT_R1="${WORKDIR}/${CURRENT_R1_NAME}"
CURRENT_R2="${WORKDIR}/${CURRENT_R2_NAME}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Rapatriement de SILVA 138.2 depuis S3..."
rclone copy "${S3_DB_PATH}" "${LOCAL_DB_DIR}" -P --transfers=4

matam_exec() {
    "${SINGULARITY}" exec --cleanenv \
        --env "LC_ALL=C" --env "LANG=C" --env "TMPDIR=/work/tmp" \
        -B "${WORKDIR}:/work" -B "${TMP_WORK}:/work/tmp" \
        "${MATAM_SIF}" bash -lc "$*"
}

# =====================================================================
# 2. NETTOYAGE PYTHON (COMPLET)
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ÉTAPE 1] Nettoyage complet des reads..."

PY_CLIP="${WORKDIR}/clip_fastq_pe.py"
cat > "${PY_CLIP}" <<'PY'
import sys, gzip

in_r1, in_r2, out_r1, out_r2, out_inter = sys.argv[1:6]
CLIP_LEN = 15
CAP = 74
MIN_Q = 34

# Lecture exhaustive : traitement de l'ensemble du fichier
with gzip.open(in_r1, "rt") as f1, gzip.open(in_r2, "rt") as f2, \
     open(out_r1, "w") as o1, open(out_r2, "w") as o2, open(out_inter, "w") as oi:
    while True:
        h1 = f1.readline(); s1 = f1.readline(); p1 = f1.readline(); q1 = f1.readline()
        h2 = f2.readline(); s2 = f2.readline(); p2 = f2.readline(); q2 = f2.readline()
        if not (s1 and s2): break 
        
        s1_cut = s1.rstrip("\r\n")[CLIP_LEN:]; q1_cut = q1.rstrip("\r\n")[CLIP_LEN:]
        s2_cut = s2.rstrip("\r\n")[CLIP_LEN:]; q2_cut = q2.rstrip("\r\n")[CLIP_LEN:]
        
        if len(s1_cut) >= 50 and len(s2_cut) >= 50:
            q1_fix = "".join(chr(max(MIN_Q, min(ord(ch), CAP))) for ch in q1_cut)
            q2_fix = "".join(chr(max(MIN_Q, min(ord(ch), CAP))) for ch in q2_cut)
            
            o1.write(f"{h1.rstrip()}\n{s1_cut}\n{p1.rstrip()}\n{q1_fix}\n")
            o2.write(f"{h2.rstrip()}\n{s2_cut}\n{p2.rstrip()}\n{q2_fix}\n")
            oi.write(f"{h1.rstrip()}\n{s1_cut}\n{p1.rstrip()}\n{q1_fix}\n")
            oi.write(f"{h2.rstrip()}\n{s2_cut}\n{p2.rstrip()}\n{q2_fix}\n")
PY

# Définition des chemins VUS DEPUIS LE CONTENEUR pour les fichiers d'entrée
CURRENT_R1_CONT="/work/${CURRENT_R1_NAME}"
CURRENT_R2_CONT="/work/${CURRENT_R2_NAME}"

# Définition des chemins VUS DEPUIS LE CONTENEUR pour les fichiers de sortie
CLIPPED_R1_CONT="/work/${SAMPLE}_clipped_R1.fastq" 
CLIPPED_R2_CONT="/work/${SAMPLE}_clipped_R2.fastq" 
CLIPPED_INTER_CONT="/work/${SAMPLE}_clipped_interleaved.fastq"

# On passe EXCLUSIVEMENT les chemins du conteneur au script Python
matam_exec "python /work/clip_fastq_pe.py '${CURRENT_R1_CONT}' '${CURRENT_R2_CONT}' '${CLIPPED_R1_CONT}' '${CLIPPED_R2_CONT}' '${CLIPPED_INTER_CONT}'" \
    || { echo "ERREUR: échec du nettoyage Python" >&2; exit 1; }

# =====================================================================
# 3. ASSEMBLAGE MATAM
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ÉTAPE 2] Assemblage MATAM Paired-End..."

MATAM_CONT_OUT="/work/02_matam"
MATAM_SCRATCH_OUT="${WORKDIR}/02_matam"
mkdir -p "${MATAM_SCRATCH_OUT}"

matam_exec "\
matam_assembly.py \
    --verbose \
    --cpu '${THREADS}' \
    --max_memory 8000 \
    --assembler SGA \
    --ref_db '/work/matam_db_138_2/${DB_PREFIX_NAME}' \
    --out_dir '${MATAM_CONT_OUT}' \
    --input_fastq '${CLIPPED_INTER_CONT}' \
" || { echo "ERREUR: échec MATAM" >&2; exit 1; }

SCAFFOLDS_FA="${MATAM_SCRATCH_OUT}/final_assembly.fa"
[[ -s "${SCAFFOLDS_FA}" ]] || { echo "ERREUR: MATAM n'a produit aucun scaffold" >&2; exit 1; }

# =====================================================================
# 4. QUANTIFICATION SALMON
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ÉTAPE 3] Quantification Salmon Paired-End..."

SALMON_INDEX="${WORKDIR}/salmon_index"
SALMON_OUT="${WORKDIR}/salmon_quant"

salmon index -t "${SCAFFOLDS_FA}" -i "${SALMON_INDEX}" -k 31 --no-version-check \
    || { echo "ERREUR: échec salmon index" >&2; exit 1; }

salmon quant -i "${SALMON_INDEX}" -l A \
    -1 "${WORKDIR}/${SAMPLE}_clipped_R1.fastq" \
    -2 "${WORKDIR}/${SAMPLE}_clipped_R2.fastq" \
    -p "${THREADS}" --validateMappings --no-version-check -o "${SALMON_OUT}" --quiet \
    || { echo "ERREUR: échec salmon quant" >&2; exit 1; }

# =====================================================================
# 5. RÉCUPÉRATION DES RÉSULTATS
# =====================================================================
echo "---------------------------------------------------------------------"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ÉTAPE 4] Récupération des résultats finaux..."
cp -L "${SCAFFOLDS_FA}" "${OUT_DIR}/${SAMPLE}_PE_scaffolds.fasta"
cp -L "${SALMON_OUT}/quant.sf" "${OUT_DIR}/${SAMPLE}_PE_abundance.tsv"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] JOB TERMINÉ AVEC SUCCÈS POUR ${SAMPLE}"
echo "====================================================================="
