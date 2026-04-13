#!/bin/bash
#SBATCH --job-name=RUN_PICRUST2
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=log/picrust2_run_%j.out
#SBATCH --error=log/picrust2_run_%j.err

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   DATE    :   16 - 03 - 2026
#   SCRIPT  :   16S_PICRUST2.slurm
#   TOOLS   :   PICRUSt2 2.5.x (via Singularity) ; Python 3.x
#
#   GOAL    :   Prédiction fonctionnelle métagénomique à partir des scaffolds MATAM.
#               - Préparation des données (FASTA + table de comptages entiers)
#               - Pipeline PICRUSt2 complet : placement phylogénétique,
#                 prédiction génomique, pathways MetaCyc, EC numbers, KO terms
#               - Extraction des scores NSTI par ASV (diagnostic qualité)
#               - Décompression et organisation des sorties pour les scripts R
#
#   INPUT   :   data/MATAM/all_matam_sequences_compiled.tsv
#
#   OUTPUTS :   results/PICRUST2/
#                 path_abun_unstrat.tsv      — pathways MetaCyc, non-stratifié
#                 path_abun_contrib.tsv      — pathways MetaCyc, stratifié par ASV
#                 EC_predicted_unstrat.tsv   — EC numbers, non-stratifié
#                 EC_predicted_contrib.tsv   — EC numbers, stratifié par ASV
#                 KO_predicted_unstrat.tsv   — KO terms, non-stratifié
#                 KO_predicted_contrib.tsv   — KO terms, stratifié par ASV
#                 picrust2_nsti_per_asv.tsv  — score NSTI par ASV
#                 picrust2_seqs.fasta        — séquences ASV (copie)
#                 picrust2_counts.tsv        — table de comptages (copie)
#######################################################################################

set -euo pipefail

# ── Chemins ───────────────────────────────────────────────────────────────────
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
CONTAINER="${PROJECT_DIR}/containers/picrust2.sif"
INPUT_TSV="${PROJECT_DIR}/results/all_matam_sequences_compiled.tsv"
OUT_FASTA="${PROJECT_DIR}/results/PICRUST2/picrust2_seqs.fasta"
OUT_COUNTS="${PROJECT_DIR}/results/PICRUST2/picrust2_counts.tsv"
PICRUST2_OUT="${PROJECT_DIR}/results/PICRUST2/PICRUST2_OUT"
RESULTS_DIR="${PROJECT_DIR}/results/PICRUST2"

mkdir -p "${PROJECT_DIR}/log"
mkdir -p "${PROJECT_DIR}/data/PICRUST2"
mkdir -p "${RESULTS_DIR}"

echo "======================================================================"
echo "  PICRUSt2 PIPELINE — $(date)"
echo "  Job ID     : ${SLURM_JOB_ID}"
echo "  Nœud       : $(hostname)"
echo "  CPUs       : ${SLURM_CPUS_PER_TASK}"
echo "  Répertoire : ${PROJECT_DIR}"
echo "======================================================================"

# ─────────────────────────────────────────────────────────────────────────────
# ÉTAPE 1 — PRÉPARATION DES DONNÉES
# Lit all_matam_sequences_compiled.tsv et produit :
#   - picrust2_seqs.fasta   : une entrée par séquence unique (ASV_1, ASV_2, ...)
#   - picrust2_counts.tsv   : matrice OTU_ID × échantillons, comptages entiers
# Note : PICRUSt2 requiert des entiers dans la table de comptages
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "=== ÉTAPE 1/4 : Préparation FASTA + matrice de comptages ==="

python3 << PYEOF
import csv, os, sys

input_file  = "${INPUT_TSV}"
out_fasta   = "${OUT_FASTA}"
out_counts  = "${OUT_COUNTS}"

seq_to_id    = {}
sample_names = set()
counts       = {}
n_rows       = 0

print(f"  Lecture : {input_file}")

with open(input_file, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        seq    = row["Sequence"]
        sample = row["Sample"]
        reads  = float(row["NumReads"])
        sample_names.add(sample)
        if seq not in seq_to_id:
            seq_to_id[seq] = f"ASV_{len(seq_to_id) + 1}"
        seq_id = seq_to_id[seq]
        if seq_id not in counts:
            counts[seq_id] = {}
        counts[seq_id][sample] = counts[seq_id].get(sample, 0) + reads
        n_rows += 1

print(f"  Lignes lues          : {n_rows}")
print(f"  ASVs uniques         : {len(seq_to_id)}")
print(f"  Échantillons uniques : {len(sample_names)}")

# FASTA
with open(out_fasta, "w") as f:
    for seq, seq_id in seq_to_id.items():
        f.write(f">{seq_id}\n{seq}\n")
print(f"  FASTA écrit : {out_fasta}")

# Matrice de comptages — entiers requis par PICRUSt2
samples = sorted(list(sample_names))
with open(out_counts, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["OTU_ID"] + samples)
    for seq, seq_id in seq_to_id.items():
        row = [seq_id] + [int(counts[seq_id].get(s, 0)) for s in samples]
        writer.writerow(row)
print(f"  Comptages écrits : {out_counts}")
PYEOF

echo "  -> Étape 1 terminée : $(date)"

# ─────────────────────────────────────────────────────────────────────────────
# ÉTAPE 2 — PIPELINE PICRUST2 COMPLET
#
# --stratified          : produit les fichiers stratifiés par ASV (contrib)
#                         nécessaire pour path_abun_contrib et EC/KO contrib
# --max_nsti 2.0        : filtre les ASVs trop divergents des références
#                         (NSTI > 2 = prédiction peu fiable)
# --min_reads 1         : conserve tous les ASVs avec ≥ 1 read
# --min_samples 1       : conserve les ASVs présents dans ≥ 1 échantillon
#
# Sorties automatiques dans PICRUST2_OUT/ :
#   pathways_out/path_abun_unstrat.tsv.gz
#   pathways_out/path_abun_strat.tsv.gz
#   EC_metagenome_out/pred_metagenome_unstrat.tsv.gz
#   EC_metagenome_out/pred_metagenome_strat.tsv.gz
#   KO_metagenome_out/pred_metagenome_unstrat.tsv.gz
#   KO_metagenome_out/pred_metagenome_strat.tsv.gz
#   intermediate/  (fichiers intermédiaires dont NSTI)
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "=== ÉTAPE 2/4 : Pipeline PICRUSt2 ==="
echo "  max_nsti : 2.0 | stratified : oui | CPUs : ${SLURM_CPUS_PER_TASK}"

# Suppression d'un éventuel run précédent
if [ -d "${PICRUST2_OUT}" ]; then
    echo "  Suppression de l'ancien dossier PICRUST2_OUT..."
    rm -rf "${PICRUST2_OUT}"
fi

singularity exec --cleanenv \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    "${CONTAINER}" \
    picrust2_pipeline.py \
        --study_fasta  "${OUT_FASTA}"            \
        --input        "${OUT_COUNTS}"           \
        --output       "${PICRUST2_OUT}"         \
        --processes    "${SLURM_CPUS_PER_TASK}"  \
        --max_nsti     2.0                       \
        --min_reads    1                         \
        --min_samples  1                         \
        --stratified                             \
        --verbose

echo "  -> Étape 2 terminée : $(date)"

# ─────────────────────────────────────────────────────────────────────────────
# ÉTAPE 3 — EXTRACTION DES SCORES NSTI PAR ASV
#
# Le NSTI (Nearest Sequenced Taxon Index) est la distance phylogénétique entre
# un ASV et son plus proche génome de référence dans PICRUSt2.
# NSTI faible  (< 0.1) = bonne référence, prédiction fiable
# NSTI moyen   (0.1–2) = prédiction acceptable
# NSTI élevé   (> 2)   = filtré par --max_nsti, prédiction non fiable
#
# Ce fichier est essentiel pour croiser avec MATAM et expliquer pourquoi
# certains genres (ex: Staphylococcus) sont sous-représentés dans PICRUSt2
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "=== ÉTAPE 3/4 : Extraction des scores NSTI ==="

# Lister les fichiers intermédiaires pour diagnostic
echo "  Fichiers intermédiaires disponibles :"
find "${PICRUST2_OUT}/intermediate/" -name "*.tsv*" 2>/dev/null | sort | head -30

python3 << PYEOF
import os, gzip, csv

picrust2_out = "${PICRUST2_OUT}"
results_dir  = "${RESULTS_DIR}"
out_nsti     = os.path.join(results_dir, "picrust2_nsti_per_asv.tsv")

# Chemins NSTI selon version PICRUSt2
# Le fichier marker_predicted_and_normalized contient une colonne metadata_NSTI
nsti_candidates = [
    os.path.join(picrust2_out, "intermediate",
                 "marker_predicted_and_normalized.tsv.gz"),
    os.path.join(picrust2_out, "intermediate",
                 "place_seqs", "query_filt.tsv"),
    os.path.join(picrust2_out, "intermediate",
                 "place_seqs", "query_filt.tsv.gz"),
]

nsti_found = None
for c in nsti_candidates:
    if os.path.exists(c):
        nsti_found = c
        print(f"  Fichier NSTI trouvé : {c}")
        break

if nsti_found is None:
    print("  ATTENTION : fichier NSTI non trouvé.")
    print("  Lancez manuellement : find PICRUST2_OUT/intermediate/ -name '*.tsv*'")
    exit(0)

rows = []
opener = gzip.open if nsti_found.endswith(".gz") else open
with opener(nsti_found, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        # Chercher la colonne NSTI selon la version
        nsti_val = None
        for col in ["metadata_NSTI", "NSTI", "nsti", "metadata_nsti"]:
            if col in row:
                nsti_val = row[col]
                break
        if nsti_val is None:
            print("  ATTENTION : colonne NSTI non trouvée dans le fichier.")
            print(f"  Colonnes disponibles : {list(row.keys())}")
            break
        asv_id = row.get("sequence", row.get("#OTU_ID", row.get("taxon", "")))
        rows.append({"ASV_ID": asv_id, "NSTI": nsti_val})

if rows:
    with open(out_nsti, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["ASV_ID","NSTI"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    vals = [float(r["NSTI"]) for r in rows if r["NSTI"] not in ("", "NA")]
    print(f"  ASVs avec NSTI       : {len(vals)}")
    print(f"  NSTI moyen           : {sum(vals)/len(vals):.4f}")
    print(f"  NSTI médian          : {sorted(vals)[len(vals)//2]:.4f}")
    print(f"  NSTI max             : {max(vals):.4f}")
    print(f"  ASVs filtrés > 2.0   : {sum(1 for v in vals if v > 2.0)}")
    print(f"  Fichier écrit        : {out_nsti}")
PYEOF

echo "  -> Étape 3 terminée : $(date)"

# ─────────────────────────────────────────────────────────────────────────────
# ÉTAPE 4 — DÉCOMPRESSION ET COPIE VERS results/PICRUST2/
#
# Fonction utilitaire : copie + décompression si .gz, copie simple sinon
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "=== ÉTAPE 4/4 : Organisation des fichiers de sortie ==="

copy_and_decompress() {
    local src="$1"
    local dst="$2"
    local label="$3"

    if [ -f "${src}.gz" ]; then
        gunzip -c "${src}.gz" > "${dst}"
        echo "  -> ${label}"
    elif [ -f "${src}" ]; then
        cp "${src}" "${dst}"
        echo "  -> ${label}"
    else
        echo "  ATTENTION : ${label} non trouvé (${src})"
    fi
}

# ── Pathways MetaCyc ─────────────────────────────────────────────────────────
copy_and_decompress \
    "${PICRUST2_OUT}/pathways_out/path_abun_unstrat.tsv" \
    "${RESULTS_DIR}/path_abun_unstrat.tsv" \
    "path_abun_unstrat.tsv"

copy_and_decompress \
    "${PICRUST2_OUT}/pathways_out/path_abun_strat.tsv" \
    "${RESULTS_DIR}/path_abun_contrib.tsv" \
    "path_abun_contrib.tsv"

# ── EC numbers ───────────────────────────────────────────────────────────────
copy_and_decompress \
    "${PICRUST2_OUT}/EC_metagenome_out/pred_metagenome_unstrat.tsv" \
    "${RESULTS_DIR}/EC_predicted_unstrat.tsv" \
    "EC_predicted_unstrat.tsv"

copy_and_decompress \
    "${PICRUST2_OUT}/EC_metagenome_out/pred_metagenome_strat.tsv" \
    "${RESULTS_DIR}/EC_predicted_contrib.tsv" \
    "EC_predicted_contrib.tsv"

# ── KO terms ─────────────────────────────────────────────────────────────────
copy_and_decompress \
    "${PICRUST2_OUT}/KO_metagenome_out/pred_metagenome_unstrat.tsv" \
    "${RESULTS_DIR}/KO_predicted_unstrat.tsv" \
    "KO_predicted_unstrat.tsv"

copy_and_decompress \
    "${PICRUST2_OUT}/KO_metagenome_out/pred_metagenome_strat.tsv" \
    "${RESULTS_DIR}/KO_predicted_contrib.tsv" \
    "KO_predicted_contrib.tsv"

# ── Fichiers de base ─────────────────────────────────────────────────────────
cp "${OUT_FASTA}"  "${RESULTS_DIR}/picrust2_seqs.fasta"
cp "${OUT_COUNTS}" "${RESULTS_DIR}/picrust2_counts.tsv"
echo "  -> picrust2_seqs.fasta"
echo "  -> picrust2_counts.tsv"

# ── Récapitulatif ─────────────────────────────────────────────────────────────
echo ""
echo "======================================================================"
echo "  FICHIERS PRODUITS DANS ${RESULTS_DIR}/"
echo "======================================================================"
ls -lh "${RESULTS_DIR}/" | awk '{print "  " $0}'

echo ""
echo "======================================================================"
echo "  PIPELINE TERMINÉ — $(date)"
echo "  Job ID : ${SLURM_JOB_ID}"
echo "======================================================================"