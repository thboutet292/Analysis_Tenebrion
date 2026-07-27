#!/bin/bash
#SBATCH --job-name=KRAKEN_TO_QIIME
#SBATCH --output=log/KRAKEN_TO_QIIME_%j.out
#SBATCH --error=log/KRAKEN_TO_QIIME_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=normal,lmge
#SBATCH --time=01:00:00

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP, LMGE)
#   DATE    :   Mars 2026
#   SCRIPT  :   kraken2_to_qiime2.slurm
#
#   GOAL    :   Conversion stricte des rapports Kraken2/Bracken en matrices QIIME2.
#               Génère des tables en counts absolus et abondances relatives,
#               sans duplication cumulative (utilisation stricte des taxon_reads).
#######################################################################################

set -euo pipefail
module purge
mkdir -p log

PROJECT_DIR="${PWD}"
IN_DIR="${PROJECT_DIR}/results"
OUT_DIR="${PROJECT_DIR}/results/qiime2_tables_addon"
mkdir -p "${OUT_DIR}"

echo "====================================================================="
echo "[$(date '+%Y-%m-%d %H:%M:%S')] START KRAKEN/BRACKEN TO QIIME2 FORMAT"
echo "IN_DIR : ${IN_DIR}"
echo "OUT_DIR: ${OUT_DIR}"
echo "====================================================================="

python3 - << 'PY' "${IN_DIR}" "${OUT_DIR}"
import os, sys, glob
from collections import defaultdict

in_dir  = sys.argv[1]
out_dir = sys.argv[2]

# -------------------------
# PARAMÈTRES DE FILTRAGE
# -------------------------
MIN_COUNT_PER_SAMPLE   = 10   
MIN_TOTAL_COUNT_GLOBAL = 10   
FILTER_ONLY_BACTERIA   = False
CONTAM_KEYWORDS        = ()

# -------------------------
# CORRESPONDANCES TAXONOMIQUES
# -------------------------
rank_map = {
    "D": "d", "K": "d",
    "P": "p", "C": "c", "O": "o", "F": "f", "G": "g", "S": "s",
}
rank_order = ["d", "p", "c", "o", "f", "g", "s"]
EXCLUDE_AS_FEATURE = {"R", "R1", "R2", "U"}

# -------------------------
# CLASSE NOEUD TAXONOMIQUE
# -------------------------
class Node:
    __slots__ = ("name","taxid","rank_code","clade_reads","taxon_reads","depth","parent","children")
    def __init__(self, name, taxid, rank_code, clade_reads, taxon_reads, depth):
        self.name        = name
        self.taxid       = taxid
        self.rank_code   = rank_code
        self.clade_reads = clade_reads
        self.taxon_reads = taxon_reads
        self.depth       = depth
        self.parent      = None
        self.children    = []

# -------------------------
# PARSING D'UN KREPORT
# -------------------------
def parse_kreport(path):
    nodes, stack = [], []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            try:
                cr = int(parts[1])
                tr = int(parts[2])
            except ValueError:
                continue
            rank_code = parts[3].strip()
            taxid     = parts[4].strip()
            raw_name  = parts[5]
            depth     = (len(raw_name) - len(raw_name.lstrip(" "))) // 2
            name      = raw_name.strip()
            n = Node(name, taxid, rank_code, cr, tr, depth)
            while stack and stack[-1].depth >= depth:
                stack.pop()
            if stack:
                n.parent = stack[-1]
                stack[-1].children.append(n)
            stack.append(n)
            nodes.append(n)
    return nodes

def get_total_reads_from_kreport(path):
    total = 0
    if not os.path.isfile(path):
        return 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            if parts[3].strip() in ("U", "R"):
                try:
                    total += int(parts[1])
                except ValueError:
                    pass
    return total

def lineage_from_node(node):
    line = {k: "" for k in rank_order}
    cur = node
    while cur:
        if cur.rank_code in rank_map:
            r = rank_map[cur.rank_code]
            if line[r] == "":
                line[r] = cur.name
        cur = cur.parent
    return line

def to_qiime_string(line):
    return ";".join(
        f"{r}__{line[r]}" if line[r] else f"{r}__" for r in rank_order
    )

def is_contaminant(feat):
    return any(k in feat.lower() for k in CONTAM_KEYWORDS)

def is_bacteria(feat):
    return feat.startswith("d__Bacteria")

# -------------------------------------------------------------------------------------
# NOUVELLE FONCTION D'EXTRACTION STRICTE (CORRECTION DU BUG DE DUPLICATION)
# -------------------------------------------------------------------------------------
def collect_terminal_nodes(nodes):
    """
    Collecte UNIQUEMENT les reads propres assignés spécifiquement à un noeud (taxon_reads).
    Cela élimine mathématiquement les doublons cumulatifs (clade_reads).
    """
    result = []
    for n in nodes:
        if n.rank_code in EXCLUDE_AS_FEATURE or n.rank_code not in rank_map:
            continue
        
        # MAGIE ICI : On ne lit que la colonne 3.
        # En Bracken, seuls les niveaux d'espèce auront des taxon_reads > 0.
        # Les niveaux intermédiaires (Metazoa, etc.) auront 0 et disparaîtront proprement.
        if n.taxon_reads > 0:
            result.append((n, float(n.taxon_reads)))
            
    return result

def apply_filters(counts, samples):
    for feat in list(counts):
        for s in samples:
            if counts[feat].get(s, 0.0) < MIN_COUNT_PER_SAMPLE:
                counts[feat][s] = 0.0
        total = sum(counts[feat].values())
        if total < MIN_TOTAL_COUNT_GLOBAL or all(v == 0 for v in counts[feat].values()):
            del counts[feat]
    return counts

def write_qiime_table(features, counts, samples, path_counts, path_rel, total_reads):
    with open(path_counts, "w", encoding="utf-8") as out:
        out.write("# Constructed from biom file\n")
        out.write("#OTU ID\t" + "\t".join(samples) + "\n")
        for feat in features:
            row = "\t".join(f"{counts[feat].get(s, 0.0):.1f}" for s in samples)
            out.write(f"{feat}\t{row}\n")

    with open(path_rel, "w", encoding="utf-8") as out:
        out.write("# Constructed from biom file\n")
        out.write("#OTU ID\t" + "\t".join(samples) + "\n")
        for feat in features:
            row = "\t".join(
                "0.00000000" if total_reads.get(s, 0) == 0
                else f"{counts[feat].get(s, 0.0) / total_reads[s]:.8f}"
                for s in samples
            )
            out.write(f"{feat}\t{row}\n")

# ===========================================================================
# DÉTECTION DES FICHIERS
# ===========================================================================
k_dir = os.path.join(in_dir, "kraken2_addon")
b_dir = os.path.join(in_dir, "bracken_addon")

kreports   = sorted(glob.glob(os.path.join(k_dir, "*.kreport")))
b_kreports = sorted(glob.glob(os.path.join(b_dir, "*_bracken_species.kreport")))

samples = sorted([os.path.basename(f).replace(".kreport", "") for f in kreports])

print(f"[INFO] Samples détectés ({len(samples)}) : {', '.join(samples)}")

if not samples:
    raise SystemExit("[ERREUR] Aucun .kreport trouvé dans results/kraken2/")

total_reads = {}
for s in samples:
    kp = os.path.join(k_dir, f"{s}.kreport")
    total_reads[s] = get_total_reads_from_kreport(kp)
    print(f"  [{s}] Total reads : {total_reads[s]:,}")

# ===========================================================================
# TABLE 1 : KRAKEN2
# ===========================================================================
print("\n[TABLE KRAKEN2] Construction depuis les .kreport bruts...")

counts_k = defaultdict(lambda: defaultdict(float))

for s in samples:
    kp = os.path.join(k_dir, f"{s}.kreport")
    if not os.path.isfile(kp):
        continue
    nodes    = parse_kreport(kp)
    terminal = collect_terminal_nodes(nodes)
    for n, val in terminal:
        feat = to_qiime_string(lineage_from_node(n))
        if FILTER_ONLY_BACTERIA and not is_bacteria(feat):
            continue
        counts_k[feat][s] += val

counts_k = apply_filters(dict(counts_k), samples)
features_k = sorted(counts_k, key=lambda f: sum(counts_k[f].values()), reverse=True)

write_qiime_table(
    features   = features_k,
    counts     = counts_k,
    samples    = samples,
    path_counts = os.path.join(out_dir, "kraken2_table_counts.tsv"),
    path_rel    = os.path.join(out_dir, "kraken2_table_relative.tsv"),
    total_reads = total_reads
)
print("  [OK] kraken2_table_counts.tsv + kraken2_table_relative.tsv")

# ===========================================================================
# TABLE 2 : BRACKEN
# ===========================================================================
print("\n[TABLE BRACKEN] Construction depuis les *_bracken_species.kreport...")

counts_b = defaultdict(lambda: defaultdict(float))

for s in samples:
    bp = os.path.join(b_dir, f"{s}_bracken_species.kreport")
    if not os.path.isfile(bp):
        continue
    nodes    = parse_kreport(bp)
    terminal = collect_terminal_nodes(nodes)
    for n, val in terminal:
        feat = to_qiime_string(lineage_from_node(n))
        if FILTER_ONLY_BACTERIA and not is_bacteria(feat):
            continue
        counts_b[feat][s] += val

counts_b = apply_filters(dict(counts_b), samples)
features_b = sorted(counts_b, key=lambda f: sum(counts_b[f].values()), reverse=True)

write_qiime_table(
    features    = features_b,
    counts      = counts_b,
    samples     = samples,
    path_counts = os.path.join(out_dir, "bracken_table_counts.tsv"),
    path_rel    = os.path.join(out_dir, "bracken_table_relative.tsv"),
    total_reads = total_reads
)
print("  [OK] bracken_table_counts.tsv + bracken_table_relative.tsv")

print("\n=====================================================")
print(f"  Samples traités     : {len(samples)}")
print(f"  Features Kraken2    : {len(features_k)}")
print(f"  Features Bracken    : {len(features_b)}")
print("=====================================================")

PY