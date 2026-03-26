#!/bin/bash
#SBATCH --job-name=install_matam_db
#SBATCH --output=log/INSTALL_MATAM_DB_%j.out
#SBATCH --error=log/INSTALL_MATAM_DB_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=debug
#SBATCH --time=02:00:00
#SBATCH --chdir=/home/thboutet/projets/16_tenebrion

#######################################################################################
#   AUTHORS :   Thomas BOUTET (IHP,LMGE)
#   DATE    :   25 - 03 - 2026
#   SCRIPT  :   pull_MATAM_sif.slurm
#######################################################################################

# sbatch /bin/pull_MATAM_sif.slurm

######################################################################################

set -euo pipefail

# Chargement des modules
module purge
module load gcc/4.8.4
module load pigz 2>/dev/null || true

# Définition de l'arborescence du projet
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
LOG_DIR="${PROJECT_DIR}/log"
CONTAINERS_DIR="${PROJECT_DIR}/containers"
RES_DIR="${PROJECT_DIR}/resources"

# Création silencieuse des dossiers nécessaires
mkdir -p "${LOG_DIR}" "${CONTAINERS_DIR}" "${RES_DIR}"

# Localisation de Singularity. Essentiel car les chemins diffèrent selon les clusters.
SINGULARITY="/usr/bin/singularity"
if [ ! -x "${SINGULARITY}" ]; then
    SINGULARITY="$(command -v singularity || true)"
fi
if [ -z "${SINGULARITY}" ] || [ ! -x "${SINGULARITY}" ]; then
    echo "ERROR: singularity/apptainer introuvable" >&2
    exit 1
fi

# Fonction de vérification des dépendances système
need_bin() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: binaire manquant: $1" >&2
        exit 1
    fi
}

# On s'assure que les outils de manipulation d'archives sont présents
need_bin awk
need_bin sed
need_bin tar
need_bin bzip2
need_bin df

# Détection de l'outil de téléchargement disponible (curl ou wget)
HAVE_CURL="FALSE"
HAVE_WGET="FALSE"
if command -v curl >/dev/null 2>&1; then
    HAVE_CURL="TRUE"
fi
if command -v wget >/dev/null 2>&1; then
    HAVE_WGET="TRUE"
fi
if [ "${HAVE_CURL}" = "FALSE" ] && [ "${HAVE_WGET}" = "FALSE" ]; then
    echo "ERROR: ni curl ni wget disponible pour télécharger la DB" >&2
    exit 1
fi

################################################################################
#   SCRATCH (TMP + CACHE) + CLEANUP robuste
################################################################################

# Définition d'un espace temporaire sur le scratch du cluster.
# Très important : Singularity a besoin de beaucoup d'espace cache pour "builder" 
# des images Docker en .sif. Si on le laisse faire dans le /home, ça sature le quota.
SCRATCH_BASE="/storage/scratch/${USER}"
SCRATCH_JOB="${SCRATCH_BASE}/install_matam_db_${SLURM_JOB_ID}"
SCRATCH_TMP="${SCRATCH_JOB}/tmp"
mkdir -p "${SCRATCH_JOB}" "${SCRATCH_TMP}"

# Redirection des variables d'environnement de Singularity vers notre scratch
export SINGULARITY_TMPDIR="${SCRATCH_TMP}"
export SINGULARITY_CACHEDIR="${SCRATCH_BASE}/singularity_cache"
mkdir -p "${SINGULARITY_CACHEDIR}" 2>/dev/null || true

# Le "trap" est une sécurité : cette fonction s'exécutera automatiquement à la fin du script,
# même s'il plante au milieu, garantissant que l'on ne laisse pas de déchets sur le scratch.
cleanup() {
    rm -rf "${SCRATCH_JOB}" 2>/dev/null || true
}
trap cleanup EXIT

################################################################################
#   LOG HEADER
################################################################################

# Traçabilité des versions et de l'environnement d'exécution
echo "START -- $(date)"
echo "HOSTNAME -- $(hostname)"
echo "PROJECT_DIR -- ${PROJECT_DIR}"
echo "SINGULARITY -- ${SINGULARITY}"
echo "SINGULARITY_VERSION -- $(${SINGULARITY} --version 2>/dev/null || true)"
echo "SINGULARITY_CACHEDIR -- ${SINGULARITY_CACHEDIR}"
echo "SINGULARITY_TMPDIR -- ${SINGULARITY_TMPDIR}"
echo "SCRATCH_JOB -- ${SCRATCH_JOB}"
echo "HOST_HAVE_CURL -- ${HAVE_CURL}"
echo "HOST_HAVE_WGET -- ${HAVE_WGET}"
df -h "${SCRATCH_BASE}" || true
echo

################################################################################
#   1) BUILD MATAM IMAGE (depuis biocontainers)
################################################################################

MATAM_VERSION="1.6.2"
MATAM_DOCKER_IMAGE="docker://quay.io/biocontainers/matam:${MATAM_VERSION}--haf24da9_0"

MATAM_SIF="${CONTAINERS_DIR}/matam.sif"
MATAM_TMP="${MATAM_SIF}.tmp"

echo "################################################################################"
echo "# 1) BUILD MATAM IMAGE"
echo "################################################################################"
echo "MATAM_VERSION -- ${MATAM_VERSION}"
echo "MATAM_DOCKER_IMAGE -- ${MATAM_DOCKER_IMAGE}"
echo "MATAM_SIF -- ${MATAM_SIF}"

# Si l'image existe déjà, on gagne du temps et on passe à la suite
if [ -s "${MATAM_SIF}" ]; then
    echo "MATAM_SIF -- déjà présent -> skip build"
    ls -lh "${MATAM_SIF}" | sed 's/^/   /'
else
    # Construction de l'image Singularity à partir de la source Docker
    rm -f "${MATAM_TMP}"
    echo "BUILD -- ${SINGULARITY} build ${MATAM_TMP} ${MATAM_DOCKER_IMAGE}"
    ${SINGULARITY} build "${MATAM_TMP}" "${MATAM_DOCKER_IMAGE}"
    mv -f "${MATAM_TMP}" "${MATAM_SIF}"
    echo "BUILD OK -- ${MATAM_SIF}"
    ls -lh "${MATAM_SIF}" | sed 's/^/   /'
fi
echo

# Vérification rapide ("sanity check") de la présence des binaires clés dans le conteneur
echo "CHECK MATAM BINARIES"
${SINGULARITY} exec --cleanenv "${MATAM_SIF}" bash -lc '
    set -e
    echo "PATH=$PATH"
    python --version 2>/dev/null || true
    for x in matam_assembly.py sga sortmerna index_default_ssu_rrna_db.py index_ref_db.py indexdb_rna; do
        echo "--- which $x"
        command -v "$x" || true
    done
'
echo

################################################################################
#   2) GET DB (download multi-sources + validate + extract)
################################################################################

MATAM_DB_DIR="${RES_DIR}/matam_db"
MATAM_DB_PREFIX="${MATAM_DB_DIR}/SILVA_128_SSURef_NR95"
MATAM_DB_SENTINEL="${MATAM_DB_PREFIX}.clustered.fasta"
MATAM_DB_ARCHIVE="${MATAM_DB_DIR}/SILVA_128_SSURef_NR95.tar.bz2"

mkdir -p "${MATAM_DB_DIR}"

echo "################################################################################"
echo "# 2) GET DB"
echo "################################################################################"
echo "MATAM_DB_DIR -- ${MATAM_DB_DIR}"
echo "MATAM_DB_SENTINEL -- ${MATAM_DB_SENTINEL}"

# Si le fichier fasta (sentinelle) est là, la DB est déjà extraite
if [ -f "${MATAM_DB_SENTINEL}" ]; then
    echo "DB -- déjà présente (sentinel) -> skip"
    ls -lh "${MATAM_DB_DIR}" | sed 's/^/    /'
    echo
else
    echo "DB -- absente -> download + validate + extract"

    # Liste de miroirs. Les liens bioinfo cassent souvent, on prévoit un fallback.
    URLS=(
        "http://bioinfo.univ-lille.fr/matam/SILVA_128_SSURef_NR95.tar.bz2"
        "https://bioinfo.univ-lille.fr/matam/SILVA_128_SSURef_NR95.tar.bz2"
    )

    download_ok="FALSE"
    rm -f "${MATAM_DB_ARCHIVE}" "${MATAM_DB_ARCHIVE}.part" 2>/dev/null || true

    # Tentative de téléchargement en bouclant sur les URL disponibles
    for url in "${URLS[@]}"; do
        echo "TRY -- ${url}"
        if [ "${HAVE_CURL}" = "TRUE" ]; then
            if curl -L --fail --retry 3 --retry-delay 5 -o "${MATAM_DB_ARCHIVE}.part" "${url}" >/dev/null 2>&1; then
                mv -f "${MATAM_DB_ARCHIVE}.part" "${MATAM_DB_ARCHIVE}"
                download_ok="TRUE"
                break
            fi
        fi
        if [ "${HAVE_WGET}" = "TRUE" ]; then
            if wget --tries=3 --waitretry=5 -O "${MATAM_DB_ARCHIVE}.part" "${url}" >/dev/null 2>&1; then
                mv -f "${MATAM_DB_ARCHIVE}.part" "${MATAM_DB_ARCHIVE}"
                download_ok="TRUE"
                break
            fi
        fi
    done

    # Si tous les miroirs échouent, on arrête proprement
    if [ "${download_ok}" != "TRUE" ]; then
        echo "ERROR: impossible de télécharger la DB MATAM (proxy/URL)." >&2
        echo "Astuce: télécharge une fois depuis un noeud qui a internet, puis relance." >&2
        exit 1
    fi

    echo "DOWNLOADED -- ${MATAM_DB_ARCHIVE}"
    ls -lh "${MATAM_DB_ARCHIVE}" | sed 's/^/    /'

    # Vérification de l'intégrité de l'archive (évite de bloquer sur une archive corrompue)
    echo "VALIDATE -- bzip2 -t"
    if ! bzip2 -t "${MATAM_DB_ARCHIVE}"; then
        echo "ERROR: archive bz2 invalide (download corrompu/proxy HTML)" >&2
        head -c 200 "${MATAM_DB_ARCHIVE}" | sed 's/^/   /'; echo >&2
        exit 1
    fi

    # Extraction
    echo "EXTRACT -- tar -xjf ${MATAM_DB_ARCHIVE}"
    tar -xjf "${MATAM_DB_ARCHIVE}" -C "${MATAM_DB_DIR}"

    # Vérification que le fichier cible est bien sorti de l'archive
    if [ ! -f "${MATAM_DB_SENTINEL}" ]; then
        echo "ERROR: extraction OK mais sentinel introuvable: ${MATAM_DB_SENTINEL}" >&2
        echo "Contenu MATAM_DB_DIR:" >&2
        ls -lh "${MATAM_DB_DIR}" >&2 || true
        exit 1
    fi

    echo "DB OK -- ${MATAM_DB_SENTINEL}"
    ls -lh "${MATAM_DB_DIR}" | sed 's/^/    /'
    echo
fi

################################################################################
#   3) INDEX DB (nécessite indexdb_rna) - multi-stratégies
################################################################################

echo "################################################################################"
echo "# 3) INDEX DB"
echo "################################################################################"

# Fonction pour vérifier si l'indexation a déjà été faite (présence des fichiers générés par SortMeRNA)
index_present() {
    ls "${MATAM_DB_PREFIX}".* 2>/dev/null | grep -Ev '\.(fasta|tab|bz2)$' >/dev/null 2>&1
}

if index_present; then
    echo "INDEX -- déjà présent -> skip"
    ls -lh "${MATAM_DB_DIR}" | sed 's/^/    /'
    echo
else
    echo "INDEX -- absent -> tentative d'indexation"
    echo "METHOD_A (MATAM image) -- indexdb_rna"
    
    # Stratégie A : On tente d'utiliser l'indexeur directement inclus dans l'image MATAM
    if ${SINGULARITY} exec --cleanenv "${MATAM_SIF}" bash -lc 'command -v indexdb_rna >/dev/null 2>&1'; then
        echo "indexdb_rna trouvé dans MATAM image -> indexation"
        ${SINGULARITY} exec --cleanenv \
            -B "${MATAM_DB_DIR}:${MATAM_DB_DIR}" \
            "${MATAM_SIF}" bash -lc "
                set -e
                indexdb_rna --ref '${MATAM_DB_SENTINEL},${MATAM_DB_PREFIX}.clustered' -v
            "
    else
        echo "METHOD_A -- indexdb_rna=NO"
    fi

    if index_present; then
        echo "INDEX OK (METHOD_A)"
        ls -lh "${MATAM_DB_DIR}" | sed 's/^/    /'
        echo
    else
        # Stratégie B : Souvent, les conteneurs MATAM n'incluent pas le binaire 'indexdb_rna' de SortMeRNA.
        # On contourne le problème en téléchargeant dynamiquement un conteneur SortMeRNA dédié pour faire l'indexation.
        echo "METHOD_B -- essayer une image sortmerna qui fournit indexdb_rna"
        
        # Liste de tags connus de SortMeRNA v2.1b (version requise par MATAM)
        SORTMERNA_TRIES=(
            "docker://quay.io/biocontainers/sortmerna:2.1b--0"
            "docker://quay.io/biocontainers/sortmerna:2.1b--h9f5acd7_0"
            "docker://quay.io/biocontainers/sortmerna:2.1b--h9f5acd7_1"
            "docker://quay.io/biocontainers/sortmerna:2.1b--h7d875b9_0"
        )

        SORT_SIF="${CONTAINERS_DIR}/sortmerna_indexer.sif"
        rm -f "${SORT_SIF}.tmp" 2>/dev/null || true

        # Boucle d'essai pour trouver un conteneur fonctionnel
        got_indexer="FALSE"
        for img in "${SORTMERNA_TRIES[@]}"; do
            echo "TRY SORTMERNA IMAGE -- ${img}"
            if ${SINGULARITY} build "${SORT_SIF}.tmp" "${img}" >/dev/null 2>&1; then
                mv -f "${SORT_SIF}.tmp" "${SORT_SIF}"
                if ${SINGULARITY} exec --cleanenv "${SORT_SIF}" bash -lc 'command -v indexdb_rna >/dev/null 2>&1'; then
                    got_indexer="TRUE"
                    echo "OK -- indexdb_rna présent dans ${SORT_SIF}"
                    break
                else
                    echo "NO -- pas de indexdb_rna dans cette image"
                    rm -f "${SORT_SIF}" 2>/dev/null || true
                fi
            else
                echo "FAIL -- build impossible (tag inexistant/proxy)"
                rm -f "${SORT_SIF}.tmp" 2>/dev/null || true
            fi
        done

        if [ "${got_indexer}" != "TRUE" ]; then
            echo "ERROR: impossible de trouver une image sortmerna fournissant indexdb_rna." >&2
            echo "SOLUTION: soit trouver le bon tag biocontainers sortmerna 2.1b qui inclut indexdb_rna," >&2
            echo "         soit on construit un conteneur indexer depuis source (sans conda)." >&2
            exit 1
        fi

        # Indexation finale avec le conteneur de secours
        echo "INDEXATION via sortmerna_indexer.sif"
        ${SINGULARITY} exec --cleanenv \
            -B "${MATAM_DB_DIR}:${MATAM_DB_DIR}" \
            "${SORT_SIF}" bash -lc "
                set -e
                indexdb_rna --ref '${MATAM_DB_SENTINEL},${MATAM_DB_PREFIX}.clustered' -v
            "

        if ! index_present; then
            echo "ERROR: indexdb_rna exécuté mais index non détecté. Vérifie les fichiers générés." >&2
            ls -lh "${MATAM_DB_DIR}" >&2 || true
            exit 1
        fi

        echo "INDEX OK (METHOD_B)"
        ls -lh "${MATAM_DB_DIR}" | sed 's/^/    /'
        echo
    fi
fi

echo "INSTALL+DB DONE -- $(date)"
echo "END -- $(date)"