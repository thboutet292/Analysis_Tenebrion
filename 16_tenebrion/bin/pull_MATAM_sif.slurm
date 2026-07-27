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
#   DATE    :   25 - 02 - 2026
#   SCRIPT  :   pull_MATAM_sif.slurm

#######################################################################################

# sbatch bin/pull_MATAM_sif.slurm

######################################################################################

set -euo pipefail

# Load modules
module purge
module load gcc/4.8.4
module load pigz 2>/dev/null || true

# Definition of the project directory tree
PROJECT_DIR="${SLURM_SUBMIT_DIR}"
LOG_DIR="${PROJECT_DIR}/log"
CONTAINERS_DIR="${PROJECT_DIR}/containers"
RES_DIR="${PROJECT_DIR}/resources"

# Silently create the necessary folders
mkdir -p "${LOG_DIR}" "${CONTAINERS_DIR}" "${RES_DIR}"

# Locate Singularity. Essential because paths differ from cluster to cluster.
SINGULARITY="/usr/bin/singularity"
if [ ! -x "${SINGULARITY}" ]; then
    SINGULARITY="$(command -v singularity || true)"
fi
if [ -z "${SINGULARITY}" ] || [ ! -x "${SINGULARITY}" ]; then
    echo "ERROR: singularity/apptainer not found" >&2
    exit 1
fi

# Function to check system dependencies
need_bin() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: missing binary: $1" >&2
        exit 1
    fi
}

# Make sure the archive-handling tools are present
need_bin awk
need_bin sed
need_bin tar
need_bin bzip2
need_bin df

# Detect the available download tool (curl or wget)
HAVE_CURL="FALSE"
HAVE_WGET="FALSE"
if command -v curl >/dev/null 2>&1; then
    HAVE_CURL="TRUE"
fi
if command -v wget >/dev/null 2>&1; then
    HAVE_WGET="TRUE"
fi
if [ "${HAVE_CURL}" = "FALSE" ] && [ "${HAVE_WGET}" = "FALSE" ]; then
    echo "ERROR: neither curl nor wget available to download the DB" >&2
    exit 1
fi

# =====================================================================
#   SCRATCH (TMP + CACHE) + robust CLEANUP
# =====================================================================

# Definition of a temporary space on the cluster's scratch storage.
# Very important: Singularity needs a lot of cache space to "build"
# Docker images into .sif files. If left to do so in /home, it fills up the quota.
SCRATCH_BASE="/storage/scratch/${USER}"
SCRATCH_JOB="${SCRATCH_BASE}/install_matam_db_${SLURM_JOB_ID}"
SCRATCH_TMP="${SCRATCH_JOB}/tmp"
mkdir -p "${SCRATCH_JOB}" "${SCRATCH_TMP}"

# Redirect Singularity's environment variables to our scratch space
export SINGULARITY_TMPDIR="${SCRATCH_TMP}"
export SINGULARITY_CACHEDIR="${SCRATCH_BASE}/singularity_cache"
mkdir -p "${SINGULARITY_CACHEDIR}" 2>/dev/null || true

# The "trap" is a safeguard: this function runs automatically at the end of the script,
# even if it crashes partway through, ensuring no leftover data is left on the scratch space.
cleanup() {
    rm -rf "${SCRATCH_JOB}" 2>/dev/null || true
}
trap cleanup EXIT

# =====================================================================
#   LOG HEADER
# =====================================================================

# Traceability of versions and the execution environment
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

# =====================================================================
#  1) BUILD MATAM IMAGE (from biocontainers)
# =====================================================================

MATAM_VERSION="1.6.2"
MATAM_DOCKER_IMAGE="docker://quay.io/biocontainers/matam:${MATAM_VERSION}--haf24da9_0"

MATAM_SIF="${CONTAINERS_DIR}/matam.sif"
MATAM_TMP="${MATAM_SIF}.tmp"

echo "---------------------------------------------------------------------"
echo "# 1) BUILD MATAM IMAGE"
echo "---------------------------------------------------------------------"
echo "MATAM_VERSION -- ${MATAM_VERSION}"
echo "MATAM_DOCKER_IMAGE -- ${MATAM_DOCKER_IMAGE}"
echo "MATAM_SIF -- ${MATAM_SIF}"

# If the image already exists, save time and move on
if [ -s "${MATAM_SIF}" ]; then
    echo "MATAM_SIF -- already present -> skip build"
    ls -lh "${MATAM_SIF}" | sed 's/^/   /'
else
    # Build the Singularity image from the Docker source
    rm -f "${MATAM_TMP}"
    echo "BUILD -- ${SINGULARITY} build ${MATAM_TMP} ${MATAM_DOCKER_IMAGE}"
    ${SINGULARITY} build "${MATAM_TMP}" "${MATAM_DOCKER_IMAGE}"
    mv -f "${MATAM_TMP}" "${MATAM_SIF}"
    echo "BUILD OK -- ${MATAM_SIF}"
    ls -lh "${MATAM_SIF}" | sed 's/^/   /'
fi
echo

# Quick "sanity check" that the key binaries are present in the container
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

# =====================================================================
# 2) GET DB (multi-source download + validate + extract)
# =====================================================================


MATAM_DB_DIR="${RES_DIR}/matam_db"
MATAM_DB_PREFIX="${MATAM_DB_DIR}/SILVA_128_SSURef_NR95"
MATAM_DB_SENTINEL="${MATAM_DB_PREFIX}.clustered.fasta"
MATAM_DB_ARCHIVE="${MATAM_DB_DIR}/SILVA_128_SSURef_NR95.tar.bz2"

mkdir -p "${MATAM_DB_DIR}"

echo "---------------------------------------------------------------------"
echo "# 2) GET DB"
echo "---------------------------------------------------------------------"
echo "MATAM_DB_DIR -- ${MATAM_DB_DIR}"
echo "MATAM_DB_SENTINEL -- ${MATAM_DB_SENTINEL}"

# If the fasta (sentinel) file is present, the DB has already been extracted
if [ -f "${MATAM_DB_SENTINEL}" ]; then
    echo "DB -- déjà présente (sentinel) -> skip"
    ls -lh "${MATAM_DB_DIR}" | sed 's/^/    /'
    echo
else
    echo "DB -- missing -> download + validate + extract"

    # List of mirrors. Bioinfo links break often, so we provide a fallback.
    URLS=(
        "http://bioinfo.univ-lille.fr/matam/SILVA_128_SSURef_NR95.tar.bz2"
        "https://bioinfo.univ-lille.fr/matam/SILVA_128_SSURef_NR95.tar.bz2"
    )

    download_ok="FALSE"
    rm -f "${MATAM_DB_ARCHIVE}" "${MATAM_DB_ARCHIVE}.part" 2>/dev/null || true

    # Attempt download by looping over the available URLs
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

    # If all mirrors fail, stop cleanly
    if [ "${download_ok}" != "TRUE" ]; then
        echo "ERROR: unable to download the MATAM DB (proxy/URL)." >&2
        exit 1
    fi

    echo "DOWNLOADED -- ${MATAM_DB_ARCHIVE}"
    ls -lh "${MATAM_DB_ARCHIVE}" | sed 's/^/    /'

    # Check the archive's integrity (avoids getting stuck on a corrupted archive)
    echo "VALIDATE -- bzip2 -t"
    if ! bzip2 -t "${MATAM_DB_ARCHIVE}"; then
        echo "ERROR: invalid bz2 archive (download corrupted/HTML proxy)" >&2
        head -c 200 "${MATAM_DB_ARCHIVE}" | sed 's/^/   /'; echo >&2
        exit 1
    fi

    # Extraction
    echo "EXTRACT -- tar -xjf ${MATAM_DB_ARCHIVE}"
    tar -xjf "${MATAM_DB_ARCHIVE}" -C "${MATAM_DB_DIR}"

    # Check that the target file was indeed extracted from the archive
    if [ ! -f "${MATAM_DB_SENTINEL}" ]; then
        echo "ERROR: extraction OK but sentinel not found: ${MATAM_DB_SENTINEL}" >&2
        echo "Contenu MATAM_DB_DIR:" >&2
        ls -lh "${MATAM_DB_DIR}" >&2 || true
        exit 1
    fi

    echo "DB OK -- ${MATAM_DB_SENTINEL}"
    ls -lh "${MATAM_DB_DIR}" | sed 's/^/    /'
    echo
fi

# =====================================================================
#  3) INDEX DB (requires indexdb_rna) - multi-strategy
# =====================================================================

echo "---------------------------------------------------------------------"
echo "# 3) INDEX DB"
echo "---------------------------------------------------------------------"

# Function to check if the indexing has already been done
index_present() {
    ls "${MATAM_DB_PREFIX}".* 2>/dev/null | grep -Ev '\.(fasta|tab|bz2)$' >/dev/null 2>&1
}

if index_present; then
    echo "INDEX -- already present -> skip"
    ls -lh "${MATAM_DB_DIR}" | sed 's/^/    /'
    echo
else
    echo "INDEX -- missing -> try to index"
    echo "METHOD_A (MATAM image) -- indexdb_rna"
    
    # METHOD A : try to use the indexer directly included in the MATAM image.
    if ${SINGULARITY} exec --cleanenv "${MATAM_SIF}" bash -lc 'command -v indexdb_rna >/dev/null 2>&1'; then
        echo "indexdb_rna found in MATAM image -> indexation"
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
        # METHOD B: Often, MATAM containers do not include the SortMeRNA 'indexdb_rna' binary.
        # Bypass the problem by dynamically downloading a dedicated SortMeRNA container for indexing.
        echo "METHOD_B -- try a sortmerna image that provides indexdb_rna"
        
        # List of known tags from SortMeRNA v2.1b (version requise par MATAM)
        SORTMERNA_TRIES=(
            "docker://quay.io/biocontainers/sortmerna:2.1b--0"
            "docker://quay.io/biocontainers/sortmerna:2.1b--h9f5acd7_0"
            "docker://quay.io/biocontainers/sortmerna:2.1b--h9f5acd7_1"
            "docker://quay.io/biocontainers/sortmerna:2.1b--h7d875b9_0"
        )

        SORT_SIF="${CONTAINERS_DIR}/sortmerna_indexer.sif"
        rm -f "${SORT_SIF}.tmp" 2>/dev/null || true

        # Test loop to find a functional container
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
                    echo "NO -- no indexdb_rna in this image"
                    rm -f "${SORT_SIF}" 2>/dev/null || true
                fi
            else
                echo "FAIL -- build impossible (tag not existing/proxy)"
                rm -f "${SORT_SIF}.tmp" 2>/dev/null || true
            fi
        done

        if [ "${got_indexer}" != "TRUE" ]; then
            echo "ERROR: Unable to find an image sortmerna providing indexdb_rna ." >&2
            exit 1
        fi

        # Final indexing with the backup container
        echo "INDEXATION with sortmerna_indexer.sif"
        ${SINGULARITY} exec --cleanenv \
            -B "${MATAM_DB_DIR}:${MATAM_DB_DIR}" \
            "${SORT_SIF}" bash -lc "
                set -e
                indexdb_rna --ref '${MATAM_DB_SENTINEL},${MATAM_DB_PREFIX}.clustered' -v
            "

        if ! index_present; then
            echo "ERROR: indexdb_rna executed but index not detected." >&2
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
