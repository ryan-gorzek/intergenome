#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   download_fastq.sh <URL> <PROJDIR> <FASTQDIR> <FOLDER> <FILE>
# Example:
#   download_fastq.sh \
#     "https://data.nemoarchive.org/.../NW_TX0090-11_S01_L003-001.fastq.tar" \
#     "/u/scratch/r/rgorzek/intergenome" \
#     "data/fastq" \
#     "U01_Lein/Opossum" \
#     "NW_TX0090-11_S01_L003-001.fastq.tar"

URL="${1}"
PROJDIR="${2}"
FASTQDIR="${3:-data/fastq}"
FOLDER="${4}"
FILE="${5}"
DATADIR="${PROJDIR}/${FASTQDIR}/${FOLDER}"

mkdir -p "${DATADIR}"
cd "${DATADIR}"

wget "${URL}" -O "${FILE}"
tar -xf "${FILE}"
rm -rf "${FILE}"

echo "Done. FASTQs in ${FOLDER}/${FILE%%.*}"

# FOR GEO, NEEDS CORRECTION TO PULL FROM SRA:
# Usage:
# scripts/download_fastq.sh GSE299387

# ACCESSION="${1:-GSE299387}"
# DATADIR="${2:-$SCRATCH/intergenome/data/fastq/${ACCESSION}}"

# mkdir -p "${DATADIR}"
# cd "${DATADIR}"

# wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=${ACCESSION}&format=file"
# tar -xf "${DATADIR}/${ACCESSION}_RAW.tar"
