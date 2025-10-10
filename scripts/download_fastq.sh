#!/usr/bin/env bash
set -euo pipefail

# Usage:
# scripts/download_fastq.sh https://data.nemoarchive.org/biccn/grant/u01_lein/lein/transcriptome/sncell/10x_v3/opossum/raw/NW_TX0090-11_S01_L003-001.fastq.tar U01_Lein/Opossum NW_TX0090-11.fastq.tar

DATADIR="$SCRATCH/intergenome/data/fastq"
URL="${1}"
FOLDER="${2}"
FILE="${3}"

mkdir -p "${DATADIR}/${FOLDER}"
cd "${DATADIR}/${FOLDER}"

wget "${URL}" -O "${FILE}"
tar -xf "${FILE}"

# FOR GEO, NEEDS CORRECTION TO PULL FROM SRA:
# Usage:
# scripts/download_fastq.sh GSE299387

# ACCESSION="${1:-GSE299387}"
# DATADIR="${2:-$SCRATCH/intergenome/data/fastq/${ACCESSION}}"

# mkdir -p "${DATADIR}"
# cd "${DATADIR}"

# wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=${ACCESSION}&format=file"
# tar -xf "${DATADIR}/${ACCESSION}_RAW.tar"
