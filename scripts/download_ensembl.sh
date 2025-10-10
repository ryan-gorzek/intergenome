#!/usr/bin/env bash

# Usage:
# scripts/download_ensembl.sh /u/scratch/r/rgorzek/intergenome ref monodelphis_domestica ASM229v1
# scripts/download_ensembl.sh /u/scratch/r/rgorzek/intergenome ref mus_musculus GRCm39

PROJDIR="${1}"
REFDIR="${2:-ref}"
SPECIES="${3:-mus_musculus}"
BUILD="${4:-GRCm39}"
OUTDIR="${PROJDIR}/${REFDIR}/${BUILD}"

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# Get the FASTA (toplevel)
curl -fsS -L "https://ftp.ensembl.org/pub/current_fasta/${SPECIES}/dna/" \
  | grep -Eo 'href="[^"]+dna\.toplevel\.fa\.gz"' \
  | sed 's/href="//; s/"$//' \
  | grep -E "${BUILD}\.dna\.toplevel\.fa\.gz$" \
  | head -n1 \
  | while read -r fa; do
      echo " -> ${fa}"
      curl -fsSLO "https://ftp.ensembl.org/pub/current_fasta/${SPECIES}/dna/${fa}"
    done
# Get the GTF
curl -fsS -L "https://ftp.ensembl.org/pub/current_gtf/${SPECIES}/" \
  | grep -Eo 'href="[^"]+\.gtf\.gz"' \
  | grep -vE '(abinitio|chr)' \
  | sed 's/href="//; s/"$//' \
  | while read gtf; do
      echo " -> ${gtf}"
      curl -fsSLO "https://ftp.ensembl.org/pub/current_gtf/${SPECIES}/${gtf}"
    done

# Unzip main files for STAR
gunzip -f *.dna.toplevel.fa.gz
gunzip -f *.gtf.gz

echo "Done. FASTA & GTF in ${OUTDIR}"
