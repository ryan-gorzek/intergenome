#!/bin/bash
# usage: subset_fastq.sh R1.fastq.gz R2.fastq.gz N out_prefix
set -euo pipefail

R1="$1"; R2="$2"; N="$3"; OUT="$4"
L=$((N*4))  # 4 lines per FASTQ record
set +o pipefail
# Name with _R1_ and _R2_ for Nextflow pipeline (input, not output)
zcat "$R1" | head -n "$L" | gzip -c > "${OUT}_R1_.fastq.gz"
zcat "$R2" | head -n "$L" | gzip -c > "${OUT}_R2_.fastq.gz"
set -o pipefail
echo "Wrote ${OUT}_R1_.fastq.gz and ${OUT}_R2_.fastq.gz with first $N reads."
