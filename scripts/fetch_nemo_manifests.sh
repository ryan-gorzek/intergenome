#!/usr/bin/env bash
#
# fetch_nemo_manifests.sh
#
# Download Cross-Species-M1 manifest TSV files from NeMO archive for specified species.
#
# Usage:
#   ./fetch_nemo_manifests.sh [SPECIES...] [OPTIONS]
#
# Arguments:
#   SPECIES  - One or more species names (default: human opossum)
#              Available: human, opossum, macaque, rat, chimpanzee, gorilla, baboon,
#              owl_monkey, squirrel_monkey, pig, rabbit, ferret, cat, coyote, etc.
#
# Options:
#   -o, --outdir DIR  - Output directory (default: manifests/nemo)
#   -l, --list        - List available species and exit
#   -h, --help        - Show this help message
#
# Output:
#   manifests/nemo/<species>_M1_manifest.tsv for each species
#
# Example:
#   ./fetch_nemo_manifests.sh human opossum rat
#   ./fetch_nemo_manifests.sh --list

set -euo pipefail

# NeMO archive base URL for U01_Lein 10x v3 data
NEMO_BASE="https://data.nemoarchive.org/biccn/grant/u01_lein/lein/transcriptome/sncell/10x_v3"

# Available species in the archive
AVAILABLE_SPECIES=(
  "arctic_ground_squirrel"
  "armadillo"
  "baboon"
  "chimpanzee"
  "common_tree_shrew"
  "coyote"
  "domestic_cat"
  "domestic_ferret"
  "green_monkey"
  "human"
  "macaque"
  "opossum"
  "owl_monkey"
  "pig"
  "pig_tailed_macaque"
  "rabbit"
  "rat"
  "rhesus_macaque"
  "small-eared_galago"
  "squirrel_monkey"
  "western_gorilla"
)

# Defaults
OUTDIR="manifests/nemo"
SPECIES_LIST=()

show_help() {
  head -35 "$0" | tail -32 | sed 's/^# \?//'
}

list_species() {
  echo "Available species in NeMO U01_Lein 10x_v3 archive:"
  printf '  %s\n' "${AVAILABLE_SPECIES[@]}"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--outdir)
      OUTDIR="$2"
      shift 2
      ;;
    -l|--list)
      list_species
      exit 0
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    -*)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
    *)
      SPECIES_LIST+=("$1")
      shift
      ;;
  esac
done

# Default species if none specified
if [[ ${#SPECIES_LIST[@]} -eq 0 ]]; then
  SPECIES_LIST=("human" "opossum")
fi

# Validate species
for sp in "${SPECIES_LIST[@]}"; do
  found=0
  for avail in "${AVAILABLE_SPECIES[@]}"; do
    if [[ "$sp" == "$avail" ]]; then
      found=1
      break
    fi
  done
  if [[ $found -eq 0 ]]; then
    echo "[ERROR] Unknown species: $sp" >&2
    echo "Run with --list to see available species" >&2
    exit 1
  fi
done

# Create output directory
mkdir -p "${OUTDIR}"

echo "[INFO] Fetching NeMO manifests for: ${SPECIES_LIST[*]}"
echo "[INFO] Output directory: ${OUTDIR}"

# Fetch manifest for each species
for species in "${SPECIES_LIST[@]}"; do
  echo ""
  echo "[INFO] Processing species: ${species}"

  raw_url="${NEMO_BASE}/${species}/raw/"

  # Get directory listing and find manifest files
  echo "  Fetching directory listing from: ${raw_url}"

  # Find Cross-Species-M1 manifest TSV files
  manifest_files=$(curl -fsSL "${raw_url}" 2>/dev/null \
    | grep -oE 'href="[^"]*Cross-Species-M1[^"]*manifest[^"]*\.tsv"' \
    | sed 's/href="//;s/"$//' \
    | sort -u) || true

  if [[ -z "$manifest_files" ]]; then
    echo "  [WARN] No Cross-Species-M1 manifest found for ${species}"
    continue
  fi

  # Download each manifest
  for manifest in $manifest_files; do
    manifest_url="${raw_url}${manifest}"
    outfile="${OUTDIR}/${species}_M1_manifest.tsv"

    echo "  Downloading: ${manifest}"
    curl -fsSL "${manifest_url}" -o "${outfile}"

    # Count samples
    n_samples=$(tail -n +2 "${outfile}" | wc -l | tr -d ' ')
    echo "  Saved: ${outfile} (${n_samples} samples)"
  done
done

echo ""
echo "[INFO] Done. Manifests saved to: ${OUTDIR}"
echo ""
echo "Next step: Run build_fastq_manifest.py to generate pipeline-ready manifest"
