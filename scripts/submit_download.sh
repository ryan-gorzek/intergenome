#!/bin/bash
#$ -cwd
#$ -o jobs/joblog.$JOB_ID
#$ -j y
#$ -l h_data=8G,h_rt=24:00:00
#$ -N intergenome_download
#$ -M rgorzek@ucla.edu
#$ -m bea

# UCLA Hoffman2 Data Transfer Job Script
# Run this on a data transfer node to download FASTQs, then use submit_nemo.sh to process

echo "=========================================="
echo "Job started on $(hostname) at $(date)"
echo "Job ID: $JOB_ID"
echo "Working directory: $(pwd)"
echo "=========================================="

# Load the job environment
. /u/local/Modules/default/init/modules.sh

# Set up environment
source load_nf.sh
conda activate rnaseq

echo ""
echo "Download-only mode: will stop after FASTQ downloads complete"
echo "After this job finishes, run submit_nemo.sh with -resume to continue processing"
echo ""

# Parse command line arguments
MANIFEST=${1:-manifests/test_nemo.tsv}

echo "Manifest file: $MANIFEST"
echo ""

# Run the pipeline in download-only mode
echo "Starting Nextflow pipeline (download only)..."
echo "=========================================="
nextflow run main.nf \
    --fastq_mode download \
    --fastq_manifest "$MANIFEST" \
    --download_only true \
    -resume

echo "=========================================="
echo "Downloads complete!"
echo "To continue processing, submit: qsub scripts/submit_nemo.sh $MANIFEST"
echo "=========================================="
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo "=========================================="
