#!/bin/bash
#$ -cwd
#$ -o jobs/joblog.$JOB_ID
#$ -j y
#$ -l h_data=48G,h_rt=12:00:00
#$ -N intergenome_nemo
#$ -M rgorzek@ucla.edu
#$ -m bea

# UCLA Hoffman2 Cluster Job Script for Intergenome Pipeline
# Modify the email address above and resource requirements as needed

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

# Display environment info
echo ""
echo "Environment Information:"
echo "----------------------------------------"
echo "Nextflow version:"
nextflow -version
echo ""
echo "STAR version:"
STAR --version
echo ""
echo "Python version:"
python --version
echo ""
echo "Aspera version:"
ascp --version | head -2
echo "----------------------------------------"
echo ""

# Parse command line arguments
MANIFEST=${1:-manifests/test_nemo.tsv}
STAR_INDEX=${2:-ref/ASM229v1/STARindex}

echo "Manifest file: $MANIFEST"
echo "STAR index: $STAR_INDEX"
echo ""

# Run the pipeline
echo "Starting Nextflow pipeline..."
echo "=========================================="
nextflow run main.nf \
    --fastq_mode download \
    --fastq_manifest "$MANIFEST" \
    --star_index "$STAR_INDEX" \
    -resume

# echo job info on joblog:
echo "=========================================="
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
echo "=========================================="
