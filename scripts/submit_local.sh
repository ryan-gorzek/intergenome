#!/bin/bash
#$ -cwd
#$ -o jobs/joblog.$JOB_ID
#$ -j y
#$ -l h_data=48G,h_rt=12:00:00
#$ -N intergenome_local
#$ -M rgorzek@ucla.edu
#$ -m bea

# UCLA Hoffman2 Cluster Job Script for Intergenome Pipeline (Local FASTQ mode)

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

# Load STAR module AFTER conda activate so it takes precedence in PATH
module load star/2.7.10a

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
echo "----------------------------------------"
echo ""

# Parse command line arguments
FASTQ_DIR=${1:-data/fastq}
STAR_INDEX=${2:-ref/ASM229v1/STARindex}

echo "FASTQ directory: $FASTQ_DIR"
echo "STAR index: $STAR_INDEX"
echo ""

# Run the pipeline
echo "Starting Nextflow pipeline..."
echo "=========================================="
nextflow run main.nf \
    --fastq_mode local \
    --fastq_dir "$FASTQ_DIR" \
    --star_index "$STAR_INDEX" \
    -resume

echo "=========================================="
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo "=========================================="
