#!/bin/bash
set -euo pipefail

module load java/jdk-17.0.12
module load anaconda3

export NXF_OPTS='-Xms256m -Xmx1g -Xss256k -XX:+UseSerialGC -XX:ActiveProcessorCount=8'
export NXF_DEFAULT_CPUS=2
export NXF_DEFAULT_MEMORY='10 GB'
mkdir -p "$SCRATCH/.nxf-conda"
export NXF_CONDA_CACHEDIR="$SCRATCH/.nxf-conda"

nextflow -version
