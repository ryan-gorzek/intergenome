#!/bin/bash
set -euo pipefail

module load java/jdk-17.0.12

export NXF_OPTS='-Xms256m -Xmx1g -Xss256k -XX:+UseSerialGC -XX:ActiveProcessorCount=2'
export NXF_DEFAULT_CPUS=2
export NXF_DEFAULT_MEMORY='2 GB'

nextflow -version
