nextflow.enable.dsl=2

workflow {
  // 1) FASTQs
  fastq_inputs = Channel
    .fromPath(params.fastq_manifest)
    .splitCsv(header: false, sep: '\t')
    .map { row -> [ (row[0] as String), (row[1] as String), (row[2] as String) ] } // [url, folder, file]

  fastqs = COPY_FASTQ(fastq_inputs)

  // 2) Reference (FASTA + GTF)
  ref = PREP_REF()

  // Expose to later workflows
  emit:
    fastqs        // Path to each downloaded/extracted folder of FASTQs
    fasta = ref.fasta
    gtf = ref.gtf
}

// Processes //

process COPY_FASTQ { // for testing, at the moment
  tag "${folder}"

  input:
    tuple val(url), val(folder), val(file) // url and file are ignored for copy
  output:
    path "${folder}/*.fastq.gz", emit: fastqs

  script:
  """
  set -euo pipefail
  mkdir -p "\$(dirname "${folder}")"
  ln -s "${projectDir}/${params.fastq_dir}/${folder}" "${folder}"
  """
}

process DOWNLOAD_FASTQ {
  tag "${folder}"
  // publishDir "${params.fastq_dir}", mode: 'copy', overwrite: true

  input:
    tuple val(url), val(folder), val(file)
  output:
    path "${folder}/*.fastq.gz", emit: fastqs

  script:
  """
  set -euo pipefail
  mkdir -p "${folder}"
  wget -O "${file}" "${url}"
  tar -C "${folder}" --strip-components=1 -xf "${file}"
  rm -f "${file}"
  """
}

process PREP_REF {
  tag "${params.build}"
  publishDir "${params.ref_dir}/${params.build}", mode: 'copy', overwrite: true

  output:
    path "genome.fa", emit: fasta
    path "genes.gtf", emit: gtf

  script:
  """
  set -euo pipefail
  ${projectDir}/scripts/download_ensembl.sh \
    "${projectDir}" \
    "${params.ref_dir}" \
    "${params.species}" \
    "${params.build}"

  DEST="${projectDir}/${params.ref_dir}/${params.build}"

  # Normalize names to 'genome.fa' and 'genes.gtf'
  # download_ensembl.sh uses .dna.toplevel.fa and non-abinitio/non-chr GTF file
  shopt -s nullglob
  fas=( "\$DEST"/*.dna.toplevel.fa )
  gts=( "\$DEST"/*.gtf )
  shopt -u nullglob

  fa="\${fas[0]-}"
  gt="\${gts[0]-}" 

  ln -sf "\$fa" "\$PWD/genome.fa"
  ln -sf "\$gt" "\$PWD/genes.gtf"

  """
}