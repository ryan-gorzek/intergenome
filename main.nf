nextflow.enable.dsl=2

def makeLocalFASTQChannel(String dir) {
  Channel
    .fromFilePairs("${params.fastq_dir}/*/*_R{1,2}_*.fastq.gz", flat: true)
    .map { output ->
      def sample = output[0]
      def r1 = output.find { it.toString().contains('_R1_') }
      def r2 = output.find { it.toString().contains('_R2_') }
      tuple(sample as String, r1, r2)
    }
}

workflow {
  // 1) FASTQs
  if (params.fastq_mode == 'local') {
    fastqs = makeLocalFASTQChannel(params.fastq_dir)
  } else if (params.fastq_mode == 'download') {
    fastq_inputs = Channel
      .fromPath(params.fastq_manifest)
      .splitCsv(header: false, sep: '\t')
      .map { row ->
             [ (row[0] as String),
               (row[1] as String), 
               (row[0].tokenize('/').last().replaceFirst(/\.fastq\.tar$/, '')) ] 
            } // url, folder, sample
    fastqs = DOWNLOAD_FASTQ(fastq_inputs)
  }

  // 2) Reference (FASTA + GTF)
  ref = PREP_REF()

  // 3) 10x Barcode Inclusion List
  Channel
    .fromPath(params.bc_inclist, checkIfExists: true)
    .set { inclist }
}

// Processes //

process DOWNLOAD_FASTQ {
  tag "${folder}"
  // publishDir "${params.fastq_dir}", mode: 'copy', overwrite: true

  input:
    tuple val(url), val(folder), val(sample)
  output:
    tuple val(sample), path("*_R1.fastq.gz"), path("*_R2.fastq.gz")

  script:

  """
  set -euo pipefail

  wget -O download.fastq.tar "${url}"
  tmpdir=\$(mktemp -d)
  tar -C "\$tmpdir" -xf download.fastq.tar

  r1=\$(find "\$tmpdir" -type f -name "*_R1_*.fastq.gz" | head -n1)
  r2=\$(find "\$tmpdir" -type f -name "*_R2_*.fastq.gz" | head -n1)

  ln -s "\$r1" "${sample}_R1.fastq.gz"
  ln -s "\$r2" "${sample}_R2.fastq.gz"
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