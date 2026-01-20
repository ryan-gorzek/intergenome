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

  // Stop here if download_only mode
  if (params.download_only) {
    return
  }

  // 2) Reference (FASTA + GTF)
  ref = PREP_REF()

  // 3) Barcode Inclusion List
  Channel
    .fromPath(params.bc_inclist, checkIfExists: true)
    .set { inclist }

  // 4) STAR index (build or specify path)
  index = params.star_index
    ? Channel.of( file(params.star_index) )
    : STAR_INDEX(ref.fasta, ref.gtf)

  // 5) STARsolo align (currently optimized for 10x v3 3')
  alignment = STARSOLO_ALIGN(fastqs, index, inclist)

  // 6) Build gene endpoint BEDs from GTF (once per reference)
  gene_beds = BUILD_GENE_BEDS(ref.gtf)

  // 7) Convert STARsolo output to AnnData (h5ad)
  h5ad = CONVERT_TO_H5AD(alignment.solo_dir)

  // 8) Extract intergenic read 3' endpoints from BAM
  intergenic = EXTRACT_INTERGENIC(
    alignment.bam_with_sample,
    gene_beds.bodies
  )

  // 9) Calculate signed distances to nearest gene endpoints
  distances = CALC_DISTANCES(
    intergenic.bed,
    gene_beds.fivep.first(),
    gene_beds.threep.first()
  )

  // 10) Generate binned histograms by gene class
  histograms = CALC_HISTOGRAMS(distances)

  // 11) Aggregate QC across all samples
  AGGREGATE_QC(histograms.collect())

}

// Processes //

process DOWNLOAD_FASTQ {
  tag "${folder}"
  publishDir "${params.out_dir}/fastq/${folder}", mode: 'copy', overwrite: true

  input:
    tuple val(url), val(folder), val(sample)
  output:
    tuple val(sample), path("*_R1.fastq.gz"), path("*_R2.fastq.gz")

  script:
  def aspera_path = url.replaceFirst('https://data.nemoarchive.org', '')

  if (params.download_method == 'aspera')
    """
    set -euo pipefail

    # Aspera download (much faster than wget)
    # Find the Aspera key (module or home install)
    ASPERA_KEY="\$(find /u/local/apps/aspera-connect -name 'asperaweb_id_dsa.openssh' 2>/dev/null | head -1)"
    if [ -z "\$ASPERA_KEY" ]; then
      ASPERA_KEY="\$HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh"
    fi

    ascp -QT -l 500m -P 33001 -k 2 \\
      -i "\$ASPERA_KEY" \\
      asp-nemo@data.nemoarchive.org:${aspera_path} \\
      download.fastq.tar

    tmpdir=\$(mktemp -d)
    tar -C "\$tmpdir" -xf download.fastq.tar

    r1=\$(find "\$tmpdir" -type f -name "*_R1_*.fastq.gz" | head -n1)
    r2=\$(find "\$tmpdir" -type f -name "*_R2_*.fastq.gz" | head -n1)

    ln -s "\$r1" "${sample}_R1.fastq.gz"
    ln -s "\$r2" "${sample}_R2.fastq.gz"
    """
  else
    """
    set -euo pipefail

    wget -c -O download.fastq.tar "${url}"
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

process STAR_INDEX {
  tag "${params.build}"
  publishDir "${params.ref_dir}/${params.build}", mode: 'copy', overwrite: true

  input:
    path fasta
    path gtf

  output:
    path "STARindex"

  script:
  """
  set -euo pipefail
  mkdir -p STARindex
  STAR \\
    --runThreadN ${task.cpus} \\
    --runMode genomeGenerate \\
    --genomeDir STARindex \\
    --genomeFastaFiles ${fasta} \\
    --sjdbGTFfile ${gtf} \\
    --sjdbOverhang \$(( ${params.read_length} - 1 ))
  """
}

process STARSOLO_ALIGN {
  tag "$sample"
  publishDir "${params.out_dir}/star", mode: 'copy', overwrite: true

  input:
    tuple val(sample), path(r1), path(r2) // FASTQ information
    path index                            // STAR index
    path inclist                          // Barcode inclusion list (.txt)

  output:
    tuple val(sample), path("${sample}/Aligned.sortedByCoord.out.bam"), emit: bam_with_sample
    path "${sample}/Aligned.sortedByCoord.out.bam", emit: bam
    path "${sample}/Log.final.out", optional: true
    tuple val(sample), path("${sample}/Solo.out"), emit: solo_dir
    path "${sample}/Solo.out/**", emit: solo

  script:
  """
  set -euo pipefail
  STAR \\
    --runThreadN ${task.cpus} \\
    --genomeDir ${index} \\
    --readFilesIn ${r2} ${r1} \\
    --readFilesCommand zcat \\
    --outFileNamePrefix ${sample}/ \\
    --outSAMtype BAM SortedByCoordinate \\
    --soloType CB_UMI_Simple \\
    --soloCBstart 1 --soloCBlen 16 \\
    --soloUMIstart 17 --soloUMIlen 12 \\
    --soloCBwhitelist ${inclist} \\
    --soloBarcodeReadLength 0 \\
    --soloFeatures GeneFull \\
    --soloMultiMappers EM \\
    --soloStrand Forward \\
    --clipAdapterType CellRanger4
  """
}

process BUILD_GENE_BEDS {
  tag "${params.build}"
  publishDir "${params.out_dir}/gene_beds", mode: 'copy', overwrite: true

  input:
    path gtf

  output:
    path "genes.bed", emit: bodies
    path "genes_5p.bed", emit: fivep
    path "genes_3p.bed", emit: threep

  script:
  // AWK script for GTF parsing - creates BED files with gene_id|biotype|3UTR_status labels
  """
  set -euo pipefail

  # Two-pass AWK: first collect biotype and 3'UTR info, then emit BED records
  awk '
    BEGIN { FS = OFS = "\\t" }

    function trim(s) {
      gsub(/^[ \\t]+|[ \\t]+\$/, "", s)
      return s
    }

    function get_attr(attr_str, key,    n, i, f) {
      n = split(attr_str, f, ";")
      for (i = 1; i <= n; i++) {
        f[i] = trim(f[i])
        if (f[i] ~ ("^" key "[ \\t]")) {
          sub(("^" key "[ \\t]+"), "", f[i])
          gsub(/"/, "", f[i])
          return f[i]
        }
      }
      return ""
    }

    function norm_biotype(bt) {
      if (bt == "protein_coding") return "protein_coding"
      if (bt == "lncRNA" || bt == "lincRNA") return "lncRNA"
      return "other"
    }

    function label(gid,    bt, utr) {
      bt = (gid in biotype ? biotype[gid] : "other")
      utr = (gid in has3utr ? "3UTR" : "no3UTR")
      return gid "|" bt "|" utr
    }

    # Pass 1: Collect biotype and 3pUTR information
    FNR == NR {
      if (\$3 == "gene") {
        gid = get_attr(\$9, "gene_id")
        bt = get_attr(\$9, "gene_biotype")
        if (bt == "") bt = get_attr(\$9, "gene_type")
        bt = norm_biotype(bt)
        if (gid != "") biotype[gid] = bt
      } else if (\$3 == "three_prime_utr") {
        gid = get_attr(\$9, "gene_id")
        if (gid != "") has3utr[gid] = 1
      }
      next
    }

    # Pass 2: Emit BED records
    \$3 == "gene" {
      gid = get_attr(\$9, "gene_id")
      if (gid == "") gid = "NA"
      nm = label(gid)

      if (mode == "GENE") {
        print \$1, \$4 - 1, \$5, nm, 0, \$7
      } else if (mode == "5P") {
        if (\$7 == "+") { s = \$4 - 1; e = \$4 } else { s = \$5 - 1; e = \$5 }
        print \$1, s, e, nm, 0, \$7
      } else if (mode == "3P") {
        if (\$7 == "+") { s = \$5 - 1; e = \$5 } else { s = \$4 - 1; e = \$4 }
        print \$1, s, e, nm, 0, \$7
      }
    }
  ' mode="GENE" "${gtf}" "${gtf}" > genes.bed

  awk '
    BEGIN { FS = OFS = "\\t" }
    function trim(s) { gsub(/^[ \\t]+|[ \\t]+\$/, "", s); return s }
    function get_attr(attr_str, key,    n, i, f) {
      n = split(attr_str, f, ";")
      for (i = 1; i <= n; i++) {
        f[i] = trim(f[i])
        if (f[i] ~ ("^" key "[ \\t]")) {
          sub(("^" key "[ \\t]+"), "", f[i]); gsub(/"/, "", f[i]); return f[i]
        }
      }
      return ""
    }
    function norm_biotype(bt) {
      if (bt == "protein_coding") return "protein_coding"
      if (bt == "lncRNA" || bt == "lincRNA") return "lncRNA"
      return "other"
    }
    function label(gid) {
      bt = (gid in biotype ? biotype[gid] : "other")
      utr = (gid in has3utr ? "3UTR" : "no3UTR")
      return gid "|" bt "|" utr
    }
    FNR == NR {
      if (\$3 == "gene") { gid = get_attr(\$9, "gene_id"); bt = get_attr(\$9, "gene_biotype"); if (bt == "") bt = get_attr(\$9, "gene_type"); biotype[gid] = norm_biotype(bt) }
      else if (\$3 == "three_prime_utr") { has3utr[get_attr(\$9, "gene_id")] = 1 }
      next
    }
    \$3 == "gene" {
      gid = get_attr(\$9, "gene_id"); nm = label(gid)
      if (\$7 == "+") { s = \$4 - 1; e = \$4 } else { s = \$5 - 1; e = \$5 }
      print \$1, s, e, nm, 0, \$7
    }
  ' "${gtf}" "${gtf}" > genes_5p.bed

  awk '
    BEGIN { FS = OFS = "\\t" }
    function trim(s) { gsub(/^[ \\t]+|[ \\t]+\$/, "", s); return s }
    function get_attr(attr_str, key,    n, i, f) {
      n = split(attr_str, f, ";")
      for (i = 1; i <= n; i++) {
        f[i] = trim(f[i])
        if (f[i] ~ ("^" key "[ \\t]")) {
          sub(("^" key "[ \\t]+"), "", f[i]); gsub(/"/, "", f[i]); return f[i]
        }
      }
      return ""
    }
    function norm_biotype(bt) {
      if (bt == "protein_coding") return "protein_coding"
      if (bt == "lncRNA" || bt == "lincRNA") return "lncRNA"
      return "other"
    }
    function label(gid) {
      bt = (gid in biotype ? biotype[gid] : "other")
      utr = (gid in has3utr ? "3UTR" : "no3UTR")
      return gid "|" bt "|" utr
    }
    FNR == NR {
      if (\$3 == "gene") { gid = get_attr(\$9, "gene_id"); bt = get_attr(\$9, "gene_biotype"); if (bt == "") bt = get_attr(\$9, "gene_type"); biotype[gid] = norm_biotype(bt) }
      else if (\$3 == "three_prime_utr") { has3utr[get_attr(\$9, "gene_id")] = 1 }
      next
    }
    \$3 == "gene" {
      gid = get_attr(\$9, "gene_id"); nm = label(gid)
      if (\$7 == "+") { s = \$5 - 1; e = \$5 } else { s = \$4 - 1; e = \$4 }
      print \$1, s, e, nm, 0, \$7
    }
  ' "${gtf}" "${gtf}" > genes_3p.bed

  # Sort BED files for bedtools -sorted operations
  for f in genes.bed genes_5p.bed genes_3p.bed; do
    LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n "\$f" -o "\$f"
  done
  """
}

process CONVERT_TO_H5AD {
  tag "$sample"
  publishDir "${params.out_dir}/h5ad", mode: 'copy', overwrite: true

  input:
    tuple val(sample), path(solo_dir)

  output:
    tuple val(sample), path("${sample}.h5ad"), emit: h5ad

  script:
  """
  #!/usr/bin/env python3
  import scanpy as sc
  import os

  # Find the GeneFull matrix directory
  matrix_dir = None
  for subdir in ['GeneFull', 'Gene', 'GeneFull_Ex50pAS']:
    candidate = os.path.join('${solo_dir}', subdir, 'raw')
    if os.path.exists(candidate):
      matrix_dir = candidate
      break

  if matrix_dir is None:
    # Try filtered directory
    for subdir in ['GeneFull', 'Gene']:
      candidate = os.path.join('${solo_dir}', subdir, 'filtered')
      if os.path.exists(candidate):
        matrix_dir = candidate
        break

  if matrix_dir is None:
    raise FileNotFoundError(f"No STARsolo matrix found in ${solo_dir}")

  # Read the matrix
  adata = sc.read_10x_mtx(matrix_dir, var_names='gene_ids', cache=False)

  # Add sample metadata
  adata.obs['sample'] = '${sample}'

  # Basic QC metrics
  adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
  sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

  # Write to h5ad
  adata.write('${sample}.h5ad')
  print(f"Wrote ${sample}.h5ad: {adata.n_obs} cells x {adata.n_vars} genes")
  """
}

process EXTRACT_INTERGENIC {
  tag "$sample"
  publishDir "${params.out_dir}/intergenic", mode: 'copy', overwrite: true

  input:
    tuple val(sample), path(bam)
    path gene_bodies_bed

  output:
    tuple val(sample), path("${sample}.intergenic.bed.gz"), emit: bed
    path "${sample}.filtered.bam", emit: filtered_bam
    path "${sample}.3p.bed", emit: read_3p_bed

  script:
  """
  set -euo pipefail

  # 1) Filter BAM: remove unmapped, secondary, supplementary, duplicates; MAPQ >= 30
  samtools view -F 0x904 -q ${params.min_mapq} -b "${bam}" > "${sample}.filtered.bam"

  # Count total good reads for normalization
  total_reads=\$(samtools view -c "${sample}.filtered.bam")
  echo "\${total_reads}" > "${sample}.total_reads.txt"

  # 2) Convert BAM to BED of read 3' endpoints (strand-aware)
  bedtools bamtobed -i "${sample}.filtered.bam" \\
    | awk 'BEGIN { FS = OFS = "\\t" }
      {
        chr = \$1; s = \$2; e = \$3; name = \$4; score = \$5; strand = \$6
        if (strand == "+") {
          print chr, e - 1, e, name, score, strand
        } else if (strand == "-") {
          print chr, s, s + 1, name, score, strand
        } else {
          mid = int((s + e) / 2)
          print chr, mid, mid + 1, name, score, "."
        }
      }' > "${sample}.3p.unsorted.bed"

  LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n "${sample}.3p.unsorted.bed" -o "${sample}.3p.bed"
  rm "${sample}.3p.unsorted.bed"

  # 3) Extract reads that do NOT overlap any gene body (intergenic)
  bedtools intersect -sorted -a "${sample}.3p.bed" -b "${gene_bodies_bed}" -v \\
    | gzip > "${sample}.intergenic.bed.gz"
  """
}

process CALC_DISTANCES {
  tag "$sample"
  publishDir "${params.out_dir}/distances", mode: 'copy', overwrite: true

  input:
    tuple val(sample), path(intergenic_bed)
    path genes_5p_bed
    path genes_3p_bed

  output:
    tuple val(sample), path("${sample}_nearest.parquet"), emit: parquet

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import subprocess
  import tempfile
  import os

  sample = '${sample}'
  intergenic_bed = '${intergenic_bed}'
  genes_5p = '${genes_5p_bed}'
  genes_3p = '${genes_3p_bed}'

  # Decompress intergenic bed if needed
  if intergenic_bed.endswith('.gz'):
    subprocess.run(['gunzip', '-k', intergenic_bed], check=True)
    intergenic_bed = intergenic_bed[:-3]

  # Sort the intergenic bed file
  sorted_bed = intergenic_bed + '.sorted'
  subprocess.run(f'LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n {intergenic_bed} -o {sorted_bed}', shell=True, check=True)

  # Find closest 5' and 3' gene endpoints
  with tempfile.NamedTemporaryFile(mode='w', suffix='.5p.txt', delete=False) as tmp5:
    tmp5_path = tmp5.name
  with tempfile.NamedTemporaryFile(mode='w', suffix='.3p.txt', delete=False) as tmp3:
    tmp3_path = tmp3.name

  subprocess.run(f'bedtools closest -sorted -a {sorted_bed} -b {genes_5p} -d > {tmp5_path}', shell=True, check=True)
  subprocess.run(f'bedtools closest -sorted -a {sorted_bed} -b {genes_3p} -d > {tmp3_path}', shell=True, check=True)

  # Read results
  cols_read = ['chrom', 'start', 'end', 'read_name', 'score', 'strand']
  cols_gene = ['g_chrom', 'g_start', 'g_end', 'g_label', 'g_score', 'g_strand', 'distance']

  df5 = pd.read_csv(tmp5_path, sep='\\t', header=None, names=cols_read + cols_gene)
  df3 = pd.read_csv(tmp3_path, sep='\\t', header=None, names=cols_read + cols_gene)

  # Clean up temp files
  os.unlink(tmp5_path)
  os.unlink(tmp3_path)
  os.unlink(sorted_bed)

  # Combine: pick nearer endpoint (ties go to 5')
  df5['region'] = '5p'
  df3['region'] = '3p'

  # Merge on read coordinates
  df = df5[['chrom', 'start', 'end', 'read_name', 'strand']].copy()
  df['dist_5p'] = df5['distance']
  df['dist_3p'] = df3['distance']
  df['label_5p'] = df5['g_label']
  df['label_3p'] = df3['g_label']

  # Choose nearer endpoint
  df['region'] = 'NA'
  df['distance'] = -1
  df['gene_label'] = 'NA'

  mask_5p = df['dist_5p'] <= df['dist_3p']
  df.loc[mask_5p, 'region'] = '5p'
  df.loc[mask_5p, 'distance'] = df.loc[mask_5p, 'dist_5p']
  df.loc[mask_5p, 'gene_label'] = df.loc[mask_5p, 'label_5p']

  df.loc[~mask_5p, 'region'] = '3p'
  df.loc[~mask_5p, 'distance'] = df.loc[~mask_5p, 'dist_3p']
  df.loc[~mask_5p, 'gene_label'] = df.loc[~mask_5p, 'label_3p']

  # Apply signed distance convention: negative for 5' (upstream), positive for 3' (downstream)
  df['signed_distance'] = df.apply(
    lambda row: -row['distance'] if row['region'] == '5p' else row['distance'],
    axis=1
  )

  # Parse gene label: gene_id|biotype|3UTR_status
  df[['gene_id', 'biotype', 'utr_status']] = df['gene_label'].str.split('|', expand=True)
  df['gene_class'] = df['biotype'] + '_' + df['utr_status']

  # Filter to relevant gene classes
  keep_classes = ['protein_coding_3UTR', 'protein_coding_no3UTR', 'lncRNA_3UTR', 'lncRNA_no3UTR']
  df = df[df['gene_class'].isin(keep_classes)]

  # Select output columns
  out_df = df[['chrom', 'start', 'end', 'read_name', 'strand', 'gene_id', 'gene_class', 'region', 'signed_distance']].copy()
  out_df['sample'] = sample

  # Write to parquet
  out_df.to_parquet(f'{sample}_nearest.parquet', index=False)
  print(f"Wrote {len(out_df)} distance records to {sample}_nearest.parquet")
  """
}

process CALC_HISTOGRAMS {
  tag "$sample"
  publishDir "${params.out_dir}/histograms", mode: 'copy', overwrite: true

  input:
    tuple val(sample), path(parquet)

  output:
    path "${sample}_distance_hist.tsv", emit: hist

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import numpy as np

  sample = '${sample}'
  bin_width = ${params.bin_width}
  min_dist = ${params.min_distance}
  max_dist = ${params.max_distance}

  # Read distance data
  df = pd.read_parquet('${parquet}')

  # Read total good reads for normalization (from intergenic step)
  total_reads_file = '${params.out_dir}/intergenic/${sample}.total_reads.txt'
  try:
    with open(total_reads_file) as f:
      total_good_reads = int(f.read().strip())
  except:
    total_good_reads = len(df)

  # Filter to distance range
  df = df[(df['signed_distance'] >= min_dist) & (df['signed_distance'] <= max_dist)]

  # Define bins
  bins = np.arange(min_dist, max_dist + bin_width, bin_width)

  # Gene classes and regions to track
  gene_classes = ['protein_coding_3UTR', 'protein_coding_no3UTR', 'lncRNA_3UTR', 'lncRNA_no3UTR']
  regions = ['5p', '3p']

  # Build histogram
  rows = []
  for gc in gene_classes:
    for region in regions:
      subset = df[(df['gene_class'] == gc) & (df['region'] == region)]

      for i in range(len(bins) - 1):
        bin_start = bins[i]
        bin_end = bins[i + 1]
        bin_mid = (bin_start + bin_end) / 2

        count = len(subset[(subset['signed_distance'] >= bin_start) & (subset['signed_distance'] < bin_end)])

        rows.append({
          'sample': sample,
          'gene_group': gc,
          'region': region,
          'bin_start_bp': int(bin_start),
          'bin_end_bp': int(bin_end),
          'bin_mid_bp': bin_mid,
          'count': count,
          'total_good_reads': total_good_reads
        })

  hist_df = pd.DataFrame(rows)
  hist_df.to_csv(f'{sample}_distance_hist.tsv', sep='\\t', index=False)
  print(f"Wrote histogram with {len(hist_df)} bins to {sample}_distance_hist.tsv")
  """
}

process AGGREGATE_QC {
  publishDir "${params.out_dir}/qc", mode: 'copy', overwrite: true

  input:
    path histograms

  output:
    path "combined_histograms.tsv"
    path "summary.json"

  script:
  """
  #!/usr/bin/env python3
  import pandas as pd
  import json
  import glob

  # Combine all histogram files
  hist_files = glob.glob('*_distance_hist.tsv')
  dfs = [pd.read_csv(f, sep='\\t') for f in hist_files]
  combined = pd.concat(dfs, ignore_index=True)
  combined.to_csv('combined_histograms.tsv', sep='\\t', index=False)

  # Generate summary statistics
  summary = {
    'n_samples': combined['sample'].nunique(),
    'samples': combined['sample'].unique().tolist(),
    'total_reads_across_samples': int(combined.groupby('sample')['total_good_reads'].first().sum()),
    'gene_classes': combined['gene_group'].unique().tolist(),
    'bin_width': int(combined['bin_end_bp'].iloc[0] - combined['bin_start_bp'].iloc[0]) if len(combined) > 0 else 0,
    'distance_range': [int(combined['bin_start_bp'].min()), int(combined['bin_end_bp'].max())] if len(combined) > 0 else [0, 0]
  }

  # Per-sample stats
  per_sample = {}
  for sample in combined['sample'].unique():
    sample_df = combined[combined['sample'] == sample]
    per_sample[sample] = {
      'total_intergenic_reads': int(sample_df['count'].sum()),
      'total_good_reads': int(sample_df['total_good_reads'].iloc[0]) if len(sample_df) > 0 else 0
    }
  summary['per_sample'] = per_sample

  with open('summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

  print(f"Combined {len(hist_files)} histogram files")
  print(f"Summary: {summary['n_samples']} samples, {summary['total_reads_across_samples']} total reads")
  """
}













