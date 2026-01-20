# intergenome

Nextflow DSL2 pipeline for STARsolo alignment and intergenic read analysis.

## Rationale

Intergenic reads are an information source for unannotated transcription, enhancer RNA, and species-specific genome organization. This pipeline produces distance distributions of intergenic reads relative to annotated genes, stratified by gene biotype (protein-coding, lncRNA) and 3'UTR status, enabling cross-species and cross-sample comparisons.

## Features

- **Reference prep**: FASTA + GTF fetch/normalize from Ensembl; STAR index build or use existing
- **FASTQ ingest**: Local files or NeMO archive manifest-driven downloads (21 species available)
- **Alignment**: STARsolo with 10x v3 CB/UMI handling; coordinate-sorted BAM output
- **AnnData export**: STARsolo matrices converted to `.h5ad` format for scanpy
- **Intergenic filtering**: Read 3' endpoints not overlapping gene bodies; exported to BED
- **Distance metrics**: Signed distances to nearest gene 5'/3' endpoints; Parquet output
- **Histograms**: Binned distance distributions stratified by gene class
- **QC summaries**: Combined histograms and JSON summary across all samples

## Quick start

```bash
# 1) Environment
module load java nextflow
conda env create -f envs/rnaseq.yml
conda activate rnaseq

# 2) Run with local FASTQs
nextflow run main.nf \
  --fastq_mode local \
  --fastq_dir /path/to/fastqs \
  --species homo_sapiens \
  --build GRCh38

# 3) Run with NeMO manifest
# First, fetch manifests for desired species:
./scripts/fetch_nemo_manifests.sh human opossum

# Convert to pipeline format:
python scripts/build_fastq_manifest.py manifests/nemo/*.tsv -o manifests/fastq_manifest.tsv

# Run pipeline:
nextflow run main.nf \
  --fastq_mode download \
  --fastq_manifest manifests/fastq_manifest.tsv \
  --species homo_sapiens \
  --build GRCh38
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--species` | `monodelphis_domestica` | Ensembl species name |
| `--build` | `ASM229v1` | Genome assembly name |
| `--fastq_mode` | `local` | `local` or `download` |
| `--fastq_dir` | `data/fastq/test` | Directory with FASTQs |
| `--fastq_manifest` | - | TSV manifest for downloads |
| `--star_index` | - | Pre-built STAR index (skip build) |
| `--min_mapq` | `30` | Minimum mapping quality |
| `--bin_width` | `500` | Histogram bin width (bp) |
| `--min_distance` | `-10000` | Min signed distance (bp) |
| `--max_distance` | `10000` | Max signed distance (bp) |

## Inputs

- **FASTQs**: 10x v3 chemistry (R1=16bp barcode + 12bp UMI, R2=read sequence)
- **Reference**: Genome FASTA + GTF (auto-fetched from Ensembl or provided)
- **Barcode whitelist**: 10x v3 inclusion list (default: 3M-february-2018.txt)

## Outputs

```
results/
├── star/{sample}/                 # STARsolo outputs
│   ├── Aligned.sortedByCoord.out.bam
│   └── Solo.out/
├── h5ad/{sample}.h5ad             # AnnData for scanpy
├── gene_beds/                     # Gene annotation BEDs
│   ├── genes.bed                  # Gene bodies with biotype|3UTR labels
│   ├── genes_5p.bed               # TSS endpoints
│   └── genes_3p.bed               # TES endpoints
├── intergenic/
│   ├── {sample}.filtered.bam      # MAPQ-filtered BAM
│   ├── {sample}.3p.bed            # Read 3' endpoints
│   └── {sample}.intergenic.bed.gz # Intergenic reads only
├── distances/
│   └── {sample}_nearest.parquet   # Per-read distances
├── histograms/
│   └── {sample}_distance_hist.tsv # Binned counts by gene class
└── qc/
    ├── combined_histograms.tsv    # All samples merged
    └── summary.json               # Cross-sample statistics
```

## Signed distance convention

For each intergenic read, we find the nearest gene endpoint (5' TSS or 3' TES):

- **Negative distance**: Read is closer to the gene's 5' end (TSS)
- **Positive distance**: Read is closer to the gene's 3' end (TES)
- **Zero/overlap**: Excluded by intergenic filter

## Gene classes

Distances are stratified by gene biotype and 3'UTR annotation status:

- `protein_coding_3UTR`: Protein-coding genes with annotated 3'UTR
- `protein_coding_no3UTR`: Protein-coding genes without annotated 3'UTR
- `lncRNA_3UTR`: lncRNA genes with annotated 3'UTR
- `lncRNA_no3UTR`: lncRNA genes without annotated 3'UTR

## NeMO archive data

One of the largest and most intriguing ways to study intergenic reads distributions is the Allen Brain Institute's Cross-Species M1 dataset. This pipeline supports bulk downloads from the BICCN NeMO archive. Available species (21 total):

```
human, opossum, macaque, rat, chimpanzee, gorilla, baboon,
owl_monkey, squirrel_monkey, pig, rabbit, ferret, cat, coyote,
arctic_ground_squirrel, armadillo, common_tree_shrew, green_monkey,
pig_tailed_macaque, rhesus_macaque, small-eared_galago
```

To download NeMO manifests for this data:
```bash
./scripts/fetch_nemo_manifests.sh --list           # See available species
./scripts/fetch_nemo_manifests.sh human opossum    # Fetch specific species
```

## License

MIT. See `LICENSE`.

