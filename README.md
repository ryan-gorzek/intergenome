# intergenome

Nextflow DSL2 pipeline for STARsolo alignment and intergenic-read analysis. Under active development. :contentReference[oaicite:0]{index=0}

## Rationale

Intergenic reads are an information source for unannotated transcription, enhancer RNA, and species-specific genome organization. The target outcome is a cross-species mapping of intergenic read distributions with strand and orientation context, benchmarked against current annotations and distance to nearest genes.

## Current capabilities

- Reference prep: FASTA + GTF fetch/normalize; STAR index build or use existing index.
- 10x/FASTQ ingest: local or manifest-driven downloads; optional subsampling for smoke tests.
- Alignment: STARsolo with CB/UMI handling; sorted BAM output.
- Intergenic filtering: retain mapped reads not overlapping gene features; export to BED.
- Distance metrics: nearest upstream/downstream gene with sign conventions (− upstream, + downstream).
- Summaries: per-chromosome and genome-wide histograms; TSVs for downstream plots.

Repo structure: `envs/`, `manifests/`, `scripts/`, `main.nf`, `nextflow.config`, `LICENSE`. :contentReference[oaicite:1]{index=1}

## Quick start

```bash
# 1) Environment (example)
module load java
# Nextflow installed or in PATH

# 2) Run: local FASTQs
nextflow run ryan-gorzek/intergenome \
  -profile local \
  --fastq_mode local \
  --fastq_dir /path/to/fastqs \
  --species mus_musculus \
  --build GRCm39

# 3) Run: manifest (URL<TAB>sample<TAB>folder)
nextflow run ryan-gorzek/intergenome \
  -profile local \
  --fastq_mode download \
  --fastq_manifest manifests/fastq_manifest.tsv \
  --species monodelphis_domestica \
  --build ASM229v1
```
# intergenome

Nextflow DSL2 pipeline for STARsolo alignment and intergenic-read analysis. Under active development. :contentReference[oaicite:0]{index=0}

## Rationale

Intergenic reads are an information source for unannotated transcription, enhancer RNA, and species-specific genome organization. The target outcome is a cross-species mapping of intergenic read distributions with strand and orientation context, benchmarked against current annotations and distance to nearest genes.

## Current capabilities

- Reference prep: FASTA + GTF fetch/normalize; STAR index build or use existing index.
- 10x/FASTQ ingest: local or manifest-driven downloads; optional subsampling for smoke tests.
- Alignment: STARsolo with CB/UMI handling; sorted BAM output.
- Intergenic filtering: retain mapped reads not overlapping gene features; export to BED.
- Distance metrics: nearest upstream/downstream gene with sign conventions (− upstream, + downstream).
- Summaries: per-chromosome and genome-wide histograms; TSVs for downstream plots.

Repo structure: `envs/`, `manifests/`, `scripts/`, `main.nf`, `nextflow.config`, `LICENSE`. :contentReference[oaicite:1]{index=1}

## Quick start

```bash
# 1) Environment (example)
module load java
# Nextflow installed or in PATH

# 2) Run: local FASTQs
nextflow run ryan-gorzek/intergenome \
  -profile local \
  --fastq_mode local \
  --fastq_dir /path/to/fastqs \
  --species mus_musculus \
  --build GRCm39

# 3) Run: manifest (URL<TAB>sample<TAB>folder)
nextflow run ryan-gorzek/intergenome \
  -profile local \
  --fastq_mode download \
  --fastq_manifest manifests/fastq_manifest.tsv \
  --species monodelphis_domestica \
  --build ASM229v1
```

Key params:

* `--species`, `--build`: Ensembl/NCBI identifiers used for reference fetching.
* `--star_index`: path to existing index to skip rebuild.
* `--n_cores`, `--mem_gb`: resource hints for heavy steps.

## Inputs

* FASTQs: 10x v3 (R1=barcode/UMI, R2=read).
* Reference: genome FASTA + GTF; provided or auto-fetched.
* Optional: inclusion list for CBs.

## Outputs

* `alignment/`: STARsolo outputs and coordinate-sorted BAMs.
* `beds/intergenic.bed.gz`: intergenic read loci.
* `nearest/nearest.tsv.gz`: read-to-gene signed distances with gene IDs and strand.
* `qc/`: basic histograms and TSV summaries.

## Signed distance convention

For a read at position x and nearest gene G:

* Negative = upstream of G’s TSS on the gene’s strand.
* Positive = downstream of G’s TES on the gene’s strand.
* Zero/overlap excluded by intergenic filter.

## Roadmap

* Export to HDF5/Parquet; modular plotting notebooks (Python/R).
* Strand-aware aggregation around regulatory landmarks (TSS/TES, enhancers).
* Aggregated species-level KDEs and permutation tests for distributional shifts.
* GC/repMask covariates and mappability normalization.
* Cross-species harmonization: liftOver/chain support where applicable.
* Optional spliced-alignment handling for long intergenic transcripts.

## License

MIT. ([GitHub][1])
