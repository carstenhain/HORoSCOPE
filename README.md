# HORoSCOPE

HORoSCOPE is a Nextflow workflow for inferring centromere structure from short-read sequencing data. It estimates two key features of centromeres in a chromosome-specific manner:

- centromeric HOR architecture
- mean alpha-satellite/HOR array length

HORoSCOPE is based on a reference atlas of complete human centromere assemblies. In the associated study, distinct centromere architectures were defined from a mean of 493 complete centromere haplotypes per chromosome, and architecture-specific k-mers were selected from these assemblies. In parallel, k-mers whose dosage correlates with alpha-satellite/HOR array length were identified for each chromosome.

For new samples, HORoSCOPE counts predefined k-mer sets, normalizes k-mer dosage to account for differences in sequencing depth, and applies pretrained chromosome-specific models to infer centromere architecture and HOR array length.

The workflow supports several input formats, including FASTQ/FASTA files, BAM/CRAM alignments, and precomputed de Bruijn graph unitigs. It also provides alternative normalization strategies to account for variable copy-number states, for example in tumor genomes.



## Overview

Main workflow:

- `horoscope.nf` (recommended)

Auxiliary/legacy workflows in this repository:

- `genotype_centromere.nf`
- `genotype_centromere_train_test.nf`
- `kmer_search.nf`

## Software Requirements

### Core

- Nextflow (DSL2)
- Java (required by Nextflow)
- Python 3

### Command-line tools used by the workflow

- `jellyfish` (k-mer counting for read/alignment inputs)
- `samtools` (BAM/CRAM to FASTA conversion)
- `gzip`/`zcat`
- `awk`
- `grep`
- `pv`

### Python packages used by scripts

- `pandas`
- `numpy`
- `polars`
- `tqdm`
- `scikit-learn` (for loading and using trained models)

### Runtime environment notes

- The provided `nextflow.config` is configured for an HPC/Slurm environment and Apptainer.
- Adapt profiles, queues, and container settings as needed for your local cluster or workstation.

## Input

The main workflow (`horoscope.nf`) expects a semicolon-separated samplesheet with this header:

```text
NAME;FILE_PATH;FILE_TYPE;CRAM_REFERENCE_PATH;MAKE_DBG;NORMALIZATION
```

Column definitions:

- `NAME`: unique sample identifier
- `FILE_PATH`: path to sample input file
- `FILE_TYPE`: one of `dbg`, `fastq`, `fasta`, `bam`, `cram`
- `CRAM_REFERENCE_PATH`: reference FASTA for CRAM decoding; use `NA` for non-CRAM samples
- `MAKE_DBG`: `true` or `false`
- `NORMALIZATION`: one of `global`, `chromosome`, `p_arm`, `q_arm`, or a numeric value

Notes:

- For `FILE_TYPE=cram`, `CRAM_REFERENCE_PATH` must be provided.
- The current `MAKE_DBG=true` process is a placeholder in `horoscope.nf`; for production runs, use precomputed DBG unitigs (`FILE_TYPE=dbg`, `MAKE_DBG=false`) or read/alignment inputs with Jellyfish.

## Usage

Run from repository root:

```bash
nextflow run horoscope.nf \
	-c nextflow.config \
	--samplesheet samplesheet.csv \
	--outdir results \
	-resume
```

Optional parameters:

- `--kmer_fasta` (default: `data/all_tagging_kmers.fasta.gz`)
- `--model_directory` (default: `models/`)
- `--help`

## Output

The workflow publishes the following files to `--outdir`:

- `final_kmer_merged.tsv`: merged k-mer count table across all samples, including k-mer cluster annotation
- `normalization_metrics.tsv`: per-sample normalization metrics for each normalization strategy
- `centromere_genotyping.tsv`: final per-sample/per-chromosome inference table with:
	- `SAMPLE`
	- `CHROM`
	- `CLUSTER_H1`
	- `CLUSTER_H2`
	- `HOR_LENGTH`

## Repository Structure

- `scripts/`: helper scripts for k-mer extraction, table merging, and genotyping
- `models/`: pretrained models used for length inference
- `data/`: tagging k-mer resources and supporting files
- `dbg_generation/`: helper script for DBG generation

## Citation

If you use this repository for a manuscript, please cite the associated publication (add DOI/reference here once available).
