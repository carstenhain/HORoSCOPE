# HORoSCOPE

HORoSCOPE (**H**igher-**O**rder **R**epeat **o**rganization and **S**ize of **C**entromeres using **O**ligonucleotide **P**rofiles for **E**stimation) is a Nextflow workflow for inferring centromere structure from short-read sequencing data. It estimates two key features of centromeres in a chromosome-specific manner:

- centromeric Higher-Order Repeat (HOR) architecture
- mean alpha-satellite/HOR array length

HORoSCOPE is based on a reference atlas of complete human centromere assemblies. In the associated study, distinct centromere architectures were defined from a mean of 493 complete centromere haplotypes per chromosome, and architecture-specific *k*-mers were selected from these assemblies. In parallel, *k*-mers whose dosage correlates with alpha-satellite/HOR array length were identified for each chromosome.

For new samples, HORoSCOPE counts predefined *k*-mer sets, normalizes *k*-mer dosage to account for differences in sequencing depth, and applies pretrained chromosome-specific models to infer centromere architecture and HOR array length.

The workflow supports several input formats, including FASTQ/FASTA files, BAM/CRAM alignments, and precomputed de Bruijn graph unitigs. It also provides alternative normalization strategies to account for variable copy-number states, for example in tumor genomes.


## Usage

Run from repository root:

```bash
nextflow run horoscope.nf \
	-c nextflow.config \
	--samplesheet samplesheet.csv \
	--outdir results
```

Optional parameters:

- `--kmer_fasta` (default: `data/all_tagging_kmers.fasta.gz`)
- `--model_directory` (default: `models/`)
- `--help`

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

- See samplesheet.csv for an example.
- For `FILE_TYPE=cram`, `CRAM_REFERENCE_PATH` must be provided.
- HORoSCOPE inference models were trained on de Bruijn graphs generated from short-read sequencing data. Therefore, predictions from precomputed de Bruijn graph unitigs may be more accurate than predictions from raw short-read data. For convenience, the workflow includes a MAKE_DBG parameter intended to automatically generate de Bruijn graph unitigs using BCALM. However, this step is not currently implemented. To generate de Bruijn graph input, please use the helper script `dbg_generation/make_dbg.sh` for now.

Normalization:

- *k*-mer dosage is normalized to account for differences in sequencing depth (and copy-number state). The normalization mode is provided in the samplesheet through the NORMALIZATION column.
- `global`: Uses all normalization *k*-mers across all chromosomes. Recommended for samples without copy-number alterations.
- `chromosome`: Uses normalization *k*-mers from the same chromosome as the architecture- and length-informative *k*-mers. Recommended for samples with whole-chromosome gains or losses.
- `p_arm` or `q_arm`: Uses normalization *k*-mers from the corresponding chromosome arm. Might be useful for cancer samples with recurrent arm-level copy-number changes.
- Number: Uses a user-defined numeric value for normalization instead of *k*-mer-derived normalization. This can be useful for cancers with widespread copy-number alterations. In the associated manuscript, all normalization *k*-mers from copy-number-neutral chromosome arms were merged, and their mean value was used as this parameter.


## Output

The workflow publishes the following files to `--outdir`:

- `final_kmer_merged.tsv`: merged *k* count table across all samples, including *k*-mer cluster annotation
- `normalization_metrics.tsv`: per-sample normalization metrics for each normalization strategy
- `centromere_genotyping.tsv`: final per-sample/per-chromosome inference table with:
	- `SAMPLE`
	- `CHROM`
	- `CLUSTER_H1`
	- `CLUSTER_H2`
	- `HOR_LENGTH`

## Repository Structure

- `scripts/`: python scripts for *k* extraction, table merging, and genotyping
- `models/`: pretrained models used for length inference
- `data/`: tagging *k*-mer resources
- `dbg_generation/`: helper script for DBG generation

## Software Requirements

### Core

- Nextflow (DSL2)
- Java (required by Nextflow)
- Python 3

### Command-line tools used by the workflow

- [`jellyfish`](https://github.com/gmarcais/jellyfish) (*k*-mer counting for read/alignment inputs)
- `samtools` (BAM/CRAM to FASTA conversion)
- `gzip`/`zcat`
- `awk`
- `grep`
- `pv`
- [`bcalm`](https://github.com/gatb/bcalm)
- [`lighter`](https://github.com/mourisl/Lighter)

### Python packages used by scripts

- `pandas`
- `numpy`
- `polars`
- `tqdm`
- `scikit-learn`

### Runtime environment notes

- The provided `nextflow.config` is configured for a Slurm environment.
- Adapt resources, queues, and settings as needed for your local cluster or workstation.

## Citation

If you use this repository, please cite the associated [publication].
