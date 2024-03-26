# pipeline-umi-amplicon

### Overview

`pipeline-umi-amplicon` is a pipeline for generating high accuracy single
molecule reads using unique molecular identifiers (UMIs) from amplicon data.
The pipeline accepts FASTQ-format sequence files as input and outputs both
aligned reads, variant calls and QC stats.

*This pipeline has been forked and modified from the original [ONT pipeline-umi-amplicon repository](https://github.com/nanoporetech/pipeline-umi-amplicon).*

### Features

The pipeline performs the following steps:
- Reads are mapped to reference genome using minimap2
- Separate into amplicons
- Extract UMI sequences for all reads
- Cluster UMI sequences per amplicon using vsearch and compute high accuracy consensus reads
- Align high accuracy consensus reads and perform simple variant calling

<br>


# Getting Started

### Requirements
The following software packages must be installed prior to running:

-  [miniconda3](https://conda.io/miniconda.html) - please refer to installation [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Installation
After installing miniconda3, install the pipeline as follows:
```bash
# Get pipeline
git clone https://github.com/rj-price/pipeline-umi-amplicon.git 
# Change to directory
cd pipeline-umi-amplicon
# Create conda environment with all dependencies
conda env create -f environment.yml
# Activate environment
conda activate pipeline-umi-amplicon
# Install python packages provided by pipeline-umi-amplicon
cd lib && pip install . && cd ..

# To test if the installation was successful run
snakemake -j 1 -pr --configfile config.yml
# Deactivate environment
conda deactivate
```

### Input

To run the pipeline the following input files are required:

| Input | Description |
|-------|-------------|
| Reference genome | FASTA file containing the reference genome (e.g. GRCh38 for human) |
| Nanopore reads | Folder containing FASTQ files or a single concatenated FASTQ file. |
| Target / Amplicon | A BED file containing the chromosome, start and end coordinate and the name of the amplicon |

### BED format
Tab separated and needs a unique name:
```
chr1    107167322       107168239       target_a_chr1_107167756_T_C
```

### Output

The main output files created by the pipeline are:

| Output | Description |
|--------|-------------|
| Aligned reads | Aligned reads in indexed and sorted BAM format |
| Variant calls | Called variants in VCF format |

After the a pipeline analysis has completed, the aligned reads can be found at `{output_folder}/{run_name}/07_align_consensus/{amplicon_name}_final.bam` e.g. `example_egfr_single_read_run/07_align_consensus/EGFR_917_final.bam`.

### Usage:

To run the pipeline with default settings invoke snakemake as follows.

```bash
$ snakemake -j 30 variants --configfile config.yml
```

`-j` specifies how many CPU cores will be used by the pipeline.

### Options

The pipeline accepts several input parameters. They can either be changed in the `config.yml` file or specified when running snakemake.

For example:
```bash
snakemake -j 30 reads --config input_fastq=data reference_fasta=data/example_egfr_reference.fasta targets_bed=data/example_egfr_amplicon.bed
```

### Required parameters

These parameters have to be specified to run the pipeline.

| Parameter | Allowed | Description |
|-----------|---------|-------------|
| sample_name | String | Name of the output folder |
| input_fastq | Absolute file path | FASTQ file or folder containing FASTQ files |
| reference_fasta | Absolute file path | FASTA file containing the reference genome |
| targets_bed | Absolute file path | BED file containing amplicon coordinates and names |

### Optional parameters

See `config.yml`

<br>

# Original Notices
### Licence and Copyright

(c) 2020 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

### References and Supporting Information
If you use this pipeline please cite:

- Li, H. (2018) Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100.
- Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25(16), 2078-2079.
- Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4:e2584.
- Koboldt, D., Zhang, Q., Larson, D., Shen, D., McLellan, M., Lin, L., Miller, C., Mardis, E., Ding, L., & Wilson, R. (2012). VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research, 22, 568-576.

### Research Release

Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools. Support for this
software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests. However much as we would like to
rectify every issue and piece of feedback users may have, the developers may
have limited resource for support of this software. Research releases may be
unstable and subject to rapid iteration by Oxford Nanopore Technologies.

<br>

******************
### TODO
- Merge pipelines, remove old env, update test env name
- Number of reads collapsed in the header of the bam file 
- Remove split by amplicon step? (only one primer set in params, so useless?)
- Add optional barcoding step
- Check that balance reads is working as its meant to
- Add handling of errors if no reads aligned to a reference
- Create Docker image

### DONE
- Add stats output and plots, number of umis, umi distribution, read alignment stats
- Skip second UMI extraction step
- Add BED file to variant calling to specify amplicon region
