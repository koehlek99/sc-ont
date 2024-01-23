# Snakemake workflow for the identification and quantification of transcripts from long-read, scRNA-seq

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32.0-brightgreen.svg)](https://snakemake.bitbucket.io)


This workflow enables an analysis of long-read, scRNA-seq data generated with 10x Genomics single-cell protocols and nanopore sequencing. It utilizes the single-cell, long-read pipeline FLAMES to identify and quantify transcripts on a single-cell level.
Integration of data from multiple samples is performed using GffCompare that compares structures of identified transcript models between samples. Based on GffCompare's output, transcript IDs are adapted and expression matrices can be used for downstream analyses, e.g. for traditional single-cell analysis with Seurat. 
Since long-read, scRNA-seq experiments introduce a lot of noise, thorough quality control of identified transcript models is recommend. Therefore, the SQANTI3 was integrated into this pipeline. It computes various QC metrics and implements two filtering methods to separate true isoforms from technical artefacts. 


## Authors

- Kristin Köhler (@koehlek99)

## Usage

Build a singularity image from the definition file to install all necessary software into the apptainer

```bash
apptainer build sc_wf.sif singularities/sc_wf.def 
```

Specify the sample names, .fastq file directories and matching CellRanger barcodes in config/samples.tsv

| Sample   | Fastq directory | CellRanger short-read barcodes
| -------  | --------------- | -------------------
| Sample1  | sample1_fastq/  | sample1_barcodes.tsv
| Sample2  | sample2_fastq/  | sample2_barcodes.tsv


Configure the config/config.yml 

```yml
##configuration for FLAMES barcode matching 
edit_distance: 1
umi_length: 12 

##reference files 
genes_gtf: "hg38_pathogens/genes/genes.gtf"
genome_fa: "hg38_pathogens/fasta/genome.fa"

##sqanti classification (rule-based or ML-based)
filteringMethod: "rules"
```

Run the snakemake workflow using the singularity apptainer
```bash
snakemake --use-singularity snakeflow.sif --use-conda -s workflow/Snakefile
```

## Workflows 

#### FLAMES
The pipeline [FLAMES](https://doi.org/10.1186/s13059-021-02525-6 (Full-Length Analysis of Mutations and Splicing in long-read RNA-seq data)
) was developed to identify and quantify transcripts from long-read, scRNA-seq experiments. It performs two major steps:

- Barcode matching: Short-read detected barcodes, identified by CellRanger, are compared against potential barcode sequences within the long reads. If sequences match, with a maximal number of mismatches specified in the config (default 1), the long-read barcode sequences are corrected and UMIs are identified as the 12bp sequence downstream to the barcode. Both sequences are integrated into the .fastq headers of the long-read data.
- Long-read, single-cell pipeline: FLAMES maps the corrected long reads against the genome. Similar reads are grouped generating consensus isoforms and assigned to a reference transcript if TSSs, TTSs and splice junctions match. Transcripts are quantified based on UMI counts and transcript x cell expression matrices are generated. Additionally, FLAMES outputs .gff files containing identified transcript annotations. 


#### GffCompare

[GffCompare](https://github.com/gpertea/gffcompare) merges the transcript annotation per sample to facilitate a comparison of identified isoforms within and between samples. Transcripts having similar structures get identical IDs while non-equivalent isoforms get unique IDs. The IDs are used to update the expression matrices. Additionally, expression counts are aggregated in case of structurally equivalent isoforms within one sample. This part outputs a merged .gtf annotation as well as the adapted transcript x cell expression matrices. 

#### SQANTI3
- Quality Control: [SQANTI3](https://github.com/ConesaLab/SQANTI3) is used for quality control of the long-read defined transcriptomes. It requires the merged .gtf file and calculates transcript-level QC features. A detailed description of each feature can be found in SQANTI's wiki. Moreover, transcript models are grouped according to SQANTI's categorization scheme.
- Filtering: Two filtering methods are implemented in SQANTI3. The rule-based approach requires manually defined rules. Those can be specified in the .json configuration file (sc-ont/resources
/SQANTI3/SqantiFilterRules.json). Alternatively, a random forest classifier can be trained. By default, it uses reference matches as true positives while novel-not-in-catalogue isoforms with non-canonical splice junctions comprise the true negatives.


## Output

This workflow outputs transcript x cell expression matrices per sample. Additionally, .gtf files containing sample-specific transcript annotations and a merged .gtf file containing transcripts identified from all samples are generated. Expression matrices can be aggregated to gene-level and analyzed using traditional single-cell analysis software such as Seurat. 
