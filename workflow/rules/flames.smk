configfile: "config/config.yaml"

import pandas as pd 
samples = pd.read_csv("config/samples.tsv", sep=",").set_index("sample", drop=False)

rule match_barcodes: 
    input:
        longread_fastq=lambda wildcards: samples["fastq_dir"][wildcards.SAMPLE],
        shortread_barcodes=lambda wildcards: samples["shortread_bcs"][wildcards.SAMPLE],
    output:
        trimmed_fastq="results/FLAMES/barcode_matching/ed{EDIST}/{SAMPLE}/{SAMPLE}_trimmed.fastq",
        barcode_stats="results/FLAMES/barcode_matching/ed{EDIST}/{SAMPLE}/{SAMPLE}_barcode.stats",
    params:
        umi_length=config["umi_length"],
    singularity: 
        "singularities/sc_wf.sif"
    log: 
        stdout="logs/FLAMES/barcode_matching/ed{EDIST}_{SAMPLE}.stdout", 
        stderr="logs/FLAMES/barcode_matching/ed{EDIST}_{SAMPLE}.stderr",
    shell:
        """
        #touch {output.trimmed_fastq}
        #touch {output.barcode_stats}
        /FLAMES/src/match_cell_barcode {input.longread_fastq} {output.barcode_stats} \
            {output.trimmed_fastq} {input.shortread_barcodes} {wildcards.EDIST} {params.umi_length} \
            > {log.stdout} 2> {log.stderr}
        """

rule sc_long_pipeline:
    input: 
        trimmed_fastq="results/FLAMES/barcode_matching/ed{EDIST}/{SAMPLE}/{SAMPLE}_trimmed.fastq",
        genes_gtf=config["genes_gtf"],
        genome_fa=config["genome_fa"],
        config="resources/FLAMES/config_sclr_nanopore_default.json",
    output:
        out_dir=directory("results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}"),
        isoforms_filtered="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/isoform_annotated.filtered.gff3",
        genome_bam="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/align2genome_{SAMPLE}.bam",
        genome_bai="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/align2genome_{SAMPLE}.bam.bai",
        transcriptome_bai="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/realign2transcript.bam.bai",
        counts="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/transcript_count.csv.gz",
    singularity: 
        "singularities/sc_wf.sif"
    log: 
        stdout="logs/FLAMES/sc_long_pipeline/ed{EDIST}_{SAMPLE}.stdout", 
        stderr="logs/FLAMES/sc_long_pipeline/ed{EDIST}_{SAMPLE}.stderr",
    shell:
        """
        ##run FLAMES single-cell, long-read pipeline
        /FLAMES/python/sc_long_pipeline.py --infq {input.trimmed_fastq} --outdir {output.out_dir} \
            --genomefa {input.genome_fa} --gff3 {input.genes_gtf} --config_file {input.config} \
            > {log.stdout} 2> {log.stderr}

        ##rename output bam files to avoid downstream problems with merging
        mv results/FLAMES/sc_long_pipeline/ed{wildcards.EDIST}/{wildcards.SAMPLE}/align2genome.bam {output.genome_bam}
        mv results/FLAMES/sc_long_pipeline/ed{wildcards.EDIST}/{wildcards.SAMPLE}/align2genome.bam.bai {output.genome_bai}
        """
