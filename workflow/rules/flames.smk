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
    conda:
        "../envs/flames.yml"
    #resources:
    #    mem_gb=32,
    #    runtime="2d"
    shell:
        """
        touch results/FLAMES/barcode_matching/ed{wildcards.EDIST}/{wildcards.SAMPLE}/{wildcards.SAMPLE}_trimmed.fastq

        #g++ -std=c++11 -lz -O2 -o resources/FLAMES/src/match_cell_barcode resources/FLAMES/src/ssw/ssw_cpp.cpp \
        #    resources/FLAMES/src/ssw/ssw.c resources/FLAMES/src/match_cell_barcode.cpp resources/FLAMES/src/kseq.h resources/FLAMES/src/edit_dist.cpp
        resources/FLAMES/src/bin/match_cell_barcode {input.longread_fastq} {output.barcode_stats} \
            {output.trimmed_fastq} {input.shortread_barcodes} {wildcards.EDIST} {params.umi_length}
        """

rule sc_long_pipeline:
    input: 
        trimmed_fastq="results/FLAMES/barcode_matching/ed{EDIST}/{SAMPLE}/{SAMPLE}_trimmed.fastq",
        genes_gtf=config["genes_gtf"],
        genome_fa=config["genome_fa"],
        config="resources/FLAMES/python/config_sclr_nanopore_default.json",
    output:
        out_dir="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}",
        isoforms_filtered="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/isoform_annotated.filtered.gff3",
        genome_bam="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/align2genome_{SAMPLE}.bam",
        genome_bai="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/align2genome_{SAMPLE}.bam.bai",
        transcriptome_bai="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/realign2transcript.bam.bai",
        counts="results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/transcript_count.csv.gz",
    conda:
        "../envs/flames.yml"
    #resources:
    #    mem_gb=32,
    #    runtime="2d"
    shell:
        """
        ##run FLAMES single-cell, long-read pipeline
        resources/FLAMES/python/sc_long_pipeline.py --infq {input.trimmed_fastq} --outdir $out \
            --genomefa {input.genome_fa} --gff3 {input.genes_gtf} --config_file {input.config}

        ##rename output bam files to avoid downstream problems with merging
        mv results/FLAMES/sc_long_pipeline/{wildcards.EDIST}/{wildcards.SAMPLE}/align2genome.bam {output.genome_bam}
        mv results/FLAMES/sc_long_pipeline/{wildcards.EDIST}/{wildcards.SAMPLE}/align2genome.bam.bai {output.genome_bai}
        """