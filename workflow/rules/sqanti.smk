rule extractCounts: 
    input:
        countList="results/gffcompare/ed{EDIST}/collapsedExpressionMatrices.rds",
    output:
        countsFile="results/SQANTI/extractCounts/ed{EDIST}/TranscriptCountsPerSample.csv",
    singularity: 
        "singularities/sc_wf.sif"
    log: 
        stdout="logs/SQANTI/extractCounts/ed{EDIST}/extractCounts.stdout", 
        stderr="logs/SQANTI/extractCounts/ed{EDIST}/extractCounts.stderr",
    script:
        """
        ../scripts/extractCountsFile.R
        """

rule sqanti_qc: 
    input:
        mergedGtf="results/gffcompare/ed{EDIST}/mergedIsoforms.combined.gtf",
        refTSS="/SQANTI3-5.2/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed",
        polyA_motif="/SQANTI3-5.2/data/polyA_motif_list.txt",
        counts="results/SQANTI/extractCounts/ed{EDIST}/TranscriptCountsPerSample.csv",
        genes_gtf=config["genes_gtf"],
        genome_fa=config["genome_fa"],
    output:
        out_dir=directory("results/SQANTI/QC/ed{EDIST}"),
        classification="results/SQANTI/QC/ed{EDIST}/mergedIsoforms.combined_classification.txt",
        report="results/SQANTI/QC/ed{EDIST}/mergedIsoforms.combined_SQANTI3_report.html",
        gtf="results/SQANTI/QC/ed{EDIST}/mergedIsoforms.combined_corrected.gtf",
    singularity: 
        "singularities/sc_wf.sif"
    log: 
        stdout="logs/SQANTI/QC/ed{EDIST}/Sqanti_QC.stdout", 
        stderr="logs/SQANTI/QC/ed{EDIST}/Sqanti_QC.stderr",
    shell:
        """
        python SQANTI3-5.2/sqanti3_qc.py --force_id_ignore --CAGE_peak {input.refTSS} \
            --fl_count {output.trimmed_fastq} --polyA_motif_list {input.polyA_motif} \
            -d {output.out_dir} -fl {input.counts} \
            {input.mergedGtf}  {input.genes_gtf} {input.genome_fa}
        """
