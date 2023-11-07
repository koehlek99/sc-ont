configfile: "config/config.yaml"

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
        "../scripts/extractCountsFile.R"

rule sqanti_qc: 
    input:
        mergedGtf="results/gffcompare/ed{EDIST}/mergedIsoforms.combined.gtf",
        refTSS="resources/SQANTI3/human.refTSS_v3.1.hg38.bed",
        polyA_motif="resources/SQANTI3/mouse_and_human.polyA_motif.txt",
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
        python /SQANTI3-5.2/sqanti3_qc.py --force_id_ignore --CAGE_peak {input.refTSS} \
            --fl_count {input.counts} --polyA_motif_list {input.polyA_motif} \
            -d {output.out_dir} -fl {input.counts} \
            {input.mergedGtf}  {input.genes_gtf} {input.genome_fa}
        """

rule sqanti_filter: 
    input:
        classificationInput="results/SQANTI/QC/ed{EDIST}/mergedIsoforms.combined_classification.txt",
        mergedGtf="results/gffcompare/ed{EDIST}/mergedIsoforms.combined.gtf",
        ruleJson="resources/SQANTI3/SqantiFilterRules.json",
        removeCols="resources/SQANTI3/removeCols.txt",
    output:
        out_dir=directory("results/SQANTI/Filtering/ed{EDIST}/{FILTER}"),
        report="results/SQANTI/Filtering/ed{EDIST}/{FILTER}/Flames_{FILTER}Filter_SQANTI3_filter_report.pdf",
        gtf="results/SQANTI/Filtering/ed{EDIST}/{FILTER}/Flames_{FILTER}Filter.filtered.gtf",
        classificationOutput="results/SQANTI/Filtering/ed{EDIST}/{FILTER}/Flames_{FILTER}Filter_{FILTER}result_classification.txt",
    singularity: 
        "singularities/sc_wf.sif"
    log: 
        stdout="logs/SQANTI/Filtering/ed{EDIST}/{FILTER}/Sqanti_Filter.stdout", 
        stderr="logs/SQANTI/Filtering/ed{EDIST}/{FILTER}/Sqanti_Filter.stderr",
    script:
        "../scripts/SqantiFilter.py"