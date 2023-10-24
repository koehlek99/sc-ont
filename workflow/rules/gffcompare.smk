configfile: "config/config.yaml"

rule writeGtfList:
    input: 
        pathList=expand("results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/isoform_annotated.filtered.gff3",
            SAMPLE=samples["sample"],
            EDIST=config["edit_distance"]
        ),
    output: 
        gtfList="results/gffcompare/gtfList.txt",
    log: 
        stdout="logs/gffcompare/writeGtfList.stdout", 
        stderr="logs/gffcompare/writeGtfList.stderr",
    run:
        with open(output.gtfList, "w") as f:
            for path in input.pathList:
                f.write(path + "\n")

rule gffcompare:
    input: 
        gtfList="results/gffcompare/gtfList.txt", 
        genes=config["genes_gtf"],
    output: 
        gtf="results/gffcompare/ed{EDIST}/mergedIsoforms.combined.gtf",
        tracking="results/gffcompare/ed{EDIST}/mergedIsoforms.tracking",
        log="results/gffcompare/ed{EDIST}/gffcompare.log",    
    log: 
        stdout="logs/gffcompare/ed{EDIST}_gffcompare.stdout", 
        stderr="logs/gffcompare/ed{EDIST}_gffcompare.stderr",
    singularity: 
        "singularities/sc_wf.sif"
    shell: 
        """
        /gffcompare/gffcompare -i {input.gtfList} -r {input.genes} -o mergedIsoforms --debug &>> {output.log}
        mv mergedIsoforms* results/gffcompare/ed{wildcards.EDIST}/
        """

rule getDuplicateInfoFromSlurm:
    input: 
        gfflog="results/gffcompare/ed{EDIST}/gffcompare.log",
    output:
        duplicate_dir=directory("results/gffcompare/ed{EDIST}/duplicates")
    log: 
        stdout="logs/gffcompare/ed{EDIST}_getDuplicates.stdout", 
        stderr="logs/gffcompare/ed{EDIST}_getDuplicates.stderr",  
    singularity: 
        "singularities/sc_wf.sif"
    script:
        "../scripts/getDuplicateInfoFromLog.py"

rule processExpressionMatrices:
    input:
        expressionCounts = expand("results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/transcript_count.csv.gz",
                                SAMPLE=samples["sample"],
                                EDIST=config["edit_distance"]),
        duplicateDict = expand("results/gffcompare/ed{EDIST}/duplicates/duplicates_{SAMPLE}.tsv",
                                SAMPLE=samples["sample"],
                                EDIST=config["edit_distance"]),
        trackingDf ="results/gffcompare/ed{EDIST}/mergedIsoforms.tracking",
    params: 
        edist=config["edit_distance"],
    output:
        collapsedExpressionMatrices="results/gffcompare/ed{EDIST}/collapsedExpressionMatrices.rds",
    singularity: 
        "singularities/sc_wf.sif"
    script:
        "../scripts/collapsExpressionMatrices.R"



