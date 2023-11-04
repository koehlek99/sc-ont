import pandas as pd 

configfile: "config/config.yaml"
samples = pd.read_csv("config/samples.tsv", sep=",").set_index("sample", drop=False)

def get_final_output():
    final_output = []
    final_output.extend(
        expand("results/FLAMES/barcode_matching/ed{EDIST}/{SAMPLE}/{SAMPLE}_trimmed.fastq",
            SAMPLE=samples["sample"],
            EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
        expand("results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/transcript_count.csv.gz",
            SAMPLE=samples["sample"],
            EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
        expand("results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/align2genome_{SAMPLE}.bam",
            SAMPLE=samples["sample"],
            EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
        expand("results/FLAMES/sc_long_pipeline/ed{EDIST}/{SAMPLE}/isoform_annotated.filtered.gff3",
            SAMPLE=samples["sample"],
            EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
         expand("results/gffcompare/gtfList.txt"),
    ),
    final_output.extend(
        expand("results/gffcompare/ed{EDIST}/mergedIsoforms.tracking",
        EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
        expand("results/gffcompare/ed{EDIST}/mergedIsoforms.combined.gtf",
            EDIST=config["edit_distance"]        
        ),
    ),
    final_output.extend(
        expand("results/gffcompare/ed{EDIST}/gffcompare.log",
            EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
        expand("results/gffcompare/ed{EDIST}/duplicates/duplicates_{SAMPLE}.tsv",
            EDIST=config["edit_distance"],
            SAMPLE=samples["sample"]
        ),
    ),
    final_output.extend(
        expand("results/gffcompare/ed{EDIST}/collapsedExpressionMatrices.rds",
            EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
        expand("results/SQANTI/extractCounts/ed{EDIST}/TranscriptCountsPerSample.csv",
            EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
        expand("results/SQANTI/QC/ed{EDIST}/mergedIsoforms.combined_SQANTI3_report.html",
            EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
        expand("results/SQANTI/Filtering/ed{EDIST}/{FILTER}/Flames_{FILTER}Filter_SQANTI3_filter_report.pdf",
            EDIST=config["edit_distance"],
            FILTER=config["filteringMethod"]
        ),
    ),
    final_output.extend(
        expand("results/SQANTI/Filtering/ed{EDIST}/{FILTER}/Flames_{FILTER}Filter_{FILTER}result_classification.txt",
            EDIST=config["edit_distance"],
            FILTER=config["filteringMethod"]
        ),
    ),
    return final_output