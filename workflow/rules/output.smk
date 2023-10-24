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
        expand("results/gffcompare/ed{EDIST}/duplicates/",
            EDIST=config["edit_distance"]
        ),
    ),
    final_output.extend(
        expand("results/gffcompare/ed{EDIST}/collapsedExpressionMatrices.rds",
            EDIST=config["edit_distance"]
        ),
    ),
    
    return final_output