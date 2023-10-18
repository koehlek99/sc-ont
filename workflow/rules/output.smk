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
    return final_output