import numpy as np
import pandas as pd 
import re
import os

samples = pd.read_csv("config/samples.tsv", sep=",")['sample'].values
os.mkdir(snakemake.output['duplicate_dir'])

def extractDuplicateInformation(slurm_path,samples):
    lines = [line.rstrip('\n') for line in list(open(slurm_path))]
    duplicatesPerSample = {}
    sample_index = -1
    for line in lines:
        if line.startswith('Loading query file'):
            sample_index = sample_index + 1
            print('processing ', samples[sample_index])
            duplicatesPerSample[samples[sample_index]] = pd.DataFrame(columns = ['dup1', 'dup2'])
            continue
        elif line.startswith('Warning'):
            continue
        elif line.startswith('\tQuery transcript'):
            dup1, dup2 = re.findall(r'transcript:(ENSG\d+_\d+_\d+_\d+|ENST\d+)', line)
            duplicatesPerSample[samples[sample_index]].loc[len(duplicatesPerSample[samples[sample_index]])] = [dup1,dup2]
    
    for sample in duplicatesPerSample.keys():
        print(sample, ': ', duplicatesPerSample[sample].shape)
        duplicatesPerSample[sample].to_csv(snakemake.output['duplicate_dir'] + '/duplicates_'+sample+'.tsv', sep='\t', index=False)

extractDuplicateInformation(snakemake.input['gfflog'], samples)