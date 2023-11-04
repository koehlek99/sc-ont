import subprocess

if snakemake.config['filteringMethod']=="Rules":
    print('running rules filter..')
    cmd = ["python", "/SQANTI3-5.2/sqanti3_filter.py", "rules", "--gtf", snakemake.input['mergedGtf'], 
           "-o", "Flames_RulesFilter", "-d", snakemake.output['out_dir'], "-j", snakemake.input['ruleJson'], 
           snakemake.input['classificationInput']]
    subprocess.run(cmd)
elif snakemake.config['filteringMethod']=='ML':
    cmd = ["python", "/SQANTI3-5.2/sqanti3_filter.py", "ml", "--gtf", snakemake.input['mergedGtf'], 
           "-o", "Flames_MLFilter", "-d", snakemake.output['out_dir'], "-r", snakemake.input['removeCols'], 
           snakemake.input['classificationInput']]
    subprocess.run(cmd)