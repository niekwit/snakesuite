'''Runs MACS2 on replicates of ChIP-seq data.
'''

import pandas as pd
import subprocess
import sys


g = snakemake.params[0]
q = snakemake.params[1]
e = snakemake.params[2]

#get sample info
csv = pd.read_csv("samples.csv")

genotypes = csv["genotype"].unique()
factor = csv["factor"].unique()
treatments = csv["treatment"].unique()

#create MACS2 command for replicates and run
if len(factor) == 1:
    
    for i in genotypes:
            
            for j in treatments:
                
                replicates = csv[(csv["genotype"] == i) & (csv["treatment"] == j)]["sample"].tolist()
                
                ip = [f"mapped/{i}_{x}_{factor[0]}_{j}_dedup.bam" for x in replicates]

                inpt = [f"mapped/{i}_{x}_input_{j}_dedup.bam" for x in replicates]
                
                macs2 = f"macs2 callpeak -t {' '.join(ip)} -c {' '.join(inpt)} -n {i}_{factor[0]}_{j} --outdir peaks/{i}_{factor[0]}_{j} -f BAMPE -g {g} -q {q} {e}"
                
                subprocess.run(macs2, shell=True)

else:
    
    print("Error: more than one factor detected. Please split data into separate analysis directories.")
    
    sys.exit(1)






