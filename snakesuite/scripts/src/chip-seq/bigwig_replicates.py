'''Creates one bigwig file for all replicates
'''

import pandas as pd
import subprocess
import sys
import os


chrom_sizes = snakemake.params[0]

#get sample info
csv = pd.read_csv("samples.csv")

genotypes = csv["genotype"].unique()
factor = csv["factor"].unique()
treatments = csv["treatment"].unique()


def log_and_run(command, log):
    '''Runs command, and writes command, stdout and stderr to log file
    '''
    
    with open(log, "w") as f:
        
        f.write(f"{command}\n")
        
        subprocess.run(command, shell=True, stdout=f, stderr=f)
            


if len(factor) == 1:
    
    for i in genotypes:
            
            for j in treatments:
                
                replicates = csv[(csv["genotype"] == i) & (csv["treatment"] == j)]["sample"].tolist()
                
                #merge bigwig files to one bedGraph file
                def bedgraph(factor):
                    
                    samples = [f"bigwig/single/{i}_{x}_{factor}_{j}.bw" for x in replicates]
                    
                    merged_bg = f"bigwig/merged/{i}_{factor}_{j}.bg"
                    
                    merge = f"bigWigMerge {' '.join(samples)} {merged_bg}" 
                                        
                    log_and_run(merge, f"logs/merged_bw/{i}_{factor}_{j}_bedgraph.log")
                    
                    
                    return merged_bg
                
                
                #create ip bedgraph file
                ip_bg = bedgraph(factor[0])
                
                #create input bedgraph file
                input_bg = bedgraph("input")
                
                
                #convert bedGraph to bigwig
                def bigwig(factor, bg):
                                        
                    merged_bw = f"bigwig/merged/{i}_{factor}_{j}.bw"
                    
                    convert = f"bedGraphToBigWig {bg} {chrom_sizes} {merged_bw}"
                    
                    log_and_run(convert, f"logs/merged_bw/{i}_{factor}_{j}_bigwig.log")  
                   
                    
                #create ip bigwig file
                bigwig(factor[0], ip_bg)
                
                #create input bigwig file
                bigwig("input", input_bg)
                
                #remove bedgraph files
                os.remove(ip_bg)
                os.remove(input_bg)
              
else:
    
    print("Error: more than one factor detected. Please split data into separate analysis directories.")
    
    sys.exit(1)







