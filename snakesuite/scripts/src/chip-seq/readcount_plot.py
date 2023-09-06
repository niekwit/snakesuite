import os
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


#get bam files
prededup = snakemake.output[0]
dedup = snakemake.output[1]

#create df to store read counts
column_names = [os.path.basename(x).replace("_dedup.bam","") for x in dedup]
df = pd.DataFrame(columns = column_names)

#count reads pre and post deduplication
for pre,post,column in zip(prededup,dedup,column_names):
    
    #count reads pre deduplication
    count = pysam.view("-@", str(snakemake.threads), "-c", "-F", "260", pre)
    df.loc[1,column] = int(count)
    
    #add condtion
    df["condition"] = "pre-deduplication"
    
    #count reads post deduplication
    count = pysam.view("-@", str(snakemake.threads), "-c", "-F", "260", post)
        
    df.loc[2,column] = int(count)

    #add condtion
    df.loc[2,"condition"] = "post-deduplication"

#create df for plotting
df_melt = pd.melt(df, id_vars = ["condition"],
                    value_vars = column_names)
df_melt["value"] = pd.to_numeric(df_melt["value"])
df_melt["value"] = df_melt["value"] / 1000000

#create plot
sns.catplot(x = 'variable', y = 'value',
            hue = 'condition',
            data = df_melt,
            kind = 'bar',
            legend_out = False,
            edgecolor = "black",)
plt.ylabel("Uniquely mapped read count (millions)")
plt.xticks(rotation = 45, ha="right")
plt.xlabel("")
plt.ylim((0, df_melt["value"].max() * 1.3))
plt.legend(title = None,
            frameon = False)
plt.tight_layout()

#save plot to file
plt.savefig(snakemake.output[0])

plt.close()



