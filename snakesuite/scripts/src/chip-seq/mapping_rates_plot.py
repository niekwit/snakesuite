from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


#create df to store alignment rates
df = pd.DataFrame(columns=["sample","alignment_rate"],index=np.arange(len(snakemake.input)))
samples = []
rates = []

#get alignment rates from log files
for i in sorted(snakemake.input):
    #get sample name from file name
    samples.append(Path(i).stem)
    
    #extract alignment rate from file
    with open(i) as f:
        lines = f.readlines()
    rate = float([x for x in lines if "overall alignment rate" in x][0].replace("% overall alignment rate\n",""))
    rates.append(rate)

#add values to empty df for plotting
df["sample"] = samples
df["alignment_rate"] = rates

#create plot
sns.set_style("white")
sns.set_style("ticks")
sns.barplot(x=list(df.keys())[0],
                y=list(df.keys())[1],
                data=df,
                color="seagreen",
                edgecolor="black",
                linewidth=1)
plt.ylabel("Overall alignment rate (%)")
plt.xticks(rotation = 'vertical')
plt.xlabel("")
plt.tight_layout()
sns.despine()
plt.savefig(snakemake.output[0])
plt.close()





