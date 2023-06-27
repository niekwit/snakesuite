suppressPackageStartupMessages({
  library(ggplot2)
library(DESeq2)
})

#create output dir
dir.create(file.path(getwd(),"plots"), showWarnings=FALSE)

print("Generating PCA plot...")

#load DESeq2 data
load(snakemake@input[[1]])

#transform data
vsd <- vst(dds, blind=FALSE)

#create PCA plot
pcaData <- plotPCA(vsd, intgroup=c("genotype", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


p <- ggplot(pcaData, aes(PC1, PC2, color=genotype, shape=condition)) +
  geom_point(size=8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw(base_size = 16)

#save PCA plot
ggsave(p, file = snakemake@output[[1]])




