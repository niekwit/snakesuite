#inspired by https://gist.github.com/slavailn/c1f926482769ee170e74848bed40f674

library(ChIPpeakAnno)
library(GenomicRanges)
library(GenomicFeatures)
library(biomaRt)
library(stringr)

#load snakemake params
genome <- snakemake@params[[1]]
regions <- snakemake@params[[2]]

#create TxDb for selected genome
if (str_detect(genome, "hg")) {
  
  organism <- "Homo Sapiens"
  dataset <- "hsapiens_gene_ensembl"
  
} else if (str_detect(genome, "mm")) {
  
  organism <- "Mus musculus"
  dataset <- "mmusculus_gene_ensembl"
  
}

EnsemblTxDb <- makeTxDbFromEnsembl(organism = organism, release = 99)

#set data format setting
if (regions == "narrow"){
  
  format <- "narrowPeak"
  
} else if (regions == "broad") {
  
  format <- "broadPeak"
  
}

#create mart object
mart <- useMart("ensembl")
mart <- useDataset(dataset, mart = mart)

#read peak file
#peaks <- read.table(snakemake@input[[1]], header=FALSE, sep = "\t")

#convert to GRanges
macs <- readLines(file)
macs <- textConnection(macs)
gr <- toGRanges(macs, format="MACS")


gr <- toGRanges(snakemake@input[[1]], format=format, header = FALSE)

#create GRanges object with annotations from TxDb database
annoData <- toGRanges(EnsemblTxDb, feature="gene")

#annotate granges with the nearest TSS
annot <- annotatePeakInBatch(gr, 
                             AnnotationData=annoData, 
                             featureType = "TSS",
                             output="nearestLocation",
                             PeakLocForDistance = "start")

#add gene information
annot <- addGeneIDs(annot, mart = mart, feature_id_type = "ensembl_gene_id",
                    IDs2Add = c("external_gene_name","chromosome_name","gene_biotype","description"))

#write to file
write.table(annot, 
            snakemake@output[[1]], 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)





