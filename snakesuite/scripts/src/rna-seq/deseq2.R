####load snakemake variables####
map.with <- snakemake@params[[1]]
genome <- snakemake@params[[2]]
gtf <- snakemake@params[[3]]


####Load packages####
library(DESeq2)
library(GenomicFeatures)
library(tximport)


####DESeq2 analysis####

#create output dir
dir.create(file.path(getwd(),"deseq2"), showWarnings=FALSE)

#load sample data from samples.csv
samples <- read.csv("samples.csv", header=TRUE)

if (map.with == "salmon"){
  ####Prepare variables required by DESeq2 (Salmon)####
  #Create txdb from GTF file if it does not exist
  txdb.filename <- paste0(gsub(".gtf","",gtf, fixed=TRUE),".txdb")
  
  if(file.exists(txdb.filename) == FALSE){
    txdb <- makeTxDbFromGFF(gtf)
    saveDb(txdb, txdb.filename)
    txdb <- loadDb(txdb.filename)
  } else {txdb <- loadDb(txdb.filename)}
  
  #Create transcript to gene file
  k <- keys(txdb,keytype="TXNAME")
  tx2gene <- select(txdb,k,"GENEID","TXNAME")
  
  ####Load experiment information (Salmon)####
  genotypes <- unique(samples$genotype)
  conditions <- unique(samples$condition)
  samples$comb <- paste0(samples$genotype,"_",samples$condition)
  
  #read Salmon quant.sf files
  library(readr) #speeds up tximport
  
  files <- snakemake@input[["salmon"]]
  txi <- tximport(files, 
                  type="salmon", 
                  tx2gene=tx2gene)

  #create DESeqDataSet
  dds <- DESeqDataSetFromTximport(txi,
                                  colData = samples,
                                  design = ~ comb)
  
  #save DESeqDataSet to file (input for other scripts)
  save(dds, file=snakemake@output[["rdata"]])
} else if (map_with == "star"){
  
}

###Following steps are common to both Salmon and STAR derived data
library(dplyr)
library(biomaRt)

#load data for gene annotation
mart <- useMart("ensembl")

if (grepl("hg",genome, fixed=TRUE)) {
  mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
} else if (grepl("mm",genome, fixed=TRUE)){
  mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
} #might add other genomes later

#load reference samples
references <- unique(samples[samples$reference == "yes" ,]$comb)

#create nested list to store all pairwise comparisons (top level:references, lower level: samples without reference)
df.list <- vector(mode="list", length=length(references))
for (i in 1:length(references)){
  df.list[[i]] <- vector(mode="list", length=(length(unique(samples$comb)) - 1 ))
}
#df.list <- vector(mode="list", length=(length(references) * (length(unique(samples$comb)) - 1)))

#create vector to store contrast names
#names <- vector(mode="character", length=(length(references) * length(unique(samples$comb)) - 1 ))

#counter <- 1

#for each reference sample, perform pairwise comparisons with all the other samples
for (r in 1:length(references)){
  print(paste0("Setting reference level: ",references[r], " (",r,"/",length(references),")"))
  
  #copy dds
  dds_relevel <- dds
  
  #for each reference sample, set it as reference in dds
  dds_relevel$comb <- relevel(dds$comb, ref = references[r])
  
  #differential expression analysis
  dds_relevel <- DESeq(dds_relevel)
  
  #get comparisons
  library(stringr)
  comparisons <- resultsNames(dds_relevel)
  comparisons <- strsplit(comparisons," ")
  comparisons[1] <- NULL

  
  #create df for each comparison
  for (c in 1:length(comparisons)){
    comparison <- comparisons[[c]]
    comparison <- str_replace(comparison,"comb_","") #get name
    print(paste0("Specifiying contrast: ",comparison, " (",c,"/",length(comparisons),")"))
    
    res <- results(dds_relevel, name=comparisons[[c]])
    
    #res <- results(dds_relevel)
    
    df <- as.data.frame(res) %>%
      mutate(ensembl_gene_id = res@rownames, .before=1) 
    
    #create empty data frame to collect all data from all reference levels and contrasts at first one
    #if (r == 1 & c == 1 ){
    #  df.all <- as.data.frame(matrix(nrow=0, ncol=ncol(df)))
    #  names(df.all) <- names(df)
    #} 
    
    #annotate df
    df$ensembl_gene_id <- gsub("\\.[0-9]*","",df$ensembl_gene_id) #tidy up gene IDs
    gene.info <- getBM(filters = "ensembl_gene_id", 
                       attributes = c("ensembl_gene_id", 
                                      "external_gene_name",
                                      "description",
                                      "gene_biotype", 
                                      "chromosome_name",
                                      "start_position",
                                      "end_position", 
                                      "percentage_gene_gc_content"), 
                       values = df$ensembl_gene_id, 
                       mart = mart)
    
    df <- left_join(df,gene.info,by="ensembl_gene_id")
    
    #remove genes with baseMean zero
    df <- df[df$baseMean != 0, ]
    
    #add normalised read counts for each sample to df  
    temp <- as.data.frame(counts(dds_relevel, normalized=TRUE))
    temp$ensembl_gene_id <- row.names(temp)
    temp$ensembl_gene_id <- gsub("\\.[0-9]*","",temp$ensembl_gene_id) #tidy up gene IDs
    names(temp)[1:length(dds_relevel@colData@listData$sample)] <- dds_relevel@colData@listData$sample
    
    
    df <- left_join(df,temp, by="ensembl_gene_id")
    
    #move some columns around
    df <- df %>%
      relocate(external_gene_name, .after=ensembl_gene_id) %>%
      relocate((ncol(df)-(length(dds_relevel@colData@listData$sample)-1)):ncol(df), .after=baseMean)
    
    #order data for padj
    df <- df[order(df$padj), ]
    
    #add column with just contrast name
    df <- df %>%
      mutate(contrast_name = comparison, .before=1)
    
    #add data from df to df.all
    #df.all <- rbind(df.all, df)
    
    #counter for index in names
    counter <- counter + 1
    
    #sheet title can be max 31 characters
    #add column with contrast name and change it to a number (counter)
    df <- df %>%
      mutate(contrast_name = comparison, .before=1)
    
   
    #save df to df.list
    df.list[[r]][[c]] <- df
    
  }
  
}

#write df.all to csv file
#write.csv(df.all,
#          file=snakemake@output[["csv"]],
#          row.names=FALSE)

#function to flatten lists (https://stackoverflow.com/questions/16300344/how-to-flatten-a-list-of-lists/41882883#41882883)
flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}


#write all data frames in df.list to sheets in Excel file
library(openxlsx)

#flatten df.list
df.list <- flattenlist(df.list)

#name data frames in list
names(df.list) <- 1: (length(references) * (length(unique(samples$comb)) - 1))

write.xlsx(df.list, 
           file = snakemake@output[["xlsx"]],
           colNames = TRUE)

