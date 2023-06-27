library(ggplot2)
library(foreach)

files <- snakemake@input

#create output dir
dir.create(file.path(getwd(),"plots"), showWarnings=FALSE)

#create df for storing mapping rates
df <- as.data.frame(matrix(ncol=2,nrow=0))
names(df) <- c("sample","mapping.rate")

#get mapping rates from log files
counter <- 1
for (x in files) {
  sample <- system(paste0("echo ",x , " | sed 's$logs/salmon/$$' | sed 's/.log//'"), intern=TRUE)
  rate <- system(paste0('grep "Mapping rate = " ', x, " | awk -F ' ' '{print $NF}' | sed 's/%//'"), intern=TRUE)
  
  df[counter,"sample"] <- sample
  df[counter,"mapping.rate"] <- rate
  
  counter <- counter + 1
}

#round values to 1 decimal
df$mapping.rate <- as.numeric(df$mapping.rate)
df$mapping.rate <- round(df$mapping.rate, digits=1)

p <- ggplot(df, aes(x=sample, y=mapping.rate)) +
  geom_bar(stat = "identity", 
           fill="aquamarine4",
           colour = "black") +
  theme_bw(base_size = 18) +
  scale_y_continuous(limits=c(0,100)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("Mapping rate (%)") +
  xlab(NULL)

ggsave(snakemake@output[[1]], p)



