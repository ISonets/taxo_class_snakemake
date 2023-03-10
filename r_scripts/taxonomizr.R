#install.packages("taxonomizr")
library(taxonomizr)
library(writexl)
args = commandArgs(trailingOnly=TRUE)[1]
if (file.size(args) == 0L) {
print("File is empty")
args2 = gsub(".txt","_anno.xlsx",args)
cat(NULL,file=args2)
quit() }
file <- read.csv(args, sep='\t',header = F)
spcs <- file[[2]]
taxaId<-getId(spcs,"~/taxonomizr_data/taxa.sql")
file[,3] <- taxaId
colnames(file) <- c('count','species','taxaid')
tax_df <- as.data.frame(getTaxonomy(taxaId,"~/taxonomizr_data/taxa.sql"))
merged_df <- merge(file,tax_df)
args2 = gsub(".txt","_anno.xlsx",args)
write_xlsx(merged_df, args2)
