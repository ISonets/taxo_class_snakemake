# making it and saving under readable_names, for all files in loop
# read files in folder
library(readxl)
library(reshape)
library(tidyverse)
library(dplyr)
library(taxonomizr)
# read paths
path1 = "/mnt/disk1/PROJECTS/SURPRISE/eshi_i_myshi/RESULTS/results/blast/contigs/xlsx"
path2 = "/mnt/disk1/PROJECTS/SURPRISE/eshi_i_myshi/RESULTS/results/catbat/xlsx"
path3 = "/mnt/disk1/PROJECTS/SURPRISE/eshi_i_myshi/RESULTS/results/kraken2/contigs/xlsx"

#read all blast, kraken2 and catbat file names
files1 <- list.files(
  path = path1, 
  recursive = TRUE,
  full.names = TRUE
)
files2 <- list.files(
  path = path2, 
  recursive = TRUE,
  full.names = TRUE
)

files3 <- list.files(
  path = path3, 
  recursive = TRUE,
  full.names = TRUE
)

l <- length(files1)
# all species found
# iterate through files; thanks to the same prefixes files are in the same order in all folders
for (i in 1:l ) {
  catbat <- read_xlsx(files3[i])
  kraken2 <- read_xlsx(files2[i])
  blast <- read_xlsx(files1[i])
  blast$taxo <- paste(blast$superkingdom,blast$phylum,blast$class,blast$order,blast$family,blast$genus,blast$species,sep=';')
  catbat$taxo <- paste(catbat$superkingdom,catbat$phylum,catbat$class,catbat$order,catbat$family,catbat$genus,catbat$species,sep=';')
  kraken2$taxo <- paste(kraken2$superkingdom,kraken2$phylum,kraken2$class,kraken2$order,kraken2$family,kraken2$genus,kraken2$species,sep=';')
  blast_4merge <- select(blast, contig_name, contig_length, coverage, taxo)
  kraken2_4merge <- select(kraken2, contig_name, taxo)
  catbat_4merge <- select(catbat, contig_name, taxo)
  data_merged <- merge(blast_4merge,kraken2_4merge,by='contig_name',suffixes=c('_blast','_kraken2'))
  data2_merged <- merge(data_merged,catbat_4merge,by='contig_name')
  names(data2_merged)[names(data2_merged) == 'taxo'] <- 'taxo_catbat'
  data2_merged$taxo_blast<- gsub(x=data2_merged$taxo_blast,pattern = ';NA',replacement = '')
  data2_merged$taxo_kraken2<- gsub(x=data2_merged$taxo_kraken2,pattern = ';NA',replacement = '')
  data2_merged$taxo_catbat<- gsub(x=data2_merged$taxo_catbat,pattern = ';NA',replacement = '')
  nn <- gsub('blastres','all_merged',files1[i])
  write_xlsx(data2_merged,nn)
}  
# viruses
  
for (i in 1:l ) {
  #read files
  catbat <- read_xlsx(files3[i])
  kraken2 <- read_xlsx(files2[i])
  blast <- read_xlsx(files1[i])
  #subsetting
  blast_viruses <- subset(blast, superkingdom=='Viruses')
  kraken2_viruses <- subset(kraken2, superkingdom=='Viruses')
  catbat_viruses <- subset(catbat, superkingdom=='Viruses')
  # merging taxo
  blast_viruses$taxo <- paste(blast_viruses$superkingdom,blast_viruses$phylum,blast_viruses$class,
                              blast_viruses$order,blast_viruses$family,blast_viruses$genus,blast_viruses$species,sep=';')
  catbat_viruses$taxo <- paste(catbat_viruses$superkingdom,catbat_viruses$phylum,catbat_viruses$class,catbat_viruses$order,
                               catbat_viruses$family,catbat_viruses$genus,catbat_viruses$species,sep=';')
  kraken2_viruses$taxo <- paste(kraken2_viruses$superkingdom,kraken2_viruses$phylum,kraken2_viruses$class,kraken2_viruses$order,
                                kraken2_viruses$family,kraken2_viruses$genus,kraken2_viruses$species,sep=';')
  # selecting only needed columns
  blast_viruses_4merge <- select(blast_viruses, contig_name, taxo)
  kraken2_viruses_4merge <- select(kraken2_viruses, contig_name, taxo)
  catbat_viruses_4merge <- select(catbat_viruses, contig_name, taxo)
  data_viruses <- merge(blast_viruses_4merge,kraken2_viruses_4merge,by='contig_name',suffixes=c('_blast','_kraken2'), all=TRUE)
  data2_viruses <- merge(data_viruses,catbat_viruses_4merge,by='contig_name', all=TRUE)
  names(data2_viruses)[names(data2_viruses) == 'taxo'] <- 'taxo_catbat'
  data2_viruses$taxo_blast<- gsub(x=data2_viruses$taxo_blast,pattern = ';NA',replacement = '')
  data2_viruses$taxo_kraken2<- gsub(x=data2_viruses$taxo_kraken2,pattern = ';NA',replacement = '')
  data2_viruses$taxo_catbat<- gsub(x=data2_viruses$taxo_catbat,pattern = ';NA',replacement = '')
  data2_viruses$contig_length <- as.numeric(lapply(strsplit(as.character(data2_viruses$contig_name),'_'),"[",4))
  data2_viruses$coverage <- as.numeric(lapply(strsplit(as.character(data2_viruses$contig_name),'_'),"[",6))
  data2_viruses <- data2_viruses %>% relocate(contig_length, .after=contig_name)
  data2_viruses<- data2_viruses %>% relocate(coverage, .after=contig_length)
  nn <- gsub('blastres','viruses_all_merged',files1[i])
  write_xlsx(data2_viruses,nn)
} 
for (i in 1:l ) {
  #read files
  catbat <- read_xlsx(files3[i])
  kraken2 <- read_xlsx(files2[i])
  blast <- read_xlsx(files1[i])
  # counter
  blast_count <- blast %>% count(taxaid)
  names(blast_count)[names(blast_count) == 'n'] <- 'n_blast'
  kr2_count <- kraken2 %>% count(taxaid)
  names(kr2_count)[names(kr2_count) == 'n'] <- 'n_kraken2'
  catbat_count <- catbat %>% count(taxaid)
  names(catbat_count)[names(catbat_count) == 'n'] <- 'n_catbat'
  df2_list <- list(blast_count, catbat_count, kr2_count)
  data_counts <- merge_recurse(df2_list,by='taxaid')
  data_counts[,5:11] <- as.data.frame(getTaxonomy(data_counts[[1]],"/mnt/disk1/DATABASES/taxonomizr_data/taxa.sql"))
  nn <- gsub('blastres','all_counts',files1[i])
  write_xlsx(data_counts,nn)
  data_counts_viruses <- subset(data_counts, superkingdom=='Viruses')
  nn1 <- gsub('blastres','viruses_counts',files1[i])
  write_xlsx(data_counts_viruses,nn1)
}
