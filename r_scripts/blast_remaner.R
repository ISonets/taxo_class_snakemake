library(readxl)
path = "/mnt/disk1/PROJECTS/SURPRISE/eshi_i_myshi/RESULTS/results/blast/contigs/xlsx"
setwd(path)
files <- list.files(
  path = path, 
  recursive = TRUE,
  full.names = TRUE
)
names <- read_xlsx('/mnt/disk1/PROJECTS/SURPRISE/eshi_i_myshi/RESULTS/mgi 2023_eji_w_sizes.xlsx')
for (f in files) {
  id <- basename(f)
  rn <- as.integer(row.names(names)[which(names$Barcode == strsplit(id,'_')[[1]][4])])
  hr_name <- names$ID[rn]
  new_name <- paste(hr_name,sep='_',id)
  file.rename(id,new_name)
}
