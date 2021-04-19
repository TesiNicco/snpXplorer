# Libraries
library(stringr)

# Change directory
setwd("/Users/home/Desktop/GWAS_catalog_sumstats/")

# Read csv with all studies
all_studies <- read.table("/Users/home/Downloads/list_gwas_summary_statistics_4_Nov_2020-2.csv", sep=",", h=T)

# First, add year of publication
date_tmp <- str_split_fixed(all_studies$Publication.date, "-", 3)
all_studies$Year <- as.numeric(date_tmp[, 1])

# Take only those from 2017
sbs_year <- all_studies[which(all_studies$Year >= 2018),]

# Main loop to download
for (i in 1:nrow(sbs_year)){
  print("Working on trait: ")
  # create a folder for the trait
  fold_name = paste0(sbs_year$PubMed.ID[i], "_", sbs_year$Trait.s.[i], "_", sbs_year$Year[i])
  fold_name_ok = str_replace_all(fold_name, pattern = " ", replacement = "_")
  system(paste("mkdir", fold_name_ok))
  
  # set wd temporarily
  setwd(paste0("/Users/home/Desktop/GWAS_catalog_sumstats/", fold_name_ok))
  
  # add files
  cmd = paste0("wget ", sbs_year$FTP.Path[i], "/*")
  system(cmd, ignore.stdout = T)
  
  # re-set working directory
  setwd("/Users/home/Desktop/GWAS_catalog_sumstats/")
}