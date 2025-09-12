### This file will create the code list used for the type 1 linkage request with CPRD.
### This must be in a tab delimited format

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Read in code list
edh.cvd.hist.icd <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_cvd_hist_icd.csv", sep = ""), 
                                     sep = ",", header = TRUE, colClasses = "character")
edh.cvd.hist.icd <- dplyr::select(edh.cvd.hist.icd, term, icd) |>
  as.data.frame()

### Write to tab delimited txt file
write(edh.cvd.hist.icd$icd, "codelists/linkage_requests/cvd_hist_icd_type1_request.txt", sep = "\t")
