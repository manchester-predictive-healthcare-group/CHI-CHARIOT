###
### This file will create an 'history of' variable, manually cycling through all the raw .txt files
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project3/")
getwd()

### We assume a cohort has been defined and established
### Read in cohort
cohort <- readRDS("data/extract/cohort.rds")
cohort$ho <- 0

### Read in codelist
codelist <- data.table::fread("codelists/analysis/edh_hypertension_medcodeid.csv", sep = ",", header = TRUE, colClasses = "character")

### Create a loop
time_in <- Sys.time()
for (number in 1:250){
  print(paste("number = ", number, Sys.time()))
  
  ### Read in raw observation file
  raw_observation <- data.table::fread(
    file = paste("data/duplicated_raw/observation/observation", number, ".txt", sep = ""), 
    sep = "\t", header = TRUE,
    colClasses = c("character","character","integer","character","character","character","character","character","character",
                   "numeric","integer","integer","numeric","numeric","character"))
  
  ### Convert to dates where relevant
  raw_observation$obsdate <- as.Date(raw_observation$obsdate, format = "%d/%m/%Y")
  
  ### Reduce to where medcodeid matches codelist
  raw_observation <- raw_observation[!is.na(fastmatch::fmatch(raw_observation$medcodeid, codelist$medcodeid))]
  
  ### Merge with cohort
  cohort_observation <- dplyr::left_join(raw_observation, cohort[,c("patid", "index_date")], by = dplyr::join_by(patid))

  ### Keep where index_date
  cohort_observation <- cohort_observation[obsdate <= index_date]
  
  ### Create a binary variable 
  cohort$ho <- pmax(as.integer(!is.na(fastmatch::fmatch(cohort$patid, cohort_observation$patid))), cohort$ho)
  
}
time_out <- Sys.time()
time_diff <- time_out - time_in
time_diff

saveRDS(time_diff, "data/extract/time_diff_manual.rds")

### Save the variable
cohort <- dplyr::select(cohort, patid, ho)
saveRDS(cohort, "data/extract/ho_manual.rds")
print("FINISHED")