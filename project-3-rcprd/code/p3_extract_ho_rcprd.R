###
### This file will create an 'history of' variable using sqlite and rcprd functionality
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project3/")
getwd()
devtools::install_github("alexpate30/rcprd")
### Load rcprd
library(rcprd)

### Read in cohort
cohort <- readRDS("data/extract/cohort.rds")

### Establish connection to sqlite database
aurum_extract <- connect_database("data/sqlite/aurum_extract.sqlite")

### Create variable
time_in <- Sys.time()
ho_var <- extract_ho(cohort,
                     codelist = "edh_hypertension_medcodeid",
                     indexdt = "index_date",
                     db_open = aurum_extract,
                     tab = "observation",
                     return_output = TRUE)
time_out <- Sys.time()
time_diff <- time_out - time_in
time_diff
saveRDS(time_diff, "data/extract/time_diff_rcprd.rds")

RSQLite::dbDisconnect(aurum_extract)

### Save the variable
ho_var <- dplyr::select(ho_var, patid, ho)
saveRDS(ho_var, "data/extract/ho_rcprd.rds")
print("FINISHED")
