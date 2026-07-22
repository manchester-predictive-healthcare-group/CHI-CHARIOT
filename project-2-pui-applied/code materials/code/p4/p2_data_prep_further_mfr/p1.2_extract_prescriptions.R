###
### Extract all prescriptions for statinsantihypertensives, and save to disk, so we can more easily access when parallelising
### functions which require access to this data
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Source functions
R.func.sources = list.files(file.path(common.data.dir, "Aurum_Jun2021_extract/R"), full.names = TRUE)
sapply(R.func.sources, source)

### Extract cohort
cohort <- readRDS("data/p4/cohort_prototype3.rds")
str(cohort)

######################################################
### Function to extract and save all prescriptions ###
######################################################
extract_pres <- function(med){
  
  ### Define code list
  codelist <- paste("uom_", med, "_prodcodeid", sep = "")
  codelist <- data.table::fread(file = file.path(common.data.dir , paste("Aurum_Jun2021_extract/codelists/analysis/", codelist, ".csv", sep = "")),
                                sep = ",", header = TRUE, colClasses = "character")
  
  ### Query database for relevant prescriptions
  db.qry <- db_query(db.filepath = file.path(common.data.dir, "Aurum_Jun2021_extract/data/sql/aurum.sqlite"),
                     tab = "drug",
                     codelist.vec = codelist$prodcodeid)
  print(paste("query complete", Sys.time()))
  str(db.qry)
  
  ### Reduce db.qry to observations in cohort and get rid of unneccesary variables
  db.qry <- db.qry[!is.na(fastmatch::fmatch(db.qry$patid, cohort$patid)), ]
  db.qry <- dplyr::select(db.qry, patid, issuedate)
  print(paste("reduced to cohort", Sys.time()))
  str(db.qry)
  
  saveRDS(db.qry, paste("data/db_qry_", med, ".rds", sep = ""))
  print(paste("saved query", Sys.time()))
}

### Run function
extract_pres("statins")
extract_pres("antihypertensives")