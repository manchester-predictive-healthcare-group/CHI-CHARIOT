###
### Extract all prescriptions for statins, antihypertensives, and records of smoking, bmi/weight/height, 
### and chol/hdl/ldl/triglycerides, and save to disk, 
### These are variables that we may be adjusting for during follow-up.
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

###############################################################################
### Function to extract and save all records of each modifiable risk factor ###
###############################################################################
extract_med <- function(med, codelist){
  
  print(paste(med,Sys.time()))
  
  ### Define code list
  codelist <- data.table::fread(file = file.path(common.data.dir , paste("Aurum_Jun2021_extract/codelists/analysis/", codelist, ".csv", sep = "")),
                                sep = ",", header = TRUE, colClasses = "character")
  
  ### Query database for relevant prescriptions
  db.qry <- db_query(db.filepath = file.path(common.data.dir, "Aurum_Jun2021_extract/data/sql/aurum.sqlite"),
                     tab = "obs",
                     codelist.vec = codelist$medcodeid)
  print(paste("query complete", Sys.time()))
  str(db.qry)
  
  ### Reduce db.qry to observations in cohort and get rid of unneccesary variables
  db.qry <- db.qry[!is.na(fastmatch::fmatch(db.qry$patid, cohort$patid)), ]
  db.qry <- dplyr::select(db.qry, patid, medcodeid, obsdate, value, numunitid)
  print(paste("reduced to cohort", Sys.time()))
  str(db.qry)
  
  saveRDS(db.qry, paste("data/p4/db_qry_", med, ".rds", sep = ""))
  print(paste("saved query", Sys.time()))
}

### Run function
extract_med("component_bmi", "edh_bmi_medcodeid")
extract_med("component_weight", "weight_medcodeid")
extract_med("component_height", "height_medcodeid")

extract_med("sbp", "edh_sbp_medcodeid")

extract_med("smoking_non", "edh_smoking_non_medcodeid")
extract_med("smoking_ex", "edh_smoking_ex_medcodeid")
extract_med("smoking_light", "edh_smoking_light_medcodeid")
extract_med("smoking_mod", "edh_smoking_mod_medcodeid")
extract_med("smoking_heavy", "edh_smoking_heavy_medcodeid")

extract_med("component_cholesterol", "muzambi_cholesterol_medcodeid")
extract_med("component_hdl", "muzambi_hdl_medcodeid")
extract_med("component_ldl", "muzambi_ldl_medcodeid")
extract_med("component_triglycerides", "uom_triglycerides_medcodeid")