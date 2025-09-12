### Clear workspace
rm(list=ls())

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

###
### Add a couple of files
###
### NB: Will eventually turn this into a function which will extract all files of a certain type and create a permanent sql database
### NB: Consider applying the codelists at this point, to avoid storing too much data
### NB: Note that in rEHR they store the database in a temporary sql database, so may want to consider that. Shouldn't affect coding though,
### as the "import" style functions would remain the same, its just whether you run them everytime, or just once
### NB: Time how long it takes to read in all files, create temporary database, to see if its feasible to do this everytime, 
### or if i should create a permanent sql database
###


### Start by loading the patids
subset.patids <- readRDS("data/extraction/cohort_exclu1_reduced.rds")
subset.patids <- subset.patids[,c("patid", "set")]

###
### Now create the sqlite databases, subsettings on the patids that met the exclusion criteria

### Observation file
print(paste("Extract Observation file into sql", Sys.time()))
cprd_extract(dbname = "aurum",
             filetype = "obs", 
             select = c("patid", "obsid", "obsdate", "enterdate", "medcodeid", "value", "numunitid", "numrangelow", "numrangehigh", "probobsid"),
             subset.patids = subset.patids,
             use.set = TRUE)
print(paste("Observation extracted", Sys.time()))

### No need to extract DrugIssue file for initial exclusion criteria
print(paste("Extract Drug Issue file into sql", Sys.time()))
cprd_extract(dbname = "aurum",
             filetype = "drug", 
             select = c("patid", "issueid", "issuedate", "enterdate", "prodcodeid", "dosageid", "quantity", "quantunitid", "duration"),
             subset.patids = subset.patids,
             use.set = TRUE)
print(paste("DrugIssue extracted", Sys.time()))

print(paste("FINISHED", Sys.time()))