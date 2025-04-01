###
### This file will create an sqlite database for querying following procedure of rcprd
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project3/")
getwd()

### Load rcprd
library(rcprd)

### Create cohort file from patient files
print(paste("create cohort"))
cohort <- extract_cohort("data/duplicated_raw/patient")

### Read in practice file
raw_prac <- rcprd:::extract_txt_prac("data/synthetic_raw/unzip/practice.txt")

### Merge the two
cohort <- merge(cohort, raw_prac, by.x = "pracid", by.y = "pracid")
head(cohort)

### Define follow-up start at max of regstartdate and uts
### NB: In synthetic data, uts is all NA, so ignore
cohort <- dplyr::mutate(cohort, fup_start = regstartdate)

### Individuals without a regstartdate are removed
cohort <- dplyr::filter(cohort, !is.na(fup_start))

### Define follow-up end at min of regenddate and practice last collection date
### If regeenddate is NA, we ignore through na.rm = TRUE, and take last collection date of practice,
### because these individuals at still actively registered
cohort <- dplyr::mutate(cohort, fup_end = pmin(regenddate, lcd, na.rm = TRUE))

### Reduce to individuals actively registered on 1st Jan 2010
cohort <- dplyr::filter(cohort, 
                        fup_start <= as.Date("01/01/2010", format = "%d/%m/%Y") & 
                          fup_end > as.Date("01/01/2010", format = "%d/%m/%Y")
)

### Add index date
cohort$index_date <- as.Date("01/01/2010", format = "%d/%m/%Y")
saveRDS(cohort, "data/extract/cohort.rds")

###
### Going to create two sqlite databases
### 1) The first just for patients that meet the inclusion/exclusion criteria
### 2) The second for all patients
###

###
### Cohort 1
###

### Create connection
print(paste("start sqlite creation"))
aurum_extract <- connect_database("data/sqlite/aurum_extract.sqlite")

### Add obs files
cprd_extract(db = aurum_extract, 
             filepath = "data/duplicated_raw/observation",
             filetype = "observation",
             subset_patids = cohort$patid)

### Disconnect
RSQLite::dbDisconnect(aurum_extract)
print("FINISHED 1")

###
### Cohort 2
###

### Create connection
print(paste("start sqlite creation"))
aurum_extract_all <- connect_database("data/sqlite/aurum_extract_all.sqlite")

### Add obs files
cprd_extract(db = aurum_extract_all, 
             filepath = "data/duplicated_raw/observation",
             filetype = "observation")

### Disconnect
RSQLite::dbDisconnect(aurum_extract_all)
print("FINISHED 2")