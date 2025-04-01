###
### This program will test the peak and total RAM used for each method
###
library(peakRAM)

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project3/")
getwd()

### Load rcprd
library(rcprd)

### This program will compare the RAM for each approach for extracting a 'history of' type variable
### Approach 1 (using rcprd, sqlite database restricted to individuals meeting criteira)
### Approach 2 (using rcprd, sqlite database of all individuals)
### Approach 3 (manually reading in)
### Approach 4 (same as approach 2, but not subsetting the data to derive the variable)
### NB: Approach 4 will not derive the variable, but provides a hard lower bound on the time and RAM used to manually
### read in all the files. This is if users think they can do the "derive the variable bit" more efficiently.

### Read in cohort
cohort <- readRDS("data/extract/cohort.rds")

###
### Define function for approach 1
###

### NB: A function doesn't need defining for this approach, as we can just use extract_ho directly

###
### Define function for approach 2
###

### NB: A function doesn't need defining for this approach, as we can just use extract_ho directly

###
### Define function for approach 3
###
approach3 <- function(cohort.in){
  
  ### Read in codelist
  codelist <- data.table::fread("codelists/analysis/edh_hypertension_medcodeid.csv", 
                                sep = ",", 
                                header = TRUE, 
                                colClasses = "character")
  
  ### Create variable
  cohort.in$ho <- 0
  
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
    cohort_observation <- dplyr::left_join(raw_observation, cohort.in[,c("patid", "index_date")], by = dplyr::join_by(patid))
    
    ### Keep where index_date
    cohort_observation <- cohort_observation[obsdate <= index_date]
    
    ### Create a binary variable 
    cohort.in$ho <- pmax(as.integer(!is.na(fastmatch::fmatch(cohort.in$patid, cohort_observation$patid))), cohort.in$ho)
    
  }
  
}

###
### Define function for approach 4
###
approach4 <- function(){
  
  for (number in 1:250){
    print(paste("number = ", number, Sys.time()))
    
    ### Read in raw observation file
    raw_observation <- data.table::fread(
      file = paste("data/duplicated_raw/observation/observation", number, ".txt", sep = ""), 
      sep = "\t", header = TRUE,
      colClasses = c("character","character","integer","character","character","character","character","character","character",
                     "numeric","integer","integer","numeric","numeric","character"))
    
  }
}


### 
### Get peakRAMs
###

### Approach 1

# Establish connection to sqlite database
aurum_extract <- connect_database("data/sqlite/aurum_extract.sqlite")

# Extract var
peakram_approach1 <- peakRAM(extract_ho(cohort,
                                        codelist = "edh_hypertension_medcodeid",
                                        indexdt = "index_date",
                                        db_open = aurum_extract,
                                        tab = "observation",
                                        return_output = TRUE))
RSQLite::dbDisconnect(aurum_extract)

print(paste("approach 1", Sys.time()))
saveRDS(peakram_approach1, "data/extract/peakram_approach1.rds")
save.image("data/extract/p6_compare_RAM.RData")

### Approach 2

# Establish connection
aurum_extract_all <- connect_database("data/sqlite/aurum_extract_all.sqlite")

# Extract var
peakram_approach2 <- peakRAM(extract_ho(cohort,
                                        codelist = "edh_hypertension_medcodeid",
                                        indexdt = "index_date",
                                        db_open = aurum_extract_all,
                                        tab = "observation",
                                        return_output = TRUE))
RSQLite::dbDisconnect(aurum_extract_all)

print(paste("approach 2", Sys.time()))
saveRDS(peakram_approach2, "data/extract/peakram_approach2.rds")
save.image("data/extract/p6_compare_RAM.RData")

### Approach 3
peakram_approach3 <- peakRAM(approach3(cohort))

print(paste("approach 3", Sys.time()))
saveRDS(peakram_approach3, "data/extract/peakram_approach3.rds")
save.image("data/extract/p6_compare_RAM.RData")

### Approach 4
peakram_approach4 <- peakRAM(approach4())

print(paste("approach 4", Sys.time()))
saveRDS(peakram_approach4, "data/extract/peakram_approach4.rds")
save.image("data/extract/p6_compare_RAM.RData")

### Create table and round
peakram_table <- rbind(peakram_approach1, peakram_approach2, peakram_approach3, peakram_approach4)
peakram_table[, 2] <- round(peakram_table[, 2], 1)
peakram_table[, 3] <- round(peakram_table[, 3], 1)
peakram_table[, 4] <- round(peakram_table[, 4], 1)

### Save as csv
write.csv(peakram_table, "data/extract/peakram_table.csv", headers = TRUE, row.names = TRUE, quotes = FALSE)
