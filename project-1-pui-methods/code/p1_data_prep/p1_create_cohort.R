### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("")
getwd()

### Define filepath to file directory system containing extracted data.
common.data.dir <- file.path("..", "..")

### Load packages
library(mice)
library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(ggplot2)

### Extract cohort
cohort <- readRDS(file.path(common.data.dir, "Aurum_Jun2021_extract/data/extraction/cohort_baseline/cohort_var.rds"))
cohort <- data.table::as.data.table(cohort)
cohort <- filter(cohort, gender != 3)

### Redefine smoking to be three levels, Never, Ex, Current
cohort$smoking <- forcats::fct_recode(cohort$smoking, Current = "Light", Current = "Moderate", Current = "Heavy")

###
### Read in exposed and unexposed variables
### Note that this requires reading in 100 different files and concatenating
###

### Write a function to do this for a variable name
read_exposure <- function(varname, lookback = NULL, exposed = FALSE, exposed.type = NULL){
  
  ### Add exposed/unexposed
  if (exposed == TRUE){
    varname <- paste(varname, "_exposed", sep = "")
  } else if (exposed == FALSE){
    varname <- paste(varname, "_unexposed", sep = "")
  }
  
  ### Add lookback to varname
  if (!is.null(lookback)){
    varname <- paste(varname, "_", lookback, sep = "")
  }
  
  ### Add exposed.type
  if (!is.null(exposed.type)){
    varname <- paste(varname, "_", exposed.type, sep = "")
  }
  
  ### Extract
  do.call("rbind", 
          lapply(1:100, function(x) {
            readRDS(file.path(common.data.dir,
                              paste("Aurum_Jun2021_extract/data/extraction/cohort_baseline/var_", varname, "_id", x, ".rds", sep = "")))
          })
  )
}

###############################
### Read in and combine SBP ###
###############################

### Read in 'all' data on unexposed and exposed variables
### Unexposed is any test value without a prescription (from specified code list) in a 180 day window prior
### Exposed is any test value with a prescription (from specified code list) in a 1-90 day window prior

### Unexposed SBP
sbp_unexposed <- read_exposure("sbp", exposed = FALSE, exposed.type = "all") |>
  dplyr::rename(sbp_unexposed = sbp_unexposed_all)

### Remove individuals with no non-missing values and group by patid (this will take a while I believe)
sbp_unexposed <- subset(sbp_unexposed, !is.na(sbp_unexposed)) |>
  group_by(patid)

### Reduce unexposed SBP to the most recent one
sbp_unexposed <- sbp_unexposed |>
  dplyr::arrange(patid, desc(obsdate)) |>
  dplyr::filter(dplyr::row_number(patid) == 1) |>
  dplyr::select(patid, sbp_unexposed)


#########################
### Merge with cohort ###
#########################

### Merge all together
cohort <- Reduce(function(x, y) merge(x, y, by.x = "patid", by.y = "patid", all.x = TRUE),
                          list(cohort, sbp_unexposed)) |>
  dplyr::select(-c(cholhdl_ratio, sbp_var, bmi_alltime, sbp_alltime, sbp_var_alltime, cholhdl_ratio_alltime))
str(cohort)

### Split into female and male cohorts
cohort_female <- subset(cohort, gender == 2)
cohort_male <- subset(cohort, gender == 1)

### Reduce to 2 million female and male
set.seed(101)
cohort_female <- cohort_female[sample(1:nrow(cohort_female), 2000000, replace = FALSE), ]
cohort_male <- cohort_male[sample(1:nrow(cohort_male), 2000000, replace = FALSE),]

patids_female <- cohort_female$patid
patids_male <- cohort_male$patid

### Reduce the cohort to just the selected females and males
cohort <- rbind(cohort_female, cohort_male)

### Save the cohort
saveRDS(cohort, "data/cohort_pui.rds")
saveRDS(cohort_female, "data/cohort_female_pui.rds")
saveRDS(cohort_male, "data/cohort_male_pui.rds")

saveRDS(patids_female, "data/patids_female.rds")
saveRDS(patids_male, "data/patids_male.rds")

### Save a complete case cohort for developing code
cohort_female_cc <- cohort_female[complete.cases(cohort_female),]
saveRDS(cohort_female_cc, "data/cohort_female_pui_cc.rds")
print(paste("IMAGE SAVED", Sys.time()))
