###
### This file will...
### 1) Double the size of the synthetic observation and patient files, creating new patids. The double 
### size files brings the file size (~800MB) more in line with the large files provided by CPRD (~1000MB)
### 2) Create 2000+ copies of the double size synthetic observation and patient files, again, 
### creating new patids
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project3/")
getwd()

###################################################################
### Double the size of the sythetic observation and patid files ###
###################################################################

### This is so the raw files are more similar in size to the files often shared by CPRD (1GB).
### Results in observation files of ~800MB.

### Read in observation file
raw_observation <- read.table("data/synthetic_raw/unzip/observation.txt", sep = "\t", header = TRUE,
                              colClasses = c("character","character","integer","character","character","character","character","character","character",
                                             "numeric","integer","integer","numeric","numeric","character"))

raw_pat <- utils::read.table("data/synthetic_raw/unzip/patient.txt", sep="\t", header = TRUE,
                  colClasses = c("character", "integer", "character", "integer", "integer", "integer", 
                                 "character", "character", "integer", "character", "integer", "character"))

str(raw_observation)
str(raw_pat)

###
### Create new patid variable from 1:n in both files
###

### Create new patid from 1:n
raw_pat <- dplyr::arrange(raw_pat, patid) |>
  dplyr::mutate(patid_new = 1:nrow(raw_pat))

### Merge with raw_observation and create new patid variable
raw_observation <- dplyr::left_join(raw_observation, raw_pat[,c("patid", "patid_new")], by = dplyr::join_by(patid)) |>
  dplyr::select(-patid) |>
  dplyr::rename(patid = patid_new) |>
  dplyr::relocate(patid) |>
  dplyr::arrange(patid)

### Change patid in patient file
raw_pat <- dplyr::select(raw_pat, -patid) |>
  dplyr::rename(patid = patid_new) |>
  dplyr::relocate(patid)

###
### Double the size of both files (i.e. double number of patients)
###

head(raw_observation)
### Create a new version of raw_observation, with new patids
raw_pat2 <- raw_pat |>
  dplyr::mutate(patid_new = (nrow(raw_pat) + 1):(2*nrow(raw_pat)))

### Create new raw_observation file
raw_observation2 <- dplyr::left_join(raw_observation, raw_pat2[,c("patid", "patid_new")], by = dplyr::join_by(patid)) |>
  dplyr::select(-patid) |>
  dplyr::rename(patid = patid_new) |>
  dplyr::relocate(patid) |>
  dplyr::arrange(patid)

### Change patid in patient file
raw_pat2 <- dplyr::select(raw_pat2, -patid) |>
  dplyr::rename(patid = patid_new) |>
  dplyr::relocate(patid)

### Concatenate
raw_observation <- rbind(raw_observation, raw_observation2)
raw_pat <- rbind(raw_pat, raw_pat2)

### Save these files
write.table(raw_pat, 
            "data/duplicated_raw/patient/patient1.txt", 
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

write.table(raw_observation, 
            "data/duplicated_raw/observation/observation1.txt", 
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")


### Read in to see if I have saved.txt file into correct format
# raw_pat_readin_original <- rcprd:::extract_txt_pat("data/synthetic_raw/unzip/patient.txt") |> 
#   dplyr::arrange(patid)
# raw_pat_readin_duplicated <- rcprd:::extract_txt_pat("data/duplicated_raw/patient/patient1.txt")
# 
# raw_obs_readin_original <- rcprd:::extract_txt_obs("data/synthetic_raw/unzip/observation.txt") |> 
#   dplyr::arrange(patid)
# raw_obs_readin_duplicated <- rcprd:::extract_txt_obs("data/duplicated_raw/observation/observation1.txt")
# 
# ### Test for equality, excluding the patid variable, which we have altered
# testthat::expect_equal(raw_pat_readin_original[1:10000, -1], raw_pat_readin_duplicated[1:10000, -1])
# testthat::expect_equal(raw_obs_readin_original[1:10000, -1], raw_obs_readin_duplicated[1:10000, -1])

### They are the same, therefore we can proceed

###########################################
### Create thousand copies of this data ###
###########################################

### Write a function to create a new copy of the data, with new patids
create_duplicate <- function(number){
  
  print(paste(number, Sys.time()))
  
  ### Create new patid from 
  raw_pat_new <- dplyr::mutate(raw_pat, patid_new =  ((number-1)*nrow(raw_pat) + 1):(number*nrow(raw_pat)))
  
  ### Merge with raw_observation and create new patid variable
  raw_observation_new <- dplyr::left_join(raw_observation, raw_pat_new[,c("patid", "patid_new")], by = dplyr::join_by(patid)) |>
    dplyr::select(-patid) |>
    dplyr::rename(patid = patid_new) |>
    dplyr::relocate(patid) |>
    dplyr::arrange(patid)
  
  ### Change patid in patient file
  raw_pat_new <- dplyr::select(raw_pat_new, -patid) |>
    dplyr::rename(patid = patid_new) |>
    dplyr::relocate(patid)
  
  ### Save these files
  write.table(raw_pat_new, 
              paste("data/duplicated_raw/patient/patient", number, ".txt", sep = ""), 
              row.names = FALSE,
              quote = FALSE,
              sep = "\t")
  
  write.table(raw_observation_new, 
              paste("data/duplicated_raw/observation/observation", number, ".txt", sep = ""),  
              row.names = FALSE,
              quote = FALSE,
              sep = "\t")
}

print("START CREATION OF FILES")
lapply(1:250, create_duplicate)