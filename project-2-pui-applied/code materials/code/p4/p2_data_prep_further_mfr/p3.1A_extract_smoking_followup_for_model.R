###
### This program will extract smoking status during follow-up
### This program is written to be parallelised, i.e. will extract status for a subgroups of patients
### We set number of groups to be 100, and task ID cycled through 1 - 100.
### They will need to be recombined after (done in program p3.1B).
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

### Extract chain seed and gender from command line
args <- commandArgs(trailingOnly = T)
num.groups <- as.numeric(args[1])
group.id <- as.numeric(args[2])
print(paste("num.groups = ", num.groups))
print(paste("group.id = ", group.id))

### Read in cohort
cohort <- readRDS("data/p4/cohort_prototype3.rds") |>
  dplyr::select(patid, cvd_time, cvd_indicator, fup_start)
# cohort <- readRDS(file.path(common.data.dir, "Aurum_Jun2021_extract/data/extraction/cohort_baseline/cohort_var.rds"))
# cohort <- dplyr::select(cohort, patid, cvd_time, cvd_indicator, fup_start)

### Split up patids
chunk.size <- nrow(cohort)/num.groups
patids.split <- split(cohort$patid, ceiling(seq_along(cohort$patid) / chunk.size))
# sum(unlist(lapply(patids.split, length))) 

### Read in the db_qry of the events of interest (contains all prescription dates of interest)
db_qry_smoking_non <- readRDS(paste("data/p4/db_qry_", "smoking_non", ".rds", sep = ""))
db_qry_smoking_ex <- readRDS(paste("data/p4/db_qry_", "smoking_ex", ".rds", sep = ""))
db_qry_smoking_light <- readRDS(paste("data/p4/db_qry_", "smoking_light", ".rds", sep = ""))
db_qry_smoking_mod <- readRDS(paste("data/p4/db_qry_", "smoking_mod", ".rds", sep = ""))
db_qry_smoking_heavy <- readRDS(paste("data/p4/db_qry_", "smoking_heavy", ".rds", sep = ""))

### Reduce to patients of interest
cohort_reduced <- cohort[!is.na(fastmatch::fmatch(cohort$patid, patids.split[[group.id]])), ]
db_qry_smoking_non <- db_qry_smoking_non[!is.na(fastmatch::fmatch(db_qry_smoking_non$patid, patids.split[[group.id]])), ]
db_qry_smoking_ex <- db_qry_smoking_ex[!is.na(fastmatch::fmatch(db_qry_smoking_ex$patid, patids.split[[group.id]])), ]
db_qry_smoking_light <- db_qry_smoking_light[!is.na(fastmatch::fmatch(db_qry_smoking_light$patid, patids.split[[group.id]])), ]
db_qry_smoking_mod <- db_qry_smoking_mod[!is.na(fastmatch::fmatch(db_qry_smoking_mod$patid, patids.split[[group.id]])), ]
db_qry_smoking_heavy <- db_qry_smoking_heavy[!is.na(fastmatch::fmatch(db_qry_smoking_heavy$patid, patids.split[[group.id]])), ]

### Add variable for smoking status (0 = never, 1 = ex, 2 = current)
db_qry_smoking_non$value <- 0
db_qry_smoking_ex$value <- 1
db_qry_smoking_light$value <- 2
db_qry_smoking_mod$value <- 2
db_qry_smoking_heavy$value <- 2

### Combine
db_qry_reduced <- rbind(db_qry_smoking_non, db_qry_smoking_ex, db_qry_smoking_light, db_qry_smoking_mod, db_qry_smoking_heavy) |>
  dplyr::arrange(patid, obsdate) |>
  dplyr::mutate(value_new = value)

### Write function to change all zeros to 1's, if it's preceded by a 1 or a 2
my_mutate <- function(df, id){
  
  for (i in 1:nrow(df)){
    if (i > 1){
      if (df$value_new[i] == 0){
        df$value_new[i] <- min(1, max(df$value[1:(i-1)]))
      }
    }
  }
  
  return(df)
  
}

### Apply function
db_qry_reduced <- dplyr::group_by(db_qry_reduced, patid) |>
  dplyr::group_modify(my_mutate) |>
  dplyr::mutate(value = value_new) |>
  dplyr::select(-c(value_new))

### Filter to 1 obs per person on the same date
db_qry_reduced <- dplyr::arrange(db_qry_reduced, patid, obsdate, desc(value)) |>
  dplyr::group_by(patid, obsdate) |>
  dplyr::filter(dplyr::row_number() == 1) |>
  as.data.frame() 

### Merge with cohort to get relative time to index date
db_qry_reduced <- dplyr::left_join(db_qry_reduced, cohort_reduced[,c("patid", "fup_start")], by = dplyr::join_by(patid)) |>
  dplyr::mutate(times_adj = obsdate - as.numeric(fup_start))
### Write a test to make sure an individual never goes back to zero after leaving it...???

### Remove observations that happened prior to index date
db_qry_reduced_after <- dplyr::filter(db_qry_reduced, times_adj > 0)

####################################################################################
### Function to split the survival times based on medication status change times ###
####################################################################################

### Currently, I have all smoking status times, I just want the times they change
### This includes changing from baseline

### Function will split survival times based on when an individuals smoking status changes
### We always assume the first row is a status change. This means we can use the same
### outcome follow-up data with every imputed dataset. If this value is the same as the imputed value,
### i.e. there is no split, this is ok, the smoking variable will just stay the same over the split.

### Writing a function that can be utilised with dplyr::group_modify and dplyr::group_map
### Data is the cohort file
### id is patient ID
### mst is the medication status times (i.e. when individual have smoking status updates)
get_cut_surv_times <- function(data, id, mst){
  
  ### Reduce mst to just patient of interest
  mstpat <- mst[!is.na(fastmatch::fmatch(mst$patid, id)), c("times_adj", "value")]
  
  ### If mstpat is non-empty, need to create split survival times
  if (nrow(mstpat) > 0){
    ### Create empty vector
    change_status_rows <- 1
    
    ### If more than 1 smoking event
    if (nrow(mstpat) > 1){
      
      for (i in 2:nrow(mstpat)){
        ### Check if next value is same as previous one, if yes do nothing
        if (mstpat$value[i] == mstpat$value[i-1]){
          ### Do nothing
          ### Else
        } else {
          ### If its different, record the row number
          change_status_rows <- append(change_status_rows, i)
        }
      }
    }
    
    ### Reduce to the rows where smoking status changes
    mstpat <- mstpat[change_status_rows, ]
    
    ### Get splits
    splits <- mstpat$times_adj
    
    ### Apply survsplit
    data.split <- survival::survSplit(data = data, cut = mstpat$times_adj, end = "cvd_time", event = "cvd_indicator")
    data.split <- data.split[,c("tstart", "cvd_time", "cvd_indicator", "fup_start")]
    
    ### Merge by tstat and times_adj
    data.split <- merge(data.split, mstpat[, c("times_adj", "value")], 
                        by.x = c("tstart"), 
                        by.y = c("times_adj"), 
                        all.x = TRUE) |>
      dplyr::rename(smoking = value) |>
      dplyr::relocate(tstart, cvd_time, cvd_indicator, smoking, fup_start)
    
  } else {
    ### If mstpat is empty, just modify data to be in correct format
    data.split <- dplyr::mutate(data, tstart = 0, smoking = NA) |>
      dplyr::relocate(tstart, cvd_time, cvd_indicator, smoking, fup_start)
  }
  
  return(tibble::tibble(data.split))
  
}

### Apply this function after grouping cohort_reduced
cohort_reduced <- dplyr::group_by(cohort_reduced, patid)
print(paste("APPLY FUNCTION", Sys.time()))
cohort_split_times <- dplyr::group_modify(.data = cohort_reduced, 
                                          .f = get_cut_surv_times, 
                                          mst = db_qry_reduced_after) |>
  as.data.frame()

### Save data
saveRDS(cohort_split_times, paste("data/p4/cohort_split_times_smoking_id", group.id, ".rds", sep = ""))
print(paste("FINISHED", Sys.time()))