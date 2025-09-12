###
### Create patids.devel and patids.valid, and the development and validation datasets
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("")
getwd()

### Define filepath to file directory system containing extracted data.
common.data.dir <- file.path("..", "..")

### Write function to create patids for the split sample for male/female cohorts
create_split_sample <- function(gender_var){
  
  ### Set seed
  set.seed(101)
  
  ### Read in imputed data
  mice_mids_imp <- readRDS(paste("data/mice_mids_prototype3_", gender_var, "_", 1, ".rds", sep = ""))
  
  ### Extract the imputed data
  df_imp <- mice::complete(mice_mids_imp, action = 1)

  ### Re-order levels in ethnicity variable
  eth_levels <- levels(df_imp$ethnicity)
  df_imp$ethnicity <- factor(df_imp$ethnicity, c(eth_levels[9],eth_levels[1:8]))
  
  ### Choose devel size
  devel_size <- 1000000
  
  ### Sample this many individuals
  patids_devel <- sample(df_imp$patid, devel_size, replace = FALSE)
  patids_valid <- df_imp$patid[is.na(fastmatch::fmatch(df_imp$patid, patids_devel))]
  
  ### Reduce
  df_imp_devel <- df_imp[!is.na(fastmatch::fmatch(df_imp$patid, patids_devel)), ]
  df_imp_valid <- df_imp[!is.na(fastmatch::fmatch(df_imp$patid, patids_valid)), ]
  
  ### Print length
  print(c("male", "female")[gender_var])
  print(paste("There are ", nrow(df_imp_devel), " individuals in the development cohort"))
  print(paste("There are ", nrow(df_imp_valid), " individuals in the validation cohort"))
  
  ### Save them
  saveRDS(patids_devel, paste("data/patids_devel_", gender_var, sep = ""))
  saveRDS(patids_valid, paste("data/patids_valid_", gender_var, sep = ""))
  saveRDS(df_imp_devel, paste("data/df_imp_devel_", gender_var, sep = ""))
  saveRDS(df_imp_valid, paste("data/df_imp_valid_", gender_var, sep = ""))
  
  print(paste("FINISHED", gender_var, Sys.time()))
  
}

### Run function
create_split_sample(1)
create_split_sample(2)

print(paste("FINISHED", Sys.time()))