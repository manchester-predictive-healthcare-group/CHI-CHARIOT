###
### Turn imputation object (mids) into complete datasets for analysis, with some minor formatting of variables
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

###
### Extract imputed data, format variables where neccesary
###
print(paste("female cohorts", Sys.time()))
lapply(1:10, function(x){
  
  ### Read in imputed data
  imp <- readRDS(paste("data/p5/mice_mids_prototype3_", 2, "_", x, ".rds", sep = ""))
  
  ### Extract the imputed data
  imp <- mice::complete(imp, action = 1)
  
  ### Re-order levels in ethnicity variable
  eth.levels <- levels(imp$ethnicity)
  imp$ethnicity <- factor(imp$ethnicity, c(eth.levels[9],eth.levels[1:8]))
  
  ### Add calendar time
  imp <- dplyr::mutate(imp, caltime = as.numeric(fup_start) - as.numeric(as.Date("01/01/2005", format = "%d/%m/%Y")))
  
  saveRDS(imp, paste0("data/p5/dfs_all_female", x, ".rds"))
  
})

print(paste("male cohorts", Sys.time()))
lapply(1:10, function(x){
  
  ### Read in imputed data
  imp <- readRDS(paste("data/p5/mice_mids_prototype3_", 1, "_", x, ".rds", sep = ""))
  
  ### Extract the imputed data
  imp <- mice::complete(imp, action = 1)
  
  ### Re-order levels in ethnicity variable
  eth.levels <- levels(imp$ethnicity)
  imp$ethnicity <- factor(imp$ethnicity, c(eth.levels[9],eth.levels[1:8]))
  
  ### Add calendar time
  imp <- dplyr::mutate(imp, caltime = as.numeric(fup_start) - as.numeric(as.Date("01/01/2005", format = "%d/%m/%Y")))
  
  saveRDS(imp, paste0("data/p5/dfs_all_male", x, ".rds"))
  
})

print(paste(Sys.time()), "FINISHED")