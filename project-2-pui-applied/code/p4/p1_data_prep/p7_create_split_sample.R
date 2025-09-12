###
### Create patids.devel and patids.valid, and the development and validation datasets
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Define filepath to file directory system containing extracted data.
common.data.dir <- file.path("..", "..")

### Extract cohort
cohort <- readRDS("data/p4/cohort_prototype3.rds")

### Write function to create patids for the split sample for male/female cohorts
create_split_sample <- function(cohort, gender.var){
  
  ### Set seed
  set.seed(101)
  
  ### Reduce to gender of interest
  cohort.reduced <- base::subset(cohort, gender == gender.var)
  
  ### Get 70% of data
  devel.size <- round(nrow(cohort.reduced)*70/100)
  
  ### Sample this many individuals
  patids.devel <- sample(cohort.reduced$patid, devel.size, replace = FALSE)
  patids.valid <- cohort.reduced$patid[is.na(fastmatch::fmatch(cohort.reduced$patid, patids.devel))]
  
  ### Print length
  print(c("male", "female")[gender.var])
  print(paste("There are ", length(patids.devel), " individuals in the development cohort"))
  print(paste("There are ", length(patids.valid), " individuals in the validation cohort"))
  
  ### Save them
  saveRDS(patids.devel, paste("data/p4/patids_devel_", gender.var, sep = ""))
  saveRDS(patids.valid, paste("data/p4/patids_valid_", gender.var, sep = ""))
  
  print(paste("FINISHED", gender.var, Sys.time()))
  
}

### Run function
create_split_sample(cohort, 1)
create_split_sample(cohort, 2)

### Read in patids.devel and patids.valid into R workspace
patids.devel.female <- readRDS(paste("data/p4/patids_devel_", 2, sep = ""))
patids.valid.female <- readRDS(paste("data/p4/patids_valid_", 2, sep = ""))

patids.devel.male <- readRDS(paste("data/p4/patids_devel_", 1, sep = ""))
patids.valid.male <- readRDS(paste("data/p4/patids_valid_", 1, sep = ""))

###
### Create development datasets
###
print(paste("female devel", Sys.time()))
devel.female <- lapply(1:10, function(x){
  
  ### Read in imputed data
  imp <- readRDS(paste("data/p4/mice_mids_prototype3_", 2, "_", x, ".rds", sep = ""))
  
  ### Extract the imputed data
  imp <- mice::complete(imp, action = 1)
  
  ### Reduce to development cohort
  imp <- imp[!is.na(fastmatch::fmatch(imp$patid, patids.devel.female)), ]
  
  ### Re-order levels in ethnicity variable
  eth.levels <- levels(imp$ethnicity)
  imp$ethnicity <- factor(imp$ethnicity, c(eth.levels[9],eth.levels[1:8]))
  
  ### Add calendar time
  imp <- dplyr::mutate(imp, caltime = as.numeric(fup_start) - as.numeric(as.Date("01/01/2005", format = "%d/%m/%Y")))
  
  return(imp)
})
## Check levels ordered correctly
print(levels(devel.female[[1]]$ethnicity))
print(levels(devel.female[[1]]$smoking))

print(paste("male devel", Sys.time()))
devel.male <- lapply(1:10, function(x){
  
  ### Read in imputed data
  imp <- readRDS(paste("data/p4/mice_mids_prototype3_", 1, "_", x, ".rds", sep = ""))
  
  ### Extract the imputed data
  imp <- mice::complete(imp, action = 1)
  
  ### Reduce to development cohort
  imp <- imp[!is.na(fastmatch::fmatch(imp$patid, patids.devel.male)), ]
  
  ### Re-order levels in ethnicity variable
  eth.levels <- levels(imp$ethnicity)
  imp$ethnicity <- factor(imp$ethnicity, c(eth.levels[9],eth.levels[1:8]))
  
  ### Add calendar time
  imp <- dplyr::mutate(imp, caltime = as.numeric(fup_start) - as.numeric(as.Date("01/01/2005", format = "%d/%m/%Y")))
  
  return(imp)
})
## Check levels ordered correctly
print(levels(devel.male[[1]]$ethnicity))
print(levels(devel.male[[1]]$smoking))

###
### Create validation datasets
###
print(paste("female valid", Sys.time()))
valid.female <- lapply(1:10, function(x){
  
  ### Read in imputed data
  imp <- readRDS(paste("data/p4/mice_mids_prototype3_", 2, "_", x, ".rds", sep = ""))
  
  ### Extract the imputed data
  imp <- mice::complete(imp, action = 1)
  
  ### Reduce to validation cohort
  imp <- imp[!is.na(fastmatch::fmatch(imp$patid, patids.valid.female)), ]
  
  ### Re-order levels in ethnicity variable
  eth.levels <- levels(imp$ethnicity)
  imp$ethnicity <- factor(imp$ethnicity, c(eth.levels[9],eth.levels[1:8]))
  
  ### Add calendar time
  imp <- dplyr::mutate(imp, caltime = as.numeric(fup_start) - as.numeric(as.Date("01/01/2005", format = "%d/%m/%Y")))
  
  return(imp)
})

print(paste("male valid", Sys.time()))
valid.male <- lapply(1:10, function(x){
  
  ### Read in imputed data
  imp <- readRDS(paste("data/p4/mice_mids_prototype3_", 1, "_", x, ".rds", sep = ""))
  
  ### Extract the imputed data
  imp <- mice::complete(imp, action = 1)
  
  ### Reduce to validation cohort
  imp <- imp[!is.na(fastmatch::fmatch(imp$patid, patids.valid.male)), ]
  
  ### Re-order levels in ethnicity variable
  eth.levels <- levels(imp$ethnicity)
  imp$ethnicity <- factor(imp$ethnicity, c(eth.levels[9],eth.levels[1:8]))
  
  ### Add calendar time
  imp <- dplyr::mutate(imp, caltime = as.numeric(fup_start) - as.numeric(as.Date("01/01/2005", format = "%d/%m/%Y")))
  
  return(imp)
})

### Print str of each object
str(devel.male)
str(devel.female)
str(valid.male)
str(valid.female)

### Save datasets
saveRDS(devel.male, "data/p4/dfs_devel_male.rds")
saveRDS(devel.female, "data/p4/dfs_devel_female.rds")
saveRDS(valid.male, "data/p4/dfs_valid_male.rds")
saveRDS(valid.female, "data/p4/dfs_valid_female.rds")
print(paste("FINISHED", Sys.time()))