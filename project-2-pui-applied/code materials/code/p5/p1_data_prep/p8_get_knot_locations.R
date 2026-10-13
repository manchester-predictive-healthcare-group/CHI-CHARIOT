###
### Get knot locations by combining all imputed datasets, and calculating relevant quantiles
### This is done for SBP, BMI and non-HDL cholesterol. Age knot locations based off age.

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Define filepath to file directory system containing extracted data.
common.data.dir <- file.path("..", "..")

### Source functions
R.func.sources = list.files("Aurum_Jun2021_extract/R", full.names = TRUE)
sapply(R.func.sources, source)

### Extract chain seed and gender from command line
for (gender in 1:2){
  
  ### Assign character gender
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  
  ### Read in imputed data
  imp.comb <-  readRDS(paste("data/p5/dfs_all_", gender_char, ".rds", sep = ""))
  
  ### Create a stacked dataset
  imp.comb <- do.call("rbind", imp.comb)
  str(imp.comb)
  
  ### Write a function to extract the knot locations for a variable
  get_knot_loc <- function(varname){
    knots <- Hmisc::rcspline.eval(imp.comb[,varname], nk = 3)
    knots.loc <- attributes(knots)$knots
    saveRDS(knots.loc, paste("data/p5/knots_loc_", varname, "_", gender, sep = ""))
    print(paste(varname, knots.loc))
  }
  
  ### Run function for SBP and nonhdl cholesterol
  print(paste("gender =", gender))
  print("sbp")
  get_knot_loc("sbp")
  print("nonhdl")
  get_knot_loc("nonhdl")
  print("bmi")
  get_knot_loc("bmi")
  
}

