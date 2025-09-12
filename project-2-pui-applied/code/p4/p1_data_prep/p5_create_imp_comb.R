###
### Combine imputed datasets into a single object
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

library(mice)
library(survival)
library(rms)

# ### Save a smaller mice_mids object fo testing code on
# for (gender in c(2,1)){
#   
#   gender_char <- c("male", "female")[gender]
#   print(paste("gender = ", gender_char))
#   
#   ### Read in the first imputed dataset
#   imp.comb <- readRDS(paste("data/p4/mice_mids_prototype3_", gender, "_", 1, ".rds", sep = ""))
#   
#   ### Create one large mids object
#   for (i in c(2:2)){
#     print(paste("COMBINE", i, Sys.time()))
#     imp.new <- readRDS(paste("data/p4/mice_mids_prototype3_", gender, "_", i, ".rds", sep = ""))
#     imp.comb <- ibind(imp.comb, imp.new)
#   }
#   
#   ### Save the combined imputation object
#   saveRDS(imp.comb, paste("data/p4/imp_comb_prototype3_small", gender, ".rds", sep = ""))
#   print(paste("FINISHED", gender, Sys.time()))
#   
# }

###
### Combine and save all the imputed datasets
###
for (gender in c(2,1)){
  
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  
  ### Read in the first imputed dataset
  imp.comb <- readRDS(paste("data/p4/mice_mids_prototype3_", gender, "_", 1, ".rds", sep = ""))
  
  ### Create one large mids object
  for (i in c(2:10)){
    print(paste("COMBINE", i, Sys.time()))
    imp.new <- readRDS(paste("data/p4/mice_mids_prototype3_", gender, "_", i, ".rds", sep = ""))
    imp.comb <- ibind(imp.comb, imp.new)
  }
  
  ### Save the combined imputation object
  saveRDS(imp.comb, paste("data/p4/imp_comb_prototype3", gender, ".rds", sep = ""))
  print(paste("FINISHED", gender, Sys.time()))
  
}

print(paste("FINISHED", Sys.time()))

# ###
# ### Imputed datasets at 20 iterations
# ###
# for (gender in c(2,1)){
#   
#   gender_char <- c("male", "female")[gender]
#   print(paste("gender = ", gender_char))
#   
#   ### Read in the first imputed dataset
#   imp.comb <- readRDS(paste("data/p4/mice_mids_prototype3_", gender, "_", 1, "_20iter.rds", sep = ""))
#   
#   ### Create one large mids object
#   for (i in c(2:10)){
#     print(paste("COMBINE", i, Sys.time()))
#     imp.new <- readRDS(paste("data/p4/mice_mids_prototype3_", gender, "_", i, "_20iter.rds", sep = ""))
#     imp.comb <- ibind(imp.comb, imp.new)
#   }
#   
#   ### Save the combined imputation object
#   saveRDS(imp.comb, paste("data/p4/imp_comb_prototype3", gender, "_20iter.rds", sep = ""))
#   print(paste("FINISHED", gender, Sys.time()))
#   
# }

print(paste("FINISHED", Sys.time()))