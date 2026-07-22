### Going to calculate sample sizes for each of the temporal validation cohorts, so I don't have to keep 
### re-running code to get these, which involves reading in of large datasets

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

###
### Get size of female and male cohorts
###
nrow_valid_female <- readRDS(paste("data/p4/dfs_valid_", "female", ".rds", sep = ""))[[1]] |> nrow()
nrow_valid_male <- readRDS(paste("data/p4/dfs_valid_", "male", ".rds", sep = ""))[[1]] |> nrow()

### Save these
saveRDS(nrow_valid_female, "data/p4/validation_cohort_N_female.rds")
saveRDS(nrow_valid_male,"data/p4/validation_cohort_N_male.rds")

###
### Get size of cohorts used for validation of initial risk estimation layer
###
N_initial_female <- do.call("c", lapply(round(365.25*(1:5)), function(x) {nrow(readRDS(paste("data/p4/df_valid_temporalv_", 2, "_", x, ".rds", sep = "")))}))
N_initial_male <- do.call("c", lapply(round(365.25*(1:5)), function(x) {nrow(readRDS(paste("data/p4/df_valid_temporalv_", 1, "_", x, ".rds", sep = "")))}))
N_initial <- data.frame(gender = c(rep(2, 5), rep(1, 5)),
                        N = c(N_initial_female, N_initial_male))
N_initial

###
### Get size of cohorts used for validation of interventional layer (non-missing data on mfr of interest)
###

### Read in the entire unimputed cohort at baseline (this is both development and validation datasets)
entire_cohort_baseline <- readRDS("data/p4/cohort_prototype3.rds")

### Write a function to get size of cohort at time t_fup, with non-missing data at baseline/visit 0,
### for a given gender and modifiable risk factor
get_N_interventional <- function(gender, mfr){
  
  N_interventional <- lapply(round(365.25*(1:5)), function(x) {
    
    ### Read in temporal validation dataset
    df_valid_temporal <- readRDS(paste("data/p4/df_valid_temporalv_", gender, "_", x, ".rds", sep = ""))
    
    ### Reduce entire_cohort_baseline to those in the validation dataset through a merge
    df_valid_baseline <- dplyr::left_join(dplyr::select(df_valid_temporal, patid), 
                                          entire_cohort_baseline, 
                                          by = dplyr::join_by(patid))
    
    ### Now reduce to those who don't have missing data at baseline on the variable of interest
    df_valid_baseline <- df_valid_baseline[!is.na(df_valid_baseline[[mfr]]), ]
    
    ### Reduce temporal validation cohort to those who had non-missing data at baseline
    df_valid_temporal <- df_valid_temporal[!is.na(fastmatch::fmatch(df_valid_temporal$patid, df_valid_baseline$patid)), ]
    
    ### Get number of individuals
    n <- nrow(df_valid_temporal)
    
    return(n)
    
  })
  
  N_interventional <- do.call("c", N_interventional)
  
  return(N_interventional)
  
}

### Apply the above function over each of the modifiable risk factors
N_interventional <- lapply(c("sbp", "bmi", "nonhdl", "smoking"), function(x) {
  # Run function female cohort
  out_female <- get_N_interventional(2, x)
  # Run function male cohort
  out_male <- get_N_interventional(1, x)
  # Combine and create data frame
  out <- data.frame(N = c(out_female, out_male)) |>
    dplyr::mutate(mfr = x,
                  gender = c(rep(2, 5), rep(1, 5))) |>
    dplyr::relocate(mfr, gender)
  return(out)
})

### Combine these into a single dataframe
N_interventional <- do.call("rbind", N_interventional)

### Save these
N_initial
N_interventional
saveRDS(N_initial, "data/p4/temporalv_N_initial.rds")
saveRDS(N_interventional,"data/p4/temporalv_N_interventional.rds")


