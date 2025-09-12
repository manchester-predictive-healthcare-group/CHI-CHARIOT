###
### Program to create validation cohorts that will be used for the temporal validation.
### Validation will be assessed in cohorts with an index date defined 1/2/3/4/5 years post start off followup
###
### We are not imputing each of the temporal validation datasets seperately (for time). We therefore have to 
### follow a process:
### 1) For predictors that have no missing data, we used the value from the cohort extracted at t_fup years.
### 2) For ethnicity and IMD (which do have missingness) we just use the value extracted (or imputed) at the 
### start of follow-up, as it cannot change.
### 3) For the other predictors that have missingness, we use the value from the cohort extracted at t_fup years 
### if it is recorded. If it is missing, then we use the value from start of follow-up (this comes from the imputed 
### dataset, and may have been either be imputed or observed at baseline). It's possible that a test value was
### recorded at baseline, and not in the temporal validation cohort, as we only look back five years for test data.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Load survival package
library(survival)
library(foreach)
library(doParallel)
library(doFuture)

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Write a function to create the dataset
create_df_temporalv <- function(gender, t_index){
  
  ### Define gender character version
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  
  ### Define the time point
  t_fup <- round(365.25*as.numeric(t_index))
  print(paste("t_fup = ", t_fup))
  
  ### Read in validation datasets
  df_imp <- readRDS(paste("data/p4/dfs_valid_", gender_char, ".rds", sep = ""))[[1]]
  
  ### Reduce and rename variables in the imputed dataset
  ### NB: These are denoted "_imp" to indicate the dataset they came from, but not all these values were imputed,
  ### however these variables were all defined at the start of follow-up
  df_imp <- dplyr::select(df_imp, -c(pracid, gender, fup_start, fup_end, 
                                     cvd_time_prim, cvd_indicator_prim, 
                                     cvd_time_hes, cvd_indicator_hes, 
                                     cvd_time_death, cvd_indicator_death, cvd_ev_prim_aj)) |>
    dplyr::relocate(patid, cvd_time, cvd_indicator)
  
  colnames(df_imp)[-c(1,2,3)] <- paste(colnames(df_imp)[-c(1,2,3)], "_imp", sep = "")
  
  ### Reduce this to individuals who aren't censored or already had an event
  df_imp <- subset(df_imp, cvd_time > t_fup)
  
  ### Create cvd_time relative to new index date
  df_imp <- dplyr::mutate(df_imp, cvd_time = cvd_time - t_fup)
  
  ### Read in cohort with data extracted at t_fup days after start of follow-up.
  df_index <- readRDS(paste("../../Aurum_Jun2021_extract/data/extraction/cohort_baseline_followup/cohort_var", t_fup, ".rds", sep = ""))
  
  ### 
  ### Combine the Light, Moderate and Heavy into a current smoker category
  ###
  
  ### Redefine smoking to be three levels, Never, Ex, Current
  df_index[,paste("smoking", "_t", t_fup, sep = "")] <- 
    forcats::fct_recode(df_index[,paste("smoking", "_t", t_fup, sep = "")], Current = "Light", Current = "Moderate", Current = "Heavy")
  
  ### Merge, with the imputed data (note, this will reduce to just those in the validation cohort)
  df_valid <- dplyr::left_join(df_imp, df_index, by = dplyr::join_by(patid))
  testthat::expect_equal(nrow(df_valid), nrow(df_imp))
  
  ###
  ### Now want to extract the appropriate values for prediction.
  ###
  
  ###
  ### For predictors that have no missing data, we used the value from the cohort extracted at t_fup years.
  
  ### Read in list of variables included in the model
  var_all_model <- readRDS(paste("data/p4/var_model", gender, ".rds"))
  
  ### Remove variables that had missingness (ethnicity, IMD, smoking status, BMI, SBP, non-HDL cholesterol).
  var_all_model_nomiss <- var_all_model[!grepl("IMD|ethnicity|bmi|sbp|nonhdl|smoking", var_all_model)]
  
  ### For each predictor, assign the appropriate value
  for (predictor in var_all_model_nomiss){
    df_valid[, predictor] <- df_valid[, paste(predictor, "_t", t_fup, sep = "")]
  }
  
  ###
  ### For ethnicity and IMD we just use the value extracted (or imputed) at the start of follow-up, as it cannot change.
  df_valid <- dplyr::mutate(df_valid, 
                            ethnicity = ethnicity_imp,
                            IMD = IMD_imp)
  
  ### For the other predictors that have missingness, we use the value from the cohort extracted at t_fup years if it is recorded. 
  ### If it is missing, then we use the value from start of follow-up (this comes from the imputed dataset, and may either be imputed or observed).
  ###
  for (predictor in c("bmi", "sbp", "nonhdl", "smoking")){
    df_valid[, predictor] <- df_valid[, paste(predictor, "_t", t_fup, sep = "")]
    df_valid[is.na(df_valid[, predictor]), predictor] <- df_valid[is.na(df_valid[, predictor]), paste(predictor, "_imp", sep = "")]
  }
  
  ### Turn appropriate variables into factors
  df_valid <- dplyr::mutate(df_valid,
                            hypertension = factor(hypertension, labels = c("Absent", "Present")),
                            ra = factor(ra, labels = c("Absent", "Present")),
                            af = factor(af, labels = c("Absent", "Present")),
                            ckd = factor(ckd, labels = c("Absent", "Present")),
                            smi = factor(smi, labels = c("Absent", "Present")),
                            fhcvd = factor(fhcvd, labels = c("Absent", "Present")),
                            migraine = factor(migraine, labels = c("Absent", "Present")),
                            sle = factor(sle, labels = c("Absent", "Present")),
                            cortico = factor(cortico, labels = c("Absent", "Present")),
                            antipsy = factor(antipsy, labels = c("Absent", "Present")),
                            copd = factor(copd, labels = c("Absent", "Present")),
                            int_dis = factor(int_dis, labels = c("Absent", "Present")),
                            downs = factor(downs, labels = c("Absent", "Present")),
                            oral_cancer = factor(oral_cancer, labels = c("Absent", "Present")),
                            lung_cancer = factor(lung_cancer, labels = c("Absent", "Present")),
                            blood_cancer = factor(blood_cancer, labels = c("Absent", "Present")),
                            brain_cancer = factor(brain_cancer, labels = c("Absent", "Present")))
  if (gender == 1){
    df_valid <- dplyr::mutate(df_valid,
                              impotence = factor(impotence, labels = c("Absent", "Present")))
  } else if (gender == 2){
    df_valid <- dplyr::mutate(df_valid,
                              pre_eclampsia = factor(pre_eclampsia, labels = c("Absent", "Present")),
                              postnatal_depression = factor(postnatal_depression, labels = c("Absent", "Present")))
  }
  
  ### Finally, add calendar time at the individuals index date
  df_valid <- dplyr::mutate(df_valid, caltime = as.numeric(fup_start) - as.numeric(as.Date("01/01/2005", format = "%d/%m/%Y")) + t_fup)
  
  ### Our data is now set up with the appropriate predictors with names and format matching that of the predictors in the development cohort
  ### We can therefore produce risk scores for the cohorts with index dates at t_fup post baseline.
  
  ### This data will also be necessary for estimating the counterfactual survival times, as we need to know smoking status at the index date
  
  ### Reduce and save
  df_valid <- dplyr::select(df_valid, dplyr::all_of(c("patid", "gender", "cvd_time", "cvd_indicator", "caltime", var_all_model)))
  saveRDS(df_valid, paste("data/p4/df_valid_temporalv_", gender, "_", t_fup, ".rds", sep = ""))
  print(str(df_valid))
  
}

### Run the function
for (gender_in in 1:2){
  for (t_index_in in 1:5){
    create_df_temporalv(gender_in, t_index_in)
  }
}
print("FINISHED")

###
### NB: XXXX
### ADD A CHECK FOR SMOKING STATUS HERE.
###