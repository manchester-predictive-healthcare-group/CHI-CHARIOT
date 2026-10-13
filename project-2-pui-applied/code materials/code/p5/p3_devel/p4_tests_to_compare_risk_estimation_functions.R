###
### Compare 10-year risks from full and redacted models.
### The full models, and function est_surv_offset_prototype3, are used in the model validation procedures.
### The redacted/reducted models, and functions est_surv_offset_prototype3_redacted/
### est_surv_offset_prototype3_redacted_mi will be exported and put into the rshiny. I therefore need to verify
### they return the same values.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Load survival package
library(survival)

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Match the ten-year horizon in internal validation.
t_eval_y <- 10L
t_eval <- round(t_eval_y * 365.25)
k <- 0L

###
### Tests for functions which estimate risk for individual models
###

### Initiate results vector
results <- vector("list", 20)

for (gender in 1:2) {
  
  ### Read the first imputed cohort and sample the same 1,000 patients for all ten fits.
  gender_char <- c("male", "female")[gender]
  df <- readRDS(paste0("data/p5/dfs_all_", gender_char, "1.rds"))
  set.seed(101)
  df_subset <- df[sample(nrow(df), 1000), ]
  rm(df)

  ### Match the zero-offset prediction strategy in the validation and instability scripts.
  df_subset$offset_statins_timevar_lnHR <- 0
  df_subset$offset_ah_timevar_lnHR <- 0
  df_subset$offset_smoking_timevar_dummy1_lnHR <- 0
  df_subset$offset_smoking_timevar_dummy2_lnHR <- 0

  ### Read reduced fits, annual hazards, and offset means for this sex.
  small_fits <- readRDS(paste0("data/p5/small_fits_model_", gender_char, ".rds"))
  small_bhaz <- readRDS(paste0("data/p5/small_bhaz_model_", gender_char, ".rds"))
  offset_means <- readRDS(paste0("data/p5/small_means_model_", gender_char, ".rds"))

  for (imp_num in c(1,4,9)) {
    
    ### Calculate risk from the original fit using est_surv_offset_prototype3.
    fit <- readRDS(paste0("data/p5/model_cox_", gender, "_imp", imp_num, ".rds"))
    bhaz <- readRDS(paste0("data/p5/model_cox_bhaz_", gender, "_imp", imp_num, ".rds"))
    risk_full <- 1 - as.numeric(est_surv_offset_prototype3(
      newdata = df_subset, fit = fit, bhaz = bhaz, time = t_eval
    ))
    
    ### Calculate risk from redacted model using est_surv_offset_prototype3_redacted
    risk_redacted <- 1 - as.numeric(est_surv_offset_prototype3_redacted(
      newdata = df_subset, fit = small_fits[[imp_num]], bhaz = small_bhaz[[imp_num]][[t_eval_y]], offset_means = offset_means[[imp_num]]
    ))
    
    ### Extract hazard, for check whether the stored annual hazard uses the same as what we are extracting in surv_offset_prototype3
    full_hazard <- bhaz$hazard[max(which(bhaz$time <= t_eval))]
    rm(fit, bhaz)

    ### Report agreement for this imputation, without averaging across models.
    k <- k + 1L
    results[[k]] <- data.frame(
      gender = gender_char,
      imputation = imp_num,
      hazard_difference = small_bhaz[[imp_num]][[t_eval_y]] - full_hazard,
      max_abs_difference = max(abs(risk_full - risk_redacted)),
      equal = isTRUE(all.equal(risk_full, risk_redacted,
                              tolerance = 1e-10, check.attributes = FALSE))
    )
  }
}

### Save
saveRDS(results, "data/p5/test_risk_estimation_redcated_models.rds")
### Print the result for every matched pair of models.
print(do.call(rbind, results), row.names = FALSE)


###
### Tests for functions which estimate risk averaged across different imputation streams
###

### Initiate results vector
results_mi <- vector("list", 2)

for (gender in 1:2) {
  
  ### Read the first imputed cohort and sample the same 1,000 patients for all ten fits.
  gender_char <- c("male", "female")[gender]
  df <- readRDS(paste0("data/p5/dfs_all_", gender_char, "1.rds"))
  set.seed(101)
  df_subset <- df[sample(nrow(df), 1000), ]
  rm(df)
  
  ### Match the zero-offset prediction strategy in the validation and instability scripts.
  df_subset$offset_statins_timevar_lnHR <- 0
  df_subset$offset_ah_timevar_lnHR <- 0
  df_subset$offset_smoking_timevar_dummy1_lnHR <- 0
  df_subset$offset_smoking_timevar_dummy2_lnHR <- 0
  
  ### Read reduced fits, annual hazards, and offset means for this sex.
  small_fits <- readRDS(paste0("data/p5/small_fits_model_", gender_char, ".rds"))
  small_bhaz <- readRDS(paste0("data/p5/small_bhaz_model_", gender_char, ".rds"))
  offset_means <- readRDS(paste0("data/p5/small_means_model_", gender_char, ".rds"))
  
  ### Read in full non-redacted models
  ### Calculate risk from the original fit using est_surv_offset_prototype3.
  fit_list <- lapply(1:10, function(x){readRDS(paste0("data/p5/model_cox_", gender, "_imp", x, ".rds"))})
  bhaz_list <- lapply(1:10, function(x){readRDS(paste0("data/p5/model_cox_bhaz_", gender, "_imp", x, ".rds"))})
    
    ### Calculate risk from the original fit using est_surv_offset_prototype3.
    ### Note, I haven't written an mi function to do this, as in the validation programs we validated each
    ### imputation stream seperately, and summarised results using robust metrics.
    ### I will therefore do the averaging manually here.
    surv_full_list <- lapply(1:10, function(x){est_surv_offset_prototype3(
      newdata = df_subset, fit = fit_list[[x]], bhaz = bhaz_list[[x]], time = t_eval)})
    
    ### Now get average on cloglog scale
    risk_full <- 1 - exp(-exp(rowMeans(do.call("cbind", lapply(surv_full_list, function(x) {log(-log(x))})))))
    
    ### Extract bhaz at 10 years from small_bhaz (which are currently vectors at years 1:10)
    small_bhaz_10y <- lapply(small_bhaz, function(x){x[[t_eval_y]]})
    
    ### Calculate risk from redacted model using est_surv_offset_prototype3_redacted
    risk_redacted <- 1 - as.numeric(est_surv_offset_prototype3_redacted_mi(
      newdata = df_subset, fit_list = small_fits, bhaz_list = small_bhaz_10y, offset_means_list = offset_means))
    
    ### Report agreement for this imputation, without averaging across models.
    results_mi[[gender]] <- data.frame(
      gender = gender_char,
      max_abs_difference = max(abs(risk_full - risk_redacted)),
      equal = isTRUE(all.equal(risk_full, risk_redacted,
                               tolerance = 1e-10, check.attributes = FALSE))
    )
}

### Save
saveRDS(results_mi, "data/p5/test_risk_estimation_redcated_models_mi.rds")
### Print the result for every matched pair of models.
print(do.call(rbind, results_mi), row.names = FALSE)