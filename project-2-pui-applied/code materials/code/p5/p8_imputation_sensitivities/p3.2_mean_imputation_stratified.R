###
### Stratified sensitivity analysis comparing age-stratified mean replacement
### with prediction using the imputed values
###
### Strata are defined by six equally spaced intervals of the imputed values
### Uses model and data from imputation stream 1 only
### Parallelised over sbp, nonhdl and bmi
###
### This is a separate program and does not alter pX_mean_imputation
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load packages
library(survival)
library(doParallel)
library(foreach)

### Define time to evaluate performance at (10-years)
t_eval <- round(10*365.25)

### Define variables for sensitivity analyses
variables <- c("sbp", "nonhdl", "bmi")

### Read each sex to define common subgroup boundaries
subgroup_data_male <- readRDS("data/p5/dfs_all_male1.rds")
subgroup_data_female <- readRDS("data/p5/dfs_all_female1.rds")

### Define identical subgroup boundaries for females and males
group_breaks <- lapply(variables, function(variable){
  
  ### Combine values across females and males
  pooled_values <- c(
    subgroup_data_male[[variable]],
    subgroup_data_female[[variable]]
  )
  
  ### Combine all non-HDL values over 6.33 into the highest group
  if (variable == "nonhdl"){
    breaks <- c(
      seq(
        min(pooled_values, na.rm = TRUE),
        6.33,
        length.out = 6
      ),
      Inf
    )
  } else {
    ### Create six equally spaced groups across the pooled observed range
    breaks <- seq(
      min(pooled_values, na.rm = TRUE),
      max(pooled_values, na.rm = TRUE),
      length.out = 7
    )
  }
  
  return(breaks)
})

### Name break vectors by variable
names(group_breaks) <- variables

### Remove data used only to define subgroup boundaries
rm(subgroup_data_male, subgroup_data_female)

###
### Run stratified analyses for one variable
###
run_stratified_analysis <- function(gender, variable, group_breaks){
  
  ### Prelim
  gender_char <- c("male", "female")[gender]
  imp_num <- 1
  print(paste("gender =", gender_char))
  print(paste("variable =", variable, Sys.time()))
  
  ### Pick the first imputed dataset
  data_valid <- readRDS(paste0("data/p5/dfs_all_", gender_char, imp_num, ".rds"))
  
  ### Read in model fit
  fit <- readRDS(paste0("data/p5/model_cox_", gender, "_imp", imp_num, ".rds"))
  bhaz <- readRDS(paste0("data/p5/model_cox_bhaz_", gender, "_imp", imp_num, ".rds"))
  
  ### Read in counterfactual survival times
  cf_surv_times <- readRDS(paste0("data/p5/cf_surv_times_", gender, "_imp", imp_num, ".rds"))
  stopifnot(!anyDuplicated(cf_surv_times$patid))
  
  ### Replace cvd_time and cvd_indicator with counterfactual survival times
  n_before <- nrow(data_valid)
  data_valid <- dplyr::select(data_valid, -cvd_time, -cvd_indicator) |>
    dplyr::left_join(cf_surv_times, by = "patid") |>
    dplyr::rename(cvd_time = cvd_time_cf,
                  cvd_indicator = cvd_indicator_cf)
  stopifnot(nrow(data_valid) == n_before)
  stopifnot(!anyNA(data_valid$cvd_time))
  stopifnot(!anyNA(data_valid$cvd_indicator))
  
  ### Save the original imputed values used to define subgroups
  data_valid$imputed_value <- data_valid[[variable]]
  
  ### Select common subgroup boundaries for this variable
  variable_group_breaks <- group_breaks[[variable]]
  
  ### Check six distinct groups have been defined
  if (length(unique(variable_group_breaks)) != 7){
    stop("Could not create six distinct intervals for ", variable)
  }
  
  ### Assign individuals to groups using the imputed values
  data_valid$imputed_value_group <- cut(
    data_valid$imputed_value,
    breaks = variable_group_breaks,
    include.lowest = TRUE,
    right = TRUE,
    ordered_result = TRUE
  )
  
  ### Check all individuals have been assigned to a group
  stopifnot(!anyNA(data_valid$imputed_value_group))
  
  ### Save subgroup definitions
  subgroup_definitions <- data.frame(
    gender = gender_char,
    variable = variable,
    group = levels(data_valid$imputed_value_group),
    lower = variable_group_breaks[1:6],
    upper = variable_group_breaks[2:7]
  ) |>
    dplyr::left_join(
      data_valid |>
        dplyr::count(imputed_value_group, name = "N") |>
        dplyr::mutate(group = as.character(imputed_value_group)) |>
        dplyr::select(group, N),
      by = "group"
    )
  
  saveRDS(
    subgroup_definitions,
    paste0("data/p5/sens_mean_imputation_stratified_groups_", gender, "_", variable, ".rds")
  )
  
  ### Define centred rolling five-year age windows
  data_valid <- dplyr::mutate(
    data_valid,
    age_band_lower = floor(age) - 2,
    age_band = paste0(age_band_lower, "-", age_band_lower + 4)
  )
  
  ### Summarise predictor values by completed year of age
  ### This avoids repeatedly scanning the entire cohort
  age_totals <- data_valid |>
    dplyr::group_by(age_year = floor(age)) |>
    dplyr::summarise(
      value_sum = sum(.data[[variable]], na.rm = TRUE),
      N = sum(!is.na(.data[[variable]])),
      .groups = "drop"
    )
  
  ### Calculate means within each overlapping five-year window
  age_stratified_means <- dplyr::bind_rows(
    lapply(sort(unique(data_valid$age_band_lower)), function(lower) {
      
      window <- age_totals[
        age_totals$age_year >= lower &
          age_totals$age_year <= lower + 4,
        ,
        drop = FALSE
      ]
      
      data.frame(
        variable = variable,
        age_band_lower = lower,
        age_band = paste0(lower, "-", lower + 4),
        mean_value = sum(window$value_sum) / sum(window$N),
        N = sum(window$N)
      )
    })
  )
  
  ### Save the rolling age-specific mean values
  saveRDS(
    age_stratified_means,
    paste0(
      "data/p5/sens_mean_imputation_values_",
      gender, "_", variable, ".rds"
    )
  )
  
  ### Add rolling mean values to validation data
  data_valid <- dplyr::left_join(
    data_valid,
    dplyr::select(age_stratified_means, age_band_lower, mean_value),
    by = "age_band_lower"
  )
  
  ### Create data using the original imputed values for prediction
  data_valid_imputed <- data_valid
  
  ### Create data using age-stratified mean values for prediction
  data_valid_mean <- data_valid
  data_valid_mean[[variable]] <- data_valid_mean$mean_value
  
  ### Remove temporary age-band variables
  data_valid_imputed <- dplyr::select(
    data_valid_imputed,
    -age_band_lower,
    -age_band,
    -mean_value
  )
  data_valid_mean <- dplyr::select(
    data_valid_mean,
    -age_band_lower,
    -age_band,
    -mean_value
  )
  
  ### Create variables for prediction using imputed values
  data_valid_imputed <- dplyr::mutate(
    data_valid_imputed,
    offset_statins_timevar_lnHR = 0,
    offset_ah_timevar_lnHR = 0,
    offset_smoking_timevar_dummy1_lnHR = 0,
    offset_smoking_timevar_dummy2_lnHR = 0
  ) |>
    dplyr::arrange(patid)
  
  ### Create variables for prediction using mean values
  data_valid_mean <- dplyr::mutate(
    data_valid_mean,
    offset_statins_timevar_lnHR = 0,
    offset_ah_timevar_lnHR = 0,
    offset_smoking_timevar_dummy1_lnHR = 0,
    offset_smoking_timevar_dummy2_lnHR = 0
  ) |>
    dplyr::arrange(patid)
  
  ### Estimate survival probabilities using imputed values
  data_valid_imputed$surv <- as.numeric(
    est_surv_offset_prototype3(
      newdata = data_valid_imputed,
      fit = fit,
      bhaz = bhaz,
      time = t_eval
    )
  )
  
  ### Estimate survival probabilities using mean values
  data_valid_mean$surv <- as.numeric(
    est_surv_offset_prototype3(
      newdata = data_valid_mean,
      fit = fit,
      bhaz = bhaz,
      time = t_eval
    )
  )
  
  ### Split both datasets using the same imputed-value subgroups
  group_levels <- levels(data_valid$imputed_value_group)
  
  data_imputed_list <- lapply(group_levels, function(x){
    data_valid_imputed[data_valid_imputed$imputed_value_group == x, ]
  })
  
  data_mean_list <- lapply(group_levels, function(x){
    data_valid_mean[data_valid_mean$imputed_value_group == x, ]
  })
  
  ### Define function to validate the six subgroups for one prediction approach
  validate_subgroups <- function(data_list, prediction_type){
    
    ### Estimate calibration within each subgroup
    calib_ph <- lapply(seq_along(data_list), function(x){
      est_calib_ph(
        data = data_list[[x]],
        fit = fit,
        bhaz = bhaz,
        time = t_eval,
        surv = data_list[[x]]$surv,
        plot = FALSE
      )
    })
    
    ### Estimate discrimination within each subgroup
    Cstat_object <- lapply(seq_along(data_list), function(x){
      survival::concordance(
        survival::Surv(cvd_time, cvd_indicator) ~ surv,
        data = data_list[[x]],
        timewt = "n",
        influence = 0
      )
    })
    
    ### Create performance table
    performance_table <- data.frame(
      gender = gender_char,
      variable = variable,
      prediction_type = prediction_type,
      group = group_levels,
      N = sapply(data_list, nrow),
      ICI = sapply(calib_ph, function(x) x$ICI),
      E50 = sapply(calib_ph, function(x) x$E50),
      E90 = sapply(calib_ph, function(x) x$E90),
      discrimination = sapply(Cstat_object, function(x) unname(x$concordance))
    )
    
    ### Return results
    return(list(
      calib_ph = calib_ph,
      Cstat_object = Cstat_object,
      performance_table = performance_table
    ))
  }
  
  ### Validate predictions using the imputed values
  results_imputed <- validate_subgroups(
    data_list = data_imputed_list,
    prediction_type = "Imputed values"
  )
  
  ### Validate predictions using age-stratified mean values
  results_mean <- validate_subgroups(
    data_list = data_mean_list,
    prediction_type = "Age-stratified mean values"
  )
  
  ### Combine performance tables
  performance_table <- dplyr::bind_rows(
    results_imputed$performance_table,
    results_mean$performance_table
  )
  
  ### Save calibration objects
  saveRDS(
    results_imputed$calib_ph,
    paste0("data/p5/sens_mean_imputation_stratified_calib_ph_imputed_", gender, "_", variable, ".rds")
  )
  saveRDS(
    results_mean$calib_ph,
    paste0("data/p5/sens_mean_imputation_stratified_calib_ph_mean_", gender, "_", variable, ".rds")
  )
  
  ### Save discrimination objects
  saveRDS(
    results_imputed$Cstat_object,
    paste0("data/p5/sens_mean_imputation_stratified_discrim_imputed_", gender, "_", variable, ".rds")
  )
  saveRDS(
    results_mean$Cstat_object,
    paste0("data/p5/sens_mean_imputation_stratified_discrim_mean_", gender, "_", variable, ".rds")
  )
  
  ### Save performance table
  saveRDS(
    performance_table,
    paste0("data/p5/sens_mean_imputation_stratified_table_", gender, "_", variable, ".rds")
  )
  
  ### Return performance table
  return(performance_table)
}

###
### Run stratified analyses for a given gender
###
run_gender_stratified <- function(gender){
  
  ### Set up parallel backend with one worker per variable
  n_cores <- min(length(variables), parallel::detectCores())
  print(paste("Parallel workers =", n_cores))
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  ### Run stratified analyses in parallel over variables
  results <- tryCatch(
    foreach::foreach(
      variable = variables,
      .packages = c("dplyr", "survival", "rms", "ggplot2", "ggExtra"),
      .export = c(
        "run_stratified_analysis",
        "est_surv",
        "est_surv_offset_prototype3",
        "est_calib_ph",
        "t_eval",
        "group_breaks"
      )
    ) %dopar% {
      run_stratified_analysis(gender, variable, group_breaks)
    },
    finally = {
      parallel::stopCluster(cl)
    }
  )
  
  ### Combine performance tables across variables
  performance_table_comb <- dplyr::bind_rows(results)
  
  ### Save combined table for this gender
  saveRDS(
    performance_table_comb,
    paste0("data/p5/sens_mean_imputation_stratified_table_", gender, ".rds")
  )
  
  print(paste("FINISHED", c("male", "female")[gender], Sys.time()))
  
  ### Return combined table
  return(performance_table_comb)
}

### Run for both genders
performance_male <- run_gender_stratified(1)
performance_female <- run_gender_stratified(2)

### Combine performance tables across genders
performance_all <- dplyr::bind_rows(performance_male, performance_female)

### Save combined performance table
saveRDS(
  performance_all,
  "data/p5/sens_mean_imputation_stratified_table.rds"
)

print(performance_all)
print(paste("FINISHED ALL", Sys.time()))
