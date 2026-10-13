###
### MNAR sensitivity analysis for internal validation
###
### Applies delta adjustments to originally missing SBP, BMI and non-HDL values
### in imputation stream 1, while keeping the developed model fixed. The raw adjusted predictor is
### passed directly to predict() for the existing prototype 3 model.
###
### This program creates and saves results only. It does not create figures and
### does not save the complete MNAR-adjusted datasets.
###

### Clear workspace
rm(list = ls())
Sys.time()

### Set working directory
setwd()
getwd()

### Source functions
R.func.sources <- list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load packages
library(survival)
library(doParallel)
library(foreach)

### Define the imputation and prediction horizon
imp_num <- 1
t_eval <- round(10 * 365.25)

### Define variables included in the MNAR sensitivity analysis
variables <- c("sbp", "bmi", "nonhdl")

### Define delta values in observed standard deviation units
### Delta = 0 is the MAR reference scenario
delta_values <- c(-0.50, -0.25, 0.25, 0.50)

### Apply one MNAR scenario and estimate performance
run_mnar_scenario <- function(gender, variable, delta){

  ### Define gender and scenario labels
  gender_char <- c("male", "female")[gender]
  delta_label <- ifelse(
    delta < 0,
    paste0("minus_", abs(delta)),
    paste0("plus_", delta)
  )
  delta_label <- gsub("\\.", "p", delta_label)

  print(paste(
    "gender =", gender_char,
    "variable =", variable,
    "delta =", delta,
    Sys.time()
  ))

  ### Read the first fully imputed dataset
  data_valid <- readRDS(
    paste0("data/p5/dfs_all_", gender_char, imp_num, ".rds")
  )

  ### Read the fixed model and baseline hazard from imputation stream 1
  fit <- readRDS(
    paste0("data/p5/model_cox_", gender, "_imp", imp_num, ".rds")
  )
  bhaz <- readRDS(
    paste0("data/p5/model_cox_bhaz_", gender, "_imp", imp_num, ".rds")
  )

  ### Read the original pre-imputation cohort
  cohort_original <- readRDS("data/p5/cohort_prototype3.rds")

  ### Restrict the original cohort to the gender being analysed
  gender_num <- gender
  cohort_original <- subset(cohort_original, gender == gender_num)

  ### Check identifiers are unique before joining missingness information
  stopifnot(!anyDuplicated(data_valid$patid))
  stopifnot(!anyDuplicated(cohort_original$patid))

  ### Create the original missingness indicator for the selected variable
  missingness_data <- data.frame(
    patid = cohort_original$patid,
    originally_missing = is.na(cohort_original[[variable]])
  )

  ### Add the original missingness indicator to the completed dataset
  n_before <- nrow(data_valid)
  data_valid <- dplyr::left_join(
    data_valid,
    missingness_data,
    by = "patid"
  )
  stopifnot(nrow(data_valid) == n_before)
  stopifnot(!anyNA(data_valid$originally_missing))

  ### Calculate the SD from values that were observed before imputation
  observed_values <- cohort_original[[variable]][
    !is.na(cohort_original[[variable]])
  ]
  observed_sd <- stats::sd(observed_values)
  stopifnot(is.finite(observed_sd), observed_sd > 0)

  ### Store values before applying the MNAR adjustment
  values_before <- data_valid[[variable]]

  ### Apply the delta only to values that were originally missing
  data_valid[[variable]][data_valid$originally_missing] <-
    data_valid[[variable]][data_valid$originally_missing] +
    delta * observed_sd

  ### Check the adjustment has not created missing or non-finite values
  stopifnot(all(is.finite(data_valid[[variable]])))

  ### No derived variables are recreated here
  ### The prototype 3 fitted model uses the adjusted raw predictor directly
  ### when predict() is called inside est_surv_offset_prototype3()

  ### Summarise the applied adjustment without saving the full dataset
  adjustment_summary <- data.frame(
    gender = gender_char,
    imputation = imp_num,
    variable = variable,
    delta_sd = delta,
    observed_sd = observed_sd,
    delta_original_units = delta * observed_sd,
    N = nrow(data_valid),
    N_originally_missing = sum(data_valid$originally_missing),
    proportion_originally_missing = mean(data_valid$originally_missing),
    mean_before_missing = mean(
      values_before[data_valid$originally_missing]
    ),
    mean_after_missing = mean(
      data_valid[[variable]][data_valid$originally_missing]
    ),
    min_after_missing = min(
      data_valid[[variable]][data_valid$originally_missing]
    ),
    max_after_missing = max(
      data_valid[[variable]][data_valid$originally_missing]
    )
  )

  ### Read the existing counterfactual validation outcomes
  cf_surv_times <- readRDS(
    paste0("data/p5/cf_surv_times_", gender, "_imp", imp_num, ".rds")
  )
  stopifnot(!anyDuplicated(cf_surv_times$patid))

  ### Replace observed outcomes with the existing counterfactual outcomes
  data_valid <- dplyr::select(
    data_valid,
    -cvd_time,
    -cvd_indicator
  ) |>
    dplyr::left_join(cf_surv_times, by = "patid") |>
    dplyr::rename(
      cvd_time = cvd_time_cf,
      cvd_indicator = cvd_indicator_cf
    )

  stopifnot(nrow(data_valid) == n_before)
  stopifnot(!anyNA(data_valid$cvd_time))
  stopifnot(!anyNA(data_valid$cvd_indicator))

  ### Set time-varying offsets to the prescribed treatment strategy
  data_valid_pred <- dplyr::mutate(
    data_valid,
    offset_statins_timevar_lnHR = 0,
    offset_ah_timevar_lnHR = 0,
    offset_smoking_timevar_dummy1_lnHR = 0,
    offset_smoking_timevar_dummy2_lnHR = 0
  )

  ### Sort data before estimating predictions
  data_valid_pred <- dplyr::arrange(data_valid_pred, patid)

  ### Estimate 10-year survival using the existing fixed model
  data_valid_pred$surv <- as.numeric(
    est_surv_offset_prototype3(
      newdata = data_valid_pred,
      fit = fit,
      bhaz = bhaz,
      time = t_eval
    )
  )

  ###################
  ### Calibration ###
  ###################

  ### Estimate PH-based calibration
  print(paste("PH", gender_char, variable, delta, Sys.time()))
  calib_ph <- est_calib_ph(
    data = data_valid_pred,
    fit = fit,
    bhaz = bhaz,
    time = t_eval,
    surv = data_valid_pred$surv,
    pred_plot_range = seq(0.0001, 0.4, length.out = 250),
    plot = FALSE
  )

  ### Estimate grouped Kaplan-Meier calibration
  print(paste("KM group", gender_char, variable, delta, Sys.time()))
  calib_km_group <- est_calib_plot_group(
    data = data_valid_pred,
    fit = fit,
    bhaz = bhaz,
    time = t_eval,
    n.groups = 100,
    surv = data_valid_pred$surv,
    single.CI = TRUE
  )

  ######################
  ### Discrimination ###
  ######################

  ### Estimate the C-statistic
  Cstat_object <- survival::concordance(
    survival::Surv(cvd_time, cvd_indicator) ~ surv,
    data = data_valid_pred,
    timewt = "n",
    influence = 0
  )

  ####################
  ### Save results ###
  ####################

  ### Define the common filename suffix
  file_suffix <- paste0(
    gender, "_", variable, "_delta_", delta_label
  )

  ### Save compact adjustment information
  saveRDS(
    adjustment_summary,
    paste0(
      "data/p5/mnar_adjustment_summary_",
      file_suffix,
      ".rds"
    )
  )

  ### Save calibration results
  saveRDS(
    calib_ph,
    paste0("data/p5/mnar_calib_ph_", file_suffix, ".rds")
  )
  saveRDS(
    calib_km_group,
    paste0("data/p5/mnar_calib_km_group_", file_suffix, ".rds")
  )

  ### Save discrimination results
  discrimination_result <- data.frame(
    gender = gender_char,
    imputation = imp_num,
    variable = variable,
    delta_sd = delta,
    Cstat = unname(Cstat_object$concordance),
    Cstat_se = sqrt(unname(Cstat_object$var))
  )
  saveRDS(
    discrimination_result,
    paste0("data/p5/mnar_discrim_", file_suffix, ".rds")
  )

  ### Create compact performance summary
  performance_summary <- data.frame(
    gender = gender_char,
    imputation = imp_num,
    variable = variable,
    delta_sd = delta,
    delta_original_units = delta * observed_sd,
    N = nrow(data_valid_pred),
    N_events = sum(data_valid_pred$cvd_indicator),
    N_originally_missing = sum(data_valid_pred$originally_missing),
    mean_predicted_risk = mean(1 - data_valid_pred$surv),
    ICI = calib_ph$ICI,
    E50 = calib_ph$E50,
    E90 = calib_ph$E90,
    Cstat = unname(Cstat_object$concordance),
    Cstat_se = sqrt(unname(Cstat_object$var))
  )

  ### Save compact performance summary
  saveRDS(
    performance_summary,
    paste0("data/p5/mnar_performance_", file_suffix, ".rds")
  )

  ### Return summaries for the combined result files
  return(list(
    adjustment = adjustment_summary,
    performance = performance_summary
  ))
}

### Run all scenarios for one gender
run_gender <- function(gender){

  ### Create all combinations of variables and delta values
  scenarios <- expand.grid(
    variable = variables,
    delta = delta_values,
    stringsAsFactors = FALSE
  )

  ### Use one worker per variable
  n_cores <- min(length(variables), parallel::detectCores())
  print(paste("Parallel workers =", n_cores))

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  ### Run MNAR scenarios in parallel
  results <- tryCatch(
    foreach::foreach(
      scenario_num = seq_len(nrow(scenarios)),
      .packages = c("dplyr", "survival", "rms", "ggplot2", "ggExtra"),
      .export = c(
        "run_mnar_scenario",
        "est_surv",
        "est_surv_offset_prototype3",
        "est_calib_ph",
        "est_calib_plot_group",
        "imp_num",
        "t_eval"
      )
    ) %dopar% {
      run_mnar_scenario(
        gender = gender,
        variable = scenarios$variable[scenario_num],
        delta = scenarios$delta[scenario_num]
      )
    },
    finally = {
      parallel::stopCluster(cl)
    }
  )

  ### Combine scenario-level summaries
  adjustment_summary <- dplyr::bind_rows(
    lapply(results, function(x) x$adjustment)
  )
  performance_summary <- dplyr::bind_rows(
    lapply(results, function(x) x$performance)
  )

  return(list(
    adjustment = adjustment_summary,
    performance = performance_summary
  ))
}

### Run the analysis for males and females
results_male <- run_gender(1)
results_female <- run_gender(2)

### Combine summaries across genders
mnar_adjustment_summary <- dplyr::bind_rows(
  results_female$adjustment,
  results_male$adjustment
)
mnar_internal_validation_summary <- dplyr::bind_rows(
  results_female$performance,
  results_male$performance
)

### Add changes from the variable-specific MAR reference scenario
mnar_internal_validation_summary <- mnar_internal_validation_summary |>
  dplyr::ungroup() |>
  dplyr::arrange(gender, variable, delta_sd)

### Save combined summaries
saveRDS(
  mnar_adjustment_summary,
  "data/p5/mnar_adjustment_summary.rds"
)
saveRDS(
  mnar_internal_validation_summary,
  "data/p5/mnar_internal_validation_summary.rds"
)

print(mnar_internal_validation_summary)
print(paste("FINISHED ALL", Sys.time()))
