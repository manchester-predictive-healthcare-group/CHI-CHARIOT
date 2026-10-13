###
### Sensitivity analysis replacing SBP, BMI and non-HDL cholesterol
### with age-stratified mean values
###
### Uses model and data from imputation stream 1 only
### Parallelised over sbp, nonhdl and bmi
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

### Load packages
library(survival)
library(doParallel)
library(foreach)

### Define time to evaluate performance at (10-years)
t_eval <- round(10*365.25)
### Define variables for sensitivity analyses
variables <- c("sbp", "nonhdl", "bmi")

###
### Run sensitivity analysis for one variable
###
run_sensitivity_analysis <- function(gender, variable){
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
  
  ### Replace selected variable with age-stratified mean
  data_valid[[variable]] <- data_valid$mean_value
  ### Remove temporary variables
  data_valid <- dplyr::select(
    data_valid,
    -age_band_lower,
    -age_band,
    -mean_value
  )
  ### Create variables for prediction under prescribed treatment strategy
  data_valid_pred <- dplyr::mutate(
    data_valid,
    offset_statins_timevar_lnHR = 0,
    offset_ah_timevar_lnHR = 0,
    offset_smoking_timevar_dummy1_lnHR = 0,
    offset_smoking_timevar_dummy2_lnHR = 0
  )
  ### Ensure data are sorted in the correct order
  data_valid_pred <- dplyr::arrange(data_valid_pred, patid)
  ### Estimate survival probabilities
  data_valid_pred$surv <- as.numeric(
    est_surv_offset_prototype3(
      newdata = data_valid_pred,
      fit = fit,
      bhaz = bhaz,
      time = t_eval
    )
  )
  ####################################
  ### Calibration in entire cohort ###
  ####################################
  ### PH regression approach
  print(paste("PH", variable, Sys.time()))
  calib_ph <- est_calib_ph(
    data = data_valid_pred,
    fit = fit,
    bhaz = bhaz,
    time = t_eval,
    surv = data_valid_pred$surv,
    pred_plot_range = seq(0.0001, 0.4, length.out = 250)
  )
  ### Grouped Kaplan-Meier approach
  print(paste("KM group", variable, Sys.time()))
  calib_km_group <- est_calib_plot_group(
    data = data_valid_pred,
    fit = fit,
    bhaz = bhaz,
    time = t_eval,
    n.groups = 100,
    surv = data_valid_pred$surv,
    single.CI = TRUE
  )
  ### Extract grouped Kaplan-Meier data
  calib_km_group_data <- calib_km_group$plot.data
  ### Save complete calibration objects
  saveRDS(calib_ph, paste0("data/p5/sens_mean_imputation_calib_ph_", gender, "_", variable, ".rds"))
  saveRDS(calib_km_group, paste0("data/p5/sens_mean_imputation_km_group_", gender, "_", variable, ".rds"))
  ######################
  ### Discrimination ###
  ######################
  ### Estimate C-statistic
  Cstat_object <- survival::concordance(
    survival::Surv(cvd_time, cvd_indicator) ~ surv,
    data = data_valid_pred,
    timewt = "n",
    influence = 0
  )
  ### Save discrimination object
  saveRDS(Cstat_object, paste0("data/p5/sens_mean_imputation_discrim_", gender, "_", variable, ".rds"))
  #########################
  ### Performance table ###
  #########################
  ### Create performance table
  performance_table <- data.frame(
    gender = gender_char,
    variable = variable,
    ICI = calib_ph$ICI,
    E50 = calib_ph$E50,
    E90 = calib_ph$E90,
    discrimination = unname(Cstat_object$concordance)
  )
  ### Save performance table
  saveRDS(performance_table,paste0("data/p5/sens_mean_imputation_table_", gender, "_", variable, ".rds"))
  ####################
  ### Create plots ###
  ####################
  ### Create PH calibration curve from underlying data
  plot_calib_ph <- ggplot2::ggplot(
    calib_ph$calib_data_plot,
    ggplot2::aes(x = pred, y = pred.obs)
  ) +
    ggplot2::geom_line(color = "red", linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::coord_equal(xlim = c(0, 0.4), ylim = c(0, 0.4)) +
    ggplot2::ggtitle("Proportional hazards calibration") +
    ggplot2::xlab("Predicted risk") +
    ggplot2::ylab("Predicted-observed risk") +
    ggplot2::theme_minimal()
  ### Create grouped Kaplan-Meier plot from underlying data
  plot_calib_km <- ggplot2::ggplot(
    calib_km_group_data,
    ggplot2::aes(x = pred, y = obs)
  ) +
    ggplot2::geom_point(color = "red") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = obs_lower, ymax = obs_upper),
      width = 0.005
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::coord_equal(xlim = c(0, 0.4), ylim = c(0, 0.4)) +
    ggplot2::ggtitle("Grouped Kaplan-Meier calibration (100 groups)") +
    ggplot2::xlab("Predicted risk") +
    ggplot2::ylab("Observed risk") +
    ggplot2::theme_minimal()
  ### Define output filename
  plot_filename <- paste0("figures/p5/sens_mean_imputation_plots_", gender, "_", variable, ".png")
  ### Check output directory exists and is writable
  if (!dir.exists("figures/p5")){
    stop("Output directory does not exist: figures/p5")
  }
  if (file.access("figures/p5", 2) != 0){
    stop("Output directory is not writable: figures/p5")
  }
  ### Save both plots in one PNG using ragg
  ragg::agg_png(
    filename = plot_filename,
    width = 3200,
    height = 1500,
    units = "px",
    res = 250
  )
  ### Ensure device closes if plotting fails
  tryCatch(
    {
      grid::grid.newpage()
      grid::pushViewport(
        grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2))
      )
      print(
        plot_calib_ph,
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
      )
      print(
        plot_calib_km,
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2)
      )
    },
    error = function(e){
      stop(
        paste0(
          "Plotting failed for gender = ", gender,
          ", variable = ", variable,
          ": ", conditionMessage(e)
        )
      )
    },
    finally = {
      grDevices::dev.off()
    }
  )
  ### Return performance table
  return(performance_table)
}
###
### Run sensitivity analyses for a given gender
###
run_gender <- function(gender){
  ### Set up parallel backend with one worker per variable
  n_cores <- min(length(variables), parallel::detectCores())
  print(paste("Parallel workers =", n_cores))
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  ### Run sensitivity analyses in parallel over variables
  results <- tryCatch(
    foreach::foreach(
      variable = variables,
      .packages = c("dplyr", "survival", "rms", "ggplot2", "ggExtra", "ragg"),
      .export = c(
        "run_sensitivity_analysis",
        "est_surv",
        "est_surv_offset_prototype3",
        "est_calib_ph",
        "est_calib_plot_group",
        "t_eval"
      )
    ) %dopar% {
      run_sensitivity_analysis(gender, variable)
    },
    finally = {
      parallel::stopCluster(cl)
    }
  )
  ### Combine performance tables across variables
  performance_table_comb <- dplyr::bind_rows(results)
  print(paste("FINISHED", c("male", "female")[gender], Sys.time()))
  ### Return combined table
  return(performance_table_comb)
}
### Run for both genders
performance_male <- run_gender(1)
performance_female <- run_gender(2)
### Combine performance tables across genders
performance_all <- dplyr::bind_rows(performance_female, performance_male)
### Save combined performance table
saveRDS(performance_all,"data/p5/sens_mean_imputation_table.rds")
print(performance_all)
print(paste("FINISHED ALL", Sys.time()))
