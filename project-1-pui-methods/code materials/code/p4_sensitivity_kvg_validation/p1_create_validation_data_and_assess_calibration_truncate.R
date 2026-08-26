###
### Estimate the IPACW's for K&VG validation approach
### Truncate the combined weights at 100
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load survival
library(survival)

### Assign burnout
burnout <- 180

### And follow up time
t_fup <- 10*365.25

### Function to assess calibration in truncated weights cohort
calibrate_truncate <- function(gender_in, truncate){
  
  print(paste("gender = ", gender_in, ", truncate = ", truncate, Sys.time()))
  
  ### Read cohort
  df_valid_ipacw <- readRDS(paste("data/df_valid_ipacw", gender_in, ".rds"))
  
  ### Truncate the weights
  df_valid_ipacw$comb_weight <- pmin(truncate, df_valid_ipacw$comb_weight)
  
  ### Create an uncensored at t_fup cohort
  df_valid_ipacw_uncensored <- subset(df_valid_ipacw, 
                                      time >= t_fup | status == 1 & time < t_fup)
  
  ###
  ### Mean calibration using Horvitz-Thompson estimator
  ###
  ht_estimate <- mean(
    (df_valid_ipacw$status == 1 & df_valid_ipacw$time <= t_fup) * 
      (df_valid_ipacw$comb_weight)
  )
  ht_estimate
  mean(df_valid_ipacw$pred)
  mean(df_valid_ipacw$pred)/ht_estimate
  # This is better, but still doesn't work!!
  
  ### Calculate calibration plot using IPCW approach in uncensored cohort
  calib_ipcw_uncensored <- est_calib_ipcw_new(data = df_valid_ipacw_uncensored, 
                                              surv = df_valid_ipacw_uncensored$surv, 
                                              t = t_fup,
                                              nk = 7,
                                              weights_in = df_valid_ipacw_uncensored$comb_weight)
  
  ### Calculate calibration plot using IPCW approach using the smoothed Horvitz-Thompson approach
  calib_smoothed_ht <- est_calib_smoothed_ht(data = df_valid_ipacw, 
                                             surv = df_valid_ipacw$surv, 
                                             t = t_fup,
                                             nk = 7,
                                             weights_in = df_valid_ipacw$comb_weight)
  
  
  ### Calculate calibration plot using weighted pseudo-value approach
  calib_pv_ipcw <- est_calib_pv_ipcw(data = df_valid_ipacw, 
                                     surv = df_valid_ipacw$surv, 
                                     t = t_fup,
                                     nk = 7,
                                     weights_in = df_valid_ipacw$comb_weight)
  
  calib_ipcw_uncensored[["ICI"]]
  calib_smoothed_ht[["ICI"]]
  calib_pv_ipcw[["ICI"]]
  calib_ipcw_uncensored[["E50"]]
  calib_smoothed_ht[["E50"]]
  calib_pv_ipcw[["E50"]]
  calib_ipcw_uncensored[["E90"]]
  calib_smoothed_ht[["E90"]]
  calib_pv_ipcw[["E90"]]
  calib_ipcw_uncensored[["plot"]]
  calib_smoothed_ht[["plot"]]
  calib_pv_ipcw[["plot"]]
  
  ### Write a function to save these plots
  save_calibration_plot <- function(calib_obj,
                                    filename,
                                    width = 7,
                                    height = 7,
                                    res = 300) {
    
    ragg::agg_png(
      filename = filename,
      width = width,
      height = height,
      units = "in",
      res = res
    )
    
    print(calib_obj$plot)
    
    dev.off()
  }
  
  ### Save the plots
  save_calibration_plot(calib_ipcw_uncensored, paste0("figures/calibration_moderate_kvg_ipcw_uncesnored", gender_in, "_trunc", truncate, ".png"))
  save_calibration_plot(calib_smoothed_ht, paste0("figures/calibration_moderate_kvg_ht", gender_in, "_trunc", truncate, ".png"))
  save_calibration_plot(calib_pv_ipcw, paste0("figures/calibration_moderate_kvg_pv", gender_in, "_trunc", truncate, ".png"))
  
  ###
  ### Function to create a summary table of calibration metrics
  ###
  create_calibration_summary <- function(...,
                                         model_names = NULL,
                                         digits = 3) {
    
    ### Store all calibration objects supplied via ...
    calib_list <- list(...)
    
    ### If model names are not supplied,
    ### create default names
    if (is.null(model_names)) {
      model_names <- paste0(
        "Model_",
        seq_along(calib_list)
      )
    }
    
    ### Check the number of names matches the number of calibration objects
    if (length(model_names) != length(calib_list)) {
      stop(
        "Length of model_names must equal the number of calibration objects."
      )
    }
    
    ### Extract ICI, E50, and E90 from each calibration object and combine into a single data frame
    summary_df <- do.call(
      rbind,
      lapply(
        seq_along(calib_list),
        function(i) {
          
          ### Create one row per calibration object
          data.frame(
            Method = model_names[i],
            
            ### Round metrics to requested number of digits
            ICI = round(
              calib_list[[i]]$ICI,
              digits = digits
            ),
            
            E50 = round(
              calib_list[[i]]$E50,
              digits = digits
            ),
            
            E90 = round(
              calib_list[[i]]$E90,
              digits = digits
            )
          )
        }
      )
    )
    
    ### Remove row names created during rbind
    rownames(summary_df) <- NULL
    
    ### Return summary table
    return(summary_df)
  }
  
  
  ###
  ### Example usage
  ###
  
  ###
  ### Create summary table
  ###
  results_table <- create_calibration_summary(
    calib_ipcw_uncensored,
    calib_smoothed_ht,
    calib_pv_ipcw,
    model_names = c(
      "BLR",
      "Smoothed HT",
      "Pseudo-value"
    ),
    digits = 3
  )
  
  ###
  ### View results
  ###
  print(results_table)
  saveRDS(results_table, paste0("data/calibration_moderate_kvg_metric_table", gender_in, "_trunc", truncate, ".rds"))
  write.csv(results_table, paste0("data/calibration_moderate_kvg_metric_table", gender_in, "_trunc", truncate, ".csv"))
  
}

calibrate_truncate(gender_in = 1, truncate = 10)
calibrate_truncate(gender_in = 2, truncate = 10)
calibrate_truncate(gender_in = 1, truncate = 100)
calibrate_truncate(gender_in = 2, truncate = 100)
print(paste("FINISHED", Sys.time()))
