###
### Read in bootstrapped models and create instability plot
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
library(ggplot2)
library(survival)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- as.numeric(180)

### Write function to create plot
create_instability_df_synthetic <- function(gender, n_pat = 300, n_model = 150, t_eval_y = 10){
  
  # gender <- 1 
  # n_pat = 1000
  # n_model = 50
  # t_eval_y = 10
  
  ### Define gender_char
  gender_char <- c("male", "female")[gender]
  
  ### Read in models, offsets and bhaz's
  fit_list <- lapply(1:n_model, function(x){readRDS(paste0("data/p5/optimism_boot_reduced_cox_", gender, "_iter", x, ".rds"))})
  bhaz_list <- lapply(1:n_model, function(x){readRDS(paste0("data/p5/optimism_boot_reduced_cox_bhaz_", gender, "_iter", x, ".rds"))})
  offset_means_list <- lapply(1:n_model, function(x){readRDS(paste0("data/p5/optimism_boot_reduced_cox_offset_means_", gender, "_iter", x, ".rds"))})
  
  ### Read in main model
  fit_list_main <- readRDS(paste0("data/p5/small_fits_model_", gender_char, ".rds"))
  bhaz_list_main <- readRDS(paste0("data/p5/small_bhaz_model_", gender_char, ".rds"))
  offset_means_list_main <- readRDS(paste0("data/p5/small_means_model_", gender_char, ".rds"))
    
  ### Read in one imputed dataset (not boostrapped) to get individuals to do the plot for
  df <- readRDS(paste("data/p5/dfs_all_", gender_char, 1, ".rds", sep = ""))
  
  ### Select rows to define each profile's age band and ethnicity
  set.seed(101)
  row_id <- sample.int(nrow(df), n_pat, replace = FALSE)
  df_subset <- df[row_id, ]
  
  ### List predictors to replace, excluding age and ethnicity
  predictor_vars <- c(
    "hypertension", "ra", "af", "ckd", "smi", "fhcvd",
    "migraine", "sle", "diabetes", "impotence", "cortico",
    "antipsy", "smoking", "statins", "antihypertensives",
    "copd", "int_dis", "downs", "oral_cancer", "brain_cancer",
    "lung_cancer", "blood_cancer", "pre_eclampsia",
    "postnatal_depression", "impotence", "bmi", "sbp", "nonhdl", "IMD",
    "cholesterol", "hdl", "ldl"
  )
  
  ### Replace each predictor using a different matched donor
  for (i in seq_len(n_pat)) {
    
    ### Find donors in the same five-year age band and ethnicity
    donor_pool <- which(
      floor(df$age / 5) == floor(df_subset$age[i] / 5) &
        df$ethnicity == df_subset$ethnicity[i]
    )
    donor_pool <- setdiff(donor_pool, row_id[i])
    
    ### Require a separate donor for every predictor
    if (length(donor_pool) < length(predictor_vars)) {
      stop("Not enough matched donors for profile ", i)
    }
    
    ### Assign a different donor to each predictor
    donor_ids <- sample(donor_pool, length(predictor_vars), replace = FALSE)
    
    for (j in seq_along(predictor_vars)) {
      var <- predictor_vars[j]
      df_subset[[var]][i] <- df[[var]][donor_ids[j]]
    }
  }
  
  ### Replace real patient identifiers
  df_subset$patid <- paste0("synthetic_", seq_len(n_pat))
  
  ### Add values for offset variables/treatment strategy
  df_subset <- 
    dplyr::mutate(df_subset,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  
  ### Estimate risks using main models
  pred_list_main <- lapply(1:10, function(x){
    
    ### Get the four offsets
    offset_means_x <- offset_means_list_main[[x]]
    
    ### Get the lp
    lp <- predict(fit_list_main[[x]], newdata = df_subset, reference = "sample") - sum(offset_means_x)
    
    ### Get the survival probability
    surv <- as.numeric(exp(-exp(lp)*as.numeric(bhaz_list_main[[x]][[t_eval_y]])))
    pred <- 1-surv
    
    return(pred)
    
  })
  
  ### Convert onto cloglog scale
  est_risks_main_cloglog <- lapply(pred_list_main, function(x){log(-log(1-x))})
  
  ### Take average (on cloglog scale)
  est_risks_main_av <- do.call("cbind", est_risks_main_cloglog) |>
    apply(MARGIN = 1, FUN = mean)
  
  ### Convert back to risk score
  est_risks_main_av <- 1 - exp(-exp(est_risks_main_av))
  
  ### Estimate risks according to the main model
  df_subset$pred_main <- est_risks_main_av
  
  ### Estimate risks using bootstrapped models
  pred_list <- lapply(1:n_model, function(x){
    
    ### Get the lp
    lp <- predict(fit_list[[x]], newdata = df_subset, reference = "sample") - sum(offset_means_list[[x]])
    
    ### Get the survival probability
    surv <- as.numeric(exp(-exp(lp)*as.numeric(bhaz_list[[x]])))
    pred <- 1-surv
    
    return(pred)
    
  })
  
  ### Combine
  pred_list <- do.call("cbind", pred_list)
  colnames(pred_list) <- paste("pred", 1:n_model, sep = "")
  
  ### Add to data frame
  df_subset <- cbind(df_subset[c("patid", "pred_main")], pred_list)
  
  ### Create data frame from the predicted probabilities
  df_subset_long <- df_subset |>
    dplyr::select(patid, pred_main, paste("pred",1:n_model, sep = "")) |>
    tidyr::pivot_longer(cols = paste("pred",1:n_model, sep = ""), names_to = "sample", values_to = "pred")
  
  ### Save this
  saveRDS(df_subset_long, paste0("data/p5/instability_plot_df_synthetic_", gender, "_", t_eval_y, "_m", n_model, ".rds"))
  
}

### Create data for instability plots
create_instability_df_synthetic(2, t_eval_y = 10, n_model = 150)
create_instability_df_synthetic(1, t_eval_y = 10, n_model = 150)

create_instability_df_synthetic(2, t_eval_y = 10, n_model = 500)
create_instability_df_synthetic(1, t_eval_y = 10, n_model = 500)

### Create ggplots
create_instability_plot_synthetic <- function(gender, t_eval_y, n_plot, n_model = 150){
  
  ### Read in data
  df_instab <- readRDS(paste0("data/p5/instability_plot_df_synthetic_", gender, "_", t_eval_y, "_m", n_model, ".rds"))
  
  ### Reduce to a random subset
  set.seed(101)
  patids_random <- sample(unique(df_instab$patid),
                          n_plot,
                          replace = FALSE)
    
  df_instab <- df_instab[df_instab$patid %in% patids_random, ]
  
  ### Create plot
  plot_out <- ggplot2::ggplot(data = df_instab) +
    ggplot2::geom_point(ggplot2::aes(x = pred_main, y = pred), alpha = 0.025, size = 0.2) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::ylab("Risks from bootstrapped models") + ggplot2::xlab("Predicted risk")
  
  ### Save plot
  ggplot2::ggsave(paste0("figures/p5/stability_plot_synthetic_", gender, "_", t_eval_y, "_m", n_model, ".png"), plot_out)
  
}

create_instability_plot_synthetic(1, 10, n_plot = 300, n_model = 150)
create_instability_plot_synthetic(2, 10, n_plot = 300, n_model = 150)

create_instability_plot_synthetic(1, 10, n_plot = 300, n_model = 500)
create_instability_plot_synthetic(2, 10, n_plot = 300, n_model = 500)

