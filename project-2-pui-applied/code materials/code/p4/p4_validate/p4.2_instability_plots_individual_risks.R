###
### Read in bootstrapped models and create instability plot
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
library(ggplot2)
library(survival)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- as.numeric(180)

### Write function to create plot
create_instability_plot <- function(gender, age_knots = FALSE, caltime = FALSE, n_pat = 3000, n_model = 500){
  
  ### Define gender_char
  gender_char <- c("male", "female")[gender]
  
  ### Extract model fits
  if (age_knots == 3){
    
    ### Read in models, offsets and bhaz's
    fit_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_stability_cox_", gender, "_samp", x, ".rds", sep = ""))})
    bhaz_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_stability_cox_bhaz_", gender, "_samp", x, ".rds", sep = ""))})
    offset_means_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_stability_offset_means_", gender, "_samp", x, ".rds", sep = ""))})
    
    ### Read in prototype3
    fit_prototype3 <- readRDS(paste("data/p4/prototype3_cox_", gender, "_imp", 1, ".rds", sep = ""))
    bhaz_prototype3 <- readRDS(paste("data/p4/prototype3_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    
  } else if (age_knots == 4) {
    
    if (caltime == TRUE){
      ### Read in models, offsets and bhaz's
      fit_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_4knots_caltime_stability_cox_", gender, "_samp", x, ".rds", sep = ""))})
      bhaz_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_4knots_caltime_stability_cox_bhaz_", gender, "_samp", x, ".rds", sep = ""))})
      offset_means_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_4knots_caltime_stability_offset_means_", gender, "_samp", x, ".rds", sep = ""))})
      
      ### Read in prototype3
      fit_prototype3 <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz_prototype3 <- readRDS(paste("data/p4/prototype3_4knots_caltime_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    }
    
    if (caltime == FALSE){
      ### Read in models, offsets and bhaz's
      fit_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_4knots_stability_cox_", gender, "_samp", x, ".rds", sep = ""))})
      bhaz_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_4knots_stability_cox_bhaz_", gender, "_samp", x, ".rds", sep = ""))})
      offset_means_list <- lapply(1:n_model, function(x){readRDS(paste("data/p4/prototype3_4knots_stability_offset_means_", gender, "_samp", x, ".rds", sep = ""))})
      
      ### Read in prototype3
      fit_prototype3 <- readRDS(paste("data/p4/prototype3_4knots_cox_", gender, "_imp", 1, ".rds", sep = ""))
      bhaz_prototype3 <- readRDS(paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", 1, ".rds", sep = ""))
    }

  }
  
  ### Read in validation dataset
  data_valid <- readRDS(paste("data/p4/dfs_valid_", gender_char, ".rds", sep = ""))[[1]]
  
  ### Pick n_pat patients at random
  set.seed(101)
  data_valid_subset <- data_valid[sample(1:nrow(data_valid), n_pat, replace = FALSE), ]
  
  ### Add values for offset variables/treatment strategy
  data_valid_subset <- 
    dplyr::mutate(data_valid_subset,
                  ### Create offset vars
                  offset_statins_timevar_lnHR = 0,
                  offset_ah_timevar_lnHR = 0,
                  offset_smoking_timevar_dummy1_lnHR = 0,
                  offset_smoking_timevar_dummy2_lnHR = 0)
  

  ### Estimate risks according to the non-bootstrapped model of prototype3
  data_valid_subset$pred_10y <- 1 - est_surv_offset_prototype3(newdata = data_valid_subset, 
                                                               fit = fit_prototype3, 
                                                               bhaz = bhaz_prototype3, 
                                                               time = 10*365.25)
  
  ### Write function to estimate risks
  estimate_risks <- function(t_eval){
    
    pred_list <- lapply(1:n_model, function(x){
      
      ### Get the lp
      lp <- predict(fit_list[[x]], newdata = data_valid_subset, reference = "sample") - sum(offset_means_list[[x]])
      
      ### Get the survival probability
      surv <- as.numeric(exp(-exp(lp)*bhaz_list[[x]]$hazard[max(which(bhaz_list[[x]]$time <= t_eval))]))
      pred <- 1-surv
      
      return(pred)
    })
    
    pred_list <- do.call("cbind", pred_list)
    colnames(pred_list) <- paste("pred", 1:n_model, sep = "")
    
    return(pred_list) 
    
  }
  
  ### Estiamte risks
  pred_10y <- estimate_risks(10*365.25)
  # pred_5y <- estimate_risks(5*365.25)
  
  ### Add to data frame
  data_valid_subset <- cbind(data_valid_subset, pred_10y)
  
  ### Create data frame from the predicted probabilities
  data_valid_subset_long <- data_valid_subset |>
    dplyr::select(patid, pred_10y, paste("pred",1:n_model[-22], sep = "")) |>
    tidyr::pivot_longer(cols = paste("pred",1:n_model[-22], sep = ""), names_to = "sample", values_to = "pred")
  
  ### Create ggplot
  plot_out <- ggplot(data = data_valid_subset_long) +
    geom_point(aes(x = pred_10y, y = pred), alpha = 0.1, size = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ylab("Risks from bootstrapped models") + xlab("Predicted risk")
  
  ### Save plot
  ggsave(paste("figures/p4/prototype3_stability_", gender, "nk", age_knots, "_caltime", as.numeric(caltime), ".png", sep = ""), plot_out)

  return(data_valid_subset_long)
  
}

### Create plots
# Female
create_instability_plot(2, age_knots = 3, caltime = FALSE)
create_instability_plot(2, age_knots = 4, caltime = FALSE)
create_instability_plot(2, age_knots = 4, caltime = TRUE)
# Male
create_instability_plot(1, age_knots = 3, caltime = FALSE)
create_instability_plot(1, age_knots = 4, caltime = FALSE)
create_instability_plot(1, age_knots = 4, caltime = TRUE)

###
### The following code will produce a Figure intended for the manuscript, although this is no longer needed if putting
### only in the supplementary material.
### 
# print("female model 4knots")
# female_plotdata <- create_instability_plot(2, knots4 = TRUE)
# print("male model  4knots")
# male_plotdata <- create_instability_plot(1, knots4 = TRUE)
# 
# ### Create combined plot
# female_plot <- ggplot(data = female_plotdata) +
#   geom_point(aes(x = pred_10y, y = pred), alpha = 0.01, size = 0.5) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   ylab("Risks from bootstrapped models") + xlab("Predicted risk") +
#   ggplot2::ggtitle("Female model")
# 
# male_plot <- ggplot(data = male_plotdata) +
#   geom_point(aes(x = pred_10y, y = pred), alpha = 0.01, size = 0.5) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   ylab("Risks from bootstrapped models") + xlab("Predicted risk") +
#   ggplot2::ggtitle("Male model")
# 
# ### Combine
# combined_plot <- ggpubr::ggarrange(female_plot, male_plot, ncol = 2, nrow = 1)
# ggplot2::ggsave(paste("figures/p4/prototype3_stability.png"), plot = combined_plot, width = 12, height = 6, dpi = 600)
