###
### Program to estimate calibration in 1 development and validation dataset (split sample)
###

### Calibration curve produced seperately for each development/validation dataset combination
### This is primarily because for each development dataset, we get a different set of counter-factual survival times in the validation dataset

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
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

### Create a combined plot for models 0 to 4
create_combined_calib_plot <- function(gender){
  
  ### Create a list to store the plots in
  output_plot_list <- vector("list", 5)
  
  ### For each model, create a plot and store it
  for (model in 0:4){
    
    ### Read in data
    df_calib_smooth <- readRDS(paste("data/calib_ph_df_smooth_", gender, "_model", model, ".rds", sep = ""))
    df_calib_grouped <- readRDS(paste("data/calib_ph_df_grouped_", gender, "_model", model, ".rds", sep = ""))
    
    ### Create a ggplot out of this
    ggplot_comb <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), 
                         data = df_calib_smooth) + 
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
      ggplot2::geom_point(data = df_calib_smooth,  
                          ggplot2::aes(x = pred, y = pred.obs), col = grDevices::rgb(0, 0, 0, alpha = 0)) +
      ggplot2::geom_point(data = df_calib_grouped, 
                          ggplot2::aes(x = pred, y = obs, col = grDevices::rgb(0, 1, 0, alpha = 1))) + 
      ggplot2::theme(legend.position = "none") + 
      ggplot2::xlim(c(0,0.6)) + ggplot2::ylim(c(0,0.6))
    
    ### Add title
    if (model == 0){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(
        paste("Model ", model, ": Non-causal model; t = 10 years; Estimand: Non-causal estimate", sep = ""), width = 46)) 
    } else if (model == 1){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(
        paste("Model ", model, ": Treatment offset model; t = 10 years; Estimand: Risk under no antihypertensives", sep = ""), width = 46))
    } else if (model == 2){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(
        paste("Model ", model, ": Unexposed mediator model; t = 10 years; Estimand: Risk under no antihypertensives", sep = ""), width = 46)) 
    } else if (model == 3){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(
        paste("Model ", model, ": Modifiable risk factor model; t = 10 years; Estimand: Risk under current intervention strategy", sep = ""), width = 46)) 
    } else if (model == 4){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(
        paste("Model ", model, ": Two-component model; t = 10 years; Estimand: Risk under current intervention strategy", sep = ""), width = 46)) 
    }
    
    ### Assign to output
    output_plot_list[[(model+1)]] <- ggplot_comb
  }
  
  ### Combine the plots into a single ggplot
  plot_out <- ggpubr::ggarrange(plotlist = output_plot_list, nrow = 3, ncol = 2, align = "hv")
  
  ### Save with generic name for writing to other R markdown files
  ragg::agg_png(paste("figures/calib_plot_grid", gender, ".png", sep = ""), 
                width = 10, height = 15, scaling = 1, unit = "in", res = 600)
  plot(plot_out)
  dev.off()
  
  ### Save with Figure number, which may change depending on peer review
  ragg::agg_png(paste("figures/Figure3_", gender, ".png", sep = ""), 
                width = 10, height = 15, scaling = 1, unit = "in", res = 600)
  plot(plot_out)
  dev.off()
  
}

### Create a combined plot for all models (including models 5 to 7, for supp material)
create_combined_calib_plot_all <- function(gender){
  
  ### Create a list to store the plots in
  output_plot_list <- vector("list", 8)
  
  ### For each model, create a plot and store it
  for (model in 0:7){
    
    ### Read in data
    df_calib_smooth <- readRDS(paste("data/calib_ph_df_smooth_", gender, "_model", model, ".rds", sep = ""))
    df_calib_grouped <- readRDS(paste("data/calib_ph_df_grouped_", gender, "_model", model, ".rds", sep = ""))
    
    ### Create a ggplot out of this
    ggplot_comb <- ggplot2::ggplot() +
      ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs), 
                         data = df_calib_smooth) + 
      ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
      ggplot2::xlab("Predicted risk") + ggplot2::ylab("Predicted-observed risk") +
      ggplot2::geom_point(data = df_calib_smooth,  
                          ggplot2::aes(x = pred, y = pred.obs), col = grDevices::rgb(0, 0, 0, alpha = 0)) +
      ggplot2::geom_point(data = df_calib_grouped, 
                          ggplot2::aes(x = pred, y = obs, col = grDevices::rgb(0, 1, 0, alpha = 1))) + 
      ggplot2::theme(legend.position = "none") + 
      ggplot2::xlim(c(0,0.6)) + ggplot2::ylim(c(0,0.6))
    
    ### Add title
    if (model == 0){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(paste("Model ", model, ": Non-causal model; t = 10 years; Estimand: Non-causal estimate", sep = ""), width = 46)) 
    } else if (model == 1){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(paste("Model ", model, ": Treatment offset model; t = 10 years; Estimand: Risk under no antihypertensives", sep = ""), width = 46)) 
    } else if (model == 2){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(paste("Model ", model, ": Unexposed mediator model; t = 10 years; Estimand: Risk under no antihypertensives", sep = ""), width = 46)) 
    } else if (model == 3){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(paste("Model ", model, ": Modifiable risk factor model; t = 10 years; Estimand: Risk under current intervention strategy", sep = ""), width = 46)) 
    } else if (model == 4){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(paste("Model ", model, ": Two-component model; Estimand: Risk under current intervention strategy", sep = ""), width = 46)) 
    } else if (model == 5){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(paste("Model ", model, ": Modifiable risk factor model (SBP only); t = 10 years; Estimand: Risk under current intervention strategy", sep = ""), width = 46)) 
    } else if (model == 6){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(paste("Model ", model, ": Modifiable risk factor model (BMI only); t = 10 years; Estimand: Risk under current intervention strategy", sep = ""), width = 46)) 
    } else if (model == 7){
      ggplot_comb <- ggplot_comb + ggplot2::ggtitle(stringr::str_wrap(paste("Model ", model, ": Modifiable risk factor model (non-HDL cholesterol only); t = 10 years; Estimand: Risk under current intervention strategy", sep = ""), width = 46))
    }
    
    ### Save plot
    ragg::agg_png(paste("figures/calib_ph_", gender, "_model", model, ".png", sep = ""), 
                  width = 1, height = 1, scaling = 1/5, unit = "in", res = 600)
    plot(ggplot_comb)
    dev.off()
    
    ### Assign to output
    output_plot_list[[(model+1)]] <- ggplot_comb
  }
  
  ### Combine the plots into a single ggplot
  plot_out <- ggpubr::ggarrange(plotlist = output_plot_list, nrow = 3, ncol = 3, align = "hv")
  
  ### Save with generic name for writing to other R markdown files
  ragg::agg_png(paste("figures/calib_plot_grid", gender, "_all.png", sep = ""), 
                width = 15, height = 17, scaling = 1, unit = "in", res = 600)
  plot(plot_out)
  dev.off()
  
}

### Run this function
for (gender_in in c(2,1)){
  
  print(paste("gender = ", gender_in))
  create_combined_calib_plot_all(gender = gender_in)
  create_combined_calib_plot(gender = gender_in)
  
}

print(paste("FINISHED", Sys.time()))

warnings()