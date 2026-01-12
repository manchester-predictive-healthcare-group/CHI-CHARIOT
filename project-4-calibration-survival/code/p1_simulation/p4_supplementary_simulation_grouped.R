###
### Supplementary simulation for the grouped approach
###
rm(list=ls())
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project4")
options(scipen = 999)

### Load functions
source("R/functions.R")
library(survival)

### Set seed
set.seed(505)

### Set m
m <- 1000000

### Define variables which are always the same
b.x1 <- 0.3
## (x4 is used to induce unmeasured confounding)
b.x4 <- 0
b.cens.x4 <- 0

### Define list elements to cycle through
crossing <- tidyr::expand_grid("DGM.y" = c(1,2,3),
                               "DGM.cens" = c(1,2,3,4),
                               "model" = c(1,2,3),
                               "splitdat" = c(2),
                               "t.cens" = c(1,2),
                               "t.eval" = c(1,2))

### Extract relevant ones
rows <- c(75, 87, 123, 135, 139)
print(crossing[rows,])

### Define variables that change from scenario to scenario
### Write a function to run the supplementary simulation
run_supplementary_simulation_grouped <- function(row, n.groups){
  
  ### Print row
  print(paste("row = ", row, Sys.time()))
  
  ### Extract t cens (this is whether stopped cox is applied or not, censoring at time t.evak)
  t.cens <- c(FALSE, TRUE)[as.numeric(crossing[row, "t.cens"])]
  print(paste("t.cens = ", t.cens))
  
  ### Extract time point
  t.eval <- c(0.25, 0.5)[as.numeric(crossing[row, "t.eval"])]
  print(paste("t.eval = ", t.eval))
  
  ### Define DGM.y
  DGM.y <- as.numeric(crossing[row, "DGM.y"])
  
  ### Define DGM.cens
  DGM.cens <- as.numeric(crossing[row, "DGM.cens"])
  
  ### Scenarios for DGM for the outcome
  b.x2 <- c(0,0.3,0.3)[as.numeric(crossing[row, "DGM.y"])]
  b.x3 <- c(0,0.3,0.3)[as.numeric(crossing[row, "DGM.y"])]
  b.x2sq <- c(0,0,0.3)[as.numeric(crossing[row, "DGM.y"])]
  b.x3log <- c(0,0,0.12)[as.numeric(crossing[row, "DGM.y"])]
  
  ### Scenarios for DGM for the censoring mechanism
  b.cens.x1 <- c(0,0.3,0.3,0.3)[as.numeric(crossing[row, "DGM.cens"])]
  b.cens.x2 <- c(0,0,0.3,0.3)[as.numeric(crossing[row, "DGM.cens"])]
  b.cens.x3 <- c(0,0,0.3,0.3)[as.numeric(crossing[row, "DGM.cens"])]
  b.cens.x2sq <- c(0,0,0,0.3)[as.numeric(crossing[row, "DGM.cens"])]
  b.cens.x3log <- c(0,0,0,0.12)[as.numeric(crossing[row, "DGM.cens"])]
  
  ### Assign model
  model <- as.numeric(crossing[row, "model"])
  
  ### Assign split
  splitdat <- as.numeric(crossing[row, "splitdat"])
  
  ### Generate baseline data
  cov <- data.frame("x1" = rnorm(m, 0, 1),
                    "x2" = rnorm(m, 0, 1),
                    "x3" = rnorm(m, 0, 1),
                    "x4" = rnorm(m, 0, 1))
  
  ### Simulate data (with and without censoring)
  print("sim dat")
  dat <- simulate_DGM(m = m, cov = cov,
                      lambda.y = 1, gamma.y = 1, lambda.cens = 1, gamma.cens = 1,
                      b.x1 = b.x1, b.x2 = b.x2, b.x3 = b.x3, b.x2sq = b.x2sq, b.x3log = b.x3log, b.x4 = b.x4,
                      b.cens.x1 = b.cens.x1, b.cens.x2 = b.cens.x2, b.cens.x3 = b.cens.x3, b.cens.x2sq = b.cens.x2sq, b.cens.x3log = b.cens.x3log, b.cens.x4 = b.cens.x4, 
                      seed = 101)
  
  dat.nocens <- simulate_DGM(m = m, cov = cov,
                             lambda.y = 1, gamma.y = 1, lambda.cens = 1, gamma.cens = 1,
                             b.x1 = b.x1, b.x2 = b.x2, b.x3 = b.x3, b.x2sq = b.x2sq, b.x3log = b.x3log, b.x4 = b.x4,
                             b.cens.x1 = b.cens.x1, b.cens.x2 = b.cens.x2, b.cens.x3 = b.cens.x3, b.cens.x2sq = b.cens.x2sq, b.cens.x3log = b.cens.x3log, b.cens.x4 = b.cens.x4, 
                             seed = 101,
                             cens = FALSE)
  
  ### Assign development and validation datasets
  if (splitdat == 1){
    ### Development and validation both have censoring
    dat.devel <- dat
    dat.valid <- dat
  } else if (splitdat == 2){
    ### Development has no censoring
    dat.devel <- dat.nocens
    dat.valid <- dat
  } else if (splitdat == 3){
    ### Validation has no censoring
    dat.devel <- dat
    dat.valid <- dat.nocens
  }
  
  ### Set cumhaz.t
  cumhaz.t = t.eval
  
  ### Define number of knots
  nk.fit <- 5
  
  ### Assign the model formula
  if (model == 1){
    formula.rhs <- "x1"
  } else if (model == 2){
    formula.rhs <- "x1 + x2 + x3"
  } else if (model == 3){
    formula.rhs <- "rms::rcs(x1,nk.fit) + rms::rcs(x2,nk.fit) + rms::rcs(x3,nk.fit)"
  } else if (model == "supplementary"){
    formula.rhs <- "x1 + x2 + x3 + x4 + x5"
  }
  
  ### Censor validation data at time t if specified
  ### Some methods have had problems if not doing this
  if (t.cens == TRUE){
    dat.valid <- dplyr::mutate(dat.valid, 
                               status = dplyr::case_when(time > t.eval ~ 0,
                                                         TRUE ~ status),
                               time = dplyr::case_when(time > t.eval ~ t.eval,
                                                       TRUE ~ time))
  }
  
  ### Fit a model
  fit1 <- coxph(as.formula(paste("Surv(time, status) ~ ", formula.rhs, sep = "")), data = dat.devel)
  bhaz1 <- basehaz(fit1, centered = TRUE)
  
  ### Estimate true calibration curve
  calib.true <- est_calib_true(data = dat.valid, 
                               fit = fit1,
                               bhaz = bhaz1,
                               t = t.eval,
                               true.lp = dat.valid$true.lp, 
                               cumhaz.t = cumhaz.t)
  
  ###
  ### Assess calibration by estimating Kaplan-Meier (observed) and mean predicted risk (predicted) within subgroups, defined by predicted risk
  ###
  calib_plot_group <- est_calib_plot_group(data = dat.valid, fit = fit1, bhaz = bhaz1, t = t.eval, n.groups = n.groups, surv = NULL)
  
  ### Add true calibration curve
  plot.data <- calib_plot_group[["plot"]]$data
  
  ### Reduce and rename true calibration data frame
  calib.true <- dplyr::slice_head(calib.true, n = 5000) |>
    dplyr::rename(obs = pred.obs) |>
    as.data.frame()
  
  ### Create plot
  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), color = "red", data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = obs), color = "blue", data = calib.true) +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle(paste("C-DGM-", 
                           c("R", "SL", "ML", "QU")[DGM.cens],
                           ", O-DGM-", 
                           c("SL", "ML", "QU")[DGM.y],
                           ", ", 
                           c("Model-SL", "Model-ML", "Model-QU")[model], 
                           sep = "")) +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk (Kaplan-Meier)")
  
  ### Save object as rds and plot file
  saveRDS(calib.true, paste("data/supplementary_simulation_grouped_c", DGM.cens, "_y", DGM.y, "model", model, "_ngroup", n.groups, "_calibtrue.rds", sep = ""))
  saveRDS(plot.data, paste("data/supplementary_simulation_grouped_c", DGM.cens, "_y", DGM.y, "model", model, "_ngroup", n.groups, "_plotdata.rds", sep = ""))
  Cairo::CairoPNG(paste("figures/supplementary_simulation_grouped_c", DGM.cens, "_y", DGM.y, "model", model, "_ngroup", n.groups, ".png", sep = ""), 
                  dpi = 600, width = 5, height = 5, unit = "in")
  plot(plot)
  dev.off()
  
}

### Run sim
lapply(rows, run_supplementary_simulation_grouped, n.groups = 10)
# lapply(rows, run_supplementary_simulation_grouped, n.groups = 50)
# lapply(rows, run_supplementary_simulation_grouped, n.groups = 100)

###
### Create Figure 6
###

### Function to read in and create plot
create_plot <- function(DGM.cens, DGM.y, model, n.groups){
  
  calib.true <- readRDS(paste("data/supplementary_simulation_grouped_c", DGM.cens, "_y", DGM.y, "model", model, "_ngroup", n.groups, "_calibtrue.rds", sep = ""))
  plot.data <- readRDS(paste("data/supplementary_simulation_grouped_c", DGM.cens, "_y", DGM.y, "model", model, "_ngroup", n.groups, "_plotdata.rds", sep = ""))
  
  ### Create plot
  plot <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = pred, y = obs), color = "red", data = plot.data) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = obs), color = "blue", data = calib.true) +
    ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
    ggplot2::ggtitle(paste("C-DGM-", 
                           c("R", "SL", "ML", "QU")[DGM.cens],
                           ", O-DGM-", 
                           c("SL", "ML", "QU")[DGM.y],
                           ", ", 
                           c("Model-SL", "Model-ML", "Model-QU")[model], 
                           sep = "")) +
    ggplot2::xlab("Predicted risk") + ggplot2::ylab("Observed risk (Kaplan-Meier)")
  
  return(plot)
  
}

### Create plots
plot1 <- create_plot(3,2,1,10)
plot2 <- create_plot(3,3,1,10)
plot3 <- create_plot(4,2,1,10)
plot4 <- create_plot(4,3,1,10)

### Add to plot
plots.comb <- gridExtra::arrangeGrob(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)

### Save png
ragg::agg_png(paste("figures/Figure6.png", sep = ""),
              units = "in", 
              res = 1000,
              width = 7, height = 7, scaling = 0.95)
grid::grid.draw(plots.comb)
dev.off()

########################
### Mean calibration ###
########################

### Define some new row
### We want scenarios where the model is perfectly specified
### Will use model-flex scenarios,  and O-DGM-poly, with increasing censoring
rows <- c(107, 119, 131, 143)

### Define variables that change from scenario to scenario
### Write a function to run the supplementary simulation
run_supplementary_simulation_grouped_mean <- function(row){
  
  ### Print row
  print(paste("row = ", row, Sys.time()))
  
  ### Extract t cens (this is whether stopped cox is applied or not, censoring at time t.evak)
  t.cens <- c(FALSE, TRUE)[as.numeric(crossing[row, "t.cens"])]
  print(paste("t.cens = ", t.cens))
  
  ### Extract time point
  t.eval <- c(0.25, 0.5)[as.numeric(crossing[row, "t.eval"])]
  print(paste("t.eval = ", t.eval))
  
  ### Define DGM.y
  DGM.y <- as.numeric(crossing[row, "DGM.y"])
  
  ### Define DGM.cens
  DGM.cens <- as.numeric(crossing[row, "DGM.cens"])
  
  ### Scenarios for DGM for the outcome
  b.x2 <- c(0,0.3,0.3)[as.numeric(crossing[row, "DGM.y"])]
  b.x3 <- c(0,0.3,0.3)[as.numeric(crossing[row, "DGM.y"])]
  b.x2sq <- c(0,0,0.3)[as.numeric(crossing[row, "DGM.y"])]
  b.x3log <- c(0,0,0.12)[as.numeric(crossing[row, "DGM.y"])]
  
  ### Scenarios for DGM for the censoring mechanism
  b.cens.x1 <- c(0,0.3,0.3,0.3)[as.numeric(crossing[row, "DGM.cens"])]
  b.cens.x2 <- c(0,0,0.3,0.3)[as.numeric(crossing[row, "DGM.cens"])]
  b.cens.x3 <- c(0,0,0.3,0.3)[as.numeric(crossing[row, "DGM.cens"])]
  b.cens.x2sq <- c(0,0,0,0.3)[as.numeric(crossing[row, "DGM.cens"])]
  b.cens.x3log <- c(0,0,0,0.12)[as.numeric(crossing[row, "DGM.cens"])]
  
  ### Assign model
  model <- as.numeric(crossing[row, "model"])
  
  ### Assign split
  splitdat <- as.numeric(crossing[row, "splitdat"])
  
  ### Generate baseline data
  cov <- data.frame("x1" = rnorm(m, 0, 1),
                    "x2" = rnorm(m, 0, 1),
                    "x3" = rnorm(m, 0, 1),
                    "x4" = rnorm(m, 0, 1))
  
  ### Simulate data (with and without censoring)
  print("sim dat")
  dat <- simulate_DGM(m = m, cov = cov,
                      lambda.y = 1, gamma.y = 1, lambda.cens = 1, gamma.cens = 1,
                      b.x1 = b.x1, b.x2 = b.x2, b.x3 = b.x3, b.x2sq = b.x2sq, b.x3log = b.x3log, b.x4 = b.x4,
                      b.cens.x1 = b.cens.x1, b.cens.x2 = b.cens.x2, b.cens.x3 = b.cens.x3, b.cens.x2sq = b.cens.x2sq, b.cens.x3log = b.cens.x3log, b.cens.x4 = b.cens.x4, 
                      seed = 101)
  
  dat.nocens <- simulate_DGM(m = m, cov = cov,
                             lambda.y = 1, gamma.y = 1, lambda.cens = 1, gamma.cens = 1,
                             b.x1 = b.x1, b.x2 = b.x2, b.x3 = b.x3, b.x2sq = b.x2sq, b.x3log = b.x3log, b.x4 = b.x4,
                             b.cens.x1 = b.cens.x1, b.cens.x2 = b.cens.x2, b.cens.x3 = b.cens.x3, b.cens.x2sq = b.cens.x2sq, b.cens.x3log = b.cens.x3log, b.cens.x4 = b.cens.x4, 
                             seed = 101,
                             cens = FALSE)
  
  ### Assign development and validation datasets
  if (splitdat == 1){
    ### Development and validation both have censoring
    dat.devel <- dat
    dat.valid <- dat
  } else if (splitdat == 2){
    ### Development has no censoring
    dat.devel <- dat.nocens
    dat.valid <- dat
  } else if (splitdat == 3){
    ### Validation has no censoring
    dat.devel <- dat
    dat.valid <- dat.nocens
  }
  
  ### Set cumhaz.t
  cumhaz.t = t.eval
  
  ### Define number of knots
  nk.fit <- 5
  
  ### Assign the model formula
  if (model == 1){
    formula.rhs <- "x1"
  } else if (model == 2){
    formula.rhs <- "x1 + x2 + x3"
  } else if (model == 3){
    formula.rhs <- "rms::rcs(x1,nk.fit) + rms::rcs(x2,nk.fit) + rms::rcs(x3,nk.fit)"
  } else if (model == "supplementary"){
    formula.rhs <- "x1 + x2 + x3 + x4 + x5"
  }
  
  ### Censor validation data at time t if specified
  ### Some methods have had problems if not doing this
  if (t.cens == TRUE){
    dat.valid <- dplyr::mutate(dat.valid, 
                               status = dplyr::case_when(time > t.eval ~ 0,
                                                         TRUE ~ status),
                               time = dplyr::case_when(time > t.eval ~ t.eval,
                                                       TRUE ~ time))
  }
  
  ### Fit a model
  fit1 <- coxph(as.formula(paste("Surv(time, status) ~ ", formula.rhs, sep = "")), data = dat.devel)
  bhaz1 <- basehaz(fit1, centered = TRUE)
  
  ### Estimate true calibration curve
  calib.true.ratio <- est_calib_true_mean(data = dat.valid, 
                                          fit = fit1,
                                          bhaz = bhaz1,
                                          t = t.eval,
                                          true.lp = dat.valid$true.lp, 
                                          cumhaz.t = cumhaz.t)
  
  ###
  ### Assess calibration by estimating Kaplan-Meier (observed) and mean predicted risk (predicted) within subgroups, defined by predicted risk
  ###
  calib.mean.ratio <- est_calib_mean_km(data = dat.valid, fit = fit1, bhaz = bhaz1, t = t.eval)
  
  ### Create output object
  output.object <- c("estimated_mean_calib" = calib.mean.ratio, "true_mean_calib" = calib.true.ratio)
  print(output.object)
  
  saveRDS(output.object, paste("data/supplementary_simulation_mean_c", DGM.cens, "_y", DGM.y, "model", model, ".rds", sep = ""))
  
}

### Calculate mean calibration
lapply(rows, run_supplementary_simulation_grouped_mean)

### Get the mean calibration data
extract_mean_calibration <- function(DGM.cens, DGM.y, model){
  
  ### Extract estimated mean calibration and true calibration ratios
  ratios <- readRDS(paste("data/supplementary_simulation_mean_c", DGM.cens, "_y", DGM.y, "model", model, ".rds", sep = ""))
  
  ### Turn into data.frame
  ratios <- t(data.frame(ratios))
  
  ### Assign names for when we create table
  rownames(ratios) <- paste("C-DGM-", 
                      c("R", "SL", "ML", "QU")[DGM.cens],
                      ", O-DGM-", 
                      c("SL", "ML", "QU")[DGM.y],
                      ", ", 
                      c("Model-SL", "Model-ML", "Model-QU")[model], 
                      sep = "")
  return(ratios)
}

### Extract them
ratio1 <- extract_mean_calibration(1,3,3)
ratio2 <- extract_mean_calibration(2,3,3)
ratio3 <- extract_mean_calibration(3,3,3)
ratio4 <- extract_mean_calibration(4,3,3)

### Combine into table
mean_calibration_table <- rbind(ratio1, ratio2, ratio3, ratio4)

saveRDS(mean_calibration_table, "data/supplementary_simulation_mean_calibration.rds")
