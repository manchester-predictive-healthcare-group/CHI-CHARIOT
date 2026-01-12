###
### Assessing calibration of model developed on the lung data
###

### Prelims
rm(list=ls())
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project4")
set.seed(505)
library(survival)

### Load functions
source("R/functions.R")
source("R/functions_boot_optimism_adjusted.R")
source("R/functions_boot_CI.R")

### Read in lung data
dat.devel <- lung

### Define t.eval
t.eval <- 200

### Add an id variable, that will be utilised when assessing calibration with pseudo-values
dat.devel$id <- 1:nrow(dat.devel)

### Change status value to be 0/1, not 1/2
dat.devel$status <- dat.devel$status - 1

### Lets create a variable for censored or not censored 
### (if an individual doesn't have an event, they are censored at end of follow up)
### This will be used in the IPCW method
dat.devel$cens_time <- dat.devel$time
dat.devel$cens_indicator <- 1 - dat.devel$status

### Get complete cases with respect to variables we will model
dat.devel <- dat.devel[complete.cases(dat.devel[,c("age", "sex", "ph.ecog", "ph.karno")]),]

### Fit model
fit1 <- coxph(Surv(time, status) ~ age + sex + ph.ecog + ph.karno, data = dat.devel)
bhaz1 <- basehaz(fit1, centered = TRUE)

###################################
### Estimate calibration curves ###
###################################

## Proportional hazards approach
calib.ph <- est_calib_ph(data = dat.devel, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 3)
## Pseudo-value approach
calib.pv <- est_calib_pv(data = dat.devel, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 3, split.n.groups = 3)
## Pseudo-value approach with ipcw (binder)
calib.pv.ipcw <-  est_calib_pv_ipcw(data = dat.devel, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 3, 
                                               ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ age + sex + ph.ecog + ph.karno"))
## IPCW approach
calib.ipcw <- est_calib_ipcw(data = dat.devel, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 3,
                             ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ age + sex + ph.ecog + ph.karno"))


#####################################################
### Estimate optimism adjusted calibration curves ###
#####################################################

### Get range of predicted risks and vector of values to plot over
p1p99 <- quantile(calib.ph[["plot"]]$data$pred, p = c(0.01, 0.99))
pred.plot.range <- seq(p1p99[1], p1p99[2], length.out = 1000)

### Estimate optimism adjusted calibration curves
calib.ph.optimism.adjusted <- est_calib_optimism_adjusted_boot_wrapper(
  calib.method = "ph", 
  data = dat.devel, fit = fit1, t = t.eval,  nk = 3, 
  pred.plot.range = pred.plot.range, 
  CI = 95, 
  R.boot = 1000)

calib.pv.optimism.adjusted <- est_calib_optimism_adjusted_boot_wrapper(
  calib.method = "pv", 
  data = dat.devel, fit = fit1, t = t.eval, nk = 3, 
  pv.split.n.groups = 3,
  pred.plot.range = pred.plot.range, 
  CI = 95, 
  R.boot = 1000)

calib.pv.ipcw.optimism.adjusted <- est_calib_optimism_adjusted_boot_wrapper(
  calib.method = "pv_ipcw", 
  data = dat.devel, fit = fit1, t = t.eval, nk = 3, 
  ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ age + sex + ph.ecog + ph.karno"),
  pred.plot.range = pred.plot.range, 
  CI = 95, 
  R.boot = 1000)

calib.ipcw.optimism.adjusted <- est_calib_optimism_adjusted_boot_wrapper(
  calib.method = "ipcw", 
  data = dat.devel, fit = fit1, t = t.eval, nk = 3, 
  ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ age + sex + ph.ecog + ph.karno"),
  pred.plot.range = pred.plot.range, 
  CI = 95, 
  R.boot = 1000)

##################################################
### Add marginal histograms and create figures ###
##################################################

### Add titles
calib.ph[["plot"]] <- calib.ph[["plot"]] + ggplot2::ggtitle("PH")
calib.pv[["plot"]] <- calib.pv[["plot"]] + ggplot2::ggtitle("PV-group-risk")
calib.pv.ipcw[["plot"]] <- calib.pv.ipcw[["plot"]] + ggplot2::ggtitle("PV-IPCW")
calib.ipcw[["plot"]] <- calib.pv.ipcw[["plot"]] + ggplot2::ggtitle("IPCW")

###
### Add marginal density plots
###

### Write a function to do so
add_marginal <- function(object){
  
  ### Extract the plot
  plot <- object[["plot"]]
  
  ### Add scatter and marginal
  plot <- plot +
    ggplot2::geom_point(ggplot2::aes(x = pred, y = pred.obs), 
                        col = grDevices::rgb(0,0,0,alpha=0)) 
  plot <- ggExtra::ggMarginal(plot, type = "density", margins = "x", size = 6)
  
  ### Add to object
  object[["plot"]] <- plot
  
  return(object)
}


### Apply the function to get marginal densities
calib.ph <- add_marginal(calib.ph)
calib.pv <- add_marginal(calib.pv)
calib.pv.ipcw <- add_marginal(calib.pv.ipcw)
calib.ipcw <- add_marginal(calib.ipcw)

### Arrange into a grid
png("figures/applied_example_internal_validation.png", res = 600, width = 7, height = 8, unit = "in")
gridExtra::grid.arrange(calib.ph[["plot"]], calib.pv[["plot"]], calib.pv.ipcw[["plot"]], calib.ipcw[["plot"]],
                        layout_matrix = base::matrix(seq_len(4),nrow = 2, ncol = 2, byrow = TRUE))
dev.off()
