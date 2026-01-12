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

### Read in df
df <- lung

### Define t.eval
t.eval <- 200

### Add an id variable, that will be utilised when assessing calibration with pseudo-values
df$id <- 1:nrow(df)

### Change status value to be 0/1, not 1/2
df$status <- df$status - 1

### Lets create a variable for censored or not censored 
### (if an individual doesn't have an event, they are censored at end of follow up)
### This will be used in the IPCW method
df$cens_time <- df$time
df$cens_indicator <- 1 - df$status

### Get complete cases with respect to variables we will model
df <- df[complete.cases(df[,c("age", "sex", "ph.ecog", "ph.karno")]),]

### Not really an external validation, as using the same dataset for development and validation
### however, splitting the dataset in half leads to erroneous results, as development and validation
### samples both too small

# dat.devel <- df[1:113, ]
# dat.valid <- df[114:226, ]
dat.devel <- df
dat.valid <- df

### Fit model
fit1 <- coxph(Surv(time, status) ~ age + sex + ph.ecog + ph.karno, data = dat.devel)
bhaz1 <- basehaz(fit1, centered = TRUE)

###################################
### Estimate calibration curves ###
###################################

## Proportional hazards approach
calib.ph <- est_calib_ph(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 3)
## Pseudo-value approach
calib.pv <- est_calib_pv(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 3, split.n.groups = 3)
## Pseudo-value approach with ipcw (binder)
calib.pv.ipcw <-  est_calib_pv_ipcw(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 3, 
                                               ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ age + sex + ph.ecog + ph.karno"))
## IPCW approach
calib.ipcw <- est_calib_ipcw(data = dat.valid, fit  = fit1, bhaz = bhaz1, t = t.eval, nk = 3,
                             ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ age + sex + ph.ecog + ph.karno"))


#####################################
### Estimate confidence intervals ###
#####################################

### Get range of predicted risks and vector of values to plot over
p1p99 <- quantile(calib.ph[["plot"]]$data$pred, p = c(0.01, 0.99))
pred.plot.range <- seq(p1p99[1], p1p99[2], length.out = 1000)

### Calculate confidence intervals
calib.ph.ci <- est_calib_CI_boot_wrapper(
  calib.method = "ph", 
  data = dat.valid, fit = fit1, bhaz = bhaz1, t = t.eval,  nk = 3, 
  pred.plot.range = pred.plot.range, 
  CI = 95, 
  R.boot = 1000)

calib.pv.ci <- est_calib_CI_boot_wrapper(
  calib.method = "pv", 
  data = dat.valid, fit = fit1, bhaz = bhaz1, t = t.eval, nk = 3, 
  pv.split.n.groups = 3,
  pred.plot.range = pred.plot.range, 
  CI = 95, 
  R.boot = 1000)

calib.pv.ipcw.ci <- est_calib_CI_boot_wrapper(
  calib.method = "pv_ipcw", 
  data = dat.valid, fit = fit1, bhaz = bhaz1, t = t.eval, nk = 3, 
  ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ age + sex + ph.ecog + ph.karno"),
  pred.plot.range = pred.plot.range, 
  CI = 95, 
  R.boot = 1000)

calib.ipcw.ci <- est_calib_CI_boot_wrapper(
  calib.method = "ipcw", 
  data = dat.valid, fit = fit1, bhaz = bhaz1, t = t.eval, nk = 3, 
  ipcw.formula = as.formula("Surv(cens_time, cens_indicator) ~ age + sex + ph.ecog + ph.karno"),
  pred.plot.range = pred.plot.range, 
  CI = 95, 
  R.boot = 1000)


#######################################
### Combine data and create Figures ###
#######################################

### Add CI's
calib.ph <- add_CI_ce(calib.ph, calib.ph.ci)
calib.pv <- add_CI_ce(calib.pv, calib.pv.ci)
calib.pv.ipcw <- add_CI_ce(calib.pv.ipcw, calib.pv.ipcw.ci)
calib.ipcw <- add_CI_ce(calib.ipcw, calib.ipcw.ci)

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
  margdata <- subset(plot$data, ltynum == "Calibration")
  
  ### Add scatter and marginal
  plot <- plot +
    ggplot2::geom_point(ggplot2::aes(x = pred, y = pred.obs), 
                        data = margdata,
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
png("figures/applied_example_external_validation.png", res = 600, width = 7, height = 8, unit = "in")
gridExtra::grid.arrange(calib.ph[["plot"]], calib.pv[["plot"]], calib.pv.ipcw[["plot"]], calib.ipcw[["plot"]],
                        layout_matrix = base::matrix(seq_len(4),nrow = 2, ncol = 2, byrow = TRUE))
dev.off()
