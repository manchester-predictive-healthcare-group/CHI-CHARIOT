###
### Produce plots for the non-linear effects of each variable interacted with age
###

### For now, just going to do this for the first imputed model

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

library(survival)

#####################
### Preliminaries ###
#####################

###
### Write function to extract development dataset
###
extract_devel <- function(gender, imp.num){

  ### Read in the first imputed dataset
  imp.list <- readRDS(paste("data/p4/dfs_devel_", c("male", "female")[gender], ".rds", sep = ""))
  
  ### Define imputed dataset
  df <- imp.list[[imp.num]]
  
  return(df)
}

###
### Extract development datasets for the first model
df.male <- extract_devel(gender = 1, imp.num = 1)
df.female <- extract_devel(gender = 2, imp.num = 1)

### Want to generate the hazard ratio in an easy way, reference to the average. 
### i.e. want hazard ratio to be 1, apart from effect of variable of interest

### The process will be different for categorical and continuous covariates
### Start with binary covariates

### Create a fake data.frame with which to generate predictions
### We assign continuous variables (non-age) to be the mean, as this is what we will evaluate the effect of age at
create_dummy <- function(df){
  dummydf <- df[1, ]
  dummydf$bmi <- mean(df$bmi)
  dummydf$sbp <- mean(df$sbp)
  dummydf$nonhdl <- mean(df$nonhdl)
  dummydf$hypertension <- "Absent"
  dummydf$ra <- "Absent"
  dummydf$af <- "Absent"
  dummydf$ckd <- "Absent"
  dummydf$smi <- "Absent"
  dummydf$fhcvd <- "Absent"
  dummydf$migraine <- "Absent"
  dummydf$sle <- "Absent"
  dummydf$impotence <- "Absent"
  dummydf$cortico <- "Absent"
  dummydf$antipsy <- "Absent"
  dummydf$copd <- "Absent"
  dummydf$int_dis <- "Absent"
  dummydf$downs <- "Absent"
  dummydf$oral_cancer <- "Absent"
  dummydf$brain_cancer <- "Absent"
  dummydf$lung_cancer <- "Absent"
  dummydf$blood_cancer <- "Absent"
  dummydf$pre_eclampsia <- "Absent"
  dummydf$postnatal_depression <- "Absent"
  dummydf$smoking[1] <- "Non-smoker"
  dummydf$ethnicity[1] <- "white"
  dummydf$diabetes[1] <- "Absent"
  dummydf$offset_statins_timevar_lnHR <- 0
  dummydf$offset_ah_timevar_lnHR <- 0
  dummydf$offset_smoking_timevar_dummy1_lnHR <- 0
  dummydf$offset_smoking_timevar_dummy2_lnHR <- 0
  return(dummydf)
}

###
### Create dummy datasets
dummydf.female <- create_dummy(df.female)
dummydf.male <- create_dummy(df.male)
str(dummydf.female)
###
### Get fits
fit.female <- readRDS(paste("data/p4/prototype3_cox_", 2, "_imp1.rds", sep = ""))
fit.male <- readRDS(paste("data/p4/prototype3_cox_", 1, "_imp1.rds", sep = ""))

#################
### Functions ###
#################

###
### Function to create plot for hazard ratio as age varies, binary variables
### This will do it for just one gender (i.e. model fit)
###
create_plot_binary_1model <- function(var, fit, dummydf, savename = "_"){

  ###
  ### Function to get hazard ratio for a given binary variable at a specific age
  get_HR_binary <- function(age, var, dummydf){
    
    ### Create dummy datasets for with/without variable
    dummydf_without <- dummydf
    dummydf_with <- dummydf
    
    ### Assign variable of interest
    dummydf_with[,var] <- "Present"
    dummydf_without[,var] <- "Absent"
    
    ### Assign age
    dummydf_with$age <- age
    dummydf_without$age <- age
    
    ### Get hazard with/without
    with <- exp(predict(fit, newdata = dummydf_with, type = "lp"))
    without <- exp(predict(fit, newdata = dummydf_without, type = "lp"))
    
    ### Return
    return(with/without)
    
  }
  
  ### Get hazard ratio across ages 18 - 100
  hr <- lapply(18:100, get_HR_binary, var = var, dummydf = dummydf)
  
  ### Plot
  png(paste("figures/p4/prototype3_assess_nonlinear", savename, var, ".png", sep = ""), width = 6, height = 6, unit = "in", res = 300)
  plot(18:100, hr, type = "l", xlab = "Age", ylab = "HR", main = var)
  dev.off()
  
}


###
### Function to create plot for hazard ratio as age varies, binary variables
### This will do it for two models (so can plot male and female on same plot)
###
create_plot_binary <- function(var, 
                               fit1 = fit.male, fit2 = fit.female, 
                               dummydf1 = dummydf.male, dummydf2 = dummydf.female, 
                               level1 = "male", level2 = "female", level.title = "gender",
                               savename = "_"){

  ###
  ### Function to get hazard ratio for a given binary variable at a specific age
  get_HR_binary <- function(age, var, fit, dummydf){
    
    ### Create dummy datasets for with/without variable
    dummydf_without <- dummydf
    dummydf_with <- dummydf
    
    ### Assign variable of interest
    dummydf_with[,var] <- "Present"
    dummydf_without[,var] <- "Absent"
    
    ### Assign age
    dummydf_with$age <- age
    dummydf_without$age <- age
    
    ### Get hazard with/without
    with <- exp(predict(fit, newdata = dummydf_with, type = "lp"))
    without <- exp(predict(fit, newdata = dummydf_without, type = "lp"))
    
    ### Return
    return(with/without)
    
  }
  
  ### Get hazard ratio across ages 18 - 100
  hr1 <- unlist(lapply(18:100, get_HR_binary, var = var, fit = fit1, dummydf = dummydf1))
  hr2 <- unlist(lapply(18:100, get_HR_binary, var = var, fit = fit2, dummydf = dummydf2))
  
  ### Create data afor ggplot
  hr.gg.data <- data.frame("age" = rep(18:100, 2), 
                           "HR" = c(hr1, hr2), 
                           "var" = c(rep(level1, length(hr1)), rep(level2, length(hr2))))
  colnames(hr.gg.data)[3] <- level.title
  
  ### Create ggplot
  hr.gg <- ggplot2::ggplot(data = hr.gg.data) +
    ggplot2::geom_line(ggplot2::aes(x = age, y = HR, color = .data[[level.title]])) +
    ggplot2::ggtitle(var)
  
  ### Plot
  png(paste("figures/p4/prototype3_assess_nonlinear", savename, var, ".png", sep = ""), width = 6, height = 6, unit = "in", res = 300)
  print(plot(hr.gg))
  dev.off()
  
}


###
### Function to create plot for hazard ratio as age varies, categorical variables
###
create_plot_cat <- function(var, fit, dummydf, savename = "_"){
  
  ###
  ### Function to get hazard ratio for a given categorical variable at a specific age
  get_HR_cat <- function(age, var, dummydf){
    
    ### Create dummy datasets for with/without variable
    dummydf_without <- dummydf
    dummydf_with <- dummydf
    
    ### Assign age
    dummydf_with$age <- age
    dummydf_without$age <- age
    
    ### Get levels of variable of interest
    levels <- levels(dummydf[,var])
    
    ### Assign without to have the reference category
    dummydf_without[,var] <- levels[1]
    
    ### Create vector to store the HRs
    hr <- vector(length = length(levels) - 1)
    
    ### Cycle through and get hr
    for (i in 2:length(levels)){
      
      ### Assign variable of interest
      dummydf_with[,var] <- levels[i]
      
      ### Get hazard with/without
      with <- exp(predict(fit, newdata = dummydf_with, type = "lp"))
      without <- exp(predict(fit, newdata = dummydf_without, type = "lp"))
      
      ### Save hr
      hr[(i-1)] <- with/without
      
    }
    
    return(hr)
    
  }
  
  ### Get hazard ratio across ages 18 - 100
  hr <- lapply(18:100, get_HR_cat, var = var, dummydf = dummydf)
  
  ### Combine into dataframes
  hr <- data.frame("age" = 18:100, do.call("rbind", hr))
  colnames(hr) <- c("age", levels(dummydf[,var])[-1])
  
  ### Melt data for ggplot
  hr.gg.data <- tidyr::pivot_longer(hr, cols = colnames(hr)[-1], names_to = "Smoking status", values_to = "HR")
  
  ### Create ggplot
  colnames(hr.gg.data)
  hr.gg <- ggplot2::ggplot(data = hr.gg.data) +
    ggplot2::geom_line(ggplot2::aes(x = age, y = HR, color = `Smoking status`)) +
    ggplot2::ggtitle(var)
  
  ### Plot
  png(paste("figures/p4/prototype3_assess_nonlinear", savename, var, ".png", sep = ""), width = 6, height = 6, unit = "in", res = 300)
  plot(hr.gg)
  dev.off()
  
}


###
### Function to create plot for hazard ratio as a continuous variable varies, for a variety of ages
###
create_plot_cont <- function(var, var.range, var.ref, ages, fit, dummydf, savename = "_"){
  
  ###
  ### Function to get hazard ratio for a given variable at a specific point, relative to a reference point
  get_HR_cont <- function(var.value, var, var.reference, dummydf){
    
    ### Create dummy datasets for with/without variable
    dummydf_without <- dummydf
    dummydf_with <- dummydf
    
    ### Assign var and var.reference
    dummydf_with[,var] <- var.value
    dummydf_without[,var] <- var.reference
    
    ### Create vector to store the HRs
    hr <- vector(length = length(ages))
    
    ### cycle through ages
    for (i in 1:length(ages)){
      
      ### Assign age
      dummydf_with$age <- ages[i]
      dummydf_without$age <- ages[i]
      
      ### Get hazard with/without
      with <- exp(predict(fit, newdata = dummydf_with, type = "lp"))
      without <- exp(predict(fit, newdata = dummydf_without, type = "lp"))
      
      ### Save hr
      hr[i] <- with/without
      
    }
    
    return(hr)
    
  }
  
  ### Get hazard ratio across ages 18 - 100
  hr <- lapply(var.range, get_HR_cont, var = var, var.reference = var.ref, dummydf = dummydf)
  
  ### Combine into dataframes
  hr <- data.frame(var = var.range, do.call("rbind", hr))
  colnames(hr) <- c(var, ages)
  
  ### Melt data for ggplot
  hr.gg.data <- tidyr::pivot_longer(hr, cols = colnames(hr)[-1], names_to = "age", values_to = "HR")
  
  ### Create ggplot
  hr.gg <- ggplot2::ggplot(data = hr.gg.data) +
    ggplot2::geom_line(ggplot2::aes(x = .data[[var]], y = HR, color = age)) +
    ggplot2::ggtitle(var)
  
  ### Plot
  png(paste("figures/p4/prototype3_assess_nonlinear", savename, var, ".png", sep = ""), width = 6, height = 6, unit = "in", res = 300)
  plot(hr.gg)
  dev.off()
  
}


###
### Function to create plot for hazard ratio of age, with other continuous variables set to mean
###
create_plot_age <- function(fit1 = fit.male, fit2 = fit.female, 
                            dummydf1 = dummydf.male, dummydf2 = dummydf.female, 
                            level1 = "male", level2 = "female", level.title = "gender", savename = "_"){
  
  ###
  ### Function to get hazard ratio for a given categorical variable at a specific age and values of the continuous variable
  get_HR_age <- function(age, age.ref, fit, dummydf){
    
    ### Create dummy datasets for with/without variable
    dummydf_without <- dummydf
    dummydf_with <- dummydf
    
    ### Assign age
    dummydf_with$age <- age
    dummydf_without$age <- age.ref
    
    ### Get hazard with/without
    with <- exp(predict(fit, newdata = dummydf_with, type = "lp"))
    without <- exp(predict(fit, newdata = dummydf_without, type = "lp"))
    
    ### Return
    return(with/without)
    
  }
  
  ### Get hazard ratio across ages 18 - 100
  hr1 <- unlist(lapply(18:100, get_HR_age, age.ref = 40, fit = fit1, dummydf = dummydf1))
  hr2 <- unlist(lapply(18:100, get_HR_age, age.ref = 40, fit = fit2, dummydf = dummydf2))
  
  ### Create data afor ggplot
  hr.gg.data <- data.frame("age" = rep(18:100, 2), 
                           "HR" = c(hr1, hr2), 
                           "var" = c(rep(level1, length(hr1)), rep(level2, length(hr2))))
  colnames(hr.gg.data)[3] <- level.title
  
  ### Create ggplot
  hr.gg <- ggplot2::ggplot(data = hr.gg.data) +
    ggplot2::geom_line(ggplot2::aes(x = age, y = HR, color = .data[[level.title]]))
  
  ### Plot
  png(paste("figures/p4/prototype3_assess_nonlinear", savename, "age.png", sep = ""), width = 6, height = 6, unit = "in", res = 300)
  print(plot(hr.gg))
  dev.off()
  
}

####################
### Create plots ###
####################

###
### Create plots for binary variables
###

### Get binary variables
var.binary <- readRDS("data/p4/var_binary.rds")

### Create plots
for (var in var.binary){

  print(paste(var, Sys.time()))
  ### Full interaction fits
  create_plot_binary(var)

}

###
### Create plots for polytomous variables
###

### Get polytomous variables
var.poly <- readRDS("data/p4/var_poly.rds")

### Create plots
for (var in var.poly){

  print(paste(var, Sys.time()))

  ### Full interaction fits
  create_plot_cat(var, fit = fit.female, dummydf = dummydf.female, savename = "_2_")
  create_plot_cat(var, fit = fit.male, dummydf = dummydf.male, savename = "_1_")

}


###
### Create plots for continuous
###

### Get continuous variables
var.cont <- readRDS("data/p4/var_cont.rds")

### Create reference value vector to refer to
refs <- c("bmi" = 25, "sbp" = 120, "nonhdl" = 3, "IMD" = 1)

### Full interaction fits

### Female
for (var in var.cont){
  
  print(paste(var, Sys.time()))
  
  ###
  ### Full interaction fits
  
  ### Female
  create_plot_cont(var = var, 
                   var.range = seq(quantile(df.female[,var], probs = c(0.025, 0.975))[1],
                                   quantile(df.female[,var], probs = c(0.025, 0.975))[2],
                                   length = 100), 
                   var.ref = refs[var],
                   ages = c(25, 37.5, 50, 62.5, 75),
                   fit = fit.female, dummydf = dummydf.female, savename = "_2_")
  
  ### Male
  create_plot_cont(var = var, 
                   var.range = seq(quantile(df.male[,var], probs = c(0.025, 0.975))[1],
                                   quantile(df.male[,var], probs = c(0.025, 0.975))[2],
                                   length = 100), 
                   var.ref = refs[var],
                   ages = c(25, 37.5, 50, 62.5, 75),
                   fit = fit.male, dummydf = dummydf.male, savename = "_1_")
  
}


###
### Create plot for age
###
create_plot_age(savename = "_")


