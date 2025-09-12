###
### Fit a standard cox model, do a sample size calculations, assess non-linear curves and interactions
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

library(mice)
library(survival)
library(rms)

### Extract chain seed and gender from command line
# args <- commandArgs(trailingOnly = T)
# gender <- as.numeric(args[1])

###
### Create plots for imputed datasets
###
for (gender in c(2,1)){
  
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  
  ### Load imp.comb
  imp.comb <- readRDS(paste("data/p4/imp_comb_prototype3", gender, ".rds", sep = ""))
  
  ### Create convergence plots
  print("convergence plots")
  png(paste("figures/p4/convergence_plot_prototype3", gender, ".png", sep = ""), width = 12, height = 15, unit = "in", res = 600)
  print(plot(imp.comb, y = c("ethnicity", "smoking", "bmi", "sbp", "nonhdl", "IMD"), layout = c(2, 6)))
  dev.off()
  
  ### Create density plots
  print("density plots")
  png(paste("figures/p4/density_plot_prototype3", gender, ".png", sep = ""), width = 6, height = 6, unit = "in", res = 600)
  print(densityplot(imp.comb, ~ IMD + bmi + sbp + nonhdl, layout = c(2, 2)))
  dev.off()
  
}

print("FINISHED")

# ###
# ### Create plots for imputed datasets at 20 iterations
# ###
# print("20 iter")
# for (gender in c(2,1)){
#   
#   gender_char <- c("male", "female")[gender]
#   print(paste("gender = ", gender_char))
#   
#   ### Load imp.comb
#   imp.comb <- readRDS(paste("data/p4/imp_comb_prototype3", gender, "_20iter.rds", sep = ""))
#   
#   ### Create convergence plots
#   print("convergence plots")
#   png(paste("figures/p4/convergence_plot_prototype3", gender, "_20iter.png", sep = ""), width = 5, height = 15, unit = "in", res = 600)
#   print(plot(imp.comb, layout = c(2, 9)))
#   dev.off()
#   
#   ### Create density plots
#   print("density plots")
#   png(paste("figures/p4/density_plot_prototype3", gender, "_20iter.png", sep = ""), width = 6, height = 9, unit = "in", res = 600)
#   print(densityplot(imp.comb, ~ IMD + bmi + sbp + nonhdl + cholesterol + hdl + ldl, layout = c(2, 4)))
#   dev.off()
#   
# }
# 
# print("FINISHED")
