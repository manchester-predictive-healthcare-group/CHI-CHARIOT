###
### Assess imputations
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("")
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

for (gender in c(2,1)){
  
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  
  ### Load imp.comb
  imp.comb <- readRDS(paste("data/mice_mids_prototype3_", gender, "_", chain.seed, ".rds", sep = ""))
  
  ### Create convergence plots
  print("convergence plots")
  png(paste("figures/convergence_plot_exposure", gender, ".png", sep = ""), width = 6, height = 9, unit = "in", res = 600)
  print(plot(imp.comb, layout = c(4, 5)))
  dev.off()
  
  ### Create density plots
  print("density plots")
  png(paste("figures/density_plot_exposure", gender, ".png", sep = ""), width = 6, height = 9, unit = "in", res = 600)
  print(densityplot(imp.comb, ~ sbp + sbp_unexposed + bmi + nonhdl + IMD, layout = c(2, 4)))
  dev.off()
  
}
