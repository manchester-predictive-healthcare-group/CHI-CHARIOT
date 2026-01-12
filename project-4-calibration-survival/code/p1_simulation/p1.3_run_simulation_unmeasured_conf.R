###
### Simulated data example
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
b.x4 <- 1
b.cens.x4 <- 1

### Define list elements to cycle through
crossing <- tidyr::expand_grid("DGM.y" = c(1,2,3),
                               "DGM.cens" = c(1,2,3,4),
                               "model" = c(1,2,3),
                               "splitdat" = c(2,3),
                               "t.cens" = c(2),
                               "t.eval" = c(1))

### Extract relevant ones
row <- as.numeric(commandArgs(trailingOnly = T)[1])
print(crossing[row,])

### Define variables that change from scenario to scenario

### Extract t cens
t.cens <- c(FALSE, TRUE)[as.numeric(crossing[row, "t.cens"])]
print(paste("t.cens = ", t.cens))

### Extract time point
t.eval <- c(0.25, 0.5)[as.numeric(crossing[row, "t.eval"])]
print(paste("t.eval = ", t.eval))

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

### Run the simulation
sim.out <- run_sim_large_sample(dat.devel, dat.valid, t.eval = t.eval, cumhaz.t = t.eval, t.cens = t.cens, model = model)

### Save object as rds file
saveRDS(sim.out, paste("data/simulation_m", nrow(dat.valid), 
                       "_b03",
                       "_y", as.numeric(crossing[row, "DGM.y"]),
                       "_c", as.numeric(crossing[row, "DGM.cens"]),
                       "_model", as.numeric(crossing[row, "model"]),
                       "_split", as.numeric(crossing[row, "splitdat"]),  
                       "_tcens", as.numeric(crossing[row, "t.cens"]),
                       "_teval", as.numeric(crossing[row, "t.eval"]), "_un.rds", sep = ""))

Cairo::CairoPNG(paste("figures/simulation_m", nrow(dat.valid),
                      "_b03", 
                      "_y", as.numeric(crossing[row, "DGM.y"]),
                      "_c", as.numeric(crossing[row, "DGM.cens"]),
                      "_model", as.numeric(crossing[row, "model"]),
                      "_split", as.numeric(crossing[row, "splitdat"]), 
                      "_tcens", as.numeric(crossing[row, "t.cens"]),
                      "_teval", as.numeric(crossing[row, "t.eval"]), "_un.png", sep = ""), 
                dpi = 600, width = 6, height = 7.5, unit = "in")
grid::grid.draw(sim.out[["plot"]])
dev.off()

print("FINISHED")
