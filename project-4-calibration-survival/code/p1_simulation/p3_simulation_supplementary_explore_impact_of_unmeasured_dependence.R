###
### Supplementary Analyses
###

### We are finding that introducing unmeasured confounding has less of an impact on scenarios with more complex
### DGMs. Want to verify this in a simulation.

### Clear workspace
rm(list=ls())
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project4")
set.seed(505)
library(survival)

### Load functions
source("R/functions.R")

### Set m
m <- 500000

### Extract time point
t.eval <- 0.25

### Generate baseline data
cov <- data.frame("x1" = rnorm(m, 0, 1),
                  "x2" = rnorm(m, 0, 1),
                  "x3" = rnorm(m, 0, 1),
                  "x4" = rnorm(m, 0, 1),
                  "x5" = rnorm(m, 0, 1),
                  "x.cens" = rnorm(m, 0, 1))

### Simulate data (with and without censoring)
print("sim dat")
dat1.devel <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.05, b.x2 = 0.05, b.x3 = 0.05, b.x4 = 0.05, b.x5 = 0.05, 
                                         b.x.cens = 1, seed = 101, cens = FALSE)
dat2.devel <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.1, b.x2 = 0.1, b.x3 = 0.1, b.x4 = 0.1, b.x5 = 0.1, 
                                         b.x.cens = 1, seed = 101, cens = FALSE)
dat3.devel <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.25, b.x2 = 0.25, b.x3 = 0.25, b.x4 = 0.25, b.x5 = 0.25, 
                                         b.x.cens = 1, seed = 101, cens = FALSE)
dat4.devel <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.5, b.x2 = 0.5, b.x3 = 0.5, b.x4 = 0.5, b.x5 = 0.5, 
                                         b.x.cens = 1, seed = 101, cens = FALSE)
dat5.devel <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.75, b.x2 = 0.75, b.x3 = 0.75, b.x4 = 0.75, b.x5 = 0.75, 
                                         b.x.cens = 1, seed = 101, cens = FALSE)
dat6.devel <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 =1, b.x2 =1, b.x3 =1, b.x4 =1, b.x5 =1, 
                                         b.x.cens = 1, seed = 101, cens = FALSE)
dat7.devel <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 3, b.x2 = 3, b.x3 = 3, b.x4 = 3, b.x5 = 3, 
                                         b.x.cens = 1, seed = 101, cens = FALSE)
dat8.devel <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 5, b.x2 = 5, b.x3 = 5, b.x4 = 5, b.x5 = 5, 
                                         b.x.cens = 1, seed = 101, cens = FALSE)

dat1.valid <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.05, b.x2 = 0.05, b.x3 = 0.05, b.x4 = 0.05, b.x5 = 0.05, 
                                         b.x.cens = 1, seed = 101)
dat2.valid <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.1, b.x2 = 0.1, b.x3 = 0.1, b.x4 = 0.1, b.x5 = 0.1, 
                                         b.x.cens = 1, seed = 101)
dat3.valid <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.25, b.x2 = 0.25, b.x3 = 0.25, b.x4 = 0.25, b.x5 = 0.25, 
                                         b.x.cens = 1, seed = 101)
dat4.valid <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.5, b.x2 = 0.5, b.x3 = 0.5, b.x4 = 0.5, b.x5 = 0.5, 
                                         b.x.cens = 1, seed = 101)
dat5.valid <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 0.75, b.x2 = 0.75, b.x3 = 0.75, b.x4 = 0.75, b.x5 = 0.75, 
                                         b.x.cens = 1, seed = 101)
dat6.valid <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 =1, b.x2 =1, b.x3 =1, b.x4 =1, b.x5 =1, 
                                         b.x.cens = 1, seed = 101)
dat7.valid <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 3, b.x2 = 3, b.x3 = 3, b.x4 = 3, b.x5 = 3, 
                                         b.x.cens = 1, seed = 101)
dat8.valid <- simulate_DGM_supplementary(m = m, cov = cov, b.x1 = 5, b.x2 = 5, b.x3 = 5, b.x4 = 5, b.x5 = 5, 
                                         b.x.cens = 1, seed = 101)
### Run the simulation
print(paste("run sim 1", Sys.time()))
sim.out1 <- run_sim_large_sample(dat1.devel, dat1.valid, t.eval = t.eval, cumhaz.t = t.eval, t.cens = TRUE, model = "supplementary")
print(paste("run sim 2", Sys.time()))
sim.out2 <- run_sim_large_sample(dat2.devel, dat2.valid, t.eval = t.eval, cumhaz.t = t.eval, t.cens = TRUE, model = "supplementary")
print(paste("run sim 3", Sys.time()))
sim.out3 <- run_sim_large_sample(dat3.devel, dat3.valid, t.eval = t.eval, cumhaz.t = t.eval, t.cens = TRUE, model = "supplementary")
print(paste("run sim 4", Sys.time()))
sim.out4 <- run_sim_large_sample(dat4.devel, dat4.valid, t.eval = t.eval, cumhaz.t = t.eval, t.cens = TRUE, model = "supplementary")
print(paste("run sim 5", Sys.time()))
sim.out5 <- run_sim_large_sample(dat5.devel, dat5.valid, t.eval = t.eval, cumhaz.t = t.eval, t.cens = TRUE, model = "supplementary")
print(paste("run sim 6", Sys.time()))
sim.out6 <- run_sim_large_sample(dat6.devel, dat6.valid, t.eval = t.eval, cumhaz.t = t.eval, t.cens = TRUE, model = "supplementary")
print(paste("run sim 7", Sys.time()))
sim.out7 <- run_sim_large_sample(dat7.devel, dat7.valid, t.eval = t.eval, cumhaz.t = t.eval, t.cens = TRUE, model = "supplementary")
print(paste("run sim 6", Sys.time()))
sim.out8 <- run_sim_large_sample(dat8.devel, dat8.valid, t.eval = t.eval, cumhaz.t = t.eval, t.cens = TRUE, model = "supplementary")


### Save object as rds file
saveRDS(sim.out1, "data/supplementary_simulation_1.rds")
saveRDS(sim.out2, "data/supplementary_simulation_2.rds")
saveRDS(sim.out3, "data/supplementary_simulation_3.rds")
saveRDS(sim.out4, "data/supplementary_simulation_4.rds")
saveRDS(sim.out5, "data/supplementary_simulation_5.rds")
saveRDS(sim.out6, "data/supplementary_simulation_6.rds")
saveRDS(sim.out7, "data/supplementary_simulation_7.rds")
saveRDS(sim.out8, "data/supplementary_simulation_8.rds")

Cairo::CairoPNG("figures/supplementary_simulation_1.png", dpi = 600, width = 6, height = 7.5, unit = "in")
grid::grid.draw(sim.out1[["plot"]])
dev.off()
Cairo::CairoPNG("figures/supplementary_simulation_2.png", dpi = 600, width = 6, height = 7.5, unit = "in")
grid::grid.draw(sim.out2[["plot"]])
dev.off()
Cairo::CairoPNG("figures/supplementary_simulation_3.png", dpi = 600, width = 6, height = 7.5, unit = "in")
grid::grid.draw(sim.out3[["plot"]])
dev.off()
Cairo::CairoPNG("figures/supplementary_simulation_4.png", dpi = 600, width = 6, height = 7.5, unit = "in")
grid::grid.draw(sim.out4[["plot"]])
dev.off()
Cairo::CairoPNG("figures/supplementary_simulation_5.png", dpi = 600, width = 6, height = 7.5, unit = "in")
grid::grid.draw(sim.out5[["plot"]])
dev.off()
Cairo::CairoPNG("figures/supplementary_simulation_6.png", dpi = 600, width = 6, height = 7.5, unit = "in")
grid::grid.draw(sim.out6[["plot"]])
dev.off()
Cairo::CairoPNG("figures/supplementary_simulation_7.png", dpi = 600, width = 6, height = 7.5, unit = "in")
grid::grid.draw(sim.out7[["plot"]])
dev.off()
Cairo::CairoPNG("figures/supplementary_simulation_8.png", dpi = 600, width = 6, height = 7.5, unit = "in")
grid::grid.draw(sim.out8[["plot"]])
dev.off()
print(paste("FINISHED", Sys.time()))