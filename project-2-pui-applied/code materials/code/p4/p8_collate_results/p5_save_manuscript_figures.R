###
### Need to create high quality images...
###


### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

library(png)
library(grid)
library(ggplot2)
library(gridExtra)

###
### Create Figure 3
###
plots <- lapply(ll <- paste("figures/p4/prototype3_calib_ph_MxM_", c(2,1), "_nk", 4, "_caltime", 1, ".png", sep = ""),
                function(x){
                  img <- as.raster(readPNG(x))
                  grid.newpage()
                  rasterGrob(img, interpolate = FALSE)
                })

ggsave("figures/p4/Figure3.png", marrangeGrob(grobs = plots, nrow = 1, ncol = 2, top = NULL), height = 3, width = 6, unit = "in", dpi = 600)
ggsave("figures/p4/Figure3.pdf", marrangeGrob(grobs = plots, nrow = 1, ncol = 2, top = NULL), height = 2, width = 4)


###
### Figure 4 is created in the file p4_validate/p3.2_prototype3_discrimination_forestplot.R
###

###
### Create Figure 5
###
plots <- lapply(ll <- paste("figures/p4/temporalv_initial_prototype3_calib_ph_d1v1_", 2, "_nk", 4, "_caltime", 1,"_tfup", round(c(1,3,5)*365.25), "_teval", 5, ".png", sep = ""),
       function(x){
         img <- as.raster(readPNG(x))
         grid.newpage()
         rasterGrob(img, interpolate = FALSE)
       })

ggsave("figures/p4/Figure5.png", marrangeGrob(grobs = plots, nrow = 1, ncol = 3, top = NULL), height = 3, width = 9, unit = "in", dpi = 600)
ggsave("figures/p4/Figure5.pdf", marrangeGrob(grobs = plots, nrow = 1, ncol = 3, top = NULL), height = 2, width = 6)


###
### Create Figure 6
###
plots <- lapply(ll <- paste("figures/p4/temporalv_intervention_prototype3_calib_ph_d1v1_", 2, "_nk", 4, "_caltime", 1,"_tfup", round(1*365.25), "_teval", 5, "_", c("sbp", "bmi", "nonhdl", "smoking"), ".png", sep = ""),
                function(x){
                  img <- as.raster(readPNG(x))
                  grid.newpage()
                  rasterGrob(img, interpolate = FALSE)
                })

ggsave("figures/p4/Figure6.png", marrangeGrob(grobs = plots, nrow = 2, ncol = 2, top = NULL), height = 6, width = 6, unit = "in", dpi = 600)
ggsave("figures/p4/Figure6.pdf", marrangeGrob(grobs = plots, nrow = 2, ncol = 2, top = NULL), height = 3, width = 3)
