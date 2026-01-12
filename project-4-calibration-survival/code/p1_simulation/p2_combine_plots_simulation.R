###
### Simulated data example
###
options(scipen = 999)
rm(list=ls())
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project4")
set.seed(505)
library(survival)

### Load functions
source("R/functions.R")

### Define m and split
m <- 1000000

### Write function to create plot
### Plots will go in a 3x3 grid, with DGM.y vs model

### Note, in the simulations "PV-ipcw" refers to grouping pseudo-value by the ipcw's.
### We are not reporting that in the paper, and are referring to the binder approach as "PV-IPCW", so
### are changing the name. Making a note of this, as to not get confused by PV-ipcw moving forwards.
combine_plot <- function(DGM.cens, b, split, t.cens, t.eval, unmeas_conf){
  
  # DGM.cens <- 4
  # b <- "03"
  # split <- 2
  # t.cens <- 1
  # t.eval <- 1
  # unmeas_conf <- FALSE
  
  ### Create object to read in output
  plotdat <- vector("list", 3)
  plotdat <- lapply(plotdat, function(x) {return(vector("list", 3))})
  
  # plotdat <- vector("list", 2)
  # plotdat <- lapply(plotdat, function(x) {return(vector("list", 1))})
  
  ### Read in output
  for (i in 1:3){
    for (j in 1:3){
      
      ### Read in plot data
      if (unmeas_conf == FALSE){
        plotdat[[i]][[j]] <- readRDS(paste("data/simulation_m", m, 
                                           "_b", b,
                                           "_y", i,
                                           "_c", DGM.cens,
                                           "_model", j,
                                           "_split", split,  
                                           "_tcens", t.cens,
                                           "_teval", t.eval, ".rds", sep = ""))[["plotdata"]]
      } else {
        plotdat[[i]][[j]] <- readRDS(paste("data/simulation_m", m, 
                                           "_b", b,
                                           "_y", i,
                                           "_c", DGM.cens,
                                           "_model", j,
                                           "_split", split,  
                                           "_tcens", t.cens,
                                           "_teval", t.eval, "_un.rds", sep = ""))[["plotdata"]]
      }
      
      ### Get upper and lower limits
      lim_min <- quantile(dplyr::pull(subset(plotdat[[i]][[j]], Calibration == "true"), pred), p = 0.01)
      lim_max <- quantile(dplyr::pull(subset(plotdat[[i]][[j]], Calibration == "true"), pred), p = 0.99)
      
      ### Reduce to 1000 individuals for plotting to reduce output file size and order the Calibration variable
      ### Change of name from "PV-binder" to "PV-IPCW"
      plotdat[[i]][[j]] <- subset(plotdat[[i]][[j]], id %in% 1:1000) |>
        subset(Calibration != "PV-ipcw") |>
        dplyr::mutate(Calibration = dplyr::case_when(Calibration == "PV" ~ "PV-group-risk", 
                                                     Calibration == "PV-binder" ~ "PV-IPCW", 
                                                     TRUE ~ Calibration)) |>
        dplyr::mutate(Calibration = factor(Calibration, levels = c("true", "PH", "PV-group-risk", "PV-IPCW",  "IPCW-x1", "IPCW-flex")))
      
      ### Replace with a plot
      plotdat[[i]][[j]] <- ggplot2::ggplot(data = plotdat[[i]][[j]]) +
        ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, color = Calibration, linetype = Calibration), linewidth = 0.7) +
        ggplot2::geom_point(ggplot2::aes(x = pred, y = pred.obs), 
                            col = grDevices::rgb(0,0,0,alpha=0) , 
                            data = subset(plotdat[[i]][[j]], Calibration == "true")) +
        ggplot2::xlim(lim_min, lim_max) + ggplot2::ylim(lim_min, lim_max) +
        ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + 
        ggplot2::theme(legend.position = "none") #+ 
      #ggplot2::theme(text = element_text(size = 3))
      
      ### Add appropriate axis labels
      if (i == 1 & j == 1){
        plotdat[[i]][[j]] <- plotdat[[i]][[j]] + 
          ggplot2::ylab("(O-DGM-SL)\nPredicted-observed risk")
      } else if (i == 2 & j == 1){
        plotdat[[i]][[j]] <- plotdat[[i]][[j]] + 
          ggplot2::ylab("(O-DGM-ML)\nPredicted-observed risk")
      } else if (i == 3 & j == 1){
        plotdat[[i]][[j]] <- plotdat[[i]][[j]] + 
          ggplot2::ylab("(O-DGM-poly)\nPredicted-observed risk") + 
          ggplot2::xlab("Predicted risk\n(Model-SL)")
      } else if (i == 3 & j == 2){
        plotdat[[i]][[j]] <- plotdat[[i]][[j]] + 
          ggplot2::xlab("Predicted risk\n(Model-ML)")
      } else if (i == 3 & j == 3){
        plotdat[[i]][[j]] <- plotdat[[i]][[j]] + 
          ggplot2::xlab("Predicted risk\n(Model-flex)")
      }
      
      ### Add marginal
      plotdat[[i]][[j]] <- ggExtra::ggMarginal(plotdat[[i]][[j]], type = "density", margins = "x", size = 6)
      
    }
  }
  
  ### Extract into a single list
  plotdat <- do.call("c", plotdat)
  
  ### Arrange into a grid
  plots.comb <- gridExtra::arrangeGrob(grobs = plotdat,
                                       layout_matrix = base::matrix(seq_len(9),nrow = 3, ncol = 3, byrow = TRUE))
  # plots.comb <- gridExtra::arrangeGrob(grobs = plotdat,
  #                                      layout_matrix = base::matrix(seq_len(2),nrow = 1, ncol = 2, byrow = TRUE))
  ### Create a plot with a legend, so we can extract and save the legend
  temp <- readRDS(paste("data/simulation_m", m, 
                        "_b03",
                        "_y", 1,
                        "_c", 1,
                        "_model", 1,
                        "_split", 2,  
                        "_tcens", 2,
                        "_teval", 1, ".rds", sep = ""))[["plotdata"]] |>
    subset(Calibration != "PV-ipcw") |>
    ### Change of name from "PV-binder" to "PV-IPCW"
    dplyr::mutate(Calibration = dplyr::case_when(Calibration == "PV" ~ "PV-group-risk", 
                                                 Calibration == "PV-binder" ~ "PV-IPCW", 
                                                 TRUE ~ Calibration)) |>
    dplyr::mutate(Calibration = factor(Calibration, levels = c("true", "PH", "PV-group-risk", "PV-IPCW",  "IPCW-x1", "IPCW-flex")))
  
  temp.plot <- ggplot2::ggplot(data = temp) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, color = Calibration, linetype = Calibration), linewidth = 0.9) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) #+
  #ggplot2::theme(text = element_text(size = 3))
  legend.saved <- ggpubr::get_legend(temp.plot)
  
  ### Add to plot
  plots.comb <- gridExtra::arrangeGrob(plots.comb, legend.saved, nrow = 2, heights = c(16,1))
  
  ### Return
  return(plots.comb)
  
}

# DGM.cens <- 4
# b <- "03"
# split <- 2
# t.cens <- 1
# t.eval <- 1
# unmeas_conf <- FALSE
# 
# ### Create plot
# plot.comb <- combine_plot(DGM.cens, b, split, t.cens, t.eval, FALSE)
# 
# ### Save pdf
# pdf(paste("figures/aaa_TEMP.pdf", sep = ""),
#     width = 9, height = 9)
# grid::grid.draw(plot.comb)
# dev.off()
# 
# ### Save png
# ragg::agg_png(paste("figures/aaa_TEMP.png", sep = ""),
#               units = "in",
#               res = 600,
#               width = 3, height = 3, scaling = 0.33)
# grid::grid.draw(plot.comb)
# dev.off()

## Create plots for coefficient = 0.3
for (DGM.cens in c(1,2,3,4)){
  for (b in c("03")){
    for (split in c(2)){
      for (t.cens in c(1,2)){
        for (t.eval in c(1,2)){
          print(paste(t.cens, t.eval, DGM.cens, b, split))
          
          ### Create plot
          plot.comb <- combine_plot(DGM.cens, b, split, t.cens, t.eval, FALSE)
          
          ### Save pdf
          pdf(paste("figures/simulation_m", m,
                    "_b", b,
                    "_split", split,
                    "_c", DGM.cens,
                    "_tcens", t.cens,
                    "_teval", t.eval,  ".pdf", sep = ""),
              width = 9, height = 9)
          grid::grid.draw(plot.comb)
          dev.off()
          
          ### Save png
          ragg::agg_png(paste("figures/simulation_m", m,
                              "_b", b,
                              "_split", split,
                              "_c", DGM.cens,
                              "_tcens", t.cens,
                              "_teval", t.eval,  ".png", sep = ""),
                        units = "in", 
                        res = 1000,
                        width = 1, height = 1, scaling = 0.1)
          grid::grid.draw(plot.comb)
          dev.off()
          
          ### Save high res images for Figures 1, 2, 3, 4
          if (b == "03" & split == 2 & t.cens == 2 & t.eval == 1){
            if (DGM.cens == 1){
              ### Save pdf
              png(paste("figures/Figure1.png", sep = ""),
                  width = 9, height = 9, res = 900, units = "in")
              grid::grid.draw(plot.comb)
              dev.off()
            } else if (DGM.cens == 2){
              ### Save pdf
              png(paste("figures/Figure2.png", sep = ""),
                  width = 9, height = 9, res = 900, units = "in")
              grid::grid.draw(plot.comb)
              dev.off()
            } else if (DGM.cens == 3){
              ### Save pdf
              png(paste("figures/Figure3.png", sep = ""),
                  width = 9, height = 9, res = 900, units = "in")
              grid::grid.draw(plot.comb)
              dev.off()
            } else if (DGM.cens == 4){
              ### Save pdf
              png(paste("figures/Figure4.png", sep = ""),
                  width = 9, height = 9, res = 900, units = "in")
              grid::grid.draw(plot.comb)
              dev.off()
            } 
          }
          
          ### close loops
        }
      }
    }
  }
}

## Create plots for coefficient = 0.5
for (DGM.cens in c(1,2,3,4)){
  for (b in c("05")){
    for (split in c(2)){
      for (t.cens in c(2)){
        for (t.eval in c(1)){
          print(paste(t.cens, t.eval, DGM.cens, b, split))
          
          ### Create plot
          plot.comb <- combine_plot(DGM.cens, b, split, t.cens, t.eval, FALSE)
          
          ### Save pdf
          pdf(paste("figures/simulation_m", m,
                    "_b", b,
                    "_split", split,
                    "_c", DGM.cens,
                    "_tcens", t.cens,
                    "_teval", t.eval,  ".pdf", sep = ""),
              width = 9, height = 9)
          grid::grid.draw(plot.comb)
          dev.off()
          
          ### Save png
          ragg::agg_png(paste("figures/simulation_m", m,
                              "_b", b,
                              "_split", split,
                              "_c", DGM.cens,
                              "_tcens", t.cens,
                              "_teval", t.eval,  ".png", sep = ""),
                        units = "in", 
                        res = 1000,
                        width = 1, height = 1, scaling = 0.1)
          grid::grid.draw(plot.comb)
          dev.off()
          
          ### Save high res images for Figures 1, 2, 3, 4
          if (b == "03" & split == 2 & t.cens == 2 & t.eval == 1){
            if (DGM.cens == 1){
              ### Save pdf
              png(paste("figures/Figure1.png", sep = ""),
                  width = 9, height = 9, res = 900, units = "in")
              grid::grid.draw(plot.comb)
              dev.off()
            } else if (DGM.cens == 2){
              ### Save pdf
              png(paste("figures/Figure2.png", sep = ""),
                  width = 9, height = 9, res = 900, units = "in")
              grid::grid.draw(plot.comb)
              dev.off()
            } else if (DGM.cens == 3){
              ### Save pdf
              png(paste("figures/Figure3.png", sep = ""),
                  width = 9, height = 9, res = 900, units = "in")
              grid::grid.draw(plot.comb)
              dev.off()
            } else if (DGM.cens == 4){
              ### Save pdf
              png(paste("figures/Figure4.png", sep = ""),
                  width = 9, height = 9, res = 900, units = "in")
              grid::grid.draw(plot.comb)
              dev.off()
            } 
          }
          
          ### close loops
        }
      }
    }
  }
}

## Create plots for unmeasured dependence
for (DGM.cens in c(1,2,3,4)){
  for (b in c("03")){
    for (split in c(2)){
      for (t.cens in c(2)){
        for (t.eval in c(1)){
          print(paste(t.cens, t.eval, DGM.cens, b, split))
          
          ### Create plot
          plot.comb <- combine_plot(DGM.cens, b, split, t.cens, t.eval, TRUE)
          
          ### Save pdf
          pdf(paste("figures/simulation_m", m,
                    "_b", b,
                    "_split", split,
                    "_c", DGM.cens,
                    "_tcens", t.cens,
                    "_teval", t.eval,  "_un.pdf", sep = ""),
              width = 9, height = 9)
          grid::grid.draw(plot.comb)
          dev.off()
          
          ### Save png
          ragg::agg_png(paste("figures/simulation_m", m,
                              "_b", b,
                              "_split", split,
                              "_c", DGM.cens,
                              "_tcens", t.cens,
                              "_teval", t.eval,  "_un.png", sep = ""),
                        units = "in",
                        res = 1000,
                        width = 1, height = 1, scaling = 0.1)
          grid::grid.draw(plot.comb)
          dev.off()
          
          ### Save high res images for Figure 5
          if (b == "03" & split == 2 & t.cens == 2 & t.eval == 1){
            if (DGM.cens == 1){
              ### Save pdf
              png(paste("figures/Figure5.png", sep = ""),
                  width = 9, height = 9, res = 900, units = "in")
              grid::grid.draw(plot.comb)
              dev.off()
            } 
          }
          
          ### close loops
          
        }
      }
    }
  }
}
print(paste("FINISHED", Sys.time()))
