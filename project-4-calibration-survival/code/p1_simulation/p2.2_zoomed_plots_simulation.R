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
combine_plot <- function(DGM.cens, unmeas_conf){
  
  # DGM.cens <- 3
  b <- "03"
  split <- 2
  t.cens <- 2
  t.eval <- 1
  DGM.y <- 3
  unmeas_conf = FALSE
  
  ### Create object to read in output
  plotdat <- vector("list", 2)

  ### Read in output
  for (model in 1:2){
  
      ### Read in plot data
      if (unmeas_conf == FALSE){
        plotdat[[model]] <- readRDS(paste("data/simulation_m", m, 
                                           "_b", b,
                                           "_y", DGM.y,
                                           "_c", DGM.cens,
                                           "_model", model,
                                           "_split", split,  
                                           "_tcens", t.cens,
                                           "_teval", t.eval, ".rds", sep = ""))[["plotdata"]]
      } else {
        plotdat[[model]] <- readRDS(paste("data/simulation_m", m, 
                                          "_b", b,
                                          "_y", DGM.y,
                                          "_c", DGM.cens,
                                          "_model", model,
                                          "_split", split,  
                                          "_tcens", t.cens,
                                          "_teval", t.eval, ".rds", sep = ""))[["plotdata"]]
      }
      
      ### Reduce to 1000 individuals for plotting to reduce output file size and order the Calibration variable
      plotdat[[model]] <- subset(plotdat[[model]], id %in% 1:1000) |>
        dplyr::mutate(Calibration = dplyr::case_when(Calibration == "PV" ~ "PV-group-risk", 
                                                     Calibration == "PV-ipcw" ~ "PV-group-ipcw",
                                                     TRUE ~ Calibration)) |>
        dplyr::mutate(Calibration = factor(Calibration, levels = c("true", "PH", "PV-group-risk", "PV-group-ipcw", "PV-binder",  "IPCW-x1", "IPCW-flex")))
      
      ### Replace with a plot
      plotdat[[model]] <- ggplot2::ggplot(data = plotdat[[model]]) +
        ggplot2::geom_abline(slope = 1, intercept = 0, lty = "dashed") +
        ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, color = Calibration, linetype = Calibration), linewidth = 0.7) +
        ggplot2::geom_point(ggplot2::aes(x = pred, y = pred.obs), 
                            col = grDevices::rgb(0,0,0,alpha=0) , 
                            data = subset(plotdat[[model]], Calibration == "true")) +
        ggplot2::xlim(c(0.3, 0.5)[model], c(0.4,0.6)[model]) + ggplot2::ylim(c(0.3, 0.5)[model], c(0.4,0.6)[model]) +
        ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + 
        ggplot2::ylab("Predicted-observed risk (O-DGM-QU)") + 
        ggplot2::xlab(c("Model = SL", "Model = ML")[model]) + 
        ggplot2::theme(legend.position = "none") #+ 
      #ggplot2::theme(text = element_text(size = 3))
      
      ### Add marginal
      plotdat[[model]] <- ggExtra::ggMarginal(plotdat[[model]], type = "density", margins = "x", size = 6)
      
    }
  
  ### Arrange into a grid
  plots.comb <- gridExtra::arrangeGrob(grobs = plotdat,
                                       layout_matrix = base::matrix(seq_len(2), nrow = 1, ncol = 2, byrow = TRUE))
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
    dplyr::mutate(Calibration = dplyr::case_when(Calibration == "PV" ~ "PV-group-risk", 
                                                 Calibration == "PV-ipcw" ~ "PV-group-ipcw",
                                                 TRUE ~ Calibration)) |>
    dplyr::mutate(Calibration = factor(Calibration, levels = c("true", "PH", "PV-group-risk", "PV-group-ipcw", "PV-binder",  "IPCW-x1", "IPCW-flex")))
  
  temp.plot <- ggplot2::ggplot(data = temp) +
    ggplot2::geom_line(ggplot2::aes(x = pred, y = pred.obs, color = Calibration, linetype = Calibration), linewidth = 0.9) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) #+
  #ggplot2::theme(text = element_text(size = 3))
  legend.saved <- ggpubr::get_legend(temp.plot)
  
  ### Add to plot
  plots.comb <- gridExtra::arrangeGrob(plots.comb, legend.saved, nrow = 2, heights = c(8,1))
  
  ### Return
  return(plots.comb)
  
}

### Create plot
plot.comb <- combine_plot(3, FALSE)

### Save pdf
png(paste("figures/Figure3_zoom.png", sep = ""),
    width = 6, height = 3, res = 900, units = "in")
grid::grid.draw(plot.comb)
dev.off()

str(plot.comb)
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

