###
### Create forest plot for discrimination
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
library(forestplot)
### Create discrimination forest plot
create_forestplot_poster <- function(gender, age_knots, caltime){
  
  ### Set gender character version
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  print(paste("age_knots = ", age_knots))
  print(paste("caltime = ", caltime))
  
  ### Read in validation dataset
  df_valid <- readRDS(paste("data/p4/dfs_valid_", gender_char, ".rds", sep = ""))[[1]]
  
  ### Discrimination in entire cohort
  discrimination_entire_cohort <- readRDS(paste("data/p4/prototype3_discrim_d1v1_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  discrimination_entire_cohort <- c("N" = nrow(df_valid), discrimination_entire_cohort)
  
  ### Discrimination by age
  discrimination_by_age <- readRDS(paste("data/p4/prototype3_discrim_table_d1v1_age_cat_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  
  ### Discrimination by region
  discrimination_by_region <- readRDS(paste("data/p4/prototype3_discrim_table_d1v1_region_", gender, "_nk", age_knots, "_caltime", as.numeric(caltime), ".rds", sep = ""))
  
  ### Combine
  discrimination_combined <- rbind(discrimination_entire_cohort[c(1,2,6,7)],
                                   discrimination_by_region[,1:4],
                                   discrimination_by_age[,1:4]) 
  
  discrimination_combined <- 
    dplyr::select(discrimination_combined, N, Cstat, Cstat_lower, Cstat_upper) |>
    dplyr::mutate(mean = Cstat, 
                  Cstat = round(Cstat, 3), 
                  grouping = c("Entire cohort", rownames(discrimination_combined)[-1])) |>
    dplyr::rename(lower = Cstat_lower, upper = Cstat_upper)
  
  plot_out <- forestplot::forestplot(discrimination_combined, 
                                     labeltext = c(grouping, N, Cstat), 
                                     boxsize = 0.2, 
                                     zero = discrimination_combined$Cstat[1],
                                     graphwidth = grid::unit(2.5, "cm")) |>
    forestplot::fp_add_header("grouping" = "Grouping",
                              "N" = "N",
                              "Cstat" = "Cstat")
  
  ### Save plot
  ragg::agg_png(paste("figures/p4/aaaa_POSTER_UKAIRS.png", sep = ""), 
                width = 1, height = 0.7, scaling = 1/8, unit = "in", res = 600)
  plot(plot_out)
  dev.off()
  
}

### Create forest plots
create_forestplot_poster(gender = 2, age_knots = 4, caltime = TRUE)




