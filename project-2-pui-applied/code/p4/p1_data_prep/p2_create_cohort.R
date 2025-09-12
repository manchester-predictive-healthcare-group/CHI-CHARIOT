### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Define filepath to file directory system containing extracted data.
common.data.dir <- file.path("..", "..")

### Load packages
library(mice)
library(survival)
library(rms)
library(dplyr)
library(tidyr)
library(ggplot2)

### Extract cohort
cohort <- readRDS(file.path(common.data.dir, "Aurum_Jun2021_extract/data/extraction/cohort_baseline/cohort_var.rds"))
cohort <- data.table::as.data.table(cohort)
cohort <- filter(cohort, gender != 3)

### Redefine smoking to be three levels, Never, Ex, Current
cohort$smoking <- forcats::fct_recode(cohort$smoking, Current = "Light", Current = "Moderate", Current = "Heavy")

### Save structure output for reproducibility
cohort_null <- cohort[0, ]
sink("data/aaa_cohort_file_structure.txt")
print(str(cohort_null))
sink()

###################################
### Read in component variables ###
###################################

### Write a function to do this for a variable name
read_component <- function(varname){
  readRDS(file.path(common.data.dir, paste("Aurum_Jun2021_extract/data/extraction/cohort_baseline/var_", varname, ".rds", sep = "")))
}

component_height <- read_component("component_height")
component_weight <- read_component("component_weight")
component_bmi <- read_component("component_bmi")
component_total_chol <- read_component("component_total_chol")
component_hdl <- read_component("component_hdl")
component_ldl <- read_component("component_ldl")
component_total_chol_hdl_ratio <- read_component("component_total_chol_hdl_ratio")
component_triglycerides <- read_component("component_triglycerides")

#########################
### Merge with cohort ###
#########################

### Merge all together
cohort <- Reduce(function(x, y) merge(x, y, by.x = "patid", by.y = "patid", all.x = TRUE),
                 list(cohort, 
                      dplyr::select(component_total_chol,-medcodeid), 
                      dplyr::select(component_hdl,-medcodeid), 
                      dplyr::select(component_ldl,-medcodeid), 
                      component_triglycerides))

### Reduce cohort
cohort_reduced <- dplyr::select(cohort, patid, gender, bmi, sbp, cholhdl_ratio, 
                                component_total_chol, component_hdl, component_ldl, component_triglycerides)
str(cohort_reduced)

### RENAME THE COMPONENT VARIABLES of interest
cohort <- dplyr::rename(cohort, hdl = component_hdl, ldl = component_ldl, cholesterol = component_total_chol, triglycerides = component_triglycerides)

### Apply limits to triglyercides
# lower_bound <- as.numeric(quantile(cohort$triglycerides, p=0.001, na.rm = TRUE))
# upper_bound <- as.numeric(quantile(cohort$triglycerides, p=0.999, na.rm = TRUE))
# sum(is.na(cohort$triglycerides))
# cohort <- dplyr::mutate(cohort, triglycerides = dplyr::case_when(triglycerides > lower_bound & triglycerides < upper_bound ~ triglycerides,
#                                                                  TRUE ~ NA))

### Remove triglycerides from cohort as decided not to work with this var
cohort <- dplyr::select(cohort, - triglycerides)

### Save the cohorts
saveRDS(cohort, "data/p4/cohort_prototype3.rds")
saveRDS(cohort_reduced, "data/p4/cohort_components.rds")
print(paste("IMAGE SAVED", Sys.time()))

### Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

### Render the Rmarkdown document
print("render to markdown")
# rmarkdown::render("code/p4/p1_data_prep/p2_table1_cohort_components.Rmd")

### Finished
print("FINISHED")