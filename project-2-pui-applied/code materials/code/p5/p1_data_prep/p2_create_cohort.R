###
### Finalise cohort creation (initial cohort already created in data extraction code base)
### Remove gender == 3, and re-factor smoking.
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("./")
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
sink("data/p5/aaa_cohort_file_structure.txt")
print(str(cohort_null))
sink()

### Save the cohorts
saveRDS(cohort, "data/p5/cohort_prototype3.rds")
print(paste("IMAGE SAVED", Sys.time()))

### Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

### Render the Rmarkdown document
print("render to markdown")
# rmarkdown::render("code/p4/p1_data_prep/p2_table1.Rmd")

### Finished
print("FINISHED")