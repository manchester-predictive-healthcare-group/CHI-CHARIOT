### Create Table 1 of data

Merge extracted variables into baseline cohort

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

print("render to markdown")

### Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

### Render the Rmarkdown document
rmarkdown::render("code/p2_cohort_baseline/p3_table1_cohort_baseline.Rmd")

### Finished
print("FINISHED")