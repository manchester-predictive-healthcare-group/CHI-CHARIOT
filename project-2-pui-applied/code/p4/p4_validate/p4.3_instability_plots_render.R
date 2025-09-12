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
rmarkdown::render("code/p4/p4_validate/p4.3_instability_plots.Rmd")
### Finished
print("FINISHED")