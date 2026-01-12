### Create Table 1 of data

### This is going to test the functions for deriving all the variables!

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project4/")

print("render to markdown")

### Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

### Render the Rmarkdown document
rmarkdown::render("code/simulation/p5_all_tables.Rmd")

### Finished
print("FINISHED")