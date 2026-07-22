### Create Table 1 of data

### This is going to test the functions for deriving all the variables!

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

print("render to markdown")

### Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

### Render the Rmarkdown document
rmarkdown::render("code/p4/p3_devel/p3_prototype3_assess_nonlinear.Rmd")
rmarkdown::render("code/p4/p3_devel/p3_prototype3_assess_nonlinear_4knots.Rmd")
rmarkdown::render("code/p4/p3_devel/p3_prototype3_assess_nonlinear_4knots_caltime.Rmd")
### Finished
print("FINISHED")