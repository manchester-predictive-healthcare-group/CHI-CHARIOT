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
rmarkdown::render("code/p5/p3_devel/p2_assess_nonlinear.Rmd")
### Finished
print("FINISHED")