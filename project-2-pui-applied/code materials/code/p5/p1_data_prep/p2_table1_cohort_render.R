### Create Table 1 of data

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("PATH/TO/ROOT")
getwd()

print("render to markdown")

### Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

### Render the Rmarkdown document
rmarkdown::render("code/p5/p1_data_prep/p2_table1.Rmd")
rmarkdown::render("code/p5/p1_data_prep/p2_table1_word.Rmd")

### Finished
print("FINISHED")