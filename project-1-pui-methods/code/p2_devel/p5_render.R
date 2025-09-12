### Create Table 1 of data

### This is going to test the functions for deriving all the variables!

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
rmarkdown::render("code/p2_devel/p5_worked_example.Rmd")
rmarkdown::render("code/p2_devel/pX_collate_results.Rmd")
rmarkdown::render("code/p2_devel/pX_collate_tables_word.Rmd")

### Finished
print("FINISHED")