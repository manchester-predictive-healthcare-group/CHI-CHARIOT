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
#rmarkdown::render("code/p4/p8_collate_results/p1_collate_tables_word.Rmd")
rmarkdown::render("code/p4/p8_collate_results/p2_collate_figures_and_tables.Rmd")
rmarkdown::render("code/p4/p8_collate_results/p3_copare_models_age_knots.Rmd")
rmarkdown::render("code/p4/p8_collate_results/p4_copare_models_caltime.Rmd")

### Finished
print("FINISHED")