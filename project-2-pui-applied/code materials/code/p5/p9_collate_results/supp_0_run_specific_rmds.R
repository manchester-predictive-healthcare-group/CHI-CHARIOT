# run_all_rmds.R
# Knits all .Rmd files in a specified directory

# Set wd
setwd()
getwd()

# ---- Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

# ---- Run Rmds ----
rmarkdown::render("code/p5/p9_collate_results/supp_3_main_validation_results.Rmd", envir = new.env())
rmarkdown::render("code/p5/p9_collate_results/supp_5_imp_diagnostics.Rmd", envir = new.env())
rmarkdown::render("code/p5/p9_collate_results/supp_6_imp_assess_mar_and_support.Rmd", envir = new.env())
rmarkdown::render("code/p5/p9_collate_results/supp_7_mean_imputation.Rmd", envir = new.env())
rmarkdown::render("code/p5/p9_collate_results/supp_8_mnar.Rmd", envir = new.env())
rmarkdown::render("code/p5/p9_collate_results/supp_4_simulation_validate_xu_extension.Rmd", envir = new.env())