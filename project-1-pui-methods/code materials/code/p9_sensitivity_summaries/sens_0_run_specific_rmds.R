# run_all_rmds.R
# Knits all .Rmd files in a specified directory

# Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

# ---- Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

# ---- Run Rmds ----
rmarkdown::render("code/p9_sensitivity_summaries/flexible_pui_supp_4_summarise_missingness_impact.Rmd", envir = new.env())
rmarkdown::render("code/p9_sensitivity_summaries/flexible_pui_supp_5_treatment_effect_modified.Rmd", envir = new.env())
rmarkdown::render("code/p9_sensitivity_summaries/flexible_pui_supp_6_validation_kvg.Rmd", envir = new.env())
rmarkdown::render("code/p9_sensitivity_summaries/flexible_pui_supp_7_statins_direct_effect.Rmd", envir = new.env())
rmarkdown::render("code/p9_sensitivity_summaries/flexible_pui_supp_8_temporal_validation.Rmd", envir = new.env())
