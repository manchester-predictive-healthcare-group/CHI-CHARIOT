# run_all_rmds.R
# Knits all .Rmd files in a specified directory

# Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

# ---- Path to folder ----
rmd_dir <- "code/p9_sensitivity_summaries"

# ---- Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

# ---- Run all Rmds ----
rmd_files <- list.files(rmd_dir, pattern = "\\.Rmd$", full.names = TRUE)

if (length(rmd_files) == 0) {
  message("No .Rmd files found in: ", rmd_dir)
} else {
  message("Found ", length(rmd_files), " .Rmd file(s). Knitting...")
  
  for (rmd in rmd_files) {
    message("\n--- Knitting: ", basename(rmd), " ---")
    tryCatch(
      rmarkdown::render(rmd, envir = new.env()),
      error = function(e) message("ERROR in ", basename(rmd), ": ", e$message)
    )
  }
  
  message("\nDone.")
}
