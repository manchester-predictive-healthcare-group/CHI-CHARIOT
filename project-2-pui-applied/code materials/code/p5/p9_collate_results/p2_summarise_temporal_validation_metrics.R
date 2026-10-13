###
### Summarise temporal validation performance for the baseline and interventional layers
###
rm(list = ls())
Sys.time()
setwd()
getwd()

results_dir <- "data/p5"
followup_years <- 1:5
mfrs <- c("sbp", "nonhdl", "bmi", "smoking")

### Read one scalar metric without silently replacing missing/invalid results
read_metric <- function(object, metric, file) {
  value <- object[[metric]]
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop("Expected one finite numeric ", metric, " in ", file)
  }
  as.numeric(value)
}

### Each gender/time/factor combination has one calibration and one discrimination object
summarise_one <- function(gender, year, t_eval_int, mfr = NA_character_) {
  interventional <- !is.na(mfr)
  stem <- if (interventional) "temporalv_interventional_" else "temporalv_"
  suffix <- paste0("_", gender, "_t", year,
                   if (interventional) paste0("_", mfr) else "", 
                   "_teval", t_eval_int,  ".rds")
  calib_file <- file.path(results_dir, paste0(stem, "calib_ph", suffix))
  discrim_file <- file.path(results_dir, paste0(stem, "discrim", suffix))
  calib <- readRDS(calib_file)
  discrim <- readRDS(discrim_file)
  x <- data.frame(
    Gender = c("Male", "Female")[gender],
    `Follow-up index (years)` = year,
    `Prediction interval` = t_eval_int,
    N = read_metric(discrim, "n", discrim_file),
    ICI = read_metric(calib, "ICI", calib_file),
    E50 = read_metric(calib, "E50", calib_file),
    E90 = read_metric(calib, "E90", calib_file),
    Discrimination = read_metric(discrim, "concordance", discrim_file),
    check.names = FALSE
  )
  if (interventional) {
    x <- data.frame(Gender = x$Gender, `Risk factor` = mfr,
                    x[, -1, drop = FALSE], check.names = FALSE)
  }
  x
}

### Keep the initial-risk and interventional results in separate data frames
baseline_results <- do.call(rbind, lapply(c(2, 1), function(gender) {
  do.call(rbind, lapply(c(5,10), function(t_eval_int){
    do.call(rbind, lapply(followup_years, function(year) summarise_one(gender, year, t_eval_int)))
  }))
}))
interventional_results <- do.call(rbind, lapply(c(2, 1), function(gender) {
  do.call(rbind, lapply(c(5,10), function(t_eval_int){
    do.call(rbind, lapply(mfrs, function(mfr) {
      do.call(rbind, lapply(followup_years, function(year) summarise_one(gender, year, t_eval_int, mfr)))
    }))
  }))
}))
rownames(baseline_results) <- NULL
rownames(interventional_results) <- NULL

### Create female and male display tables, rounded to three decimal places
factor_labels <- c(sbp = "SBP", nonhdl = "Non-HDL cholesterol",
                   bmi = "BMI", smoking = "Smoking")
make_tables <- function(results, interventional = FALSE) {
  tables <- list()
  for (gender in c("Female", "Male")) {
    x <- results[results$Gender == gender, , drop = FALSE]
    x$Gender <- NULL
    if (interventional) x[["Risk factor"]] <- unname(factor_labels[x[["Risk factor"]]])
    for (metric in c("ICI", "E50", "E90", "Discrimination")) {
      x[[metric]] <- round(x[[metric]], 3)
    }
    rownames(x) <- NULL
    tables[[gender]] <- x
  }
  tables
}
baseline_tables_3dp <- make_tables(baseline_results)
interventional_tables_3dp <- make_tables(interventional_results, interventional = TRUE)

### Save separate RDS files, each with full-precision results and female/male tables
saveRDS(list(results = baseline_results, tables_3dp = baseline_tables_3dp),
        file.path(results_dir, "temporal_validation_performance_metrics.rds"))
saveRDS(list(results = interventional_results, tables_3dp = interventional_tables_3dp),
        file.path(results_dir, "temporal_validation_interventional_performance_metrics.rds"))

