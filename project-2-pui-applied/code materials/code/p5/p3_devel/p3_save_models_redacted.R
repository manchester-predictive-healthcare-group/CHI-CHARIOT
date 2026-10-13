###
### Reduce and redact saved p5 Cox models for prediction with newdata.
### These models can be used in rshiny app, and called more efficient on worked example, stability plots, etc
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Function to shrink the output
shrink_fit <- function(fit, knot_values) {
  stopifnot(inherits(fit, "coxph"), !is.null(fit$terms))

  ### If the saved terms still contain knot values, check against the files.
  ### Non-parallel fits may instead point to .GlobalEnv, where they are absent.
  old_env <- attr(fit$terms, ".Environment")
  if (is.environment(old_env)) {
    for (nm in names(knot_values)) {
      if (exists(nm, envir = old_env, inherits = FALSE) &&
          !identical(get(nm, envir = old_env, inherits = FALSE), knot_values[[nm]])) {
        stop("Saved knot values differ from the knot file: ", nm)
      }
    }
  }

  ### Keep only the external knot vectors in the formula environment.
  prediction_env <- list2env(knot_values, parent = asNamespace("stats"))
  attr(fit$terms, ".Environment") <- prediction_env
  if (!is.null(fit$formula)) environment(fit$formula) <- prediction_env

  ### Remove development-row data; retain coefficients, terms, means,
  ### contrasts and other small fit components needed by predict.coxph.
  fit[c("model", "y", "linear.predictors", "residuals",
        "x", "weights", "na.action")] <- NULL
  fit
}

for (gender in c(1, 2)) {
  gender_char <- c("male", "female")[gender]
  message("Processing ", gender_char, " models")

  ### Read the same sex-specific knots used by either fitting script.
  knot_values <- list(
    knot_loc_sbp = readRDS(paste0("data/p5/knots_loc_sbp_", gender)),
    knot_loc_nonhdl = readRDS(paste0("data/p5/knots_loc_nonhdl_", gender)),
    knot_loc_bmi = readRDS(paste0("data/p5/knots_loc_bmi_", gender))
  )

  small_fits <- vector("list", 10)
  small_bhaz <- vector("list", 10)
  offset_names <- c("statins", "ah", "smoking1", "smoking2")
  model_columns <- c(
    "offset(offset_statins_timevar_lnHR)",
    "offset(offset_ah_timevar_lnHR)",
    "offset(offset_smoking_timevar_dummy1_lnHR)",
    "offset(offset_smoking_timevar_dummy2_lnHR)"
  )
  offset_means <- vector("list", 10)

  for (imp_num in seq_len(10)) {
    fit_path <- paste0("data/p5/model_cox_", gender, "_imp", imp_num, ".rds")
    bhaz_path <- paste0("data/p5/model_cox_bhaz_", gender, "_imp", imp_num, ".rds")
    fit <- readRDS(fit_path)
    bhaz <- readRDS(bhaz_path)

    ### Check full model contains appropriate offset columns in development data.
    if (is.null(fit$model) || !all(model_columns %in% names(fit$model))) {
      stop("Missing model frame or offset columns in ", fit_path)
    }
    
    ### Extract the means
    offset_means[[imp_num]] <- setNames(
      vapply(model_columns, function(column) {mean(fit$model[[column]])}, numeric(1)), 
      offset_names)
    
    ### Extract cumulative baseline hazards
    stopifnot(is.data.frame(bhaz), all(c("time", "hazard") %in% names(bhaz)))
    small_bhaz[[imp_num]] <- lapply(seq_len(10), function(year) {
      eligible <- which(bhaz$time <= year * 365.25)
      if (length(eligible)) bhaz$hazard[max(eligible)] else 0
    })

    ### Shrink/redact the model fit
    small <- shrink_fit(fit, knot_values)
    stopifnot(!any(c("model", "y", "linear.predictors", "residuals",
                     "x", "weights", "na.action") %in% names(small)))

    ### Assign redacted fit to list
    small_fits[[imp_num]] <- small
    
    ### Output message
    message(gender_char, " imputation ", imp_num, ": original ",
            round(file.info(fit_path)$size / 1024^2, 1), " MiB; reduced in memory ",
            round(as.numeric(object.size(small)) / 1024^2, 1), " MiB")
    rm(fit, bhaz, small)
    gc(verbose = FALSE)
  }

  ### Save
  saveRDS(small_fits, paste0("data/p5/small_fits_model_", gender_char, ".rds"))
  saveRDS(offset_means, paste0("data/p5/small_means_model_", gender_char, ".rds"))
  saveRDS(small_bhaz, paste0("data/p5/small_bhaz_model_", gender_char, ".rds"))

  message("Saved ", gender_char, " outputs; fits file size: ",
          round(file.info(paste0("data/p5/small_fits_model_", gender_char, ".rds"))$size /
                  1024^2, 1), " MiB")
  rm(small_fits, small_bhaz, offset_means, knot_values)
  gc(verbose = FALSE)
}

message("FINISHED ", Sys.time())
