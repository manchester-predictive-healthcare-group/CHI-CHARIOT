### Assess the plausibility of the MAR assumption and observed-data support
### Parallelise over missing variables, separately for females and males

### Clear workspace
rm(list = ls())
Sys.time()

### Set working directory
setwd()
getwd()

### Load packages
library(dplyr)
library(gtsummary)
library(rms)
library(survival)
library(doParallel)
library(foreach)

### Variables containing missing data
missing_variables <- c("ethnicity", "smoking", "bmi", "sbp", "nonhdl", "IMD")

### Variables shown in descriptive tables
summary_variables <- c(
  "age", "ethnicity", "smoking", "bmi", "sbp", "nonhdl", "IMD",
  "cvd_indicator", "cvd_ev_prim_aj"
)

### Fit a model for missingness and return performance and predictions
assess_missingness <- function(data, variable, model_matrix) {
  print(paste("assess missingness", variable, Sys.time()))

  ### Create missingness indicator: 1 = missing, 0 = observed
  missing_indicator <- as.integer(is.na(data[[variable]]))

  ### Fit logistic regression using the pre-calculated design matrix
  model <- stats::glm.fit(
    x = model_matrix,
    y = missing_indicator,
    family = stats::binomial()
  )

  ### Extract predicted probabilities before discarding the full fitted model
  predicted_probability <- model$fitted.values
  rm(model)

  ### Calculate C-statistic
  c_statistic <- survival::concordance(
    missing_indicator ~ predicted_probability
  )$concordance

  ### Calculate model and null log-likelihoods directly
  loglik_model <- sum(stats::dbinom(
    missing_indicator,
    size = 1,
    prob = predicted_probability,
    log = TRUE
  ))
  missing_probability <- mean(missing_indicator)
  loglik_null <- sum(stats::dbinom(
    missing_indicator,
    size = 1,
    prob = missing_probability,
    log = TRUE
  ))

  ### Calculate Nagelkerke R-squared
  n <- length(missing_indicator)
  r_squared <- (1 - exp((2 / n) * (loglik_null - loglik_model))) /
    (1 - exp((2 / n) * loglik_null))

  ### Return only the results needed for reporting
  list(
    performance = data.frame(
      variable = variable,
      missing_proportion = mean(missing_indicator),
      c_statistic = c_statistic,
      nagelkerke_r_squared = r_squared
    ),
    predictions = data.frame(
      missing_indicator = missing_indicator,
      predicted_missingness = predicted_probability
    )
  )
}

### Create a compact cohort summary stratified by missingness
summarise_by_missingness <- function(data, variable) {
  ### Create labelled missingness indicator
  table_data <- data |>
    dplyr::select(dplyr::all_of(summary_variables)) |>
    dplyr::mutate(
      missing_indicator = factor(
        is.na(data[[variable]]),
        levels = c(FALSE, TRUE),
        labels = c("Observed", "Missing")
      )
    )

  ### Calculate the same summary and overall column, then retain only displayed cells
  gtsummary::tbl_summary(
    data = table_data,
    by = missing_indicator,
    include = dplyr::all_of(summary_variables),
    missing = "ifany",
    statistic = list(
      gtsummary::all_continuous() ~ "{mean} ({sd})",
      gtsummary::all_categorical() ~ "{n} ({p}%)"
    )
  ) |>
    gtsummary::add_overall() |>
    tibble::as_tibble(col_labels = TRUE) |>
    as.data.frame()
}

### Summarise support across groups of a continuous predictor
summarise_missingness_by_group <- function(data, variable, predictor, groups = 20) {
  data |>
    dplyr::mutate(
      predictor_group = dplyr::ntile(.data[[predictor]], groups),
      missing_indicator = is.na(.data[[variable]])
    ) |>
    dplyr::group_by(predictor_group) |>
    dplyr::summarise(
      predictor_min = min(.data[[predictor]], na.rm = TRUE),
      predictor_max = max(.data[[predictor]], na.rm = TRUE),
      n = dplyr::n(),
      proportion_missing = mean(missing_indicator),
      n_observed = sum(!missing_indicator),
      .groups = "drop"
    ) |>
    dplyr::mutate(variable = variable, predictor = predictor, .before = 1)
}

### Run analyses separately for females and males
for (gender in c(2, 1)) {
  gender_char <- c("male", "female")[gender]
  print(paste("gender =", gender_char))

  ### Load and restrict the cohort
  cohort <- readRDS("data/p5/cohort_prototype3.rds")
  cohort <- cohort[cohort$gender == gender, ]
  # cohort <- cohort[sample(1:nrow(cohort), 100000), ] # xx

  ### Match variables removed before imputation
  cohort <- dplyr::select(
    cohort,
    -c(cholhdl_ratio, sbp_var, bmi_alltime, sbp_alltime,
       sbp_var_alltime, cholhdl_ratio_alltime)
  )

  ### Add the Nelson-Aalen cumulative hazard used in imputation
  survfit_object <- survival::survfit(
    survival::Surv(cvd_time, cvd_indicator) ~ 1,
    data = cohort
  )

  ### Match each follow-up time to the latest cumulative hazard
  cumhaz_index <- findInterval(cohort$cvd_time, survfit_object$time)
  cohort$cvd_ev_prim_aj <- 0
  has_cumhaz <- cumhaz_index > 0
  cohort$cvd_ev_prim_aj[has_cumhaz] <-
    survfit_object$cumhaz[cumhaz_index[has_cumhaz]]

  ### Recreate predictor terms used in the imputation models
  inter_age_rcs <- readRDS(
    paste0("data/p5/var_imp_inter_age_rcs", gender, ".rds")
  )
  predictor_terms <- c(
    "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(bmi, 4)",
    "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(sbp, 4)",
    "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(nonhdl, 4)",
    "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(IMD, c(1, 10, 20))",
    paste0("rms::rcs(age, c(25, 40, 57.5, 75))*", inter_age_rcs),
    "cvd_indicator",
    "cvd_ev_prim_aj"
  )

  ### Remove terms containing any incomplete variable
  incomplete_pattern <- paste(missing_variables, collapse = "|")
  predictor_terms <- predictor_terms[!grepl(incomplete_pattern, predictor_terms)]

  ### Create the spline and interaction design matrix once
  predictor_formula <- stats::as.formula(
    paste("~", paste(predictor_terms, collapse = " + "))
  )
  model_matrix <- stats::model.matrix(predictor_formula, data = cohort)

  ### Set up parallel backend with one worker per variable
  n_cores <- min(length(missing_variables), parallel::detectCores())
  print(paste("Parallel workers =", n_cores))
  cl <- parallel::makeCluster(
    n_cores,
    outfile = paste0("data/p5/sens_imputation_missingness_workers_", gender, ".log")
  )
  doParallel::registerDoParallel(cl)

  ### Fit missingness models in parallel over variables
  model_results <- tryCatch(
    foreach::foreach(
      variable = missing_variables,
      .packages = "survival",
      .export = c("assess_missingness", "cohort", "model_matrix")
    ) %dopar% {
      assess_missingness(data = cohort, variable = variable, model_matrix = model_matrix)
    },
    finally = {
      parallel::stopCluster(cl)
    }
  )
  names(model_results) <- missing_variables
  print(paste("Missingness models finished for", gender_char, Sys.time()))

  ### Combine model-performance results
  performance <- dplyr::bind_rows(
    lapply(model_results, function(x) x$performance)
  ) |>
    dplyr::mutate(gender = gender_char, .before = 1)

  ### Summarise predicted missingness and high-probability support
  ### (distribution of probability of being missing, for observed and missing data)
  prediction_summary_overlap <- dplyr::bind_rows(
    lapply(names(model_results), function(variable) {
      predictions <- model_results[[variable]]$predictions
      predictions |>
        dplyr::group_by(missing_indicator) |>
        dplyr::summarise(
          q01 = stats::quantile(predicted_missingness, 0.01),
          q10 = stats::quantile(predicted_missingness, 0.10),
          median = stats::median(predicted_missingness),
          q90 = stats::quantile(predicted_missingness, 0.90),
          q99 = stats::quantile(predicted_missingness, 0.99),
          .groups = "drop"
        ) |>
        dplyr::mutate(variable = variable, .before = 1)
    })
  ) |>
    dplyr::mutate(gender = gender_char, .before = 1)

  ### Save model outputs before creating the descriptive tables
  saveRDS(performance, paste0("data/p5/sens_imputation_missingness_model_performance_", gender, ".rds"))
  saveRDS(prediction_summary_overlap, paste0("data/p5/sens_imputation_prediction_summary_overlap_", gender, ".rds"))
  saveRDS(model_results, paste0("data/p5/sens_imputation_missingness_models_data_", gender, ".rds"))
  rm(model_results)
  print(paste("Missingness model outputs saved for", gender_char, Sys.time()))

  ### Create compact descriptive tables stratified by missingness
  missingness_tables <- stats::setNames(
    lapply(
      missing_variables,
      summarise_by_missingness,
      data = cohort
    ),
    missing_variables
  )
  print(paste("Missingness tables finished for", gender_char, Sys.time()))

  # ### Assess observed-data support across age twentieths
  # age_support <- dplyr::bind_rows(
  #   lapply(
  #     missing_variables,
  #     summarise_missingness_by_group,
  #     data = cohort,
  #     predictor = "age"
  #   )
  # ) |>
  #   dplyr::mutate(gender = gender_char, .before = 1)

  ### Save compact sex-specific descriptive tables under the original filename
  saveRDS(missingness_tables, paste0("data/p5/sens_imputation_missingness_summary_tables_", gender, ".rds"))
  # write.csv(
  #   age_support,
  #   paste0("data/p5/sens_imputation_missingness_by_age_", gender, ".csv"),
  #   row.names = FALSE
  # )
}

### Mark completion of both gender-specific analyses
print(paste("FINISHED", Sys.time()))
