###
### Summarise variation in test performance across bootstrap iterations
###

### Clear workspace
rm(list = ls())
Sys.time()

### Set working directory
setwd()
getwd()

### Define inputs and output locations
results_dir <- "data/p5"
figures_dir <- "figures/p5"
iterations <- 1:500
metrics <- c("ICI", "E50", "E90", "Discrimination")

### Read test-set performance metrics from each bootstrap iteration
### These scalar objects are saved directly by p2_optimism_loop.R
performance <- do.call(
  rbind,
  lapply(c(2, 1), function(gender) {
    do.call(
      rbind,
      lapply(iterations, function(iter) {
        data.frame(
          Gender = c("Male", "Female")[gender],
          Iteration = iter,
          ICI = readRDS(file.path(results_dir,paste0( "optimism_boot_calib_ph_ICI_", gender, "_iter", iter, "_dforig.rds"))),
          E50 = readRDS(file.path(results_dir,paste0( "optimism_boot_calib_ph_E50_", gender, "_iter", iter, "_dforig.rds"))),
          E90 = readRDS(file.path(results_dir,paste0( "optimism_boot_calib_ph_E90_", gender, "_iter", iter, "_dforig.rds"))),
          Discrimination = readRDS(file.path(results_dir,paste0( "optimism_boot_discrim_",gender,"_iter", iter,"_dforig.rds"))),
          check.names = FALSE
        )
      })
    )
  })
)
rownames(performance) <- NULL

### Reshape to one row per sex and metric
performance_long <- tidyr::pivot_longer(
  performance,
  cols = dplyr::all_of(metrics),
  names_to = "Metric",
  values_to = "Test performance"
)

### Calculate percentile summaries across bootstrap iterations
performance_percentiles <- performance_long |>
  dplyr::group_by(Gender, Metric) |>
  dplyr::summarise(
    `2.5th percentile` = as.numeric(stats::quantile(
      `Test performance`,
      probs = 0.025,
      names = FALSE
    )),
    `50th percentile` = as.numeric(stats::quantile(
      `Test performance`,
      probs = 0.50,
      names = FALSE
    )),
    `97.5th percentile` = as.numeric(stats::quantile(
      `Test performance`,
      probs = 0.975,
      names = FALSE
    )),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Gender = factor(Gender, levels = c("Female", "Male")),
    Metric = factor(Metric, levels = metrics)
  ) |>
  dplyr::arrange(Gender, Metric) |>
  dplyr::mutate(
    Gender = as.character(Gender),
    Metric = as.character(Metric)
  )

### Create a three-decimal version for presentation
performance_percentiles_3dp <- performance_percentiles
percentile_columns <- c(
  "2.5th percentile",
  "50th percentile",
  "97.5th percentile"
)
performance_percentiles_3dp[percentile_columns] <- lapply(
  performance_percentiles_3dp[percentile_columns],
  function(x) sprintf("%.3f", x)
)

### Save full-precision input results and percentile summaries
saveRDS(
  list(
    test_performance = performance,
    percentiles = performance_percentiles,
    percentiles_3dp = performance_percentiles_3dp
  ),
  file.path(results_dir, "variance_bootstrap_test_metrics.rds")
)

### Print the final table
print(performance_percentiles_3dp)
print(paste("FINISHED", Sys.time()))
