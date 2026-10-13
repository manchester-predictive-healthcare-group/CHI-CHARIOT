###
### Summarise overall and subgroup optimism-adjusted performance
###
rm(list = ls())
Sys.time()
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Settings
results_dir <- "data/p5"
metrics <- c("ICI", "E50", "E90", "Discrimination")

### Choose iterations
i <- 1
j <- 500

### Stop and report missing requested iterations
check_files <- function(files, description) {
  missing <- which(!file.exists(files))
  if (length(missing) > 0) {
    stop(paste0("Missing ", description, " for iteration(s): ", paste(missing, collapse = ", ")))
  }
}

### Summarise apparent, paired optimism and test performance
summarise_values <- function(apparent, boot, test, gender, analysis, subgroup, N) {
  optimism <- Map(function(x, y) x - y, boot, test)
  data.frame(
    Gender = gender,
    Analysis = analysis,
    Subgroup = subgroup,
    N = rep(N, length(apparent)),
    Metric = names(apparent),
    `Apparent performance (median)` = as.numeric(apparent),
    `Optimism-adjusted performance` = as.numeric(apparent) - vapply(optimism, mean, numeric(1)),
    `Optimism mean` = vapply(optimism, mean, numeric(1)),
    `Optimism 2.5th percentile` = vapply(optimism, quantile, numeric(1), probs = 0.025, names = FALSE),
    `Optimism 50th percentile` = vapply(optimism, median, numeric(1)),
    `Optimism 97.5th percentile` = vapply(optimism, quantile, numeric(1), probs = 0.975, names = FALSE),
    `Test performance 2.5th percentile` = vapply(test, quantile, numeric(1), probs = 0.025, names = FALSE),
    `Test performance 50th percentile` = vapply(test, median, numeric(1)),
    `Test performance 97.5th percentile` = vapply(test, quantile, numeric(1), probs = 0.975, names = FALSE),
    check.names = FALSE
  )
}

### Overall results for one gender
summarise_overall <- function(gender) {
  apparent <- setNames(numeric(4), metrics)
  boot <- setNames(vector("list", 4), metrics)
  test <- setNames(vector("list", 4), metrics)
  for (metric in metrics) {
    if (metric == "Discrimination") {
      apparent[metric] <- median(vapply(1:10, function(i) {
        readRDS(file.path(results_dir, paste0("discrim_", gender, "_", i, ".rds")))$concordance
      }, numeric(1)))
      stem <- "optimism_boot_discrim_"
    } else {
      apparent[metric] <- median(vapply(1:10, function(i) {
        readRDS(file.path(results_dir, paste0("calib_ph_", gender, "_", i, ".rds")))[[metric]]
      }, numeric(1)))
      stem <- paste0("optimism_boot_calib_ph_", metric, "_")
    }
    boot_files <- file.path(results_dir, paste0(stem, gender, "_iter", i:j, "_dfboot.rds"))
    test_files <- file.path(results_dir, paste0(stem, gender, "_iter", i:j, "_dforig.rds"))
    check_files(boot_files, paste("overall", metric, "boot gender", gender))
    check_files(test_files, paste("overall", metric, "test gender", gender))
    boot[[metric]] <- vapply(boot_files, readRDS, numeric(1))
    test[[metric]] <- vapply(test_files, readRDS, numeric(1))
  }
  ### Overall N from the first imputed dataset for this gender
  gender_char <- c("male", "female")[gender]
  N <- nrow(readRDS(file.path(results_dir, paste0("dfs_all_", gender_char, "1.rds"))))
  ### Calculate optimism
  summarise_values(apparent, boot, test, c("Male", "Female")[gender], "Overall", "Overall", N)
}

### Subgroup results for one gender and grouping variable
summarise_subgroups <- function(gender, byvar, file_label) {
  # Read in internal calibration results
  calib_internal <- lapply(1:10, function(i) {
    readRDS(file.path(results_dir, paste0("calib_ph_", gender, "_", i, "_", byvar, ".rds")))
  })
  discrim_internal <- lapply(1:10, function(i) {
    readRDS(file.path(results_dir, paste0("discrim_", gender, "_", i, "_", byvar, ".rds")))
  })
  levels <- as.character(discrim_internal[[1]]$levels)

  ### Read in apparent (dfboot) and test (dforig) performance from boostrap procedure
  ### (Create filepaths, check exist, read in)
  calib_boot_files <- file.path(results_dir, paste0("optimism_boot_calib_ph_by_", file_label, "_", gender, "_iter", i:j, "_dfboot.rds"))
  calib_test_files <- file.path(results_dir, paste0("optimism_boot_calib_ph_by_", file_label, "_", gender, "_iter", i:j, "_dforig.rds"))
  discrim_boot_files <- file.path(results_dir, paste0("optimism_boot_discrim_by_", file_label, "_", gender, "_iter", i:j, "_dfboot.rds"))
  discrim_test_files <- file.path(results_dir, paste0("optimism_boot_discrim_by_", file_label, "_", gender, "_iter", i:j, "_dforig.rds"))
  check_files(calib_boot_files, paste(byvar, "calibration boot gender", gender))
  check_files(calib_test_files, paste(byvar, "calibration test gender", gender))
  check_files(discrim_boot_files, paste(byvar, "discrimination boot gender", gender))
  check_files(discrim_test_files, paste(byvar, "discrimination test gender", gender))

  calib_boot <- lapply(calib_boot_files, readRDS)
  calib_test <- lapply(calib_test_files, readRDS)
  discrim_boot <- lapply(discrim_boot_files, readRDS)
  discrim_test <- lapply(discrim_test_files, readRDS)

  ### Add N for subgroups (for ethnicity we take mean across imputations)
  do.call(rbind, lapply(seq_along(levels), function(j) {
    level <- levels[j]
    ### Use the N saved for each validation subgroup, matched by label
    N_by_imputation <- vapply(discrim_internal, function(x) {
      k <- match(level, as.character(x$levels))
      if (is.na(k)) stop("Subgroup missing from discrimination results: ", byvar, ": ", level)
      as.numeric(x$N[k])
    }, numeric(1))
    if (anyNA(N_by_imputation)) stop("Missing N for ", byvar, ": ", level)

    ### Age and region N should not vary across imputations
    if (byvar != "ethnicity" && length(unique(N_by_imputation)) != 1L) {
      stop("N differs across imputations for ", byvar, ": ", level)
    }

    ### Ethnicity N is the rounded mean across the ten imputations
    if (byvar == "ethnicity") {
      N <- round(mean(N_by_imputation))
    } else {
      N <- N_by_imputation[1]
    }

    ### Bootstrap calibration labels capitalise the first character of each level
    calib_level <- paste0(toupper(substr(level, 1, 1)), substring(level, 2))
    ### Prepare data to calculate optimism
    apparent <- c(
      ICI = median(vapply(calib_internal, function(x) as.numeric(x[[j]]$ICI), numeric(1))),
      E50 = median(vapply(calib_internal, function(x) as.numeric(x[[j]]$E50), numeric(1))),
      E90 = median(vapply(calib_internal, function(x) as.numeric(x[[j]]$E90), numeric(1))),
      Discrimination = median(vapply(discrim_internal, function(x) x$Cstat[match(level, x$levels)], numeric(1)))
    )
    boot <- list(
      ICI = vapply(calib_boot, function(x) as.numeric(x$ICI[[calib_level]]), numeric(1)),
      E50 = vapply(calib_boot, function(x) as.numeric(x$E50[[calib_level]]), numeric(1)),
      E90 = vapply(calib_boot, function(x) as.numeric(x$E90[[calib_level]]), numeric(1)),
      Discrimination = vapply(discrim_boot, function(x) x[level, "Cstat"], numeric(1))
    )
    test <- list(
      ICI = vapply(calib_test, function(x) as.numeric(x$ICI[[calib_level]]), numeric(1)),
      E50 = vapply(calib_test, function(x) as.numeric(x$E50[[calib_level]]), numeric(1)),
      E90 = vapply(calib_test, function(x) as.numeric(x$E90[[calib_level]]), numeric(1)),
      Discrimination = vapply(discrim_test, function(x) x[level, "Cstat"], numeric(1))
    )
    ### Calculate optimism
    summarise_values(apparent, boot, test, c("Male", "Female")[gender],
                     tools::toTitleCase(gsub("_", " ", byvar)), level, N)
  }))
}

### Create all overall and subgroup results
optimism_results <- do.call(rbind, lapply(1:2, function(gender) {
  rbind(
    summarise_overall(gender),
    summarise_subgroups(gender, "ethnicity", "ethnicity"),
    summarise_subgroups(gender, "age_cat", "age"),
    summarise_subgroups(gender, "region", "region")
  )
}))
rownames(optimism_results) <- NULL

### Create a 3dp version
optimism_tables_3dp <- optimism_results
optimism_tables_3dp$N <- round(optimism_tables_3dp$N)
num <- vapply(optimism_tables_3dp, is.numeric, logical(1))
optimism_tables_3dp[num] <- lapply(optimism_tables_3dp[num], round, 3)
optimism_tables <- split(optimism_tables_3dp,
                         list(optimism_tables_3dp$Gender, optimism_tables_3dp$Analysis), drop = TRUE)

### Save full-precision results and three-decimal tables for later use
saveRDS(list(results = optimism_results, tables_3dp = optimism_tables),
        "data/p5/optimism_adjusted_performance_metrics.rds")
print(paste("FINISHED", Sys.time()))