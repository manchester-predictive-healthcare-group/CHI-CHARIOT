###
### Assess imputation convergence and distributions
###

### Clear workspace
rm(list = ls())
Sys.time()

### Set working directory
setwd()
getwd()

### Source functions
R.func.sources <- list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

library(mice)

### Create a table of category percentages for the observed data and imputations
create_percentage_table <- function(imp, variable) {

  ### Extract the observed data and each completed imputation
  completed_data <- lapply(0:imp$m, function(imputation) {
    mice::complete(imp, action = imputation)[[variable]]
  })

  ### Retain a consistent set and order of categories across all rows
  if (is.factor(imp$data[[variable]])) {
    categories <- levels(imp$data[[variable]])
  } else {
    categories <- sort(unique(unlist(completed_data)))
    categories <- categories[!is.na(categories)]
  }

  ### Calculate the percentage in each category, excluding missing values
  percentages <- lapply(completed_data, function(values) {
    counts <- table(factor(values, levels = categories), useNA = "no")
    100 * counts / sum(counts)
  })

  ### Return observed data first, followed by each imputation
  percentage_table <- as.data.frame(do.call(rbind, percentages))
  percentage_table <- data.frame(
    imputation = c("Observed", paste0("Imputation ", seq_len(imp$m))),
    percentage_table,
    check.names = FALSE
  )

  rownames(percentage_table) <- NULL
  percentage_table
}

###
### Create diagnostics for imputed datasets
###

for (gender in c(2, 1)) {

  gender_char <- c("male", "female")[gender]
  print(paste("gender =", gender_char))

  ### Load combined imputation object
  imp.comb <- readRDS(
    paste0("data/p5/imp_comb_prototype3", gender, ".rds")
  )

  ### Get prop missing data
  variables <- c("ethnicity", "smoking", "bmi", "sbp", "nonhdl", "IMD")
  proportion_missing <- colMeans(is.na(imp.comb$data[variables]))
  proportion_missing <- data.frame(variable = variables, proportion_missing = proportion_missing)
  saveRDS(proportion_missing, paste0("data/p5/imputation_diagnostics_missingness_proportions_", gender, ".rds"))
  
  ### Create convergence plots
  print("convergence plots")

  ragg::agg_png(
    paste0("figures/p5/convergence_plot_prototype3", gender, ".png"),
    width = 12,
    height = 15,
    units = "in",
    res = 120
  )

  print(
    plot(
      imp.comb,
      y = c("ethnicity", "smoking", "bmi", "sbp", "nonhdl", "IMD"),
      layout = c(2, 6)
    )
  )

  dev.off()

  # ### Randomly select individuals for density plots for code testing
  # set.seed(123)
  # 
  # plot_rows <- sample.int(
  #   nrow(imp.comb$data),
  #   min(100000, nrow(imp.comb$data))
  # )
  # 
  # imp.comb <- dplyr::filter(
  #   imp.comb,
  #   seq_len(nrow(imp.comb$data)) %in% plot_rows
  # )

  ### Create density plots for the continuous variables in one call
  print("density plots")

  ragg::agg_png(
    paste0("figures/p5/density_plot_prototype3", gender, ".png"),
    width = 12,
    height = 5,
    units = "in",
    res = 120
  )

  print(
    mice::densityplot(
      imp.comb,
      ~ sbp + nonhdl + bmi,
      plot.points = FALSE,
      layout = c(3, 1)
    )
  )

  dev.off()

  ### Create category percentage tables using the full imputation object
  imd_table <- create_percentage_table(imp.comb, "IMD")
  ethnicity_table <- create_percentage_table(imp.comb, "ethnicity")
  smoking_table <- create_percentage_table(imp.comb, "smoking")

  ### Save one table per categorical variable
  saveRDS(imd_table, paste0("data/p5/imputation_diagnostics_IMD_percentage_table", gender, ".rds"))

  saveRDS(ethnicity_table, paste0("data/p5/imputation_diagnostics_ethnicity_percentage_table", gender, ".rds"))

  saveRDS(smoking_table, paste0("data/p5/imputation_diagnostics_smoking_percentage_table", gender, ".rds"))
}

Sys.time()
