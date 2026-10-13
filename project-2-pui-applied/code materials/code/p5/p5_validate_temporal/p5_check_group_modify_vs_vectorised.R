###
### Check equivalence of the original group_modify() adjustment and the new
### vectorised adjustment using the first 500 individuals.
###

### Clear workspace
rm(list = ls())
print(Sys.time())

### Set working directory
setwd()
print(getwd())

### Define filepath to extracted data and shared functions
common.data.dir <- file.path("..", "..")

### Source functions
R.func.sources <- list.files(
  file.path(common.data.dir, "Aurum_Jun2021_extract/R"),
  full.names = TRUE
)
invisible(lapply(R.func.sources, source))

### Extract scenario from command line
args <- commandArgs(trailingOnly = TRUE)
t_fup <- round(365.25 * as.numeric(args[1]))
burnout <- 180

if (length(args) < 1L || is.na(t_fup)) {
  stop("Supply the follow-up time in years as the first command-line argument.")
}

### Read in split survival times
input_file <- paste0(
  "data/p5/cohort_split_times_smoking_statins_antihypertensives_burnout",
  burnout,
  ".rds"
)
cohort_split_times <- readRDS(input_file)

### Apply the preprocessing used in the full program
cohort_split_times <- cohort_split_times |>
  dplyr::mutate(
    tstart = tstart - t_fup,
    cvd_time = cvd_time - t_fup
  ) |>
  dplyr::filter(cvd_time > 0) |>
  dplyr::mutate(tstart = pmax(tstart, 0)) |>
  dplyr::arrange(patid, tstart)

### Retain all rows belonging to the first 500 individuals after preprocessing
first_500_patids <- unique(cohort_split_times$patid)[seq_len(
  min(500L, dplyr::n_distinct(cohort_split_times$patid))
)]
cohort_500 <- dplyr::filter(cohort_split_times, patid %in% first_500_patids)

print(paste(
  "Number of individuals included in check =",
  dplyr::n_distinct(cohort_500$patid)
))

### Original patient-level function
my_mutate <- function(df, id) {
  if (df[1, "med_status_adj_statins"] == 1) {
    df <- dplyr::mutate(
      df,
      med_status_adj_statins = med_status_adj_statins - 1
    )
  } else if (df[1, "med_status_adj_statins"] == -1) {
    df <- dplyr::mutate(
      df,
      med_status_adj_statins = med_status_adj_statins + 1
    )
  }

  if (df[1, "med_status_adj_ah"] == 1) {
    df <- dplyr::mutate(df, med_status_adj_ah = med_status_adj_ah - 1)
  } else if (df[1, "med_status_adj_ah"] == -1) {
    df <- dplyr::mutate(df, med_status_adj_ah = med_status_adj_ah + 1)
  }

  df
}

### Run original process
print(paste("Start original group_modify process", Sys.time()))
result_old <- cohort_500 |>
  dplyr::group_by(patid) |>
  dplyr::group_modify(.f = my_mutate) |>
  dplyr::ungroup() |>
  dplyr::arrange(patid, tstart) |>
  as.data.frame()
print(paste("Finish original group_modify process", Sys.time()))

### Run vectorised process
print(paste("Start vectorised process", Sys.time()))
result_new <- cohort_500 |>
  dplyr::group_by(patid) |>
  dplyr::mutate(
    med_status_adj_statins = med_status_adj_statins - dplyr::first(med_status_adj_statins),
    med_status_adj_ah = med_status_adj_ah - dplyr::first(med_status_adj_ah)
  ) |>
  dplyr::ungroup() |>
  dplyr::arrange(patid, tstart) |>
  as.data.frame()
print(paste("Finish vectorised process", Sys.time()))

### Check the complete outputs and each adjusted variable
comparison_all <- all.equal(
  result_old,
  result_new,
  check.attributes = TRUE
)
comparison_statins <- identical(
  result_old$med_status_adj_statins,
  result_new$med_status_adj_statins
)
comparison_ah <- identical(
  result_old$med_status_adj_ah,
  result_new$med_status_adj_ah
)

print(comparison_all)
print(paste("Statin adjustment identical:", comparison_statins))
print(paste("Antihypertensive adjustment identical:", comparison_ah))

if (!isTRUE(comparison_all)) {
  stop("The original and vectorised outputs are not equivalent.")
}

print("SUCCESS: original and vectorised processes are equivalent for the first 500 individuals.")
