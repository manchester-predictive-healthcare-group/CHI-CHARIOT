### Test to check baseline smoking against outcome split-survival files.
### Individuals that are smoker or ex-smoker at baseline, should never move to non-smoker during follow-up.

### Preliminaries
rm(list = ls())
library(data.table)

### Locate the input files
project_dir <- "./"
data_dir <- file.path(project_dir, "data/p5")
results <- list()

### Record a test without stopping execution
record_test <- function(file, test, failures = NA_real_, detail = "", status = NULL) {
  if (is.null(status)) {
    status <- if (is.na(failures)) "NOT TESTED" else if (failures == 0) "PASS" else "FAIL"
  }
  results[[length(results) + 1L]] <<- data.table(
    file = file, test = test, status = status,
    failures = as.numeric(failures), detail = detail
  )
  cat(sprintf("[%s] %s: %s; failures=%s %s\n",
              status, file, test, as.character(failures), detail))
}

### Read and validate baseline data
baseline <- tryCatch({
  x <- readRDS(file.path(data_dir, "cohort_prototype3.rds"))
  setDT(x)
  if (!all(c("patid", "smoking") %in% names(x))) {
    stop("Required columns: patid, smoking")
  }
  x <- x[, .(patid, baseline_smoking = as.character(smoking))]
  x[baseline_smoking %in% c("Light", "Moderate", "Heavy"),
    baseline_smoking := "Current"]
  record_test("baseline", "Missing patient IDs", sum(is.na(x$patid)))
  record_test("baseline", "Duplicate patient IDs", sum(duplicated(x$patid)))
  record_test("baseline", "Unexpected smoking labels",
              sum(!is.na(x$baseline_smoking) &
                    !x$baseline_smoking %in% c("Non-smoker", "Ex-smoker", "Current")))
  cat("Baseline smoking counts:\n")
  print(x[, .N, by = baseline_smoking])
  x
}, error = function(e) {
  record_test("baseline", "Read/prepare baseline", detail = conditionMessage(e), status = "ERROR")
  NULL
})

### Check one split file at a time to limit peak memory
check_split <- function(filename, status_column) {
  label <- filename
  cat("\nChecking", label, "at", as.character(Sys.time()), "\n")
  x <- tryCatch({
    z <- readRDS(file.path(data_dir, filename))
    setDT(z)
    needed <- c("patid", "tstart", "cvd_time", status_column)
    if (!all(needed %in% names(z))) {
      stop("Missing required columns: ", paste(setdiff(needed, names(z)), collapse = ", "))
    }
    z <- z[, needed, with = FALSE]
    setnames(z, status_column, "smoking_status")
    z
  }, error = function(e) {
    record_test(label, "Read/prepare split file", detail = conditionMessage(e), status = "ERROR")
    NULL
  })
  if (is.null(x)) return(invisible(NULL))

  ### Isolate test errors so later checks still run
  run_test <- function(name, fun) {
    tryCatch(fun(), error = function(e) {
      record_test(label, name, detail = conditionMessage(e), status = "ERROR")
    })
  }

  run_test("Status coding", function() {
    if (!is.numeric(x$smoking_status)) stop("Smoking status must be numeric: 0=never, 1=ex, 2=current")
    record_test(label, "Unexpected non-missing status codes",
                sum(!is.na(x$smoking_status) & !x$smoking_status %in% c(0, 1, 2)))
    cat("Split smoking counts (intervals, not patients):\n")
    print(x[, .N, by = smoking_status])
  })

  run_test("Interval structure", function() {
    record_test(label, "Missing patient IDs", sum(is.na(x$patid)))
    record_test(label, "Invalid interval bounds",
                sum(!is.finite(x$tstart) | !is.finite(x$cvd_time) |
                      x$tstart < 0 | x$cvd_time <= x$tstart))
    record_test(label, "Duplicate patient/start-time pairs",
                sum(duplicated(x, by = c("patid", "tstart"))))
  })

  ### Compare unique IDs, without joining baseline onto every interval
  run_test("Patient coverage", function() {
    if (is.null(baseline)) stop("Baseline unavailable")
    ids <- unique(x$patid)
    record_test(label, "Baseline patients absent from split file",
                sum(!baseline$patid %in% ids))
    record_test(label, "Split patients absent from baseline",
                sum(!ids %in% baseline$patid))
  })

  ### Main test: observed baseline ever-smokers must never have explicit status zero
  run_test("Baseline versus follow-up smoking", function() {
    if (is.null(baseline)) stop("Baseline unavailable")
    if (!is.numeric(x$smoking_status)) stop("Unexpected smoking status type")
    zero_ids <- x[!is.na(smoking_status) & smoking_status == 0, patid]
    zero_counts <- data.table(patid = zero_ids)[, .(never_intervals = .N), by = patid]
    ever <- baseline[baseline_smoking %in% c("Current", "Ex-smoker")]
    bad <- zero_counts[ever, on = "patid", nomatch = 0L]
    record_test(label, "Baseline Current/Ex-smoker with any never-smoker interval",
                nrow(bad), "Failure count is patients, not intervals")
    cat("Affected patients and intervals by observed baseline status:\n")
    print(bad[, .(patients = .N, intervals = sum(never_intervals)), by = baseline_smoking])
    rm(zero_ids, zero_counts, ever, bad)
  })

  ### Missing split status is not automatically an error: initial intervals may be unresolved
  run_test("Missing-status diagnostics", function() {
    cat("Missing split status by initial versus later interval:\n")
    print(x[, .(intervals = .N, missing_status = sum(is.na(smoking_status))),
            by = .(initial_interval = tstart == 0)])
    record_test(label, "Missing split smoking status (informational)",
                sum(is.na(x$smoking_status)),
                "Initial NA statuses may be filled from baseline/imputation downstream", status = "INFO")
  })

  ### Check longitudinal reversions without group_modify or a full row-level join
  run_test("Within-file smoking history", function() {
    if (!is.numeric(x$smoking_status)) stop("Unexpected smoking status type")
    setorder(x, patid, tstart, cvd_time)
    failures <- x[, {
      previous_ever <- shift(cummax(as.integer(smoking_status %in% c(1, 2))), fill = 0L)
      n_bad <- sum(smoking_status == 0 & previous_ever == 1L, na.rm = TRUE)
      .(bad_intervals = n_bad)
    }, by = patid]
    record_test(label, "Patients reverting to never-smoker after Ex/Current",
                sum(failures$bad_intervals > 0))
    record_test(label, "Intervals reverting to never-smoker after Ex/Current",
                sum(failures$bad_intervals))
  })
  invisible(NULL)
}

### Check both files even if one is missing or fails
files <- c("cohort_split_times_smoking.rds",
           "cohort_split_times_smoking_statins_antihypertensives_burnout180.rds")
columns <- c("smoking", "med_status_adj_smoking")
for (i in seq_along(files)) {
  tryCatch(check_split(files[i], columns[i]), error = function(e) {
    record_test(files[i], "Unexpected file-level error", detail = conditionMessage(e), status = "ERROR")
  })
  gc()
}

### Print every result, without stopping on failures
summary_results <- rbindlist(results, fill = TRUE)
cat("\nFINAL TEST SUMMARY\n")
print(as.data.frame(summary_results), row.names = FALSE)
cat("\nCompleted all available checks at", as.character(Sys.time()), "\n")
