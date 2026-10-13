############################################################
### Benchmark implementations of Harrell's C-statistic
############################################################

### Check that required packages are installed
required_packages <- c(
  "survival",
  "intsurv",
  "Hmisc",
  "yardstick",
  "testthat",
  "bench"
)

packages_missing <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(packages_missing) > 0L) {
  stop(
    paste(
      "Install the following packages before running:",
      paste(packages_missing, collapse = ", ")
    )
  )
}

############################################################
### Simulate survival data
############################################################

set.seed(20260917)

n <- 100000L

### Simulate two predictors
x1 <- stats::rnorm(n)
x2 <- stats::rbinom(
  n = n,
  size = 1L,
  prob = 0.4
)

### Define the true linear predictor
risk_score <- 0.8 * x1 + 0.5 * x2

### Simulate event times from a proportional hazards model
baseline_hazard <- 0.05

event_time <- stats::rexp(
  n = n,
  rate = baseline_hazard * exp(risk_score)
)

### Simulate independent censoring times
censoring_time <- stats::rexp(
  n = n,
  rate = 0.03
)

### Create observed follow-up time and event indicator
time <- pmin(
  event_time,
  censoring_time
)

event <- as.integer(
  event_time <= censoring_time
)

### Create survival outcome once
surv_object <- survival::Surv(
  time = time,
  event = event
)

### Check the simulated data
data.frame(
  n = n,
  number_events = sum(event),
  event_proportion = mean(event),
  median_follow_up = stats::median(time)
)

############################################################
### Define one wrapper for each implementation
############################################################

### Current implementation used in your script
calculate_c_intsurv <- function() {

  result <- intsurv::cIndex(
    time = time,
    event = event,
    risk_score = risk_score
  )

  unname(result["index"])
}

### survival package implementation
calculate_c_survival <- function() {

  result <- survival::concordance(
    surv_object ~ risk_score,
    reverse = TRUE,
    timewt = "n",
    influence = 0
  )

  unname(result$concordance)
}

### Hmisc implementation
###
### rcorr.cens() expects larger predictions to indicate longer
### survival. Our risk_score has the opposite direction, so its
### sign is reversed.
calculate_c_hmisc <- function() {

  result <- Hmisc::rcorr.cens(
    x = -risk_score,
    S = surv_object,
    outx = FALSE
  )

  unname(result["C Index"])
}

### yardstick implementation
###
### Larger estimates are documented as indicating greater risk,
### so risk_score is supplied without changing its sign.
calculate_c_yardstick <- function() {

  unname(
    yardstick::concordance_survival_vec(
      truth = surv_object,
      estimate = risk_score
    )
  )
}

############################################################
### Calculate each estimate once
############################################################

c_estimates <- c(
  intsurv = calculate_c_intsurv(),
  survival = calculate_c_survival(),
  Hmisc = calculate_c_hmisc(),
  yardstick = calculate_c_yardstick()
)

print(c_estimates)

############################################################
### Test numerical agreement
############################################################

### Show differences relative to intsurv
differences_from_intsurv <-
  c_estimates - c_estimates["intsurv"]

print(differences_from_intsurv)

### Require agreement to numerical tolerance
testthat::expect_equal(
  unname(c_estimates["survival"]),
  unname(c_estimates["intsurv"]),
  tolerance = 1e-10
)

testthat::expect_equal(
  unname(c_estimates["Hmisc"]),
  unname(c_estimates["intsurv"]),
  tolerance = 1e-10
)

testthat::expect_equal(
  unname(c_estimates["yardstick"]),
  unname(c_estimates["intsurv"]),
  tolerance = 1e-10
)

message(
  "All implementations agree within the specified tolerance."
)

############################################################
### Benchmark runtime and memory allocation
############################################################

benchmark_results <- bench::mark(
  intsurv = calculate_c_intsurv(),
  survival = calculate_c_survival(),
  Hmisc = calculate_c_hmisc(),
  yardstick = calculate_c_yardstick(),
  iterations = 5L,
  check = FALSE,
  memory = FALSE,
  filter_gc = TRUE
)

print(
  benchmark_results[
    ,
    c(
      "expression",
      "min",
      "median",
      "itr/sec",
      "n_gc"
    )
  ]
)