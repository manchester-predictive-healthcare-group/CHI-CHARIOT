### Clear workspace
rm(list=ls())
# Sys.time()

### Set wd
setwd()
# getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
# sapply(R.func.sources, source)

### Read in cohort
cohort <- readRDS("data/extraction/cohort_baseline/cohort_var.rds")

### Create an age variable
cohort <- dplyr::mutate(cohort, age.group = cut(age, breaks = c(0,seq(25, 80, 5), 86)))
cohort.male <- subset(cohort, gender == 1)
cohort.female <- subset(cohort, gender == 2)

###
### Calculate incidence rates for outcome, by age group
###

### All
outcome.table <- cohort |>
  dplyr::group_by(age.group) |>
  dplyr::summarise(sum_t = sum(cvd_time/365.25), sum_i = sum(cvd_indicator)) |>
  dplyr::mutate(rate = 1000*sum_i/sum_t)
colnames(outcome.table) <- c("Age", "total follow up (years)", "n events", "rate (per 1000 years)")
outcome.table

### Female
outcome.table.female <- cohort.female |>
  dplyr::group_by(age.group) |>
  dplyr::summarise(sum_t = sum(cvd_time/365.25), sum_i = sum(cvd_indicator)) |>
  dplyr::mutate(rate = 1000*sum_i/sum_t)
colnames(outcome.table.female) <- c("Age", "total follow up (years)", "n events", "rate (per 1000 years)")
outcome.table.female

### Male
outcome.table.male <- cohort.male |>
  dplyr::group_by(age.group) |>
  dplyr::summarise(sum_t = sum(cvd_time/365.25), sum_i = sum(cvd_indicator)) |>
  dplyr::mutate(rate = 1000*sum_i/sum_t)
colnames(outcome.table.male) <- c("Age", "total follow up (years)", "n events", "rate (per 1000 years)")
outcome.table.male

### Save these these
saveRDS(outcome.table, "data/extraction/cohort_baseline/incidence_by_age.rds")
saveRDS(outcome.table.female, "data/extraction/cohort_baseline/incidence_by_age_female.rds")
saveRDS(outcome.table.male, "data/extraction/cohort_baseline/incidence_by_age_male.rds")
