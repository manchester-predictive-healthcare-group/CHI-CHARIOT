###
### Create incidence by year tables
###
rm(list = ls())
Sys.time()
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Read current full cohort
cohort <- readRDS("data/p5/cohort_prototype3.rds")

### Ensure follow-up start is stored as a Date
cohort <- dplyr::mutate(
  cohort,
  fup_start = as.Date(fup_start),
  followup_end = fup_start + cvd_time,
  event_date = dplyr::if_else(
    cvd_indicator == 1,
    followup_end,
    as.Date(NA)
  ),
  Gender = factor(
    gender,
    levels = c(1, 2),
    labels = c("Male", "Female")
  )
)

### Define calendar years covered by follow-up
calendar_years <- seq(
  as.integer(format(min(cohort$fup_start, na.rm = TRUE), "%Y")),
  as.integer(format(max(cohort$followup_end, na.rm = TRUE), "%Y"))
)

### Calculate person-time and events occurring within each calendar year
incidence_by_year <- do.call(
  rbind,
  lapply(calendar_years, function(year){
    
    ### Define year start and end
    year_start <- as.Date(paste0(year, "-01-01"))
    next_year_start <- as.Date(paste0(year + 1, "-01-01"))
    
    ### Start and end of each person's contribution within this year
    time_start <- pmax(
      cohort$fup_start,
      year_start,
      na.rm = TRUE
    )
    
    time_end <- pmin(
      cohort$followup_end,
      next_year_start,
      na.rm = TRUE
    )
    
    ### A person contributes only if under follow-up during this year
    contributes <- (
      !is.na(time_start) &
        !is.na(time_end) &
        time_end > time_start
    )
    
    ### Identify events occurring within this calendar year
    event_in_year <- (
      cohort$cvd_indicator == 1 &
        !is.na(cohort$event_date) &
        cohort$event_date >= year_start &
        cohort$event_date < next_year_start
    )
    
    ### Output data frame
    data.frame(
      Gender = cohort$Gender[contributes],
      Year = year,
      person_time_days = as.numeric(
        time_end[contributes] - time_start[contributes]
      ),
      event = as.integer(event_in_year[contributes])
    )
  })
) |>
  ### Group by gender and year, and summarise
  dplyr::group_by(Gender, Year) |>
  dplyr::summarise(
    `Person-years` = sum(person_time_days) / 365.25,
    `Number of events` = sum(event),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    `Incidence rate per 1,000 person-years` =
      1000 * `Number of events` / `Person-years`
  ) |>
  dplyr::arrange(dplyr::desc(Gender), Year)

### Save
saveRDS(incidence_by_year, "data/p5/outcome_incidence_by_year.rds")
print(paste("FINISHED BY YEAR", Sys.time()))

###
### Also do break down of outcome type (primary care, secondary care, ONS)
###

###
### Calculate incidence rates for outcome
###

### Primary + Secondary + ONS
outcome.table.all <- cohort |>
  dplyr::filter(gender %in% c(1,2)) |>
  dplyr::group_by(Gender) |>
  dplyr::summarise(sum_t = sum(cvd_time/365.25), sum_i = sum(cvd_indicator)) |>
  dplyr::mutate(rate = 1000*sum_i/sum_t) |>
  dplyr::arrange(dplyr::desc(Gender))
colnames(outcome.table.all) <- c("Gender", "total follow up (years)", "n events", "rate (per 1000 person years)")

### Primary
outcome.table.primary <- cohort |>
  dplyr::filter(gender %in% c(1,2)) |>
  dplyr::group_by(Gender) |>
  dplyr::summarise(sum_t = sum(cvd_time_prim/365.25), sum_i = sum(cvd_indicator_prim)) |>
  dplyr::mutate(rate = 1000*sum_i/sum_t) |>
  dplyr::arrange(dplyr::desc(Gender))
colnames(outcome.table.primary) <- c("Gender", "total follow up (years)", "n events", "rate (per 1000 person years)")

### Secondary
outcome.table.secondary <- cohort |>
  dplyr::filter(gender %in% c(1,2)) |>
  dplyr::group_by(Gender) |>
  dplyr::summarise(sum_t = sum(cvd_time_hes/365.25), sum_i = sum(cvd_indicator_hes)) |>
  dplyr::mutate(rate = 1000*sum_i/sum_t) |>
  dplyr::arrange(dplyr::desc(Gender))
colnames(outcome.table.secondary) <- c("Gender", "total follow up (years)", "n events", "rate (per 1000 person years)")

### ONS
outcome.table.ons <- cohort |>
  dplyr::filter(gender %in% c(1,2)) |>
  dplyr::group_by(Gender) |>
  dplyr::summarise(sum_t = sum(cvd_time_death/365.25), sum_i = sum(cvd_indicator_death)) |>
  dplyr::mutate(rate = 1000*sum_i/sum_t) |>
  dplyr::arrange(dplyr::desc(Gender))
colnames(outcome.table.ons) <- c("Gender", "total follow up (years)", "n events", "rate (per 1000 person years)")

saveRDS(outcome.table.all, "data/p5/outcome_incidence_all.rds")
saveRDS(outcome.table.primary, "data/p5/outcome_incidence_primary.rds")
saveRDS(outcome.table.secondary, "data/p5/outcome_incidence_secondary.rds")
saveRDS(outcome.table.ons, "data/p5/outcome_incidence_ons.rds")
print(paste("FINISHED", Sys.time()))
