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
cohort.dates <- readRDS("data/extraction/cohort_baseline/cohort_var.rds")

###
### Calculate secular trend in incidence rates for outcome
###
cohort.dates <- dplyr::select(cohort.dates, patid, gender, fup_start, cvd_time, cvd_indicator)
# set.seed(101)
# cohort.dates <- cohort.dates[sample(1:nrow(cohort.dates), 5000000, replace = FALSE), ]

### Create a survival variable based on fup_start and cvd_time 
###
cohort.dates <- dplyr::mutate(cohort.dates, 
                              cvd_time_calendar = as.numeric(fup_start + cvd_time),
                              fup_start = as.numeric(fup_start))

### Create cut-offs for data
splits <- as.numeric(c(seq(as.Date("01/01/2005", format = "%d/%m/%Y"), as.Date("01/01/2019", format = "%d/%m/%Y"), by = "year"), 
                       as.Date("01/03/2020", format = "%d/%m/%Y")))

### Split data
cohort.dates <- survival::survSplit(cohort.dates, cut = splits, end = "cvd_time_calendar", event = "cvd_indicator", start = "fup_start")

### Add an indicator variable for the period
cohort.dates$year <- sapply(cohort.dates$fup_start, function(x){max(which(splits <= x))})

### Now create incidence table, by year
calc_incidence_year <- function(x, gender.var){
  
  ### Primary + Secondary + ONS
  outcome.table <- cohort.dates |>
    dplyr::filter(year == x) |>
    dplyr::filter(gender == gender.var) |>
    dplyr::summarise(sum_t = sum((cvd_time_calendar - fup_start)/365.25), sum_i = sum(cvd_indicator)) |>
    dplyr::mutate(rate = 1000*sum_i/sum_t,
                  year = as.Date(splits[x]))
  colnames(outcome.table) <- c("total follow up (years)", "n events", "rate (per 1000 years)", "year")
  
  return(outcome.table)
}

### Calculate incidence rate tables
incidence_by_year_male <- do.call("rbind", lapply(1:max(cohort.dates$year), calc_incidence_year, gender.var = 1))
incidence_by_year_female <- do.call("rbind", lapply(1:max(cohort.dates$year), calc_incidence_year, gender.var = 2))

### Print these
print("Incidence by year female")
incidence_by_year_female

print("Incidence by year male")
incidence_by_year_male

saveRDS(incidence_by_year_male, "data/extraction/cohort_baseline/incidence_by_year_male.rds")
saveRDS(incidence_by_year_female, "data/extraction/cohort_baseline/incidence_by_year_female.rds")
