### Info for copilot

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load survival package
library(survival)

### define gender and burnout
gender <- 1
burnout <- 180

### Read in imputed dataset
### Load imp.comb
#imp.comb <- readRDS(paste("data/mice_mids_prototype3_", gender, "_", 1, ".rds", sep = ""))

### Read in raw cohort with missingness
if (gender == 1){
  cohort <- readRDS("data/cohort_male_pui.rds")
} else if (gender == 2){
  cohort <- readRDS("data/cohort_female_pui.rds")
}

### Read in interval censored outcome data, which indicates treatment use at baseline
cohort_split_times <- readRDS(paste("data/cohort_split_times_burnout", burnout, ".rds", sep = "")) |>
  dplyr::filter(tstart == 0) |>
  dplyr::select(patid, med_status_statins, med_status_ah)

### Merge
df <- merge(cohort, cohort_split_times,
            by.x = "patid",
            by.y = "patid", 
            all.x = TRUE)

#--------------------------------------------
# Identify variables with any missingness
#--------------------------------------------
vars_with_na <- names(df)[
  base::colSums(base::is.na(df)) > 0
]

# Remove treatment variables themselves if present
vars_with_na <- base::setdiff(
  vars_with_na,
  c("med_status_statins", "med_status_ah")
)

#--------------------------------------------
# Function to summarise missingness by group
#--------------------------------------------
#--------------------------------------------
# Function to summarise missingness by group
# (no pivot_longer -> no type conflict)
#--------------------------------------------
summarise_missing_by_group <- function(data, treatment_var, vars) {
  
  data |>
    dplyr::group_by(!!rlang::sym(treatment_var)) |>
    
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(vars),
        ~ 100 * base::mean(base::is.na(.x)),
        .names = "pct_missing_{.col}"
      ),
      n_total = dplyr::n(),
      .groups = "drop"
    ) |>
    
    # Long format
    tidyr::pivot_longer(
      cols = dplyr::starts_with("pct_missing_"),
      names_to = "variable",
      values_to = "pct_missing"
    ) |>
    
    # Clean variable names and relabel treatment info
    dplyr::mutate(
      variable = base::sub("^pct_missing_", "", variable),
      
      # Map treatment variable names
      treatment = dplyr::case_when(
        treatment_var == "med_status_statins" ~ "statins",
        treatment_var == "med_status_ah" ~ "antihypertensives"
      ),
      
      # Map treatment group values
      receiving = dplyr::case_when(
        !!rlang::sym(treatment_var) == 1 ~ "yes",
        !!rlang::sym(treatment_var) == 0 ~ "no",
        TRUE ~ NA_character_
      )
    ) |>
    
    # Rename columns as requested
    dplyr::rename(
      N = n_total
    ) |>
    
    # Drop original numeric treatment column
    dplyr::select(
      treatment,
      receiving,
      variable,
      N,
      pct_missing
    )
}

#--------------------------------------------
# Run separately for each treatment variable
#--------------------------------------------
missing_statins <- summarise_missing_by_group(
  data = df,
  treatment_var = "med_status_statins",
  vars = vars_with_na
)

missing_ah <- summarise_missing_by_group(
  data = df,
  treatment_var = "med_status_ah",
  vars = vars_with_na
)

#--------------------------------------------
# Optional: combine results
#--------------------------------------------
missing_summary <- dplyr::bind_rows(
  missing_statins,
  missing_ah
)

#--------------------------------------------
# Optional: arrange for readability
#--------------------------------------------
missing_summary <- missing_summary |>
  dplyr::arrange(
    treatment,
    variable,
    receiving
  ) |> 
  as.data.frame() |>
  dplyr::filter(treatment == antihypertensives)

# Print results
missing_summary
