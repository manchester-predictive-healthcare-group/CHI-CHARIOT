### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("")
getwd()

### Define filepath to file directory system containing extracted data.
common.data.dir <- file.path("..", "..")

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

library(mice)
library(survival)
library(rms)

### Extract chain seed and gender from command line
args <- commandArgs(trailingOnly = T)
chain.seed <- as.numeric(args[1])
gender <- as.numeric(args[2])
# chain.seed <- as.numeric(1)
# gender <- as.numeric(1)
gender_char <- c("male", "female")[gender]
print(paste("chain seed = ", chain.seed))
print(paste("gender = ", gender_char))

### Extract cohort
if (gender == 1){
  cohort <- readRDS("data/cohort_male_pui.rds")
} else if (gender == 2){
  cohort <- readRDS("data/cohort_female_pui.rds")
}
cohort <- data.table::as.data.table(cohort)
colnames(cohort)

### Get column names for those with missing data, to help visualise
var_miss <- colnames(cohort)[apply(cohort[1:50000, ], 2, function(x) {any(is.na(x))})]
print("VAR MISS")
var_miss

### Check missing pattern
# missing_pattern <- mice::md.pattern(cohort[, ..var_miss])
# str(missing_pattern)

###
### Add Nelson-Aalen estimator of outcome, this will be used in imuptation instead of the time until event
###

### Fit survival model to obtain Nelson-Aalen estimator
survfit_obj <- survfit(Surv(cvd_time, cvd_indicator) ~ 1, data = cohort)

### Write function to extract cumhaz at time point of interest
extract_cumhaz <- function(t){survfit_obj$cumhaz[max(which(survfit_obj$time <= t))]}

### Obtain this for all individuals
cohort$cvd_ev_prim_aj <- sapply(cohort$cvd_time, extract_cumhaz)
str(cohort)

###
### Define the formulas argument
###

### Create first attempt using .make formulas function
### This is to just get the object structure, and will be completely changed
mice_formulas <- make.formulas(cohort)
str(mice_formulas)

### Read in variables that will be interacted with the spline for age
inter_age_rcs <- readRDS(paste("data/var_imp_inter_age_rcs", gender, ".rds", sep = ""))

### Create vector with all terms for the formula
### Note at this point, we add in cvd_indicator and cvd_ev_prim_aj, and all the continuous variables
full_formula_vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))",
                      "age*rms::rcs(bmi, 4)",
                      "age*rms::rcs(sbp, 4)",
                      "age*rms::rcs(sbp_unexposed, 4)",
                      "age*rms::rcs(nonhdl, 4)",
                      "age*rms::rcs(IMD, c(1,10,20))", 
                      paste("age*", inter_age_rcs, sep = ""), "cvd_indicator", "cvd_ev_prim_aj")

### Remove term containing variable being imputed, also remove nonhdl, as this isn't used to impute any other variables
### This is because its effectively derived from the other variables
ethnicity_formula_vec <- full_formula_vec[!grepl("ethnicity", full_formula_vec)]
bmi_formula_vec <- full_formula_vec[!grepl("bmi", full_formula_vec)]
smoking_formula_vec <- full_formula_vec[!grepl("smoking", full_formula_vec)]
IMD_formula_vec <- full_formula_vec[!grepl("IMD", full_formula_vec)]
sbp_formula_vec <- full_formula_vec[!grepl("sbp", full_formula_vec)]
sbp_unexposed_formula_vec <- full_formula_vec[!grepl("sbp_unexposed", full_formula_vec)]
nonhdl_formula_vec <- full_formula_vec[!grepl("nonhdl", full_formula_vec)]

### Assign formulas
print("ethnicity")
mice_formulas[["ethnicity"]] <- as.formula(paste("ethnicity ~ ", paste(ethnicity_formula_vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("smoking")
mice_formulas[["smoking"]] <- as.formula(paste("smoking ~ ", paste(smoking_formula_vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("bmi")
mice_formulas[["bmi"]] <- as.formula(paste("bmi ~ ", paste(bmi_formula_vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("sbp")
mice_formulas[["sbp"]] <- as.formula(paste("sbp ~ ", paste(sbp_formula_vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("sbp_unexposed")
mice_formulas[["sbp_unexposed"]] <- as.formula(paste("sbp_unexposed ~ ", paste(sbp_unexposed_formula_vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("nonhdl")
mice_formulas[["nonhdl"]] <- as.formula(paste("nonhdl ~ ", paste(nonhdl_formula_vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("IMD")
mice_formulas[["IMD"]] <- as.formula(paste("IMD ~ ",paste(IMD_formula_vec, sep = "", collapse = "+"), sep = "", collapse = ""))

### Print the final formulas
mice_formulas

###
### Dry-run
###
print(paste("DRY RUN 1", Sys.time()))
dry <- mice(cohort, formulas = mice_formulas, maxit = 0)
str(dry)
names(dry)

### Choose methods for imputation
### Use predictive mean matching for all variables(including discrete)
### (see Austin and Stef van Buuren 2023 paper for some support for this approach)
meth <- dry$method
meth["smoking"] <- "pmm"
meth["ethnicity"] <- "pmm"
meth

### Real run
print(paste("REAL RUN 1", Sys.time()))
time.in <- Sys.time()
set.seed(chain.seed)
imp <- mice(cohort, maxit = 20, m = 1, formulas = mice_formulas, method = meth)
time.out <- Sys.time()
time.out - time.in

saveRDS(imp, paste("data/mice_mids_prototype3_", gender, "_", chain.seed, ".rds", sep = ""))