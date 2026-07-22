### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
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
cohort <- readRDS("data/p4/cohort_prototype3.rds")
cohort <- data.table::as.data.table(cohort)

### Reduce to gender of interest
cohort <- subset(cohort, gender == get("gender", pos = 1))
# cohort <- cohort[1:50000, ]

### Remove continuous variables we don't want
cohort <- dplyr::select(cohort, -c(cholhdl_ratio, sbp_var, bmi_alltime, sbp_alltime, sbp_var_alltime, cholhdl_ratio_alltime))
colnames(cohort)

### Get column names for those with missing data, to help visualise
var.miss <- colnames(cohort)[apply(cohort[1:50000, ], 2, function(x) {any(is.na(x))})]
print("VAR MISS")
var.miss

### Check missing pattern
# missing.pattern <- mice::md.pattern(cohort[, ..var.miss])
# str(missing.pattern)

###
### Add Nelson-Aalen estimator of outcome, this will be used in imuptation instead of the time until event
###

### Fit survival model to obtain Nelson-Aalen estimator
survfit.obj <- survfit(Surv(cvd_time, cvd_indicator) ~ 1, data = cohort)

### Write function to extract cumhaz at time point of interest
extract_cumhaz <- function(t){survfit.obj$cumhaz[max(which(survfit.obj$time <= t))]}

### Obtain this for all individuals
cohort$cvd_ev_prim_aj <- sapply(cohort$cvd_time, extract_cumhaz)
str(cohort)

###
### Define the formulas argument
###

### Create first attempt using .make formulas function
### This is to just get the object structure, and will be completely changed
mice.formulas <- make.formulas(cohort)
str(mice.formulas)

### Read in variables that will be interacted with the spline for age
inter.age.rcs <- readRDS(paste("data/p4/var_imp_inter_age_rcs", gender, ".rds", sep = ""))

### Create vector with all terms for the formula
### Note at this point, we add in cvd_indicator and cvd_ev_prim_aj, and all the continuous variables
full.formula.vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(bmi, 4)",
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(sbp, 4)",
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(nonhdl, 4)",
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(IMD, c(1,10,20))", 
                      paste("rms::rcs(age, c(25, 40, 57.5, 75))*", inter.age.rcs, sep = ""), "cvd_indicator", "cvd_ev_prim_aj")

### Remove term containing variable being imputed, also remove nonhdl, as this isn't used to impute any other variables
### This is because its effectively derived from the other variables
ethnicity.formula.vec <- full.formula.vec[!grepl("ethnicity", full.formula.vec)]
bmi.formula.vec <- full.formula.vec[!grepl("bmi", full.formula.vec)]
smoking.formula.vec <- full.formula.vec[!grepl("smoking", full.formula.vec)]
IMD.formula.vec <- full.formula.vec[!grepl("IMD", full.formula.vec)]
sbp.formula.vec <- full.formula.vec[!grepl("sbp", full.formula.vec)]
nonhdl.formula.vec <- full.formula.vec[!grepl("nonhdl", full.formula.vec)]

### Assign formulas
print("ethnicity")
mice.formulas[["ethnicity"]] <- as.formula(paste("ethnicity ~ ", paste(ethnicity.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("smoking")
mice.formulas[["smoking"]] <- as.formula(paste("smoking ~ ", paste(smoking.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("bmi")
mice.formulas[["bmi"]] <- as.formula(paste("bmi ~ ", paste(bmi.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("sbp")
mice.formulas[["sbp"]] <- as.formula(paste("sbp ~ ", paste(sbp.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("nonhdl")
mice.formulas[["nonhdl"]] <- as.formula(paste("nonhdl ~ ", paste(nonhdl.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("IMD")
mice.formulas[["IMD"]] <- as.formula(paste("IMD ~ ",paste(IMD.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

### Print the final formulas
mice.formulas

###
### Dry-run
###
print(paste("DRY RUN 1", Sys.time()))
dry <- mice(cohort, formulas = mice.formulas, maxit = 0)
str(dry)
names(dry)

### Choose methods for imputation
### Use predictive mean matching for all variables(including discrete)
### (see Austin and Stef van Buuren 2023 paper for some support for this approach)
meth <- dry$method
meth["smoking"] <- "pmm"
meth["ethnicity"] <- "pmm"
meth["hdl"] <- ""
meth["ldl"] <- ""
meth["cholesterol"] <- ""
meth

### Dry-run 2
print(paste("DRY RUN 2", Sys.time()))
dry <- mice(cohort, maxit = 0, formulas = mice.formulas, method = meth)

### Dry-run 3, with some actual iterations
print(paste("REAL RUN 1", Sys.time()))
time.in <- Sys.time()
set.seed(chain.seed)
imp <- mice(cohort, maxit = 40, m = 1, formulas = mice.formulas, method = meth)
time.out <- Sys.time()
time.out - time.in

saveRDS(imp, paste("data/p4/mice_mids_prototype3_", gender, "_", chain.seed, ".rds", sep = ""))