###
### This program will run the multiple imputation for cohort_baseline
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
getwd()

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
gender_char <- c("male", "female")[gender]
print(paste("chain seed = ", chain.seed))
print(paste("gender = ", gender_char))
# chain.seed <- 1
# gender <- 1

### Extract cohort
cohort <- readRDS("data/extraction/cohort_baseline/cohort_var.rds")
cohort <- data.table::as.data.table(cohort)

### Reduce to gender of interest
cohort <- subset(cohort, gender == get("gender", pos = 1))
# cohort <- cohort[1:100000, ]

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

### Retain only formulas for missing variables
mice.formulas <- mice.formulas[var.miss]
str(mice.formulas)

### Read in variables that will be interacted with the spline for age
inter.age.rcs <- readRDS(paste("data/extraction/cohort_baseline/inter_age_rcs_imp", gender, ".rds", sep = ""))

### Create vector with all terms for the formula
### Note at this point, we add in cvd_indicator and cvd_ev_prim_aj
full.formula.vec <- c("rms::rcs(age, c(40, 57.5, 75))*rms::rcs(bmi, 4)",
                      "rms::rcs(age, c(40, 57.5, 75))*rms::rcs(sbp, 4)",
                      "rms::rcs(age, c(40, 57.5, 75))*rms::rcs(sbp_var, 4)",
                      "rms::rcs(age, c(40, 57.5, 75))*rms::rcs(cholhdl_ratio, 4)", 
                      "rms::rcs(age, c(40, 57.5, 75))*rms::rcs(IMD, c(1,10,20))", 
                      paste("rms::rcs(age, c(40, 57.5, 75))*", var_imp_inter_age_rcs, sep = ""), "cvd_indicator", "cvd_ev_prim_aj")

### Remove term containing variable being imputed
ethnicity.formula.vec <- full.formula.vec[!grepl("ethnicity", full.formula.vec)]
bmi.formula.vec <- full.formula.vec[!grepl("bmi", full.formula.vec)]
sbp.formula.vec <- full.formula.vec[!grepl("sbp,", full.formula.vec)]
sbp_var.formula.vec <- full.formula.vec[!grepl("sbp_var", full.formula.vec)]
cholhdl_ratio.formula.vec <- full.formula.vec[!grepl("cholhdl_ratio", full.formula.vec)]
smoking.formula.vec <- full.formula.vec[!grepl("smoking", full.formula.vec)]
IMD.formula.vec <- full.formula.vec[!grepl("IMD", full.formula.vec)]

### Assign formulas
print("ethnicity")
mice.formulas[["ethnicity"]] <- as.formula(paste("ethnicity ~ ",paste(ethnicity.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("smoking")
mice.formulas[["smoking"]] <- as.formula(paste("smoking ~ ",paste(smoking.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("bmi")
mice.formulas[["bmi"]] <- as.formula(paste("bmi ~ ",paste(bmi.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("sbp")
mice.formulas[["sbp"]] <- as.formula(paste("sbp ~ ",paste(sbp.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("sbp_var")
mice.formulas[["sbp_var"]] <- as.formula(paste("sbp_var ~ ",paste(sbp_var.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("cholhdl_ratio")
mice.formulas[["cholhdl_ratio"]] <- as.formula(paste("cholhdl_ratio ~ ",paste(cholhdl_ratio.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

print("IMD")
mice.formulas[["IMD"]] <- as.formula(paste("IMD ~ ",paste(IMD.formula.vec, sep = "", collapse = "+"), sep = "", collapse = ""))

# ### Keeping this here but commented out in case I want code again in future. An attempt to create the formulas manually.
# ### Add interaction terms between age and relevant predictor varaibles
# mice.formulas <- lapply(1:length(mice.formulas), function(x){
#   
#   ## Get vector of variable names to be interacted with age
#   # Not cvd_indicator or anything containing age
#   inter.age <- var.all[!(var.all %in% c("cvd_indicator", "age", paste("age.rcs", 1:ncol(age.rcs), sep = "")))]
#   # Not the variable we are imputing
#   inter.age <- inter.age[!(inter.age == names[x])]
#   # Not the splines of the varaible we are imputing
#   inter.age <- inter.age[!grepl(paste(names[x], ".rcs", sep = ""), inter.age)]
#   
#   ## Get vector of variable names to be interacted with age splines (just the other continuous rcs terms for spline x spline interaction)
#   # Not cvd_indicator or anything containing age
#   inter.age.rcs <- var.all[!(var.all %in% c("cvd_indicator", "age", paste("age.rcs", 1:ncol(age.rcs), sep = "")))]
#   # Nothing without an rcs in it
#   inter.age.rcs <- inter.age.rcs[grepl("rcs", inter.age.rcs)]
#   # Not the splines of the varaible we are imputing
#   inter.age.rcs <- inter.age.rcs[!grepl(paste(names[x], ".rcs", sep = ""), inter.age.rcs)]
#   
#   ## Create interaction terms
#   inter1 <- paste(paste(inter.age, "age", sep = "*"), collapse = "+")
#   inter2 <- lapply(1:3, function(y) {paste(paste(paste(inter.age.rcs, "age.rcs", sep = "*"), y, sep = ""), collapse = "+")})
#   
#   ## Combine (choose whether or not to include spline x spline interactions)
#   inter.comb <- paste(inter1, paste(inter2, collapse = "+"), sep = "+")
#   
#   ## Add to formulas
#   #out <- as.formula(paste(mice.formulas[x], inter1, sep = "+"))
#   out <- as.formula(paste(mice.formulas[x], inter.comb, sep = "+"))
#   
#   ## Return
#   return(out)})
# 
# ### Assign names
# names(mice.formulas) <- names
# str(mice.formulas)
# print(mice.formulas)

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
meth

### Dry-run 2
print(paste("DRY RUN 2", Sys.time()))
dry <- mice(cohort, maxit = 0, formulas = mice.formulas, meth = meth)

### Dry-run 3, with some actual iterations
print(paste("REAL RUN 1", Sys.time()))
time.in <- Sys.time()
set.seed(chain.seed)
imp <- mice(cohort, maxit = 20, m = 1, formulas = mice.formulas, meth = meth)
time.out <- Sys.time()
time.out - time.in

saveRDS(imp, paste("data/extraction/cohort_baseline/mice_mids_", gender, "_", chain.seed, ".rds", sep = ""))
