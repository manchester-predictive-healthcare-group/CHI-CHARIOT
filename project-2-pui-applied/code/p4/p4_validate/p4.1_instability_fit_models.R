###
### Fit models for instability plots
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define filepath to file directory system containing extracted data, and functions for extracting.
common.data.dir <- file.path("..", "..")

### Define burnout
burnout <- as.numeric(180)

### Extract relevant scenario from command line
args <- as.numeric(commandArgs(trailingOnly = T))
seed <- as.numeric(args[1])
gender <- as.numeric(args[2])
gender_char <- c("male", "female")[gender]
print(paste("seed = ", seed))
print(paste("gender = ", gender_char))
print(paste("burnout = ", burnout))

### Set seed
set.seed(seed)

###
### Define formula
###

### Read in variables that will be interacted with the spline for age
inter.age.rcs <- readRDS(paste("data/p4/var_model_inter_age_rcs", gender, ".rds", sep = ""))

### Create vector with all terms for the formula
full.formula.vec <- c("rms::rcs(age, c(40, 57.5, 75))*rms::rcs(IMD, c(1,10,20))",
                      "rms::rcs(age, c(40, 57.5, 75))*rms::rcs(sbp, knot_loc_sbp)",
                      "rms::rcs(age, c(40, 57.5, 75))*rms::rcs(nonhdl, knot_loc_nonhdl)",
                      "rms::rcs(age, c(40, 57.5, 75))*rms::rcs(bmi, knot_loc_bmi)",
                      "rms::rcs(age, c(40, 57.5, 75))*smoking",
                      paste("rms::rcs(age, c(40, 57.5, 75))*", inter.age.rcs, sep = ""),
                      "offset(offset_statins_timevar_lnHR)",
                      "offset(offset_ah_timevar_lnHR)",
                      "offset(offset_smoking_timevar_dummy1_lnHR)",
                      "offset(offset_smoking_timevar_dummy2_lnHR)")

### Create formula
model.formula.full <- as.formula(
  paste("survival::Surv(tstart, cvd_time, cvd_indicator) ~ ", paste(full.formula.vec, sep = "", collapse = "+"), sep = "", collapse = "")
)
print(model.formula.full)

### Define HR for offsets
lnHR_statins <- readRDS("data/p4/offsets_lnHR_statins.rds")
lnHR_ah <- readRDS("data/p4/offsets_lnHR_ah.rds")
lnHR_smoking_dummy1_total <- readRDS("data/p4/offsets_lnHR_smoking_dummy1_total.rds")
lnHR_smoking_dummy2_total <- readRDS("data/p4/offsets_lnHR_smoking_dummy2_total.rds")
lnHR_smoking_dummy1_direct <- readRDS("data/p4/offsets_lnHR_smoking_dummy1_direct.rds")
lnHR_smoking_dummy2_direct <- readRDS("data/p4/offsets_lnHR_smoking_dummy2_direct.rds")
lnHR_sbp <- readRDS("data/p4/offsets_lnHR_sbp.rds")
lnHR_bmi <- readRDS("data/p4/offsets_lnHR_bmi.rds")
lnHR_nonhdl <- readRDS("data/p4/offsets_lnHR_nonhdl.rds")

### Read in the outcome data
cohort.split.times <- readRDS(paste("data/p4/cohort_split_times_smoking_statins_antihypertensives_burnout", burnout, ".rds", sep = ""))

### Get the knot locations for continuous variables
knot_loc_sbp <- readRDS(paste("data/p4/knots_loc_sbp_", gender, sep = ""))
knot_loc_sbp
knot_loc_nonhdl <- readRDS(paste("data/p4/knots_loc_nonhdl_", gender, sep = ""))
knot_loc_nonhdl
knot_loc_bmi <- readRDS(paste("data/p4/knots_loc_bmi_", gender, sep = ""))
knot_loc_bmi

### Read in the imputed development datasets
imp.list <- readRDS(paste("data/p4/dfs_devel_", gender_char, ".rds", sep = ""))

### Define imputed dataset
imp <- imp.list[[1]]

### Create bootstrap sample
imp <- imp[sample(1:nrow(imp), nrow(imp), replace = TRUE), ]

### Create a new patid
imp$patid2 <- 1:nrow(imp)

### Merge imputed dataset with cohort.split.times
data.devel <- dplyr::left_join(dplyr::select(imp, -c("cvd_time", "cvd_indicator",
                                                     "cvd_time_prim", "cvd_time_hes", "cvd_time_death",
                                                     "cvd_indicator_prim", "cvd_indicator_hes", "cvd_indicator_death")),
                               cohort.split.times,
                               by = dplyr::join_by("patid"),
                               relationship = "many-to-many")

print("Is there correct number of individuals in data.devel")
print(length(unique(data.devel$patid2))==nrow(imp))

###
### Create dummy variables for smoking (at baseline, and during follow-up)
###

### The first dummy variable = 0 if never smoker, = 1 if current or past smoker.
### NB: This dummy effect (HR1) is the effect of current smoker vs never smoker
###
### The second dummy variable = 0 if never or current smoker, = 1 if past smoker
### NB: This dummy effect (HR2) is the effect of smoking cessation (stopping smoking, past smoker vs current smoker)
###
### These dummies are only appropriate when applied in tandem.
### A never smoker will have: LP*1*1
### A current smoker will have: LP*HR1*1
### A past smoker will have: LP*HR1*HR2, i.e., the effect of initiating smoking, then the effect of stopping
data.devel <- create_smoking_dummies(data.devel)

###
### Define adjusted versions of SBP, nonhdl and BMI to be relative to some baseline, to apply the offsets
###

### Sort by patid and tstart, just for when looking at the dataset
data.devel <- dplyr::arrange(data.devel, patid, tstart)

### Create appropriate offsets for timevarying variables (interventions applied during follow-up)
data.devel$offset_statins_timevar_lnHR <- lnHR_statins*data.devel$med_status_adj_statins
data.devel$offset_ah_timevar_lnHR <- lnHR_ah*data.devel$med_status_adj_ah
data.devel$offset_smoking_timevar_dummy1_lnHR <- lnHR_smoking_dummy1_total*data.devel$med_status_adj_smoking_dummy1
data.devel$offset_smoking_timevar_dummy2_lnHR <- lnHR_smoking_dummy2_total*data.devel$med_status_adj_smoking_dummy2

### Fit a cox model and get bhaz
print(paste("fit model", 1, Sys.time()))
fit <- survival::coxph(model.formula.full, data = data.devel, model = TRUE)
print(paste("bhaz fit model", 1, Sys.time()))
bhaz <- survival::basehaz(fit, centered = TRUE)

### Calculate means
offset_means <- c("statins" = mean(fit$model[,"offset(offset_statins_timevar_lnHR)"]),
                  "ah" = mean(fit$model[,"offset(offset_ah_timevar_lnHR)"]),
                  "smoking1" = mean(fit$model[,"offset(offset_smoking_timevar_dummy1_lnHR)"]),
                  "smoking2" = mean(fit$model[,"offset(offset_smoking_timevar_dummy2_lnHR)"]))

### Get rid of residuals and linear predictors to reduce file size
fit[names(fit) == "model"] <- "DELETED"
fit[names(fit) == "y"] <- "DELETED"
fit[names(fit) == "linear.predictors"] <- "DELETED"
fit[names(fit) == "residuals"] <- "DELETED"

### Save fit
saveRDS(fit, paste("data/p4/prototype3_stability_cox_", gender, "_samp", seed, ".rds", sep = ""))
saveRDS(bhaz, paste("data/p4/prototype3_stability_cox_bhaz_", gender, "_samp", seed, ".rds", sep = ""))
saveRDS(offset_means, paste("data/p4/prototype3_stability_offset_means_", gender, "_samp", seed, ".rds", sep = ""))

print("FINISHED")





