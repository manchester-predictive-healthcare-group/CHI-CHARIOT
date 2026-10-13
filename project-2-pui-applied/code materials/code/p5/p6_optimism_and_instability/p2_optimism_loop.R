###
### This file will run the entire procedure of:
### - Take bootstrap resample of data and impute
### - Fit model in bootstrapped dataset
### - Generate risk scores in bootstrapped and M original datasets
### - Estimate counterfactual outcomes in bootstrapped and M original datasets
### - Estimate performance of model in boostrapped (bootstrap performance) and M original (test performance) datasets
### - Calculate optimism
###

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

### Define burnout
burnout <- as.numeric(180)

library(mice)
library(survival)
library(rms)

### Define iteration
args <- commandArgs(trailingOnly = T)
gender <- as.numeric(args[1])
iter <- as.numeric(args[2])
# gender <- as.numeric(2)
# iter <- as.numeric(1)
print(paste("iter = ", iter))
gender_char <- c("male", "female")[gender]
print(paste("gender = ", gender_char))

### Set seed
set.seed(iter)

#####################################################
#####################################################
### 1. Take bootstrap resample of data and impute ###
#####################################################
#####################################################

### Extract cohort
cohort <- readRDS("data/p5/cohort_prototype3.rds")
cohort <- data.table::as.data.table(cohort)

### Reduce to gender of interest
cohort <- subset(cohort, gender == get("gender", pos = 1))
# cohort <- cohort[1:50000, ]

### Resample with replacement
cohort <- cohort[sample(nrow(cohort), replace = TRUE), ]
# cohort <- cohort[sample(200000, replace = TRUE), ] xx
cohort <- dplyr::mutate(cohort, patid_boot = 1:nrow(cohort))
saveRDS(cohort, paste("data/p5/optimism_boot_df_", gender, "_", iter, ".rds", sep = ""))

### Remove continuous variables we don't want
cohort <- dplyr::select(cohort, -c(cholhdl_ratio, sbp_var, bmi_alltime, sbp_alltime, sbp_var_alltime, cholhdl_ratio_alltime))
colnames(cohort)

### Get column names for those with missing data, to help visualise
var.miss <- colnames(cohort)[apply(cohort[1:50000, ], 2, function(x) {any(is.na(x))})]
print("VAR MISS")
var.miss

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
inter.age.rcs <- readRDS(paste("data/p5/var_imp_inter_age_rcs", gender, ".rds", sep = ""))

### Create vector with all terms for the formula
### Note at this point, we add in cvd_indicator and cvd_ev_prim_aj, and all the continuous variables
full.formula.vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(bmi, 4)",
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(sbp, 4)",
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(nonhdl, 4)",
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(IMD, c(1,10,20))", 
                      paste("rms::rcs(age, c(25, 40, 57.5, 75))*", inter.age.rcs, sep = ""), "cvd_indicator", "cvd_ev_prim_aj")

### Remove term containing variable being imputed
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

### Real run, with some actual iterations
print(paste("REAL RUN 1", Sys.time()))
time.in <- Sys.time()
imp <- mice(cohort, maxit = 10, m = 1, formulas = mice.formulas, method = meth)
# imp <- mice(cohort, maxit = 3, m = 1, formulas = mice.formulas, method = meth) xx
time.out <- Sys.time()
time.out - time.in

### Save imputation object
saveRDS(imp, paste("data/p5/optimism_boot_mice_mids_", gender, "_", iter, ".rds", sep = ""))

### Extract the imputed data
imp <- mice::complete(imp, action = 1)

### Re-order levels in ethnicity variable
eth.levels <- levels(imp$ethnicity)
imp$ethnicity <- factor(imp$ethnicity, c(eth.levels[9],eth.levels[1:8]))

### Add calendar time
imp <- dplyr::mutate(imp, caltime = as.numeric(fup_start) - as.numeric(as.Date("01/01/2005", format = "%d/%m/%Y")))

### Save imputed dataset
saveRDS(imp, paste("data/p5/optimism_boot_df_imp_", gender, "_", iter, ".rds", sep = ""))
print(paste("Step 1 complete", Sys.time()))

#########################################
#########################################
### 2. Fit model in bootstrapped data ###
#########################################
#########################################

###
### Define formula
###

### Read in variables that will be interacted with the spline for age
inter.age.rcs <- readRDS(paste("data/p5/var_model_inter_age_rcs", gender, ".rds", sep = ""))

### Create vector with all terms for the formula
full.formula.vec <- c("rms::rcs(caltime, 4)", 
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(IMD, c(1,10,20))", 
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(sbp, knot_loc_sbp)",
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(nonhdl, knot_loc_nonhdl)", 
                      "rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(bmi, knot_loc_bmi)", 
                      "rms::rcs(age, c(25, 40, 57.5, 75))*smoking", 
                      paste("rms::rcs(age, c(25, 40, 57.5, 75))*", inter.age.rcs, sep = ""), 
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
lnHR_statins <- readRDS("data/p5/offsets_total_lnHR_statins.rds")
lnHR_ah <- readRDS("data/p5/offsets_total_lnHR_ah.rds")
lnHR_smoking_dummy1_total <- readRDS("data/p5/offsets_total_lnHR_smok_init.rds")
lnHR_smoking_dummy2_total <- readRDS("data/p5/offsets_total_lnHR_smok_cess.rds")

### Read in the outcome data
cohort_split_times <- readRDS(paste("data/p5/cohort_split_times_smoking_statins_antihypertensives_burnout", burnout, ".rds", sep = ""))

### Estimate the knot locations for continuous variables
knots <- Hmisc::rcspline.eval(imp[,"sbp"], nk = 3)
knot_loc_sbp <- attributes(knots)$knots
print(paste("sbp knots", knot_loc_sbp))

knots <- Hmisc::rcspline.eval(imp[,"bmi"], nk = 3)
knot_loc_bmi <- attributes(knots)$knots
print(paste("bmi knots", knot_loc_bmi))

knots <- Hmisc::rcspline.eval(imp[,"nonhdl"], nk = 3)
knot_loc_nonhdl <- attributes(knots)$knots
print(paste("nonhdl knots", knot_loc_nonhdl))

### Merge bootstrapped dataset with cohort_split_times
data_devel <- dplyr::left_join(dplyr::select(imp, -c("cvd_time", "cvd_indicator",
                                                     "cvd_time_prim", "cvd_time_hes", "cvd_time_death",
                                                     "cvd_indicator_prim", "cvd_indicator_hes", "cvd_indicator_death")),
                               cohort_split_times,
                               by = dplyr::join_by("patid"),
                               relationship = "many-to-many") |>
  dplyr::arrange(patid_boot, tstart)

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
data_devel <- create_smoking_dummies(data_devel)

### Sort by patid and tstart, just for when looking at the dataset
data_devel <- dplyr::arrange(data_devel, patid, tstart)

### Create appropriate offsets for timevarying variables (interventions applied during follow-up)
data_devel$offset_statins_timevar_lnHR <- lnHR_statins*data_devel$med_status_adj_statins
data_devel$offset_ah_timevar_lnHR <- lnHR_ah*data_devel$med_status_adj_ah
data_devel$offset_smoking_timevar_dummy1_lnHR <- lnHR_smoking_dummy1_total*data_devel$med_status_adj_smoking_dummy1
data_devel$offset_smoking_timevar_dummy2_lnHR <- lnHR_smoking_dummy2_total*data_devel$med_status_adj_smoking_dummy2

### Fit a cox model and get bhaz
print(paste("fit model", Sys.time()))
fit <- survival::coxph(model.formula.full, data = data_devel, model = TRUE)
print(paste("bhaz fit model", Sys.time()))
bhaz <- survival::basehaz(fit, centered = TRUE)
bhaz.uncent <- survival::basehaz(fit, centered = FALSE)

### Save fit
saveRDS(fit, paste("data/p5/optimism_boot_cox_", gender, "_iter", iter, ".rds", sep = ""))
saveRDS(bhaz, paste("data/p5/optimism_boot_cox_bhaz_", gender, "_iter", iter, ".rds", sep = ""))
saveRDS(bhaz.uncent, paste("data/p5/optimism_boot_cox_bhaz_uncent_", gender, "_iter", iter, ".rds", sep = ""))

### Finished
print(paste("Step 2. complete", Sys.time()))

######################################################
######################################################
### Step 3. Estimate counterfactual survival times ###
######################################################
######################################################

###
### Want to estimate counterfactual survival times in the bootstrapped dataset, and one of the original datasets
###

### Define datasets
# Bootstrapped dataset (defined above)
df_boot <- readRDS(paste("data/p5/optimism_boot_df_imp_", gender, "_", iter, ".rds", sep = ""))
# Original dataset (iterate systematically through the different imputed original datasets)
if (gender == 2){
  df_orig <- readRDS(paste0("data/p5/dfs_all_female", (((iter - 1) %% 10) + 1), ".rds"))
} else if (gender == 1){
  df_orig <- readRDS(paste0("data/p5/dfs_all_male", (((iter - 1) %% 10) + 1), ".rds"))
}
# df_orig <- df_orig[1:200000, ] xx

### Read in the outcome data, which contains the split survival times
cohort_split_times <- readRDS(paste("data/p5/cohort_split_times_smoking_statins_antihypertensives_burnout", burnout, ".rds", sep = ""))

### Estimate counterfactual survival times in both datasets
cf_surv_times_boot <- estimate_cf_surv_times(df_boot, cohort_split_times, bhaz, type = "boot")
cf_surv_times_orig <- estimate_cf_surv_times(df_orig, cohort_split_times, bhaz, type = "orig")

saveRDS(cf_surv_times_boot, paste("data/p5/optimism_boot_cf_surv_times_", gender, "_iter", iter, "_dfboot", ".rds", sep = ""))
saveRDS(cf_surv_times_orig, paste("data/p5/optimism_boot_cf_surv_times_", gender, "_iter", iter, "_dforig", ".rds", sep = ""))
print(paste("Step 3. complete", Sys.time()))

##########################################################################
##########################################################################
### Step 4. Estimate calibration in bootstrapped and original datasets ###
##########################################################################
##########################################################################

### 
### Create variables for assessing performance within subgroups
###

### Create age variable categorised
df_boot$age_cat <- cut(df_boot$age, breaks = c(17,30,40,50,60,70,80,Inf))
df_orig$age_cat <- cut(df_orig$age, breaks = c(17,30,40,50,60,70,80,Inf))

###
### Load the region information
region_df <- readRDS("../../Aurum_Jun2021_extract/data/extraction/cohort_exclu3.rds") |>
  dplyr::select(patid, region)

### Read in lookup
region_lookup <- read.table("../../Aurum_Jun2021_extract/data/unzip/zLookups/region.txt", header = TRUE, sep = "\t") |>
  dplyr::rename(region = regionid, region_desc = Description)

region_df <- dplyr::left_join(region_df, region_lookup, by = dplyr::join_by(region)) |>
  dplyr::select(-region) |>
  dplyr::rename(region = region_desc) |>
  dplyr::mutate(region = as.factor(region))

### Create region variable
df_boot <- dplyr::left_join(df_boot, region_df, by = dplyr::join_by(patid))
df_orig <- dplyr::left_join(df_orig, region_df, by = dplyr::join_by(patid))

### Create variable with forms that can be used for prediction, and agrees with 'prescribed treatment strategy'
### which is what the individual is doing
df_boot <- 
  dplyr::mutate(df_boot,
                ### Create offset vars
                offset_statins_timevar_lnHR = 0,
                offset_ah_timevar_lnHR = 0,
                offset_smoking_timevar_dummy1_lnHR = 0,
                offset_smoking_timevar_dummy2_lnHR = 0)

df_orig <- 
  dplyr::mutate(df_orig,
                ### Create offset vars
                offset_statins_timevar_lnHR = 0,
                offset_ah_timevar_lnHR = 0,
                offset_smoking_timevar_dummy1_lnHR = 0,
                offset_smoking_timevar_dummy2_lnHR = 0)

### Estimate survival probabilities in each cohort
df_boot$surv <-as.numeric(est_surv_offset_prototype3(newdata = df_boot, 
                                                     fit = fit, 
                                                     bhaz = bhaz, 
                                                     time = round(10*365.25)))

df_orig$surv <-as.numeric(est_surv_offset_prototype3(newdata = df_orig, 
                                                     fit = fit, 
                                                     bhaz = bhaz, 
                                                     time = round(10*365.25)))

### Merge counter factual survival times with data frames
df_boot_cf <- dplyr::left_join(dplyr::select(df_boot, -cvd_time, -cvd_indicator), cf_surv_times_boot, by = dplyr::join_by(patid_boot)) |>
  dplyr::rename(cvd_time = cvd_time_cf, cvd_indicator = cvd_indicator_cf)

df_orig_cf <- dplyr::left_join(dplyr::select(df_orig, -cvd_time, -cvd_indicator), cf_surv_times_orig, by = dplyr::join_by(patid)) |>
  dplyr::rename(cvd_time = cvd_time_cf, cvd_indicator = cvd_indicator_cf)

###
### Function that will produce a calibration curve, and estimate ICI, E50 and E90
### Note, this is a wrapper for est_calib_ph, which is the underpinning function that applies the methodology
### Variable type defines whether validation occurs in bootstrapped dataset, or original dataset
###
est_calib_ph_wrapper <- function(type){
  
  ### Choose validation data frame and survival probabilities
  if (type == "boot"){
    data_valid_pred <- df_boot_cf
  } else if (type == "orig"){
    data_valid_pred <- df_orig_cf
  }
  
  ##########################
  ### Assess calibration ###
  ##########################
  
  ##############################
  ### PH regression approach ###
  ##############################
  print(paste("PH", Sys.time()))
  calib.ph <- est_calib_ph(data = data_valid_pred, 
                           fit = fit, 
                           bhaz = bhaz, 
                           time = round(10*365.25),
                           surv = data_valid_pred$surv,
                           pred_plot_range = seq(0.0001, 0.4, length.out = 250),
                           plot = FALSE)
  
  ### Create output.object
  output.object <- list("calib_data_plot" = calib.ph[["calib_data_plot"]],
                        "ICI" = as.numeric(calib.ph[["ICI"]]),
                        "E50" = as.numeric(calib.ph[["E50"]]),
                        "E90" =  as.numeric(calib.ph[["E90"]]))
  
  ### Return output object
  return(output.object)
  
}

###
### Write a function that will estimate calibration metrics (ICI, E50, E90) within subgroups defined by a variable (must be a factor)
### Variable type defines whether validation occurs in bootstrapped dataset, or original dataset
###
est_calib_ph_wrapper_byvar <- function(type,
                                       byvar){
  
  ### Choose validation data frame and survival probabilities
  if (type == "boot"){
    data_valid_pred <- df_boot_cf
  } else if (type == "orig"){
    data_valid_pred <- df_orig_cf
  }
  
  ### Remove people with missing region variable
  if (byvar == "region"){
    data_valid_pred <- subset(data_valid_pred, !is.na(region))
  }
  
  ### Subset within subgroups of interest
  bylevels <- levels(data_valid_pred[[byvar]])
  data_valid_pred_list <- lapply(bylevels, function(x) {data_valid_pred[data_valid_pred[,byvar] == x, ]})
  
  ##########################
  ### Assess calibration ###
  ##########################
  
  ##############################
  ### PH regression approach ###
  ##############################
  print(paste("PH", Sys.time()))
  calib.ph <- lapply(1:length(data_valid_pred_list), function(x) {
    est_calib_ph(data = data_valid_pred_list[[x]], 
                 fit = fit, 
                 bhaz = bhaz, 
                 time = round(10*365.25),
                 surv = data_valid_pred_list[[x]]$surv,
                 plot = FALSE)
  })
  
  ####################################################
  ### Function to extract list of ICI, E50 and E90 ###
  ####################################################
  
  ### Get metrics
  ICI_out <- lapply(1:length(data_valid_pred_list), function(x){calib.ph[[x]][["ICI"]]})
  E50_out <- lapply(1:length(data_valid_pred_list), function(x){calib.ph[[x]][["E50"]]})
  E90_out <- lapply(1:length(data_valid_pred_list), function(x){calib.ph[[x]][["E90"]]})
  names(ICI_out) <- paste0(toupper(substr(bylevels, 1, 1)), substring(bylevels, 2))
  names(E50_out) <- paste0(toupper(substr(bylevels, 1, 1)), substring(bylevels, 2))
  names(E90_out) <- paste0(toupper(substr(bylevels, 1, 1)), substring(bylevels, 2))
  
  ## Get N
  N_out <- sapply(data_valid_pred_list, nrow)
  names(N_out) <- paste0(toupper(substr(bylevels, 1, 1)), substring(bylevels, 2))
  
  ## Get calibration plot data
  calib_data_plot <- lapply(1:length(data_valid_pred_list), function(x){calib.ph[[x]][["calib_data_plot"]]}) 
  names(calib_data_plot) <- paste0(toupper(substr(bylevels, 1, 1)), substring(bylevels, 2))
  
  ### Create output.object
  output.object <- list("ICI" = ICI_out,
                        "E50" = E50_out,
                        "E90" = E90_out,
                        "N" = N_out,
                        "calib_data_plot" = calib_data_plot)
  
  ### Return output object
  return(output.object)
  
}

##########################
### Assess calibration ###
##########################

############################
### Within entire cohort ###
############################

calib_entire_boot <- est_calib_ph_wrapper(type = "boot")
calib_entire_orig <- est_calib_ph_wrapper(type = "orig")

### Save output
saveRDS(calib_entire_boot[["ICI"]], paste0("data/p5/optimism_boot_calib_ph_ICI_", gender, "_iter", iter, "_dfboot", ".rds"))
saveRDS(calib_entire_boot[["E50"]], paste0("data/p5/optimism_boot_calib_ph_E50_", gender, "_iter", iter, "_dfboot", ".rds"))
saveRDS(calib_entire_boot[["E90"]], paste0("data/p5/optimism_boot_calib_ph_E90_", gender, "_iter", iter, "_dfboot", ".rds"))
saveRDS(calib_entire_boot[["calib_data_plot"]], paste0("data/p5/optimism_boot_calib_ph_data_range_", gender, "_iter", iter, "_dfboot", ".rds"))

saveRDS(calib_entire_orig[["ICI"]], paste0("data/p5/optimism_boot_calib_ph_ICI_", gender, "_iter", iter, "_dforig", ".rds"))
saveRDS(calib_entire_orig[["E50"]], paste0("data/p5/optimism_boot_calib_ph_E50_", gender, "_iter", iter, "_dforig", ".rds"))
saveRDS(calib_entire_orig[["E90"]], paste0("data/p5/optimism_boot_calib_ph_E90_", gender, "_iter", iter, "_dforig", ".rds"))
saveRDS(calib_entire_orig[["calib_data_plot"]], paste0("data/p5/optimism_boot_calib_ph_data_range_", gender, "_iter", iter, "_dforig", ".rds"))

########################
### Within subgroups ###
########################

calib_data_by_region_boot <- est_calib_ph_wrapper_byvar(type = "boot", byvar = "region")
calib_data_by_region_orig <- est_calib_ph_wrapper_byvar(type = "orig", byvar = "region")

calib_data_by_age_boot <- est_calib_ph_wrapper_byvar(type = "boot", byvar = "age_cat")
calib_data_by_age_orig <- est_calib_ph_wrapper_byvar(type = "orig", byvar = "age_cat")

calib_data_by_ethnicity_boot <- est_calib_ph_wrapper_byvar(type = "boot", byvar = "ethnicity")
calib_data_by_ethnicity_orig <- est_calib_ph_wrapper_byvar(type = "orig", byvar = "ethnicity")

### Save output
saveRDS(calib_data_by_region_boot, paste0("data/p5/optimism_boot_calib_ph_by_region_", gender, "_iter", iter, "_dfboot", ".rds"))
saveRDS(calib_data_by_age_boot, paste0("data/p5/optimism_boot_calib_ph_by_age_", gender, "_iter", iter, "_dfboot", ".rds"))
saveRDS(calib_data_by_ethnicity_boot, paste0("data/p5/optimism_boot_calib_ph_by_ethnicity_", gender, "_iter", iter, "_dfboot", ".rds"))

saveRDS(calib_data_by_region_orig, paste0("data/p5/optimism_boot_calib_ph_by_region_", gender, "_iter", iter, "_dforig", ".rds"))
saveRDS(calib_data_by_age_orig, paste0("data/p5/optimism_boot_calib_ph_by_age_", gender, "_iter", iter, "_dforig", ".rds"))
saveRDS(calib_data_by_ethnicity_orig, paste0("data/p5/optimism_boot_calib_ph_by_ethnicity_", gender, "_iter", iter, "_dforig", ".rds"))

print(paste("Step 4. complete", Sys.time()))

#############################################################################
#############################################################################
### Step 5. Estimate discrimination in bootstrapped and original datasets ###
#############################################################################
#############################################################################

###
### Write a function that will estimate discrimination.
### Variable type defines whether validation occurs in bootstrapped dataset, or original dataset
###
est_discrim_wrapper <- function(type){
  
  ### Choose validation data frame and survival probabilities
  if (type == "boot"){
    data_valid_pred <- df_boot_cf
  } else if (type == "orig"){
    data_valid_pred <- df_orig_cf
  }
  
  #############################
  ### Assess discrimination ###
  #############################
  
  ### Estimate C-statistics
  Cstat.object <- survival::concordance(survival::Surv(cvd_time, cvd_indicator) ~ surv, 
                                        data = data_valid_pred,
                                        timewt = "n",
                                        influence = 0)
  
  ### Extract C-stat
  Cstat <- unname(Cstat.object$concordance)

  ### Return
  return(Cstat)
  
}

###
### Write a function that will estimate discrimination within subgroups defined by a variable (must be a factor)
### Variable type defines whether validation occurs in bootstrapped dataset, or original dataset
###
est_discrim_wrapper_byvar <- function(type,
                                      byvar){
  
  ### Choose validation data frame and survival probabilities
  if (type == "boot"){
    data_valid_pred <- df_boot_cf
  } else if (type == "orig"){
    data_valid_pred <- df_orig_cf
  }
  
  ### Remove people with missing region variable
  if (byvar == "region"){
    data_valid_pred <- subset(data_valid_pred, !is.na(region))
  }
  
  ### Subset within subgroups of interest
  bylevels <- levels(data_valid_pred[[byvar]])
  data_valid_pred_list <- lapply(bylevels, function(x) {data_valid_pred[data_valid_pred[,byvar] == x, ]})
  
  #############################
  ### Assess discrimination ###
  #############################
  discrim_out <- lapply(1:length(data_valid_pred_list), function(x) {
    
    ### Estimate C-statistics
    Cstat.object <- survival::concordance(survival::Surv(cvd_time, cvd_indicator) ~ surv, 
                                          data = data_valid_pred_list[[x]],
                                          timewt = "n",
                                          influence = 0)
    
    ### Extract C-stat
    Cstat <- unname(Cstat.object$concordance)
    
    ### Return
    return(Cstat)
    
  })
  
  ### Extract data in vector
  Cstat <- unlist(discrim_out, use.names = FALSE)
  
  ### Get N
  N_out <- sapply(data_valid_pred_list, nrow)
  
  ### Create table
  discrim_data_table <- data.frame("N" = N_out,
                                   "Cstat" = Cstat)
  rownames(discrim_data_table) <- bylevels
  
  ### Return output object
  return(discrim_data_table)
  
}

#############################
### Assess discrimination ###
#############################

############################
### Within entire cohort ###
############################

discrim_entire_boot <- est_discrim_wrapper(type = "boot")
discrim_entire_orig <- est_discrim_wrapper(type = "orig")

saveRDS(discrim_entire_boot, paste0("data/p5/optimism_boot_discrim_", gender, "_iter", iter, "_dfboot", ".rds"))
saveRDS(discrim_entire_orig, paste0("data/p5/optimism_boot_discrim_", gender, "_iter", iter, "_dforig", ".rds"))


########################
### Within subgroups ###
########################

discrim_by_region_boot <- est_discrim_wrapper_byvar(type = "boot", byvar = "region")
discrim_by_region_orig <- est_discrim_wrapper_byvar(type = "orig", byvar = "region")

discrim_by_age_boot <- est_discrim_wrapper_byvar(type = "boot", byvar = "age_cat")
discrim_by_age_orig <- est_discrim_wrapper_byvar(type = "orig", byvar = "age_cat")

discrim_by_ethnicity_boot <- est_discrim_wrapper_byvar(type = "boot", byvar = "ethnicity")
discrim_by_ethnicity_orig <- est_discrim_wrapper_byvar(type = "orig", byvar = "ethnicity")

### Save
saveRDS(discrim_by_region_boot, paste0("data/p5/optimism_boot_discrim_by_region_", gender, "_iter", iter, "_dfboot", ".rds"))
saveRDS(discrim_by_age_boot, paste0("data/p5/optimism_boot_discrim_by_age_", gender, "_iter", iter, "_dfboot", ".rds"))
saveRDS(discrim_by_ethnicity_boot, paste0("data/p5/optimism_boot_discrim_by_ethnicity_", gender, "_iter", iter, "_dfboot", ".rds"))

saveRDS(discrim_by_region_orig, paste0("data/p5/optimism_boot_discrim_by_region_", gender, "_iter", iter, "_dforig", ".rds"))
saveRDS(discrim_by_age_orig, paste0("data/p5/optimism_boot_discrim_by_age_", gender, "_iter", iter, "_dforig", ".rds"))
saveRDS(discrim_by_ethnicity_orig, paste0("data/p5/optimism_boot_discrim_by_ethnicity_", gender, "_iter", iter, "_dforig", ".rds"))

print(paste("Step 5. complete", Sys.time()))
