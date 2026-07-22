############################################################
############################################################
### The goal of programs p5.1 is to estimate relative risk ratios
### for sustaining statin/antihypertensives use for a period of time.
###
### This will be done using artificial censoring. We create cohorts
### that are on or off treatment at baseline, and artificially censor
### as soon as someone deviates from this strategy.
### 
### We will estimate inverse probability of artificial censoring weights,
### inverse probability of censoring weights, and inverse probability of treatment
### weights. These combined weights should target an estimate in a pseudo-population
### where everybody received the assigned treatment at baseline and sustained for the
### period of follow up.
###
### The goal of this program is to estimate the weights for this analysis. We 
### know we will run into positivity issues, so we then create two overlap cohorts.
###
### The first overlap cohort only considers individuals with at least 1% probability
### of receiving either treatment at baseline.
###
### The second overlap cohort reduces to individuals with at least 1% probability of
### receiving and sustaining treatment.
###
### In both overlap cohorts, the weights must be re-estimated.
### and possible violation of the positivity assumption.
###
### Everything will be done seperately for men/women, statins/antihypertensives,
### interactions or no interactions in the models for estimating weights, and 1/3/5/10 year follow-up.
###
### We do the following:
### - Step 1: Create validation dataset with appropriate variable for all analyses
### - Step 2: Estimate the various weights
### - Step 3: Create an overlap cohorts based
### - Step 4: Re-estimate weights in overlap cohort
###
### NB: This program is slightly inefficient in terms of disk space usage.
### We create a new dataset of the entire validation cohort(df_valid_ipacw) dependent on t_fup, 
### as this changes the definition of the outcome variables, and weight estimation variables. We
### save these as seperate datasets, which is wasteful, as the baseline predictors variables are the same.
### Do not think this is a major issue, but flagging in case changes needed in the future.
############################################################

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load survival
library(survival)

### Define medication, gender and read in dataset
### Extract gender and medication 
args <- commandArgs(trailingOnly = T)
gender_in <- as.numeric(args[1])
med <- args[2]
t_fup_int <- as.numeric(args[3])
t_fup <- t_fup_int*365.25
model_flex <- args[4]
print(paste("gender = ", gender_in))
print(paste("medication = ", med))
print(paste("t_fup = ", t_fup))
print(paste("model flexibility = ", model_flex))

### Define burnout
burnout <- 180

### Tuning parameter (can be 0.01, 0.05, 0.1, etc.)
alpha <- 0.01

### Read in variables that were interacted with the spline for age in prediction model
inter_age_rcs <- readRDS(paste("data/var_model_inter_age_rcs", gender_in, ".rds", sep = ""))

##################################################################################
##################################################################################
### Step 1: Start by reading in the validation dataset, merging with outcome data, and
### defining outcome variables for estimating the weights
##################################################################################
##################################################################################

### Read in validation dataset
df_valid <- readRDS(paste("../project2/data/p4/mice_mids_prototype3_", 2, "_", 1, ".rds", sep = ""))
df_valid <- mice::complete(df_valid, action = 1)

### Create indicators for estimating probability of being censored (ipcw's)
df_valid <- dplyr::mutate(df_valid, 
                          ipcw_indicator = dplyr::case_when(cvd_time <= t_fup ~ 1 - cvd_indicator,
                                                            cvd_time > t_fup ~ 0),
                          ipcw_time = dplyr::case_when(cvd_time <= t_fup ~ cvd_time,
                                                       cvd_time > t_fup ~ t_fup)
)

###
### Read in times at which individuals change medication status, and the time-varying variable for medication status
### NB: We don't use the augmented medication times (used for fitting the prediction model), where
### individuals medication always takes value 0 at zero, and changes when individual deviates from medication
### status at baseline. Instead, we use the original medication status times, as we want to stratify dependent
### on treatment status at baseline when estimating the ipacw weights.
medication_times <- readRDS(paste("../project2/data/cohort_split_times_", med, "_burnout", burnout, ".rds", sep = ""))
colnames(medication_times)

### Create variables for original cvd_time and cvd_indcator in validation dataset
df_valid <- dplyr::mutate(df_valid, time_original = cvd_time, status_original = cvd_indicator)

### Create data.tables for merge
data.table::setDT(df_valid)
data.table::setDT(medication_times)

### Setkeys
data.table::setkey(medication_times, patid)

### Cols to drop
cols_to_drop <- c("cvd_time", "cvd_indicator", "cvd_time_prim", "cvd_indicator_prim", 
                  "cvd_time_hes", "cvd_indicator_hes", "cvd_time_death", "cvd_indicator_death", 
                  "cvd_ev_prim_aj")

### Inner join
df_valid_ipacw <- medication_times[
  df_valid[ , !cols_to_drop, with = FALSE],
  on = "patid",
  nomatch = NULL
][order(gender, patid, tstart)]
colnames(medication_times)

### Create ipacw indicator and time, and reduce dataset to first row for each individual, 
### which is the time they deviate from baseline treatment strategy
df_valid_ipacw <- df_valid_ipacw[
  # Create ipacw event and time
  , `:=`(
    ipacw_indicator = as.integer(.N > 1),
    ipacw_time = cvd_time[1]
  ),
  by = patid][tstart == 0] |>
  as.data.frame()

### View these
# View(df_valid_ipacw[1:20,])

### Note that in this dataset, 'cvd_time' and 'cvd_indicator' variables are already correctly defined for
### assessing calibration in the artificially censored cohort. cvd_indicator = 1 if the individual had
### an event before being censored by either mechanism, and = 0 otherwise. cvd_time = time until either
### event, or censored by either mechanism. This is because they were taken from the 'medication_times'
### dataset when merging above, not from the df_valid dataset.

### Just need to rename these to 'time' and 'status' to align with functions used later.
df_valid_ipacw <- dplyr::mutate(df_valid_ipacw, time = cvd_time, status = cvd_indicator)

### Create offset variables with appropriate names and values, for predictions under treatment
### strategy of 'do nothing'.
df_valid_ipacw <-
  dplyr::mutate(df_valid_ipacw,
                offset_statins_timevar_lnHR = 0,
                offset_ah_timevar_lnHR = 0)

###
### NB:
### ipacw_time = time until min of artificial censoring, censoring, event
### ipacw_indicator = indicator = 1 when artificial censoring occurs first, 0 otherwise
###
### ipcw_time = time until min of censoring, event
### ipcw_indicator = indicator = 1 when censoring occurs first, 0 otherwise
###
### cvd_time == time = time until min of artificial censoring, censoring, event
### cvd_indicator == status = 1 when event occurs first, 0 otherwise
###
### time_original = time until min of censoring, event
### status_original = 1 when event occurs first, 0 otherwise
###

###########################################################
###########################################################
### Step 2: Estimate weights in the full validation cohort
###########################################################
###########################################################

### Create a vector of terms
if (model_flex == "flex"){
  formula_vec <-
    c("rms::rcs(age, c(25, 40, 57.5, 75))",
      paste("age*", inter_age_rcs, sep = ""),
      "age*rms::rcs(IMD, c(1,10,20))",
      "age*rms::rcs(sbp, 4)",
      "age*rms::rcs(bmi, 4)",
      "age*rms::rcs(nonhdl, 4)")
} else if (model_flex == "linear"){
  formula_vec <-
    c("rms::rcs(age, c(25, 40, 57.5, 75))",
      paste(inter_age_rcs, sep = ""),
      "rms::rcs(IMD, c(1,10,20))",
      "rms::rcs(sbp, 4)",
      "rms::rcs(bmi, 4)",
      "rms::rcs(nonhdl, 4)")
}

### Function to estimate IPCWs
estimate_ipcws <- function(df){
  
  ### Define formula for estimating ipcw's
  formula_ipcw <- stats::as.formula(
    paste("survival::Surv(ipcw_time, ipcw_indicator) ~ ", paste(formula_vec, collapse = " + "))
  )
  
  ### Fit model to estimate ipcws
  fit_est_ipcw <- survival::coxph(formula_ipcw, data = df, model = TRUE)
  bhaz_est_ipcw <- survival::basehaz(fit_est_ipcw, centered = TRUE)
  
  ### Estimate the IPCWs
  ### These are estimated at the minimum of event time or censoring in the artificially censored cohort (ipacw_time)
  ### NB: Estimating at ipcw_time does not actually impact analysis, as only individuals for which ipcw_time != ipacw_time
  ### are artificially censored individuals, who contribute a zero*weight to the outcome
  df$ipcw_weight <- 1/est_surv_eventtime_or_t(newdata = df,
                                              fit = fit_est_ipcw,
                                              bhaz = bhaz_est_ipcw,
                                              eventtime = "ipacw_time",
                                              t = t_fup)
  
  ### Add survival probabilities at t_fup (used to define overlap cohort)
  df$ipcw_surv <- est_surv(
    newdata = df,
    fit     = fit_est_ipcw,
    bhaz    = bhaz_est_ipcw,
    t       = t_fup
  )
  
  return(df)
  
}

### Function to estimate IPACWs
estimate_ipacws <- function(df){
  
  ### We are going to build stratified models for estimating the IPACW's, given the probability of remaining
  ### off treatment (if off treatment at baseline), will be distributed differently than the probability of
  ### remaining on treatment (if on treatment at baseline).
  
  ### Create formula for predicting ipacws
  formula_ipacw <- stats::as.formula(
    paste("survival::Surv(ipacw_time, ipacw_indicator) ~ ", paste(formula_vec, collapse = " + "))
  )
  
  ### Create cohorts
  df_stratified_list <- 
    lapply(c(0,1), function(baseline_group){
      df_out <- df
      df_out <- subset(df_out, med_status == baseline_group)
    })
  
  ### Check they sum to 1,000,000
  testthat::expect_equal(nrow(df), nrow(do.call("rbind", df_stratified_list)))
  ### And print sizes of each cohort
  lapply(1:2, function(baseline_group){
    nrow(df_stratified_list[[baseline_group]])
  })
  
  ### Fit models for estimating the IPACWs and estimate IPACWs
  print("FIT MODELS")
  df_stratified_list <- lapply(0:1, function(baseline_group){
    
    print(paste("gender = ", gender_in, ", baseline_group = ", baseline_group, ",", Sys.time(), sep = ""))
    
    ### Extract relevant development dataset
    df_stratified <- df_stratified_list[[(baseline_group+1)]]
    
    ### Fit model
    fit_est_ipacw <- survival::coxph(formula_ipacw, data = df_stratified, model = TRUE)
    bhaz_est_ipacw <- survival::basehaz(fit_est_ipacw, centered = TRUE)
    
    ### The output from est_surv_eventtime_or_t is a survival probability, so the probability of having 
    ### not been artificially censored by min of (t, eventtime)
    ipacw_temp <- est_surv_eventtime_or_t(newdata = df_stratified,
                                          fit = fit_est_ipacw,
                                          bhaz = bhaz_est_ipacw,
                                          eventtime = "ipacw_time",
                                          t = t_fup)
    
    ### The weights are 1/surv
    ### Assign these to the dataset
    df_stratified$ipacw_weight <- 1/ipacw_temp
    
    ### The output from est_surv is a survival probability, 
    ### so the probability of having not been artificially censored by time t (used for defining overlap cohort)
    ipacw_temp_surv <- est_surv(newdata = df_stratified,
                                fit = fit_est_ipacw,
                                bhaz = bhaz_est_ipacw,
                                t = t_fup)
    
    ### Assign survival probabilities for defining overlap cohort
    df_stratified$ipacw_surv <- ipacw_temp_surv
    
    return(df_stratified)
    
  })
  
  ### Recombine into a single dataset
  df <- do.call("rbind", df_stratified_list) |> as.data.frame()
  
  return(df)
  
}

### Function to estimate IPTWs
estimate_iptws <- function(df){
  
  ### Define formula for estimating IPTWs
  formula_iptw <- stats::as.formula(
    paste("med_status ~ ", paste(formula_vec, collapse = " + "))
  )
  
  ### Fit a logistic regression model to predict this outcome
  fit_iptw <-
    glm(formula_iptw,
        family = binomial(), 
        data = df)
  
  ### Get predictions
  df$iptw_p <- predict(fit_iptw, 
                       newdata = df, 
                       type = "response")
  
  ### Assign weight for status at baseline
  ### I use stabilised weights, but shouldn't matter for this particular estimator
  pA <- mean(df$med_status == 1)
  df <- dplyr::mutate(df,
                      iptw_weight =
                        dplyr::case_when(
                          med_status == 0 ~ (1-pA)/(1-iptw_p),
                          med_status == 1 ~ pA/iptw_p
                        ))
  
  return(df)
  
}

#################################################
### Estimate survival probabilities for IPCWs ###
#################################################

### Estimate IPCWs IPACWs and IPTWs
df_valid_ipacw <- estimate_ipcws(df_valid_ipacw)
print(paste("Estimated IPCWs in entire cohort", Sys.time()))
df_valid_ipacw <- estimate_ipacws(df_valid_ipacw)
print(paste("Estimated IPACWs in entire cohort", Sys.time()))
df_valid_ipacw <- estimate_iptws(df_valid_ipacw)
print(paste("Estimated IPTWs in entire cohort", Sys.time()))

###
### Function to calculate combined weights
###
create_combined_weights <- function(df){
  
  #### Combined censoring weights
  df$comb_weight_censoring <- df$ipacw_weight * df$ipcw_weight
  
  ### Create combined weights to include treatment at baseline (probability of receiving assigned 
  ### treatment and sustaining it)
  df$comb_weight_baseline <- 
    df$iptw_weight*df$comb_weight_censoring
  
  return(df)
  
}

### Apply function
df_valid_ipacw <- create_combined_weights(df_valid_ipacw)

### Save cohort
saveRDS(df_valid_ipacw, paste("data/df_valid_ipacw_", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
print(paste("Saved df_valid_ipacw", Sys.time()))

##################################################
##################################################
### Step 2: Create overlap cohort based on IPTWs
##################################################
##################################################

### This overlap cohort is individuals with a probability > 1% of being assigned either treatment

### Remove weights to avoid possible mis-use of these variables within the overlap cohort, as they will be recalculated
df_overlap <- dplyr::select(df_valid_ipacw, -c(ipacw_weight, ipcw_weight, iptw_weight, 
                                               comb_weight_censoring, comb_weight_baseline))
print(paste("colnames = ", colnames(df_overlap)))

### For IPTWs we use lower and upper bound
lower_cut <- alpha
upper_cut <- 1-alpha

df_overlap <- df_overlap[
  df_overlap$iptw_p >= lower_cut &
    df_overlap$iptw_p <= upper_cut,
]

### Treatment prevalence shift
mean(df_valid_ipacw$med_status)
mean(df_overlap$med_status)

cat("Original N:", nrow(df_valid_ipacw), "\n")
cat("Overlap N:", nrow(df_overlap), "\n\n")

cat("Treatment distribution (overlap):\n")
print(table(df_overlap$med_status))

cat("\nIPACW summary (overlap):\n")
print(summary(df_overlap$ipacw_surv))

### Save cohort
saveRDS(df_overlap, paste("data/df_valid_overlap_cohort", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
print(paste("Saved overlap cohort", Sys.time()))


##################################################################
##################################################################
### Step 3: Create overlap cohort based on IPTWs, IPACWs and IPCWs
##################################################################
##################################################################

### This overlap cohort is individuals with a probability > 1% of being assigned either treatment, and
### a > 1% probability of sustaining the treatment they did receieve, and a > 1% probability of being observed
### until t_fup (i.e. not being censored).
df_overlap_comb <- df_overlap

### For IPACWs, one-sided lower bound only (no positivity issue with people with high probability of
### sustaining treatment)
df_overlap_comb <- df_overlap_comb[df_overlap_comb$ipacw_surv >= alpha, ]

### For IPCWs, one-sided lower bound only
df_overlap_comb <- df_overlap_comb[df_overlap_comb$ipcw_surv >= alpha, ]

### Note, we used survival probabilities at 10-years (ipacw_surv and ipcw_surv), rather than the weights,
### which are evaluated at the minimum of censoring time and t_fup.

### Treatment prevalence shift
mean(df_valid_ipacw$med_status)
mean(df_overlap_comb$med_status)

cat("Original N:", nrow(df_valid_ipacw), "\n")
cat("Overlap N:", nrow(df_overlap_comb), "\n\n")

cat("Treatment distribution (overlap):\n")
print(table(df_overlap_comb$med_status))

cat("\nIPACW summary (overlap):\n")
print(summary(df_overlap_comb$ipacw_surv))

### Save cohort
saveRDS(df_overlap_comb, paste("data/df_valid_overlap_comb_cohort", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
print(paste("Saved overlap cohort comb", Sys.time()))

### Note, we could try and define the overlap cohort based on a combined probability of all three,
### such as the following code, but it's a bit unclear what this is. Inparticular,
### when it's the probability of treatment you were not assigned to, multiplied by the probability of sustaining 
### treatment you were assigned to.

# ### Create combined survival probabilities for being on assigned treatment at baseline and sustaining
# df_valid_ipacw$comb_surv_on <- 
#   df_valid_ipacw$ipcw_surv*df_valid_ipacw$ipacw_surv*df_valid_ipacw$iptw_p
# 
# ### Create combined survival probabilities for being off assigned treatment at baseline and sustaining
# df_valid_ipacw$comb_surv_off <- 
#   df_valid_ipacw$ipcw_surv*df_valid_ipacw$ipacw_surv*(1 - df_valid_ipacw$iptw_p)
# 
# ### Apply overlap cohort
# df_overlap_comb <- df_overlap[
#   df_valid_ipacw$comb_surv_on >= lower_cut &
#     df_valid_ipacw$comb_surv_on >= lower_cut,
# ]

############################################################
############################################################
### Step 4: Estimate weights in overlap cohorts 
############################################################
############################################################

print(paste("Estimate weights in overlap cohort", Sys.time()))

########################
### Estimate IPCW's  ###
########################

### Apply function
df_overlap <- estimate_ipcws(df_overlap)
df_overlap_comb <- estimate_ipcws(df_overlap_comb)

#########################
### Estimate IPACW's  ###
#########################

### Apply function
df_overlap <- estimate_ipacws(df_overlap)
df_overlap_comb <- estimate_ipacws(df_overlap_comb)

######################
### Estimate IPTWs ###
######################

### Apply function
df_overlap <- estimate_iptws(df_overlap)
df_overlap_comb <- estimate_iptws(df_overlap_comb)

##################################
### Calculate combined weights ###
##################################

### Apply function
df_overlap <- create_combined_weights(df_overlap)
df_overlap_comb <- create_combined_weights(df_overlap_comb)

### Save overlap cohorts
saveRDS(df_overlap, paste("data/df_valid_overlap_cohort", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
saveRDS(df_overlap_comb, paste("data/df_valid_overlap_comb_cohort", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
print(paste("Saved overlap cohorts with estimated weights", Sys.time()))
