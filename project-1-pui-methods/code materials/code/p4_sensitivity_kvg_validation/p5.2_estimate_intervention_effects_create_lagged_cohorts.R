###
### Want to do a sensitivity analysis where we do not consider statin/antihypertensives
### in the 1-year prior to a CVD event
### 
### This will change the artificial censoring times and the outcome in the artificially censored cohort.
### These will therefore need to be re-calculated in each cohort.
###
### It will also change the make up of the overlap cohort which is dependent on surv_ipacw, so this cohort will
### need to be re-defined (df_overlap_comb).
###

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
# gender_in <- 1
# med <- "statins"
# t_fup_int <- as.numeric(5)
# t_fup <- t_fup_int*365.25
# model_flex <- "linear"
gender_in <- as.numeric(args[1])
med <- args[2]
t_fup_int <- as.numeric(args[3])
t_fup <- t_fup_int*365.25
model_flex <- args[4]

### Tuning parameter (can be 0.01, 0.05, 0.1, etc.)
alpha <- 0.01

### Read in variables that were interacted with the spline for age in prediction model
inter_age_rcs <- readRDS(paste("data/var_model_inter_age_rcs", gender_in, ".rds", sep = ""))

### Read in the three cohorts
# Entire cohort
df_valid_ipacw <- readRDS(paste("data/df_valid_ipacw_", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
# Overlap cohort
df_overlap <- readRDS(paste("data/df_valid_overlap_cohort", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))

### For individuals that were artificially censored (changed treatment) within 1-year of CVD event, 
### don't artificially censor them

### Create a function to create variables for estimating artificial censoring weights, and CVD outcome,
### where we don't consider a change in prescription if it happens within 1-year prior to a CVD event
### Call these 'lagged'
create_lag_df <- function(df){
  dplyr::mutate(df,
                ipacw_indicator_lag = dplyr::case_when(ipacw_time + 365 <= time_original ~ ipacw_indicator,
                                                       TRUE ~ 0),
                ipacw_time_lag = dplyr::case_when(ipacw_time + 365 <= time_original ~ ipacw_time,
                                                  TRUE ~ time_original),
                status_lag = dplyr::case_when(ipacw_time + 365 <= time_original ~ status,
                                              TRUE ~ status_original),
                time_lag = dplyr::case_when(ipacw_time + 365 <= time_original ~ time,
                                            TRUE ~ time_original))
}

### Apply the function
df_valid_ipacw <- create_lag_df(df_valid_ipacw)
df_overlap <- create_lag_df(df_overlap)

##########################################################################################
### Define functions for estimating weights. For IPACW, it's on the lagged variables.
### For IPCW and IPTW, its same as program 5.1. For entire cohort and first overlap cohort,
### don't need to recalculate these as they will be the same. For the second overlap cohort,
### which is dependent on the ipacw's, they will need to be recalculated
###########################################################################################

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

### Function to estimate ipacws on lagged variables
reestimate_ipacws_lag <- function(df){
  
  ### Create formula for predicting ipacws
  formula_ipacw <- stats::as.formula(
    paste("survival::Surv(ipacw_time_lag, ipacw_indicator_lag) ~ ", paste(formula_vec, collapse = " + "))
  )
  
  ### Create cohorts for estimating weights of staying on/off treatment
  df_stratified_list <- 
    lapply(c(0,1), function(baseline_group){
      df_out <- df
      df_out <- subset(df_out, med_status == baseline_group)
    })
  
  ### Check they sum
  testthat::expect_equal(nrow(df), nrow(do.call("rbind", df_stratified_list)))
  
  ### Fit models for estimating the IPACWs and estimate IPACWs
  print("FIT MODELS")
  df_stratified_list <- lapply(0:1, function(baseline_group){
    
    print(paste("gender = ", gender_in, ", baseline_group = ", baseline_group, ",", Sys.time(), sep = ""))
    
    ### Extract relevant development dataset
    df_stratified <- df_stratified_list[[(baseline_group+1)]]
    
    ### Fit model
    fit_est_ipacw <- survival::coxph(formula_ipacw, data = df_stratified, model = TRUE)
    bhaz_est_ipacw <- survival::basehaz(fit_est_ipacw, centered = TRUE)
    
    ### The output from est_surv is a survival probability, so the probability of having not been artificially censored by time t
    ipacw_temp_surv <- est_surv(newdata = df_stratified,
                                fit = fit_est_ipacw,
                                bhaz = bhaz_est_ipacw,
                                t = t_fup)
    
    ### Assign survival probabilities for defining overlap cohort
    df_stratified$ipacw_surv_lag <- ipacw_temp_surv
    
    ### The output from est_surv_eventtime_or_t is a survival probability, so the probability of having 
    ### not been artificially censored by min of (t, eventtime)
    ipacw_temp <- est_surv_eventtime_or_t(newdata = df_stratified,
                                          fit = fit_est_ipacw,
                                          bhaz = bhaz_est_ipacw,
                                          eventtime = "ipacw_time_lag",
                                          t = t_fup)
    
    ### The weights are 1/surv
    ### Assign these to the dataset
    df_stratified$ipacw_weight_lag <- 1/ipacw_temp
    
    return(df_stratified)
    
  })
  
  ### Recombine into a single dataset
  df <- do.call("rbind", df_stratified_list) |> as.data.frame()
  
  return(df)
  
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
                                              eventtime = "ipcw_time",
                                              t = t_fup)
  
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

### Re-estimate IPACWs in overall cohort and original cohort
df_valid_ipacw <- reestimate_ipacws_lag(df_valid_ipacw)
df_overlap <- reestimate_ipacws_lag(df_overlap)

### Function to calculate combined weights
create_combined_weights <- function(df){
  
  #### Combined censoring weights
  df$comb_weight_censoring_lag <- df$ipacw_weight_lag * df$ipcw_weight
  
  ### Create combined weights to include treatment at baseline (probability of receiving assigned 
  ### treatment and sustaining it)
  df$comb_weight_baseline_lag <- 
    df$iptw_weight*df$comb_weight_censoring_lag
  
  return(df)
  
}

### Calculate combined weights
df_valid_ipacw <- create_combined_weights(df_valid_ipacw)
df_overlap <- create_combined_weights(df_overlap)

### Save cohort
### NB: We overwrite the files created in p5.1, because we are just creating additional variables. The cohort is the same
### group of people
saveRDS(df_valid_ipacw, paste("data/df_valid_ipacw_", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
saveRDS(df_overlap, paste("data/df_valid_overlap_cohort", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
print(paste("Added lagged outcome and ipacw weights to entire cohort and overlap cohort"))


###
### Create new overlap cohort based on combined weights
###

### Tuning parameter (can be 0.01, 0.05, 0.1, etc.)
alpha <- 0.01

### This overlap cohort is individuals with a probability > 1% of being assigned either treatment, and
### a > 1% probability of sustaining the treatment they did receieve, and a > 1% probability of being observed
### until t_fup (i.e. not being censored).
df_overlap_comb <- df_overlap

### For IPACWs, one-sided lower bound only (no positivity issue with people with high probability of
### sustaining treatment)
### 
### NB: We are using updated ipacw_surv_lag to define this cohort now
df_overlap_comb <- df_overlap_comb[df_overlap_comb$ipacw_surv_lag >= alpha, ]

### For IPCWs, one-sided lower bound only ()
df_overlap_comb <- df_overlap_comb[df_overlap_comb$ipcw_surv >= alpha, ]

### Remove weights that have just been calculated, 
### to avoid possible mis-use of these variables as they will be recalculated
df_overlap_comb <- dplyr::select(df_overlap_comb, -c(ipacw_weight_lag, ipacw_weight, ipcw_weight, iptw_weight, 
                                               comb_weight_censoring_lag, comb_weight_baseline_lag))
print(paste("colnames = ", colnames(df_overlap_comb)))

### Create lag outcome variables for re-estimating weights
df_overlap_comb <- create_lag_df(df_overlap_comb)

### Estimate lagged IPACWs
df_overlap_comb <- reestimate_ipacws_lag(df_overlap_comb)
### Estimate lagged IPCWs
df_overlap_comb <- estimate_ipcws(df_overlap_comb)
### Estimate lagged IPTWs
df_overlap_comb <- estimate_iptws(df_overlap_comb)
### Create combined weights
df_overlap_comb <- create_combined_weights(df_overlap_comb)

### Save cohort
### NB: We create a new cohort files compared to the one created in p5.1, 
### because this overlap cohort now represents a different group of people, and we have re-calculated a number of
### pre-existing variables
saveRDS(df_overlap_comb, paste("data/df_valid_overlap_comb_lag_cohort", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
print(paste("Saved overlap cohort comb", Sys.time()))
