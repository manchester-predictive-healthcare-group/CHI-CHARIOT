###
### Program to fit causal model
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

for (gender in c(2,1)){
  
  gender_char <- c("male", "female")[gender]
  print(paste("gender = ", gender_char))
  print(paste("burnout = ", burnout))
  
  ###
  ### Define formula
  ###
  
  ### Read in variables that will be interacted with the spline for age
  inter.age.rcs <- readRDS(paste("data/p4/var_model_inter_age_rcs", gender, ".rds", sep = ""))
  
  ### Create vector with all terms for the formula
  full.formula.vec <- c("rms::rcs(age, c(25, 40, 57.5, 75))*rms::rcs(IMD, c(1,10,20))", 
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
  
  ###
  ### Write a function to fit model and save it for a given imputed dataset
  ###
  fit_model <- function(imp.num){
    
    ### Define imputed dataset
    imp <- imp.list[[imp.num]]
    
    # # ### REDUCE### TO REMOVE
    # imp <- imp[!is.na(fastmatch::fmatch(imp$patid, imp$patid[1:5000])), ]
    # cohort.split.times <- cohort.split.times[!is.na(fastmatch::fmatch(cohort.split.times$patid, imp$patid[1:5000])), ]
    # # 
    
    ### Merge imputed dataset with cohort.split.times
    data.devel <- dplyr::left_join(dplyr::select(imp, -c("cvd_time", "cvd_indicator",
                                                         "cvd_time_prim", "cvd_time_hes", "cvd_time_death",
                                                         "cvd_indicator_prim", "cvd_indicator_hes", "cvd_indicator_death")),
                                   cohort.split.times,
                                   by = dplyr::join_by("patid"))
    
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
    
    ### SBP relative to 120 (no benefit for being lower than 120)
    ### BMI relative to 25 (no benefit for being lower than 25)
    ### nonhdl relative to 4
    data.devel <- dplyr::mutate(data.devel,
                                sbp_adj = dplyr::case_when(sbp < 120 ~ 0,
                                                           TRUE ~ sbp - 120),
                                bmi_adj = dplyr::case_when(bmi < 25 ~ 0,
                                                           TRUE ~ bmi - 25),
                                nonhdl_adj = nonhdl - 4
    )
    
    ### Sort by patid and tstart, just for when looking at the dataset
    data.devel <- dplyr::arrange(data.devel, patid, tstart)
    
    ### Create appropriate offsets for timevarying variables (interventions applied during follow-up)
    data.devel$offset_statins_timevar_lnHR <- lnHR_statins*data.devel$med_status_adj_statins
    data.devel$offset_ah_timevar_lnHR <- lnHR_ah*data.devel$med_status_adj_ah
    data.devel$offset_smoking_timevar_dummy1_lnHR <- lnHR_smoking_dummy1_total*data.devel$med_status_adj_smoking_dummy1
    data.devel$offset_smoking_timevar_dummy2_lnHR <- lnHR_smoking_dummy2_total*data.devel$med_status_adj_smoking_dummy2
    
    ### Create appropriate offsets for baseline variables (i.e. for dropping in effect of variables at baseline)
    data.devel$offset_smoking_dummy1_lnHR <- lnHR_smoking_dummy1_direct*data.devel$smoking_baseline_dummy1
    data.devel$offset_smoking_dummy2_lnHR <- lnHR_smoking_dummy2_direct*data.devel$smoking_baseline_dummy2
    data.devel$offset_sbp_lnHR <- lnHR_sbp*data.devel$sbp_adj
    data.devel$offset_bmi_lnHR <- lnHR_bmi*data.devel$bmi_adj
    data.devel$offset_nonhdl_lnHR <- lnHR_nonhdl*data.devel$nonhdl_adj
    
    ### Fit a cox model and get bhaz
    print(paste("fit model", imp.num, Sys.time()))
    fit <- survival::coxph(model.formula.full, data = data.devel, model = TRUE)
    print(paste("bhaz fit model", imp.num, Sys.time()))
    bhaz <- survival::basehaz(fit, centered = TRUE)
    bhaz.uncent <- survival::basehaz(fit, centered = FALSE)
    
    ### Save fit
    saveRDS(fit, paste("data/p4/prototype3_4knots_cox_", gender, "_imp", imp.num, ".rds", sep = ""))
    saveRDS(bhaz, paste("data/p4/prototype3_4knots_cox_bhaz_", gender, "_imp", imp.num, ".rds", sep = ""))
    saveRDS(bhaz.uncent, paste("data/p4/prototype3_4knots_cox_bhaz_uncent_", gender, "_imp", imp.num, ".rds", sep = ""))
    
  }
  
  ### Apply function
  for (imp.num in 1:10){
    print(paste("Model = ", imp.num, Sys.time()))
    fit_model(imp.num)
  }
  
  ### Finished
  print(paste("FINISHED", Sys.time()))
  
}
