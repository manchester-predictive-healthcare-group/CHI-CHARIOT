###
### Program to fit causal model
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
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
  
  ### Get the knot locations for continuous variables
  knot_loc_sbp <- readRDS(paste("data/p5/knots_loc_sbp_", gender, sep = ""))
  knot_loc_sbp
  knot_loc_nonhdl <- readRDS(paste("data/p5/knots_loc_nonhdl_", gender, sep = ""))
  knot_loc_nonhdl
  knot_loc_bmi <- readRDS(paste("data/p5/knots_loc_bmi_", gender, sep = ""))
  knot_loc_bmi
  
  ### Read in the imputed development datasets
  imp_list <- readRDS(paste("data/p5/dfs_all_", gender_char, ".rds", sep = ""))
  
  ###
  ### Write a function to fit model and save it for a given imputed dataset
  ###
  fit_model <- function(imp_num){
    
    ### Define imputed dataset
    imp <- imp_list[[imp_num]]
    
    ### Merge imputed dataset with cohort_split_times
    data_devel <- dplyr::left_join(dplyr::select(imp, -c("cvd_time", "cvd_indicator",
                                                         "cvd_time_prim", "cvd_time_hes", "cvd_time_death",
                                                         "cvd_indicator_prim", "cvd_indicator_hes", "cvd_indicator_death")),
                                   cohort_split_times,
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
    data_devel <- create_smoking_dummies(data_devel)
    
    ### Sort by patid and tstart, just for when looking at the dataset
    data_devel <- dplyr::arrange(data_devel, patid, tstart)
    
    ### Create appropriate offsets for time varying variables (interventions applied during follow-up)
    data_devel$offset_statins_timevar_lnHR <- lnHR_statins*data_devel$med_status_adj_statins
    data_devel$offset_ah_timevar_lnHR <- lnHR_ah*data_devel$med_status_adj_ah
    data_devel$offset_smoking_timevar_dummy1_lnHR <- lnHR_smoking_dummy1_total*data_devel$med_status_adj_smoking_dummy1
    data_devel$offset_smoking_timevar_dummy2_lnHR <- lnHR_smoking_dummy2_total*data_devel$med_status_adj_smoking_dummy2
    
    ### Fit a cox model and get bhaz
    print(paste("fit model", imp_num, Sys.time()))
    fit <- survival::coxph(model.formula.full, data = data_devel, model = TRUE)
    print(paste("bhaz fit model", imp_num, Sys.time()))
    bhaz <- survival::basehaz(fit, centered = TRUE)
    bhaz.uncent <- survival::basehaz(fit, centered = FALSE)
    
    ### Save fit
    saveRDS(fit, paste("data/p5/model_cox_", gender, "_imp", imp_num, ".rds", sep = ""))
    saveRDS(bhaz, paste("data/p5/model_cox_bhaz_", gender, "_imp", imp_num, ".rds", sep = ""))
    saveRDS(bhaz.uncent, paste("data/p5/model_cox_bhaz_uncent_", gender, "_imp", imp_num, ".rds", sep = ""))
    
  }
  
  ### Apply function
  for (imp_num in 1:10){
    print(paste("Model = ", imp_num, Sys.time()))
    fit_model(imp_num)
  }
  
  ### Finished
  print(paste("FINISHED", Sys.time()))
  
}
