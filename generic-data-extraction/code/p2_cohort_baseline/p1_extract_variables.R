### This program will be used to extract variables at baseline.
### It will call on a variety of functions extract_... (stored in /R/)

### BASIC CONCEPTS

### READ IN WHICH VARIABLE TO EXTRACT
### READ IN A RANGE OF INPUT CRITERIA TO THE EXTRACT FUNCTIONS
### - cohort
### - indexdt name
### - t values
### - upper and low bounds for tests data

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Create variable list
variable.list <- c(
  ### Demographics
  "age"#1
  ,"ethnicity"#2
  ### Medical diagnoses
  ,"hypertension"#3
  ,"ra"#4
  ,"af"#5
  ,"ckd"#6
  ,"smi"#7
  ,"fhcvd"#8
  ,"migraine"#9
  ,"sle"#10
  ,"diabetes"#11
  ### Medical or drug
  ,"impotence"#12
  ### Drug Products
  ,"cortico"#13
  ,"antipsy"#14
  ### Test data
  ,"bmi"#15
  ,"sbp"#16
  ,"sbp_var"#17
  ,"cholhdl_ratio"#18
  ,"smoking"#19
  ### Outcomes
  ,"time_until_cvd"#20
  ### Treatments at baseline
  ,"antihypertensives"#21
  ,"statins"#22
  ### Other medical diagnoses
  ,"copd"#23
  ,"int_dis"#24
  ,"downs"#25
  ,"oral_cancer"#26
  ,"brain_cancer"#27
  ,"lung_cancer"#28
  ,"blood_cancer"#29
  ,"pre_eclampsia"#30
  ,"postnatal_depression"#31
  
  ###
  ### Exposed/unexposed 180-day window for exposure status
  ### 
  
  ### Unexposed
  ,"sbp_unexposed"#32
  ,"cholhdl_ratio_unexposed"#33
  ,"sbp_var_unexposed"#34
  ### SBP var exposed to prescriptions (no notion of mr/np/md for this variable, as not picking one test observation)
  ,"sbp_var_exposed"#35
  ### Test data exposed to prescriptions, most recent
  ,"sbp_exposed_mr"#36
  ,"cholhdl_ratio_exposed_mr"#37
  ### Test data exposed to prescriptions, highest number of prescriptions in window
  ,"sbp_exposed_np"#38
  ,"cholhdl_ratio_exposed_np"#39
  ### Test data exposed to prescriptions, closest to precription
  ,"sbp_exposed_md"#40
  ,"cholhdl_ratio_exposed_md"#41
  
  ###
  ### Exposed/unexposed 90-day window for exposure status
  ###
  
  ### Unexposed
  ,"sbp_unexposed_90"#42
  ,"cholhdl_ratio_unexposed_90"#43
  ,"sbp_var_unexposed_90"#44
  ### SBP var exposed to prescriptions, 90 day window (no notion of mr/np/md for this variable, as not picking one test observation)
  ,"sbp_var_exposed_90"#45
  ### Test data exposed to prescriptions, 90 day window, most recent
  ,"sbp_exposed_90_mr"#46
  ,"cholhdl_ratio_exposed_90_mr"#47
  ### Test data exposed to prescriptions, 90 day window, highest number of prescriptions in window
  ,"sbp_exposed_90_np"#48
  ,"cholhdl_ratio_exposed_90_np"#49
  ### Test data exposed to prescriptions, 90 day window, closest to prescription
  ,"sbp_exposed_90_md"#50
  ,"cholhdl_ratio_exposed_90_md"#51
  
  ###
  ### Exposed 90-day window for exposure status, require at least 2 prescriptions in window, 
  ###
  
  ### Test data exposed to prescriptions, 90 day window, most recent
  ,"sbp_exposed_90_mr_rn2"#52
  ,"cholhdl_ratio_exposed_90_mr_rn2"#53
  ### Test data exposed to prescriptions, 90 day window, highest number of prescriptions in window
  ,"sbp_exposed_90_np_rn2"#54
  ,"cholhdl_ratio_exposed_90_np_rn2"#55
  ### Test data exposed to prescriptions, 90 day window, closest to prescription
  ,"sbp_exposed_90_md_rn2"#56
  ,"cholhdl_ratio_exposed_90_md_rn2"#57
  
  ###
  ### Technically I shouldn't need any of the above, I am repeating lots of computation here.
  ### I should be extract ALL exposed values, that have at least 1 prescription in the X days before, and extract
  ### observation date, number in window, min diff, etc. Then I can pay around with window, exposure type, etc.
  ###
  
  ### Test data unexposed to prescriptions, in 180 day window
  ,"sbp_unexposed_all"#58
  ,"sbp_expanded_unexposed_all"#59 (using expanded definition of antihypertensives)
  ,"cholhdl_ratio_unexposed_all"#60
  
  ### Test data exposed to prescriptions, in 180 day window
  ,"sbp_exposed_all"#61  
  ,"sbp_expanded_exposed_all"#62 (using expanded definition of antihypertensives)
  ,"cholhdl_ratio_exposed_all"#63
  
  ### Test data exposed to prescriptions, 90 day window
  ,"sbp_exposed_90_all"#64  
  ,"sbp_expanded_exposed_90_all"#65 (using expanded definition of antihypertensives)
  ,"cholhdl_ratio_exposed_90_all"#66
  
  ### Nonhdl
  ,"nonhdl"#67
  
)


### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
variable.name <- variable.list[as.numeric(args[1])]
cohort.base.name <- args[2]
cohort.base <- readRDS(paste("data/extraction/", cohort.base.name ,".rds", sep = ""))
cohort.base <- data.table::as.data.table(cohort.base)
indexdt.name <- args[3]
out.subdir <- args[4]
t <- as.numeric(args[5])
if (t == 0){t <- NULL}
try(task.id <- as.numeric(args[6]))

print(paste("variable = ", variable.name, sep = ""))
print(paste("cohort.base.name = ", cohort.base.name, sep = ""))
print(paste("indexdt.name = ", indexdt.name, sep = ""))
print(paste("out.subdir = ", out.subdir, sep = ""))
print(paste("t = ", t, sep = ""))
if (exists("task.id")){print(paste("task.id = ", task.id, sep = ""))}

### Apply extraction
if (variable.name == "age"){
  print(paste("age", Sys.time()))
  var_age <- extract_age(cohort = cohort.base, 
                         indexdt = indexdt.name, 
                         out.save.disk = TRUE, 
                         out.filepath = NULL, 
                         out.subdir = out.subdir, 
                         return.output = TRUE)
  quantile(var_age$age)
} else if (variable.name == "ethnicity"){
  print(paste("ethnicity", Sys.time()))
  var_ethnicity <- extract_ethnicity(cohort = cohort.base, 
                                     db = "aurum",
                                     db.linked = "aurum_linked",
                                     out.save.disk = TRUE, 
                                     out.filepath = NULL, 
                                     out.subdir = out.subdir, 
                                     return.output = TRUE)
  table(var_ethnicity$ethnicity)
} else if (variable.name == "hypertension"){
  print(paste("hypertension", Sys.time()))
  var_hypertension <- extract_ho(cohort = cohort.base, 
                                 varname = "hypertension", 
                                 codelist = "edh_hypertension_medcodeid", 
                                 indexdt = indexdt.name, 
                                 db = "aurum", 
                                 tab = "obs",
                                 out.save.disk = TRUE, 
                                 out.filepath = NULL, 
                                 out.subdir = out.subdir, 
                                 return.output = TRUE)
  table(var_hypertension$hypertension)
} else if (variable.name == "ra"){
  print(paste("ra", Sys.time()))
  var_ra <- extract_ho(cohort = cohort.base, 
                       varname = "ra", 
                       codelist = "edh_ra_medcodeid", 
                       indexdt = indexdt.name, 
                       db = "aurum", 
                       tab = "obs",
                       out.save.disk = TRUE, 
                       out.filepath = NULL, 
                       out.subdir = out.subdir, 
                       return.output = TRUE)
  table(var_ra$ra)
} else if (variable.name == "af"){
  print(paste("af", Sys.time()))
  var_af <- extract_ho(cohort = cohort.base, 
                       varname = "af", 
                       codelist = "edh_af_medcodeid", 
                       indexdt = indexdt.name, 
                       db = "aurum", 
                       tab = "obs",
                       out.save.disk = TRUE, 
                       out.filepath = NULL, 
                       out.subdir = out.subdir, 
                       return.output = TRUE)
  table(var_af$af)
} else if (variable.name == "ckd"){
  print(paste("ckd", Sys.time()))
  var_ckd <- extract_ho(cohort = cohort.base, 
                        varname = "ckd", 
                        codelist = "edh_ckd_medcodeid", 
                        indexdt = indexdt.name, 
                        db = "aurum", 
                        tab = "obs",
                        out.save.disk = TRUE, 
                        out.filepath = NULL, 
                        out.subdir = out.subdir, 
                        return.output = TRUE)
  table(var_ckd$ckd)
} else if (variable.name == "smi"){
  print(paste("smi", Sys.time()))
  var_smi <- extract_ho(cohort = cohort.base, 
                        varname = "smi", 
                        codelist = "edh_smi_medcodeid", 
                        indexdt = indexdt.name, 
                        db = "aurum", 
                        tab = "obs",
                        out.save.disk = TRUE, 
                        out.filepath = NULL, 
                        out.subdir = out.subdir, 
                        return.output = TRUE)
  table(var_smi$smi)
} else if (variable.name == "fhcvd"){
  print(paste("fhcvd", Sys.time()))
  var_fhcvd <- extract_ho(cohort = cohort.base, 
                          varname = "fhcvd", 
                          codelist = "edh_fhcvd_medcodeid", 
                          indexdt = indexdt.name, 
                          db = "aurum", 
                          tab = "obs",
                          out.save.disk = TRUE, 
                          out.filepath = NULL, 
                          out.subdir = out.subdir, 
                          return.output = TRUE)
  table(var_fhcvd$fhcvd)
} else if (variable.name == "migraine"){
  print(paste("migraine", Sys.time()))
  var_migraine <- extract_ho(cohort = cohort.base, 
                             varname = "migraine", 
                             codelist = "edh_migraine_medcodeid", 
                             indexdt = indexdt.name, 
                             db = "aurum", 
                             tab = "obs",
                             out.save.disk = TRUE, 
                             out.filepath = NULL, 
                             out.subdir = out.subdir, 
                             return.output = TRUE)
  table(var_migraine$migraine)
} else if (variable.name == "sle"){
  print(paste("sle", Sys.time()))
  var_sle <- extract_ho(cohort = cohort.base, 
                        varname = "sle", 
                        codelist = "edh_sle_medcodeid",  
                        indexdt = indexdt.name, 
                        db = "aurum", 
                        tab = "obs",
                        out.save.disk = TRUE, 
                        out.filepath = NULL, 
                        out.subdir = out.subdir, 
                        return.output = TRUE)
  table(var_sle$sle)
} else if (variable.name == "diabetes"){
  print(paste("diabetes", Sys.time()))
  var_diabetes <- extract_diabetes(cohort = cohort.base, 
                                   codelist.type1 = "edh_t1dia_medcodeid", 
                                   codelist.type2 = "edh_t2dia_medcodeid",
                                   indexdt = indexdt.name, 
                                   db = "aurum", 
                                   out.save.disk = TRUE, 
                                   out.filepath = NULL, 
                                   out.subdir = out.subdir, 
                                   return.output = TRUE)
  table(var_diabetes$diabetes)
} else if (variable.name == "impotence"){
  print(paste("impotence", Sys.time()))
  var_impotence <- extract_impotence(cohort = cohort.base, 
                                     codelist.med = "edh_impotence_medcodeid", 
                                     codelist.drug = "uom_impotence_prodcodeid",
                                     indexdt = indexdt.name, 
                                     db = "aurum", 
                                     out.save.disk = TRUE, 
                                     out.filepath = NULL, 
                                     out.subdir = out.subdir, 
                                     return.output = TRUE)
  table(var_impotence$impotence)
} else if (variable.name == "cortico"){
  print(paste("cortico", Sys.time()))
  var_cortico <- extract_ho(varname = "cortico", 
                            codelist = "uom_oral_corticosteroids_prodcodeid", 
                            cohort = cohort.base, 
                            indexdt = indexdt.name, 
                            time.prev = round(365.25*0.5),
                            numobs = 2,
                            db = "aurum", 
                            tab = "drug",
                            out.save.disk = TRUE, 
                            out.filepath = NULL, 
                            out.subdir = out.subdir, 
                            return.output = TRUE)
  table(var_cortico$cortico)
} else if (variable.name == "antipsy"){
  print(paste("antipsy", Sys.time()))
  var_antipsy <- extract_ho(varname = "antipsy", 
                            codelist = "uom_antipsychotics_prodcodeid", 
                            cohort = cohort.base, 
                            indexdt = indexdt.name, 
                            time.prev = round(365.25*0.5),
                            numobs = 2,
                            db = "aurum", 
                            tab = "drug",
                            out.save.disk = TRUE, 
                            out.filepath = NULL, 
                            out.subdir = out.subdir, 
                            return.output = TRUE)
  table(var_antipsy$antipsy)
} else if (variable.name == "bmi"){
  print(paste("bmi", Sys.time()))
  var_bmi <- extract_bmi(cohort = cohort.base, 
                         codelist.bmi = "edh_bmi_medcodeid",
                         codelist.weight = "weight_medcodeid",
                         codelist.height = "height_medcodeid",
                         indexdt = indexdt.name, 
                         lower.bound = 18, 
                         upper.bound = 47,
                         db = "aurum", 
                         out.save.disk = TRUE, 
                         out.filepath = NULL, 
                         out.subdir = out.subdir, 
                         return.output = TRUE)
  quantile(var_bmi$bmi, na.rm = TRUE)
} else if (variable.name == "sbp"){
  print(paste("sbp", Sys.time()))
  var_sbp <- extract_sbp(cohort = cohort.base, 
                         codelist = "edh_sbp_medcodeid",
                         indexdt = indexdt.name, 
                         lower.bound = 70, 
                         upper.bound = 210,
                         db = "aurum", 
                         out.save.disk = TRUE, 
                         out.filepath = NULL, 
                         out.subdir = out.subdir, 
                         return.output = TRUE)
  quantile(var_sbp$sbp, na.rm = TRUE)
} else if (variable.name == "sbp_var"){
  print(paste("sbp var", Sys.time()))
  var_sbp_var <- extract_sbp_var(cohort = cohort.base, 
                                 codelist = "edh_sbp_medcodeid",
                                 indexdt = indexdt.name, 
                                 lower.bound = 70, 
                                 upper.bound = 210,
                                 db = "aurum", 
                                 out.save.disk = TRUE, 
                                 out.filepath = NULL, 
                                 out.subdir = out.subdir, 
                                 return.output = TRUE)
  quantile(var_sbp_var$sbp_var, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio"){
  print(paste("cholhdl_ratio", Sys.time()))
  var_cholhdl_ratio <- extract_cholhdl_ratio(cohort = cohort.base, 
                                             codelist.ratio = "muzambi_ratio_medcodeid",
                                             codelist.chol = "muzambi_cholesterol_medcodeid",
                                             codelist.hdl = "muzambi_hdl_medcodeid",
                                             indexdt = indexdt.name, 
                                             lower.bound = 1,
                                             upper.bound = 12,
                                             db = "aurum",
                                             out.save.disk = TRUE, 
                                             out.filepath = NULL, 
                                             out.subdir = out.subdir, 
                                             return.output = TRUE)
  quantile(var_cholhdl_ratio$cholhdl_ratio, na.rm = TRUE)
} else if (variable.name == "smoking"){
  print(paste("smoking", Sys.time()))
  var_smoking <- extract_smoking(cohort = cohort.base, 
                                 codelist.non = "edh_smoking_non_medcodeid",
                                 codelist.ex = "edh_smoking_ex_medcodeid",
                                 codelist.light = "edh_smoking_light_medcodeid",
                                 codelist.mod = "edh_smoking_mod_medcodeid",
                                 codelist.heavy = "edh_smoking_heavy_medcodeid",
                                 indexdt = indexdt.name, 
                                 db = "aurum",
                                 out.save.disk = TRUE, 
                                 out.filepath = NULL, 
                                 out.subdir = out.subdir, 
                                 return.output = TRUE)
  table(var_smoking$smoking)
  str(var_smoking)
} else if (variable.name == "time_until_cvd"){
  print(paste("outcome", Sys.time()))
  var_time_until_cvd <- extract_time_until_cvd(cohort.base, 
                                               varname.time = NULL, 
                                               varname.indicator = NULL,
                                               codelist.primary = "edh_cvd_medcodeid",
                                               codelist.hes = "uom_cvd_icd",
                                               indexdt = indexdt.name, 
                                               censdt = "fup_end",
                                               t = NULL, 
                                               t.varname = TRUE,
                                               db.primary = "aurum",
                                               db.hes = "aurum_linked",
                                               db.filepath = NULL,
                                               out.save.disk = TRUE, 
                                               out.filepath = NULL, 
                                               out.subdir = out.subdir,  
                                               return.output = TRUE)
  str(var_time_until_cvd)
  table(var_time_until_cvd$cvd_indicator_prim)
  table(var_time_until_cvd$cvd_indicator_hes)
  table(var_time_until_cvd$cvd_indicator_death)
  table(var_time_until_cvd$cvd_indicator)
} else if (variable.name == "antihypertensives"){
  print(paste("antihypertensives", Sys.time()))
  var_antihypertensives <- extract_ho(cohort = cohort.base, 
                                      varname = "antihypertensives", 
                                      codelist = "uom_antihypertensives_prodcodeid", 
                                      indexdt = indexdt.name, 
                                      time.prev = round(365.25*0.5),
                                      numobs = 2,
                                      db = "aurum", 
                                      tab = "drug",
                                      out.save.disk = TRUE, 
                                      out.filepath = NULL, 
                                      out.subdir = out.subdir, 
                                      return.output = TRUE)
  table(var_antihypertensives$antihypertensives)
} else if (variable.name == "statins"){
  print(paste("statins", Sys.time()))
  var_statins <- extract_ho(cohort = cohort.base, 
                            varname = "statins", 
                            codelist = "uom_statins_prodcodeid", 
                            indexdt = indexdt.name, 
                            time.prev = round(365.25*0.5),
                            numobs = 2,
                            db = "aurum", 
                            tab = "drug",
                            out.save.disk = TRUE, 
                            out.filepath = NULL, 
                            out.subdir = out.subdir, 
                            return.output = TRUE)
  table(var_statins$statins)
} else if (variable.name == "copd"){
  print(paste("copd", Sys.time()))
  var_copd <- extract_ho(cohort = cohort.base, 
                         varname = "copd", 
                         codelist = "ah_copd",  
                         indexdt = indexdt.name, 
                         db = "aurum", 
                         tab = "obs",
                         out.save.disk = TRUE, 
                         out.filepath = NULL, 
                         out.subdir = out.subdir, 
                         return.output = TRUE)
  table(var_copd$copd)
} else if (variable.name == "int_dis"){
  print(paste("int_dis", Sys.time()))
  var_int_dis <- extract_ho(cohort = cohort.base, 
                            varname = "int_dis", 
                            codelist = "ah_intellectual_disability",  
                            indexdt = indexdt.name, 
                            db = "aurum", 
                            tab = "obs",
                            out.save.disk = TRUE, 
                            out.filepath = NULL, 
                            out.subdir = out.subdir, 
                            return.output = TRUE)
  table(var_int_dis$int_dis)
} else if (variable.name == "downs"){
  print(paste("downs", Sys.time()))
  var_downs <- extract_ho(cohort = cohort.base, 
                          varname = "downs", 
                          codelist = "ah_downs_syndrome",  
                          indexdt = indexdt.name, 
                          db = "aurum", 
                          tab = "obs",
                          out.save.disk = TRUE, 
                          out.filepath = NULL, 
                          out.subdir = out.subdir, 
                          return.output = TRUE)
  table(var_downs$downs)
} else if (variable.name == "brain_cancer"){
  print(paste("brain_cancer", Sys.time()))
  var_brain_cancer <- extract_ho(cohort = cohort.base, 
                                 varname = "brain_cancer", 
                                 codelist = "ah_brain_cancer",  
                                 indexdt = indexdt.name, 
                                 db = "aurum", 
                                 tab = "obs",
                                 out.save.disk = TRUE, 
                                 out.filepath = NULL, 
                                 out.subdir = out.subdir, 
                                 return.output = TRUE)
  table(var_brain_cancer$brain_cancer)
} else if (variable.name == "lung_cancer"){
  print(paste("lung_cancer", Sys.time()))
  var_lung_cancer <- extract_ho(cohort = cohort.base, 
                                varname = "lung_cancer", 
                                codelist = "ah_lung_cancer",  
                                indexdt = indexdt.name, 
                                db = "aurum", 
                                tab = "obs",
                                out.save.disk = TRUE, 
                                out.filepath = NULL, 
                                out.subdir = out.subdir, 
                                return.output = TRUE)
  table(var_lung_cancer$lung_cancer)
} else if (variable.name == "oral_cancer"){
  print(paste("oral_cancer", Sys.time()))
  var_oral_cancer <- extract_ho(cohort = cohort.base, 
                                varname = "oral_cancer", 
                                codelist = "ah_oral_cancer",  
                                indexdt = indexdt.name, 
                                db = "aurum", 
                                tab = "obs",
                                out.save.disk = TRUE, 
                                out.filepath = NULL, 
                                out.subdir = out.subdir, 
                                return.output = TRUE)
  table(var_oral_cancer$oral_cancer)
} else if (variable.name == "blood_cancer"){
  print(paste("blood_cancer", Sys.time()))
  var_blood_cancer <- extract_ho(cohort = cohort.base, 
                                 varname = "blood_cancer", 
                                 codelist = "ah_blood_cancer",  
                                 indexdt = indexdt.name, 
                                 db = "aurum", 
                                 tab = "obs",
                                 out.save.disk = TRUE, 
                                 out.filepath = NULL, 
                                 out.subdir = out.subdir, 
                                 return.output = TRUE)
  table(var_blood_cancer$blood_cancer)
} else if (variable.name == "pre_eclampsia"){
  print(paste("pre_eclampsia", Sys.time()))
  var_pre_eclampsia <- extract_ho(cohort = cohort.base, 
                                  varname = "pre_eclampsia", 
                                  codelist = "uom_pre_eclampsia_medcodeid",  
                                  indexdt = indexdt.name, 
                                  db = "aurum", 
                                  tab = "obs",
                                  out.save.disk = TRUE, 
                                  out.filepath = NULL, 
                                  out.subdir = out.subdir, 
                                  return.output = TRUE)
  table(var_pre_eclampsia$pre_eclampsia)
} else if (variable.name == "postnatal_depression"){
  print(paste("postnatal_depression", Sys.time()))
  var_postnatal_depression <- extract_postnatal_depression(cohort = cohort.base, 
                                                           varname = "postnatal_depression", 
                                                           codelist.med = "uom_postnatal_depression_medcodeid",  
                                                           codelist.test = "uom_postnatal_depression_score_medcodeid", 
                                                           indexdt = indexdt.name, 
                                                           db = "aurum", 
                                                           out.save.disk = TRUE, 
                                                           out.filepath = NULL, 
                                                           out.subdir = out.subdir, 
                                                           return.output = TRUE)
  table(var_postnatal_depression$postnatal_depression)
  ###
  ### Exposed and unexposed, 180-day window for exposure status
  ###
} else if (variable.name == "sbp_unexposed"){
  print(paste("sbp_unexposed", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_unexposed <- extract_test_unexposed(cohort = cohort.base, 
                                              varname = "sbp_unexposed",
                                              codelist.med = "edh_sbp_medcodeid",
                                              codelist.drug = "uom_antihypertensives_prodcodeid",
                                              query.med = query.med,
                                              query.drug = query.drug,
                                              task.id = task.id,
                                              task.id.n = 100,
                                              range.max = 180,
                                              range.min = 1,
                                              indexdt = indexdt.name, 
                                              lower.bound = 70, 
                                              upper.bound = 210,
                                              db = "aurum", 
                                              out.save.disk = TRUE, 
                                              out.filepath = NULL, 
                                              out.subdir = out.subdir, 
                                              return.output = TRUE)
  quantile(var_sbp_unexposed$sbp_unexposed, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_unexposed"){
  print(paste("cholhdl_ratio_unexposed", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_unexposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                                 varname = "cholhdl_ratio_unexposed",
                                                                 codelist.ratio = "muzambi_ratio_medcodeid",
                                                                 codelist.chol = "muzambi_cholesterol_medcodeid",
                                                                 codelist.hdl = "muzambi_hdl_medcodeid",
                                                                 codelist.drug = "uom_statins_prodcodeid",
                                                                 query.ratio = query.cholhdl,
                                                                 query.chol = query.chol,
                                                                 query.hdl = query.hdl,
                                                                 query.drug = query.drug,
                                                                 task.id = task.id,
                                                                 task.id.n = 100,
                                                                 range.max = 180,
                                                                 range.min = 1,
                                                                 indexdt = indexdt.name, 
                                                                 lower.bound = 1, 
                                                                 upper.bound = 12,
                                                                 db = "aurum", 
                                                                 out.save.disk = TRUE, 
                                                                 out.filepath = NULL, 
                                                                 out.subdir = out.subdir, 
                                                                 return.output = TRUE)
  quantile(var_cholhdl_ratio_unexposed$cholhdl_ratio_unexposed, na.rm = TRUE)
} else if (variable.name == "sbp_var_unexposed"){
  print(paste("sbp var", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_var_unexposed <-  extract_test_var_unexposed(cohort = cohort.base, 
                                                       varname = "sbp_var_unexposed",
                                                       codelist.med = "edh_sbp_medcodeid",
                                                       codelist.drug = "uom_antihypertensives_prodcodeid",
                                                       query.med = query.med,
                                                       query.drug = query.drug,
                                                       task.id = task.id,
                                                       task.id.n = 100,
                                                       range.max = 180,
                                                       range.min = 1,
                                                       indexdt = indexdt.name, 
                                                       lower.bound = 70, 
                                                       upper.bound = 210,
                                                       db = "aurum", 
                                                       out.save.disk = TRUE, 
                                                       out.filepath = NULL, 
                                                       out.subdir = out.subdir, 
                                                       return.output = TRUE)
  quantile(var_sbp_var_unexposed$sbp_var_unexposed, na.rm = TRUE)
} else if (variable.name == "sbp_var_exposed"){
  print(paste("sbp var", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_var_exposed <-  extract_test_var_unexposed(cohort = cohort.base, 
                                                     varname = "sbp_var_exposed",
                                                     codelist.med = "edh_sbp_medcodeid",
                                                     codelist.drug = "uom_antihypertensives_prodcodeid",
                                                     query.med = query.med,
                                                     query.drug = query.drug,
                                                     task.id = task.id,
                                                     task.id.n = 100,
                                                     range.max = 180,
                                                     range.min = 14,
                                                     keep.exposed = TRUE,
                                                     indexdt = indexdt.name, 
                                                     lower.bound = 70, 
                                                     upper.bound = 210,
                                                     db = "aurum", 
                                                     out.save.disk = TRUE, 
                                                     out.filepath = NULL, 
                                                     out.subdir = out.subdir, 
                                                     return.output = TRUE)
  quantile(var_sbp_var_exposed$sbp_var_exposed, na.rm = TRUE)
} else if (variable.name == "sbp_exposed_mr"){
  print(paste("sbp_exposed_mr", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_mr",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 180,
                                            range.min = 14,
                                            keep.exposed = TRUE,
                                            exposed.type = "mr",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_mr, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_mr"){
  print(paste("cholhdl_ratio_exposed_mr", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_mr",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 180,
                                                               range.min = 14,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "mr",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_mr, na.rm = TRUE)
} else if (variable.name == "sbp_exposed_np"){
  print(paste("sbp_exposed_np", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_np",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 180,
                                            range.min = 14,
                                            keep.exposed = TRUE,
                                            exposed.type = "np",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_np, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_np"){
  print(paste("cholhdl_ratio_exposed_np", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_np",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 180,
                                                               range.min = 14,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "np",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_np, na.rm = TRUE)
}  else if (variable.name == "sbp_exposed_md"){
  print(paste("sbp_exposed_md", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_md",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 180,
                                            range.min = 14,
                                            keep.exposed = TRUE,
                                            exposed.type = "md",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_md, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_md"){
  print(paste("cholhdl_ratio_exposed_md", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_md",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 180,
                                                               range.min = 14,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "md",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_md, na.rm = TRUE)
  ###
  ### Exposed and unexposed 90-day window for exposure status
  ###
} else if (variable.name == "sbp_unexposed_90"){
  print(paste("sbp_unexposed_90", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_unexposed <- extract_test_unexposed(cohort = cohort.base, 
                                              varname = "sbp_unexposed_90",
                                              codelist.med = "edh_sbp_medcodeid",
                                              codelist.drug = "uom_antihypertensives_prodcodeid",
                                              query.med = query.med,
                                              query.drug = query.drug,
                                              task.id = task.id,
                                              task.id.n = 100,
                                              range.max = 90,
                                              range.min = 1,
                                              indexdt = indexdt.name, 
                                              lower.bound = 70, 
                                              upper.bound = 210,
                                              db = "aurum", 
                                              out.save.disk = TRUE, 
                                              out.filepath = NULL, 
                                              out.subdir = out.subdir, 
                                              return.output = TRUE)
  quantile(var_sbp_unexposed$sbp_unexposed_90, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_unexposed_90"){
  print(paste("cholhdl_ratio_unexposed_90", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_unexposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                                 varname = "cholhdl_ratio_unexposed_90",
                                                                 codelist.ratio = "muzambi_ratio_medcodeid",
                                                                 codelist.chol = "muzambi_cholesterol_medcodeid",
                                                                 codelist.hdl = "muzambi_hdl_medcodeid",
                                                                 codelist.drug = "uom_statins_prodcodeid",
                                                                 query.ratio = query.cholhdl,
                                                                 query.chol = query.chol,
                                                                 query.hdl = query.hdl,
                                                                 query.drug = query.drug,
                                                                 task.id = task.id,
                                                                 task.id.n = 100,
                                                                 range.max = 90,
                                                                 range.min = 1,
                                                                 indexdt = indexdt.name, 
                                                                 lower.bound = 1, 
                                                                 upper.bound = 12,
                                                                 db = "aurum", 
                                                                 out.save.disk = TRUE, 
                                                                 out.filepath = NULL, 
                                                                 out.subdir = out.subdir, 
                                                                 return.output = TRUE)
  quantile(var_cholhdl_ratio_unexposed$cholhdl_ratio_unexposed_90, na.rm = TRUE)
} else if (variable.name == "sbp_var_unexposed_90"){
  print(paste("sbp var", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_var_unexposed <-  extract_test_var_unexposed(cohort = cohort.base, 
                                                       varname = "sbp_var_unexposed_90",
                                                       codelist.med = "edh_sbp_medcodeid",
                                                       codelist.drug = "uom_antihypertensives_prodcodeid",
                                                       query.med = query.med,
                                                       query.drug = query.drug,
                                                       task.id = task.id,
                                                       task.id.n = 100,
                                                       range.max = 90,
                                                       range.min = 1,
                                                       indexdt = indexdt.name, 
                                                       lower.bound = 70, 
                                                       upper.bound = 210,
                                                       db = "aurum", 
                                                       out.save.disk = TRUE, 
                                                       out.filepath = NULL, 
                                                       out.subdir = out.subdir, 
                                                       return.output = TRUE)
  quantile(var_sbp_var_unexposed$sbp_var_unexposed, na.rm = TRUE)
} else if (variable.name == "sbp_var_exposed_90"){
  print(paste("sbp var", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_var_exposed <-  extract_test_var_unexposed(cohort = cohort.base, 
                                                     varname = "sbp_var_exposed_90",
                                                     codelist.med = "edh_sbp_medcodeid",
                                                     codelist.drug = "uom_antihypertensives_prodcodeid",
                                                     query.med = query.med,
                                                     query.drug = query.drug,
                                                     task.id = task.id,
                                                     task.id.n = 100,
                                                     range.max = 90,
                                                     range.min = 14,
                                                     keep.exposed = TRUE,
                                                     indexdt = indexdt.name, 
                                                     lower.bound = 70, 
                                                     upper.bound = 210,
                                                     db = "aurum", 
                                                     out.save.disk = TRUE, 
                                                     out.filepath = NULL, 
                                                     out.subdir = out.subdir, 
                                                     return.output = TRUE)
  quantile(var_sbp_var_exposed$sbp_var_exposed_90, na.rm = TRUE)
} else if (variable.name == "sbp_exposed_90_mr"){
  print(paste("sbp_exposed_90_mr", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_90_mr",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 90,
                                            range.min = 14,
                                            keep.exposed = TRUE,
                                            exposed.type = "mr",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_90_mr, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_90_mr"){
  print(paste("cholhdl_ratio_exposed_90_mr", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_90_mr",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 90,
                                                               range.min = 14,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "mr",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_90_mr, na.rm = TRUE)
} else if (variable.name == "sbp_exposed_90_np"){
  print(paste("sbp_exposed_90_np", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_90_np",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 90,
                                            range.min = 14,
                                            keep.exposed = TRUE,
                                            exposed.type = "np",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_90_np, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_90_np"){
  print(paste("cholhdl_ratio_exposed_90_np", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_90_np",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 90,
                                                               range.min = 14,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "np",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_90_np, na.rm = TRUE)
}  else if (variable.name == "sbp_exposed_90_md"){
  print(paste("sbp_exposed_90_md", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_90_md",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 90,
                                            range.min = 14,
                                            keep.exposed = TRUE,
                                            exposed.type = "md",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_90_md, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_90_md"){
  print(paste("cholhdl_ratio_exposed_90_md", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_90_md",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 90,
                                                               range.min = 14,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "md",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_90_md, na.rm = TRUE)
  ###
  ### Exposed and unexposed 90-day window for exposure status, require at least 2 prescriptions in window, 
  ### only look backwards in 3-year period for SBP value.
  ###
} else if (variable.name == "sbp_unexposed_90"){
  print(paste("sbp_unexposed_90", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_unexposed <- extract_test_unexposed(cohort = cohort.base, 
                                              varname = "sbp_unexposed_90",
                                              codelist.med = "edh_sbp_medcodeid",
                                              codelist.drug = "uom_antihypertensives_prodcodeid",
                                              query.med = query.med,
                                              query.drug = query.drug,
                                              task.id = task.id,
                                              task.id.n = 100,
                                              range.max = 90,
                                              range.min = 1,
                                              indexdt = indexdt.name, 
                                              lower.bound = 70, 
                                              upper.bound = 210,
                                              db = "aurum", 
                                              out.save.disk = TRUE, 
                                              out.filepath = NULL, 
                                              out.subdir = out.subdir, 
                                              return.output = TRUE)
  quantile(var_sbp_unexposed$sbp_unexposed_90, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_unexposed_90"){
  print(paste("cholhdl_ratio_unexposed_90", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_unexposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                                 varname = "cholhdl_ratio_unexposed_90",
                                                                 codelist.ratio = "muzambi_ratio_medcodeid",
                                                                 codelist.chol = "muzambi_cholesterol_medcodeid",
                                                                 codelist.hdl = "muzambi_hdl_medcodeid",
                                                                 codelist.drug = "uom_statins_prodcodeid",
                                                                 query.ratio = query.cholhdl,
                                                                 query.chol = query.chol,
                                                                 query.hdl = query.hdl,
                                                                 query.drug = query.drug,
                                                                 task.id = task.id,
                                                                 task.id.n = 100,
                                                                 range.max = 90,
                                                                 range.min = 1,
                                                                 indexdt = indexdt.name, 
                                                                 lower.bound = 1, 
                                                                 upper.bound = 12,
                                                                 db = "aurum", 
                                                                 out.save.disk = TRUE, 
                                                                 out.filepath = NULL, 
                                                                 out.subdir = out.subdir, 
                                                                 return.output = TRUE)
  quantile(var_cholhdl_ratio_unexposed$cholhdl_ratio_unexposed_90, na.rm = TRUE)
} else if (variable.name == "sbp_var_unexposed_90"){
  print(paste("sbp var", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_var_unexposed <-  extract_test_var_unexposed(cohort = cohort.base, 
                                                       varname = "sbp_var_unexposed_90",
                                                       codelist.med = "edh_sbp_medcodeid",
                                                       codelist.drug = "uom_antihypertensives_prodcodeid",
                                                       query.med = query.med,
                                                       query.drug = query.drug,
                                                       task.id = task.id,
                                                       task.id.n = 100,
                                                       range.max = 90,
                                                       range.min = 1,
                                                       indexdt = indexdt.name, 
                                                       lower.bound = 70, 
                                                       upper.bound = 210,
                                                       db = "aurum", 
                                                       out.save.disk = TRUE, 
                                                       out.filepath = NULL, 
                                                       out.subdir = out.subdir, 
                                                       return.output = TRUE)
  quantile(var_sbp_var_unexposed$sbp_var_unexposed, na.rm = TRUE)
} else if (variable.name == "sbp_var_exposed_90"){
  print(paste("sbp var", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_var_exposed <-  extract_test_var_unexposed(cohort = cohort.base, 
                                                     varname = "sbp_var_exposed_90",
                                                     codelist.med = "edh_sbp_medcodeid",
                                                     codelist.drug = "uom_antihypertensives_prodcodeid",
                                                     query.med = query.med,
                                                     query.drug = query.drug,
                                                     task.id = task.id,
                                                     task.id.n = 100,
                                                     range.max = 90,
                                                     range.min = 14,
                                                     range.num = 2,
                                                     keep.exposed = TRUE,
                                                     indexdt = indexdt.name, 
                                                     lower.bound = 70, 
                                                     upper.bound = 210,
                                                     db = "aurum", 
                                                     out.save.disk = TRUE, 
                                                     out.filepath = NULL, 
                                                     out.subdir = out.subdir, 
                                                     return.output = TRUE)
  quantile(var_sbp_var_exposed$sbp_var_exposed_90, na.rm = TRUE)
} else if (variable.name == "sbp_exposed_90_mr_rn2"){
  print(paste("sbp_exposed_90_mr_rn2", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_90_mr_rn2",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 90,
                                            range.min = 14,
                                            range.num = 2,
                                            keep.exposed = TRUE,
                                            exposed.type = "mr",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_90_mr_rn2, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_90_mr_rn2"){
  print(paste("cholhdl_ratio_exposed_90_mr_rn2", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_90_mr_rn2",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 90,
                                                               range.min = 14,
                                                               range.num = 2,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "mr",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_90_mr_rn2, na.rm = TRUE)
} else if (variable.name == "sbp_exposed_90_np_rn2"){
  print(paste("sbp_exposed_90_np_rn2", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_90_np_rn2",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 90,
                                            range.min = 14,
                                            range.num = 2,
                                            keep.exposed = TRUE,
                                            exposed.type = "np",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_90_np_rn2, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_90_np_rn2"){
  print(paste("cholhdl_ratio_exposed_90_np_rn2", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_90_np_rn2",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 90,
                                                               range.min = 14,
                                                               range.num = 2,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "np",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_90_np_rn2, na.rm = TRUE)
} else if (variable.name == "sbp_exposed_90_md_rn2"){
  print(paste("sbp_exposed_90_md_rn2", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_90_md_rn2",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 90,
                                            range.min = 14,
                                            range.num = 2,
                                            keep.exposed = TRUE,
                                            exposed.type = "md",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_90_md_rn2, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_90_md_rn2"){
  print(paste("cholhdl_ratio_exposed_90_md_rn2", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_90_md_rn2",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 90,
                                                               range.min = 14,
                                                               range.num = 2,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "md",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_90_md_rn2, na.rm = TRUE)
  ###
  ### Extract ALL exposed/unexposed values that meet the criteria (at least 1 prescription in the X days before, or no 
  ### prescriptions in the X days before), and extract observation date, number in window, min diff, etc. 
  ### Then we can change the time window for searching and exposure type,without having to re-extract data
  ###
  ### Start with unexposed data (no prescription in 180 days prior)
} else if (variable.name == "sbp_unexposed_all"){
  print(paste("sbp_unexposed_all", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_unexposed <- extract_test_unexposed(cohort = cohort.base, 
                                              varname = "sbp_unexposed_all",
                                              codelist.med = "edh_sbp_medcodeid",
                                              codelist.drug = "uom_antihypertensives_prodcodeid",
                                              query.med = query.med,
                                              query.drug = query.drug,
                                              task.id = task.id,
                                              task.id.n = 100,
                                              range.max = 180,
                                              range.min = 1,
                                              range.num = 1,
                                              keep.exposed = FALSE,
                                              exposed.type = "all",
                                              indexdt = indexdt.name, 
                                              lower.bound = 70, 
                                              upper.bound = 210,
                                              db = "aurum", 
                                              out.save.disk = TRUE, 
                                              out.filepath = NULL, 
                                              out.subdir = out.subdir, 
                                              return.output = TRUE)
  quantile(var_sbp_unexposed$sbp_unexposed_all, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_unexposed_all"){
  print(paste("cholhdl_ratio_unexposed_all", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_unexposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                                 varname = "cholhdl_ratio_unexposed_all",
                                                                 codelist.ratio = "muzambi_ratio_medcodeid",
                                                                 codelist.chol = "muzambi_cholesterol_medcodeid",
                                                                 codelist.hdl = "muzambi_hdl_medcodeid",
                                                                 codelist.drug = "uom_statins_prodcodeid",
                                                                 query.ratio = query.cholhdl,
                                                                 query.chol = query.chol,
                                                                 query.hdl = query.hdl,
                                                                 query.drug = query.drug,
                                                                 task.id = task.id,
                                                                 task.id.n = 100,
                                                                 range.max = 180,
                                                                 range.min = 1,
                                                                 range.num = 1,
                                                                 keep.exposed = FALSE,
                                                                 exposed.type = "all",
                                                                 indexdt = indexdt.name, 
                                                                 lower.bound = 1, 
                                                                 upper.bound = 12,
                                                                 db = "aurum", 
                                                                 out.save.disk = TRUE, 
                                                                 out.filepath = NULL, 
                                                                 out.subdir = out.subdir, 
                                                                 return.output = TRUE)
  quantile(var_cholhdl_ratio_unexposed$cholhdl_ratio_unexposed_all, na.rm = TRUE)
  ###
  ### Exposed data, 180-day window
  ###
} else if (variable.name == "sbp_exposed_all"){
  print(paste("sbp_exposed_all", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_all",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 180,
                                            range.min = 1,
                                            range.num = 1,
                                            keep.exposed = TRUE,
                                            exposed.type = "all",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_all, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_all"){
  print(paste("cholhdl_ratio_exposed_all", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_all",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 180,
                                                               range.min = 1,
                                                               range.num = 1,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "all",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_all, na.rm = TRUE)
  ###
  ### 90-day window
  ###
} else if (variable.name == "sbp_exposed_90_all"){
  print(paste("sbp_exposed_90_all", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_exposed_90_all",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 90,
                                            range.min = 1,
                                            range.num = 1,
                                            keep.exposed = TRUE,
                                            exposed.type = "all",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_exposed_90_all, na.rm = TRUE)
} else if (variable.name == "cholhdl_ratio_exposed_90_all"){
  print(paste("cholhdl_ratio_exposed_90_all", Sys.time()))
  # Read in queries then extract variable
  query.cholhdl <- readRDS("data/extraction/cohort_baseline/query_cholhdl.rds")
  query.chol <- readRDS("data/extraction/cohort_baseline/query_chol.rds")
  query.hdl <- readRDS("data/extraction/cohort_baseline/query_hdl.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_statins.rds")
  var_cholhdl_ratio_exposed <- extract_cholhdl_ratio_unexposed(cohort = cohort.base, 
                                                               varname = "cholhdl_ratio_exposed_90_all",
                                                               codelist.ratio = "muzambi_ratio_medcodeid",
                                                               codelist.chol = "muzambi_cholesterol_medcodeid",
                                                               codelist.hdl = "muzambi_hdl_medcodeid",
                                                               codelist.drug = "uom_statins_prodcodeid",
                                                               query.ratio = query.cholhdl,
                                                               query.chol = query.chol,
                                                               query.hdl = query.hdl,
                                                               query.drug = query.drug,
                                                               task.id = task.id,
                                                               task.id.n = 100,
                                                               range.max = 90,
                                                               range.min = 1,
                                                               range.num = 1,
                                                               keep.exposed = TRUE,
                                                               exposed.type = "all",
                                                               indexdt = indexdt.name, 
                                                               lower.bound = 1, 
                                                               upper.bound = 12,
                                                               db = "aurum", 
                                                               out.save.disk = TRUE, 
                                                               out.filepath = NULL, 
                                                               out.subdir = out.subdir, 
                                                               return.output = TRUE)
  quantile(var_cholhdl_ratio_exposed$cholhdl_ratio_exposed_90_all, na.rm = TRUE)
} else if (variable.name == "sbp_expanded_unexposed_all"){
  print(paste("sbp_expanded_unexposed_all", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_unexposed <- extract_test_unexposed(cohort = cohort.base, 
                                              varname = "sbp_expanded_unexposed_all",
                                              codelist.med = "edh_sbp_medcodeid",
                                              codelist.drug = "uom_antihypertensives_expanded_prodcodeid",
                                              query.med = query.med,
                                              query.drug = query.drug,
                                              task.id = task.id,
                                              task.id.n = 100,
                                              range.max = 180,
                                              range.min = 1,
                                              range.num = 1,
                                              keep.exposed = FALSE,
                                              exposed.type = "all",
                                              indexdt = indexdt.name, 
                                              lower.bound = 70, 
                                              upper.bound = 210,
                                              db = "aurum", 
                                              out.save.disk = TRUE, 
                                              out.filepath = NULL, 
                                              out.subdir = out.subdir, 
                                              return.output = TRUE)
  quantile(var_sbp_unexposed$sbp_expanded_unexposed_all, na.rm = TRUE)
} else if (variable.name == "sbp_expanded_exposed_90_all"){
  print(paste("sbp_expanded_exposed_90_all", Sys.time()))
  # Read in queries then extract variable
  query.med <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
  query.drug <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
  var_sbp_exposed <- extract_test_unexposed(cohort = cohort.base, 
                                            varname = "sbp_expanded_exposed_90_all",
                                            codelist.med = "edh_sbp_medcodeid",
                                            codelist.drug = "uom_antihypertensives_expanded_prodcodeid",
                                            query.med = query.med,
                                            query.drug = query.drug,
                                            task.id = task.id,
                                            task.id.n = 100,
                                            range.max = 90,
                                            range.min = 1,
                                            range.num = 1,
                                            keep.exposed = TRUE,
                                            exposed.type = "all",
                                            indexdt = indexdt.name, 
                                            lower.bound = 70, 
                                            upper.bound = 210,
                                            db = "aurum", 
                                            out.save.disk = TRUE, 
                                            out.filepath = NULL, 
                                            out.subdir = out.subdir, 
                                            return.output = TRUE)
  quantile(var_sbp_exposed$sbp_expanded_exposed_90_all, na.rm = TRUE)
} else if (variable.name == "nonhdl"){
  print(paste(variable.name, Sys.time()))
  extract_nonhdl(cohort = cohort.base, 
                 varname = "nonhdl",
                 codelist.chol = "muzambi_cholesterol_medcodeid",
                 codelist.hdl = "muzambi_hdl_medcodeid",
                 codelist.ldl = "muzambi_ldl_medcodeid",
                 indexdt = indexdt.name, 
                 lower.bound = 0.4, 
                 upper.bound = 20.7,
                 db = "aurum", 
                 out.save.disk = TRUE, 
                 out.filepath = NULL, 
                 out.subdir = out.subdir, 
                 return.output = FALSE)
  print("FINISHED")
} 