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
)


### Extract scenario from command line
args <- commandArgs(trailingOnly = T)

cohort.base.name <- "cohort_exclu3"
cohort.base <- readRDS(paste("data/extraction/", cohort.base.name ,".rds", sep = ""))
cohort.base <- data.table::as.data.table(cohort.base)
indexdt.name <- "fup_start"
out.subdir <- "cohort_baseline_followup"
t_fup <- round(365.25*as.numeric(args[1]))

print(paste("cohort.base.name = ", cohort.base.name, sep = ""))
print(paste("indexdt.name = ", indexdt.name, sep = ""))
print(paste("out.subdir = ", out.subdir, sep = ""))
print(paste("t_fup = ", t_fup, sep = ""))

# ### Apply extraction
# print(paste("age", Sys.time()))
# var_age <- extract_age(cohort = cohort.base,
#                        indexdt = indexdt.name,
#                        t = t_fup,
#                        out.save.disk = TRUE,
#                        out.filepath = NULL,
#                        out.subdir = out.subdir,
#                        return.output = TRUE)
# quantile(var_age$age)
# 
# print(paste("hypertension", Sys.time()))
# var_hypertension <- extract_ho(cohort = cohort.base,
#                                varname = "hypertension",
#                                codelist = "edh_hypertension_medcodeid",
#                                indexdt = indexdt.name,
#                                t = t_fup,
#                                db = "aurum",
#                                tab = "obs",
#                                out.save.disk = TRUE,
#                                out.filepath = NULL,
#                                out.subdir = out.subdir,
#                                return.output = TRUE)
# table(var_hypertension$hypertension)
# 
# print(paste("ra", Sys.time()))
# var_ra <- extract_ho(cohort = cohort.base,
#                      varname = "ra",
#                      codelist = "edh_ra_medcodeid",
#                      indexdt = indexdt.name,
#                      t = t_fup,
#                      db = "aurum",
#                      tab = "obs",
#                      out.save.disk = TRUE,
#                      out.filepath = NULL,
#                      out.subdir = out.subdir,
#                      return.output = TRUE)
# table(var_ra$ra)
# 
# print(paste("af", Sys.time()))
# var_af <- extract_ho(cohort = cohort.base,
#                      varname = "af",
#                      codelist = "edh_af_medcodeid",
#                      indexdt = indexdt.name,
#                      t = t_fup,
#                      db = "aurum",
#                      tab = "obs",
#                      out.save.disk = TRUE,
#                      out.filepath = NULL,
#                      out.subdir = out.subdir,
#                      return.output = TRUE)
# table(var_af$af)
# 
# print(paste("ckd", Sys.time()))
# var_ckd <- extract_ho(cohort = cohort.base,
#                       varname = "ckd",
#                       codelist = "edh_ckd_medcodeid",
#                       indexdt = indexdt.name,
#                       t = t_fup,
#                       db = "aurum",
#                       tab = "obs",
#                       out.save.disk = TRUE,
#                       out.filepath = NULL,
#                       out.subdir = out.subdir,
#                       return.output = TRUE)
# table(var_ckd$ckd)
# 
# print(paste("smi", Sys.time()))
# var_smi <- extract_ho(cohort = cohort.base,
#                       varname = "smi",
#                       codelist = "edh_smi_medcodeid",
#                       indexdt = indexdt.name,
#                       t = t_fup,
#                       db = "aurum",
#                       tab = "obs",
#                       out.save.disk = TRUE,
#                       out.filepath = NULL,
#                       out.subdir = out.subdir,
#                       return.output = TRUE)
# table(var_smi$smi)
# 
# print(paste("fhcvd", Sys.time()))
# var_fhcvd <- extract_ho(cohort = cohort.base,
#                         varname = "fhcvd",
#                         codelist = "edh_fhcvd_medcodeid",
#                         indexdt = indexdt.name,
#                         t = t_fup,
#                         db = "aurum",
#                         tab = "obs",
#                         out.save.disk = TRUE,
#                         out.filepath = NULL,
#                         out.subdir = out.subdir,
#                         return.output = TRUE)
# table(var_fhcvd$fhcvd)
# 
# print(paste("migraine", Sys.time()))
# var_migraine <- extract_ho(cohort = cohort.base,
#                            varname = "migraine",
#                            codelist = "edh_migraine_medcodeid",
#                            indexdt = indexdt.name,
#                            t = t_fup,
#                            db = "aurum",
#                            tab = "obs",
#                            out.save.disk = TRUE,
#                            out.filepath = NULL,
#                            out.subdir = out.subdir,
#                            return.output = TRUE)
# table(var_migraine$migraine)
# 
# print(paste("sle", Sys.time()))
# var_sle <- extract_ho(cohort = cohort.base,
#                       varname = "sle",
#                       codelist = "edh_sle_medcodeid",
#                       indexdt = indexdt.name,
#                       t = t_fup,
#                       db = "aurum",
#                       tab = "obs",
#                       out.save.disk = TRUE,
#                       out.filepath = NULL,
#                       out.subdir = out.subdir,
#                       return.output = TRUE)
# table(var_sle$sle)
# 
# print(paste("diabetes", Sys.time()))
# var_diabetes <- extract_diabetes(cohort = cohort.base,
#                                  codelist.type1 = "edh_t1dia_medcodeid",
#                                  codelist.type2 = "edh_t2dia_medcodeid",
#                                  indexdt = indexdt.name,
#                                  t = t_fup,
#                                  db = "aurum",
#                                  out.save.disk = TRUE,
#                                  out.filepath = NULL,
#                                  out.subdir = out.subdir,
#                                  return.output = TRUE)
# table(var_diabetes$diabetes)
# 
# print(paste("impotence", Sys.time()))
# var_impotence <- extract_impotence(cohort = cohort.base,
#                                    codelist.med = "edh_impotence_medcodeid",
#                                    codelist.drug = "uom_impotence_prodcodeid",
#                                    indexdt = indexdt.name,
#                                    t = t_fup,
#                                    db = "aurum",
#                                    out.save.disk = TRUE,
#                                    out.filepath = NULL,
#                                    out.subdir = out.subdir,
#                                    return.output = TRUE)
# table(var_impotence$impotence)
#
# print(paste("cortico", Sys.time()))
# var_cortico <- extract_ho(varname = "cortico",
#                           codelist = "uom_oral_corticosteroids_prodcodeid",
#                           cohort = cohort.base,
#                           indexdt = indexdt.name,
#                           t = t_fup,
#                           time.prev = round(365.25*0.5),
#                           numobs = 2,
#                           db = "aurum",
#                           tab = "drug",
#                           out.save.disk = TRUE,
#                           out.filepath = NULL,
#                           out.subdir = out.subdir,
#                           return.output = TRUE)
# table(var_cortico$cortico)
# 
# print(paste("antipsy", Sys.time()))
# var_antipsy <- extract_ho(varname = "antipsy",
#                           codelist = "uom_antipsychotics_prodcodeid",
#                           cohort = cohort.base,
#                           indexdt = indexdt.name,
#                           t = t_fup,
#                           time.prev = round(365.25*0.5),
#                           numobs = 2,
#                           db = "aurum",
#                           tab = "drug",
#                           out.save.disk = TRUE,
#                           out.filepath = NULL,
#                           out.subdir = out.subdir,
#                           return.output = TRUE)
# table(var_antipsy$antipsy)
# 
# print(paste("bmi", Sys.time()))
# var_bmi <- extract_bmi(cohort = cohort.base,
#                        codelist.bmi = "edh_bmi_medcodeid",
#                        codelist.weight = "weight_medcodeid",
#                        codelist.height = "height_medcodeid",
#                        indexdt = indexdt.name,
#                        t = t_fup,
#                        lower.bound = 18,
#                        upper.bound = 47,
#                        db = "aurum",
#                        out.save.disk = TRUE,
#                        out.filepath = NULL,
#                        out.subdir = out.subdir,
#                        return.output = TRUE)
# quantile(var_bmi$bmi, na.rm = TRUE)
# 
# print(paste("sbp", Sys.time()))
# var_sbp <- extract_sbp(cohort = cohort.base,
#                        codelist = "edh_sbp_medcodeid",
#                        indexdt = indexdt.name,
#                        t = t_fup,
#                        lower.bound = 70,
#                        upper.bound = 210,
#                        db = "aurum",
#                        out.save.disk = TRUE,
#                        out.filepath = NULL,
#                        out.subdir = out.subdir,
#                        return.output = TRUE)
# quantile(var_sbp$sbp, na.rm = TRUE)
# 
# print(paste("sbp var", Sys.time()))
# var_sbp_var <- extract_sbp_var(cohort = cohort.base,
#                                codelist = "edh_sbp_medcodeid",
#                                indexdt = indexdt.name,
#                                t = t_fup,
#                                lower.bound = 70,
#                                upper.bound = 210,
#                                db = "aurum",
#                                out.save.disk = TRUE,
#                                out.filepath = NULL,
#                                out.subdir = out.subdir,
#                                return.output = TRUE)
# quantile(var_sbp_var$sbp_var, na.rm = TRUE)
# 
print(paste("smoking", Sys.time()))
var_smoking <- extract_smoking(cohort = cohort.base,
                               codelist.non = "edh_smoking_non_medcodeid",
                               codelist.ex = "edh_smoking_ex_medcodeid",
                               codelist.light = "edh_smoking_light_medcodeid",
                               codelist.mod = "edh_smoking_mod_medcodeid",
                               codelist.heavy = "edh_smoking_heavy_medcodeid",
                               indexdt = indexdt.name,
                               t = t_fup,
                               db = "aurum",
                               out.save.disk = TRUE,
                               out.filepath = NULL,
                               out.subdir = out.subdir,
                               return.output = TRUE)
table(var_smoking$smoking)
str(var_smoking)
# 
# print(paste("outcome", Sys.time()))
# var_time_until_cvd <- extract_time_until_cvd(cohort.base,
#                                              varname.time = NULL,
#                                              varname.indicator = NULL,
#                                              codelist.primary = "edh_cvd_medcodeid",
#                                              codelist.hes = "uom_cvd_icd",
#                                              indexdt = indexdt.name,
#                                              censdt = "fup_end",
#                                              t = t_fup,
#                                              t.varname = TRUE,
#                                              db.primary = "aurum",
#                                              db.hes = "aurum_linked",
#                                              db.filepath = NULL,
#                                              out.save.disk = TRUE,
#                                              out.filepath = NULL,
#                                              out.subdir = out.subdir,
#                                              return.output = TRUE)
# str(var_time_until_cvd)
# table(var_time_until_cvd$cvd_indicator_prim)
# table(var_time_until_cvd$cvd_indicator_hes)
# table(var_time_until_cvd$cvd_indicator_death)
# table(var_time_until_cvd$cvd_indicator)
# 
# print(paste("antihypertensives", Sys.time()))
# var_antihypertensives <- extract_ho(cohort = cohort.base,
#                                     varname = "antihypertensives",
#                                     codelist = "uom_antihypertensives_prodcodeid",
#                                     indexdt = indexdt.name,
#                                     t = t_fup,
#                                     time.prev = round(365.25*0.5),
#                                     numobs = 2,
#                                     db = "aurum",
#                                     tab = "drug",
#                                     out.save.disk = TRUE,
#                                     out.filepath = NULL,
#                                     out.subdir = out.subdir,
#                                     return.output = TRUE)
# table(var_antihypertensives$antihypertensives)
# 
# print(paste("statins", Sys.time()))
# var_statins <- extract_ho(cohort = cohort.base,
#                           varname = "statins",
#                           codelist = "uom_statins_prodcodeid",
#                           indexdt = indexdt.name,
#                           t = t_fup,
#                           time.prev = round(365.25*0.5),
#                           numobs = 2,
#                           db = "aurum",
#                           tab = "drug",
#                           out.save.disk = TRUE,
#                           out.filepath = NULL,
#                           out.subdir = out.subdir,
#                           return.output = TRUE)
# table(var_statins$statins)
# 
# print(paste("copd", Sys.time()))
# var_copd <- extract_ho(cohort = cohort.base,
#                        varname = "copd",
#                        codelist = "ah_copd",
#                        indexdt = indexdt.name,
#                        t = t_fup,
#                        db = "aurum",
#                        tab = "obs",
#                        out.save.disk = TRUE,
#                        out.filepath = NULL,
#                        out.subdir = out.subdir,
#                        return.output = TRUE)
# table(var_copd$copd)
# 
# print(paste("int_dis", Sys.time()))
# var_int_dis <- extract_ho(cohort = cohort.base,
#                           varname = "int_dis",
#                           codelist = "ah_intellectual_disability",
#                           indexdt = indexdt.name,
#                           t = t_fup,
#                           db = "aurum",
#                           tab = "obs",
#                           out.save.disk = TRUE,
#                           out.filepath = NULL,
#                           out.subdir = out.subdir,
#                           return.output = TRUE)
# table(var_int_dis$int_dis)
# 
# print(paste("downs", Sys.time()))
# var_downs <- extract_ho(cohort = cohort.base,
#                         varname = "downs",
#                         codelist = "ah_downs_syndrome",
#                         indexdt = indexdt.name,
#                         t = t_fup,
#                         db = "aurum",
#                         tab = "obs",
#                         out.save.disk = TRUE,
#                         out.filepath = NULL,
#                         out.subdir = out.subdir,
#                         return.output = TRUE)
# table(var_downs$downs)
# 
# print(paste("brain_cancer", Sys.time()))
# var_brain_cancer <- extract_ho(cohort = cohort.base,
#                                varname = "brain_cancer",
#                                codelist = "ah_brain_cancer",
#                                indexdt = indexdt.name,
#                                t = t_fup,
#                                db = "aurum",
#                                tab = "obs",
#                                out.save.disk = TRUE,
#                                out.filepath = NULL,
#                                out.subdir = out.subdir,
#                                return.output = TRUE)
# table(var_brain_cancer$brain_cancer)
# 
# print(paste("lung_cancer", Sys.time()))
# var_lung_cancer <- extract_ho(cohort = cohort.base,
#                               varname = "lung_cancer",
#                               codelist = "ah_lung_cancer",
#                               indexdt = indexdt.name,
#                               t = t_fup,
#                               db = "aurum",
#                               tab = "obs",
#                               out.save.disk = TRUE,
#                               out.filepath = NULL,
#                               out.subdir = out.subdir,
#                               return.output = TRUE)
# table(var_lung_cancer$lung_cancer)
# 
# print(paste("oral_cancer", Sys.time()))
# var_oral_cancer <- extract_ho(cohort = cohort.base,
#                               varname = "oral_cancer",
#                               codelist = "ah_oral_cancer",
#                               indexdt = indexdt.name,
#                               t = t_fup,
#                               db = "aurum",
#                               tab = "obs",
#                               out.save.disk = TRUE,
#                               out.filepath = NULL,
#                               out.subdir = out.subdir,
#                               return.output = TRUE)
# table(var_oral_cancer$oral_cancer)
# 
# print(paste("blood_cancer", Sys.time()))
# var_blood_cancer <- extract_ho(cohort = cohort.base,
#                                varname = "blood_cancer",
#                                codelist = "ah_blood_cancer",
#                                indexdt = indexdt.name,
#                                t = t_fup,
#                                db = "aurum",
#                                tab = "obs",
#                                out.save.disk = TRUE,
#                                out.filepath = NULL,
#                                out.subdir = out.subdir,
#                                return.output = TRUE)
# table(var_blood_cancer$blood_cancer)
# 
# print(paste("pre_eclampsia", Sys.time()))
# var_pre_eclampsia <- extract_ho(cohort = cohort.base,
#                                 varname = "pre_eclampsia",
#                                 codelist = "uom_pre_eclampsia_medcodeid",
#                                 indexdt = indexdt.name,
#                                 t = t_fup,
#                                 db = "aurum",
#                                 tab = "obs",
#                                 out.save.disk = TRUE,
#                                 out.filepath = NULL,
#                                 out.subdir = out.subdir,
#                                 return.output = TRUE)
# table(var_pre_eclampsia$pre_eclampsia)
# 
# print(paste("postnatal_depression", Sys.time()))
# var_postnatal_depression <- extract_postnatal_depression(cohort = cohort.base,
#                                                          varname = "postnatal_depression",
#                                                          codelist.med = "uom_postnatal_depression_medcodeid",
#                                                          codelist.test = "uom_postnatal_depression_score_medcodeid",
#                                                          indexdt = indexdt.name,
#                                                          t = t_fup,
#                                                          db = "aurum",
#                                                          out.save.disk = TRUE,
#                                                          out.filepath = NULL,
#                                                          out.subdir = out.subdir,
#                                                          return.output = TRUE)
# table(var_postnatal_depression$postnatal_depression)
# 
# print(paste("cholhdl_ratio", Sys.time()))
# var_cholhdl_ratio <- extract_cholhdl_ratio(cohort = cohort.base, 
#                                            codelist.ratio = "muzambi_cholesterol_medcodeid",
#                                            codelist.chol = "muzambi_hdl_medcodeid",
#                                            codelist.hdl = "muzambi_ldl_medcodeid",
#                                            indexdt = indexdt.name, 
#                                            t = t_fup,
#                                            lower.bound = 1,
#                                            upper.bound = 12,
#                                            db = "aurum",
#                                            out.save.disk = TRUE, 
#                                            out.filepath = NULL, 
#                                            out.subdir = out.subdir, 
#                                            return.output = TRUE)
# quantile(var_cholhdl_ratio$cholhdl_ratio, na.rm = TRUE)
# 
# print(paste("nonhdl", Sys.time()))
# var_nonhdl <- extract_nonhdl(cohort = cohort.base, 
#                              varname = "nonhdl",
#                              codelist.chol = "muzambi_cholesterol_medcodeid",
#                              codelist.hdl = "muzambi_hdl_medcodeid",
#                              codelist.ldl = "muzambi_ldl_medcodeid",
#                              indexdt = indexdt.name, 
#                              t = t_fup,
#                              lower.bound = 0.4, 
#                              upper.bound = 20.7,
#                              db = "aurum", 
#                              out.save.disk = TRUE, 
#                              out.filepath = NULL, 
#                              out.subdir = out.subdir, 
#                              return.output = FALSE)
# quantile(var_nonhdl$nonhdl, na.rm = TRUE)