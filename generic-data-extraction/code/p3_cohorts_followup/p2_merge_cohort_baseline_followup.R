###
### Merge variables extracted at a specific time (1-5 years) after index date into a single dataset
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

###
### Load the cohort
### 
cohort <- readRDS("data/extraction/cohort_exclu3.rds")
cohort <- data.table::as.data.table(cohort)

### Initialize empty list
variable.list <- vector("list", 0)

### Extract time point from command line
args <- commandArgs(trailingOnly = T)
t_fup <- round(365.25*as.numeric(args[1]))

### Put all variables and the cohort into a list object
variable.list[["cohort"]] <- dplyr::select(cohort, patid, pracid, gender, fup_start, fup_end)
for (variable in c("age"
                   ,"hypertension"
                   ,"ra"
                   ,"af"
                   ,"ckd"
                   ,"smi"
                   ,"fhcvd"
                   ,"migraine"
                   ,"sle"#10
                   ,"diabetes" 
                   ,"impotence"
                   ,"cortico"
                   ,"antipsy"
                   ,"smoking"
                   ,"statins"
                   ,"antihypertensives"
                   ,"copd"
                   ,"int_dis"
                   ,"downs"
                   ,"oral_cancer"
                   ,"brain_cancer"
                   ,"lung_cancer"
                   ,"blood_cancer"
                   ,"pre_eclampsia"
                   ,"postnatal_depression"
                   ,"bmi"
                   ,"sbp"
                   ,"nonhdl"
)){
  variable.list[[variable]] <- readRDS(paste(getwd(), "/data/extraction/cohort_baseline_followup/var_", variable, "_t", t_fup, ".rds", sep = ""))
}

### Check sizes of datasets
for (i in 1:length(variable.list)){
  print(i)
  str(variable.list[[i]])
}

### Merge the vaiables and cohort
cohort.var <- Reduce(function(x, y) merge(x, y, by.x = "patid", by.y = "patid"), variable.list)
str(cohort.var)

### Turn appropriate variables into factors
# cohort.var <- dplyr::mutate(cohort.var,
#                             hypertension = factor(hypertension, labels = c("Absent", "Present")),
#                             ra = factor(ra, labels = c("Absent", "Present")),
#                             af = factor(af, labels = c("Absent", "Present")),
#                             ckd = factor(ckd, labels = c("Absent", "Present")),
#                             smi = factor(smi, labels = c("Absent", "Present")),
#                             fhcvd = factor(fhcvd, labels = c("Absent", "Present")),
#                             migraine = factor(migraine, labels = c("Absent", "Present")),
#                             sle = factor(sle, labels = c("Absent", "Present")),
#                             impotence = factor(impotence, labels = c("Absent", "Present")),
#                             cortico = factor(cortico, labels = c("Absent", "Present")),
#                             antipsy = factor(antipsy, labels = c("Absent", "Present")),
#                             statins = factor(statins, labels = c("Absent", "Present")),
#                             antihypertensives = factor(antihypertensives, labels = c("Absent", "Present")),
#                             copd = factor(copd, labels = c("Absent", "Present")),
#                             int_dis = factor(int_dis, labels = c("Absent", "Present")),
#                             downs = factor(downs, labels = c("Absent", "Present")),
#                             oral_cancer = factor(oral_cancer, labels = c("Absent", "Present")),
#                             lung_cancer = factor(lung_cancer, labels = c("Absent", "Present")),
#                             blood_cancer = factor(blood_cancer, labels = c("Absent", "Present")),
#                             brain_cancer = factor(brain_cancer, labels = c("Absent", "Present")),
#                             pre_eclampsia = factor(pre_eclampsia, labels = c("Absent", "Present")),
#                             postnatal_depression = factor(postnatal_depression, labels = c("Absent", "Present")))
print("STR AFTER MERGE")
str(cohort.var)

### Remove variables that we don't want in final dataset (e.g. ethnicity derived only from AURUM or HES)
# cohort.var <- dplyr::select(cohort.var, -c(ethnicity.prim, ethnicity.hes))

### Save the cohort
saveRDS(cohort.var, paste("data/extraction/cohort_baseline_followup/cohort_var", t_fup, ".rds", sep = ""))
print(paste("IMAGE SAVED", Sys.time()))