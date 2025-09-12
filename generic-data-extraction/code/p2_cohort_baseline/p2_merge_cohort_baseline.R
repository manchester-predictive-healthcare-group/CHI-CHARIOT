### Create Table 1 of data

### This is going to test the functions for deriving all the variables!

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
### XXX TBC
cohort <- readRDS("data/extraction/cohort_exclu3.rds")
cohort <- data.table::as.data.table(cohort)

### Initialize empty list
variable.list <- vector("list", 0)

### Put all variables and the cohort into a list object
variable.list[["cohort"]] <- dplyr::select(cohort, patid, pracid, gender, fup_start, fup_end)
for (variable in c("age"
                   ,"ethnicity"
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
                   ,"time_until_cvd"
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
                   ,"sbp_var"
                   ,"cholhdl_ratio"
                   ,"nonhdl"
                   ,"bmi_alltime"
                   ,"sbp_alltime"
                   ,"sbp_var_alltime"
                   ,"cholhdl_ratio_alltime"
)){
  variable.list[[variable]] <- readRDS(paste(getwd(), "/data/extraction/cohort_baseline/var_", variable, ".rds", sep = ""))
}

### Check sizes of datasets
for (i in 1:length(variable.list)){
  print(i)
  str(variable.list[[i]])
}

### Merge the vaiables and cohort
cohort.var <- Reduce(function(x, y) merge(x, y, by.x = "patid", by.y = "patid"), variable.list)
str(cohort.var)

### Merge with IMD
### NB: Previously this was read directly in from the raw data, which has now been deleted.
### In p7_save_IMD.R I saved the IMD scores as an .rds object.
### When running this code subsequently, I read in the processed .rds file with the IMD scores,
### not the raw data (which is the commented out code)
###
# ## Get filename for IMD data
# filename <- "data/unzip/Linkage_type2/22_002333_Type2/Aurum_linked/Final/patient_2019_imd_22_002333.txt"
# 
# ## Read in IMD data
# imd.raw <- read.table(filename, sep = "\t", colClasses = c("character", "integer", "numeric"), header = TRUE)
# 
# ### Merge with cohort
# cohort.var <- merge(cohort.var, imd.raw[,c("patid", "e2019_imd_20")], by.x = "patid", by.y = "patid") |>
#   dplyr::rename(IMD = e2019_imd_20)

### Read in from saved file
imd.raw <- readRDS("data/extraction/cohort_baseline/cohort_IMD.rds")

### Merge with cohort
cohort.var <- merge(cohort.var, imd.raw, by.x = "patid", by.y = "patid")

### Turn appropriate variables into factors
cohort.var <- dplyr::mutate(cohort.var,
                            hypertension = factor(hypertension, labels = c("Absent", "Present")),
                            ra = factor(ra, labels = c("Absent", "Present")),
                            af = factor(af, labels = c("Absent", "Present")),
                            ckd = factor(ckd, labels = c("Absent", "Present")),
                            smi = factor(smi, labels = c("Absent", "Present")),
                            fhcvd = factor(fhcvd, labels = c("Absent", "Present")),
                            migraine = factor(migraine, labels = c("Absent", "Present")),
                            sle = factor(sle, labels = c("Absent", "Present")),
                            impotence = factor(impotence, labels = c("Absent", "Present")),
                            cortico = factor(cortico, labels = c("Absent", "Present")),
                            antipsy = factor(antipsy, labels = c("Absent", "Present")),
                            statins = factor(statins, labels = c("Absent", "Present")),
                            antihypertensives = factor(antihypertensives, labels = c("Absent", "Present")),
                            copd = factor(copd, labels = c("Absent", "Present")),
                            int_dis = factor(int_dis, labels = c("Absent", "Present")),
                            downs = factor(downs, labels = c("Absent", "Present")),
                            oral_cancer = factor(oral_cancer, labels = c("Absent", "Present")),
                            lung_cancer = factor(lung_cancer, labels = c("Absent", "Present")),
                            blood_cancer = factor(blood_cancer, labels = c("Absent", "Present")),
                            brain_cancer = factor(brain_cancer, labels = c("Absent", "Present")),
                            pre_eclampsia = factor(pre_eclampsia, labels = c("Absent", "Present")),
                            postnatal_depression = factor(postnatal_depression, labels = c("Absent", "Present")))
print("STR AFTER MERGE")
str(cohort.var)

### Remove variables that we don't want in final dataset (e.g. ethnicity derived only from AURUM or HES)
cohort.var <- dplyr::select(cohort.var, -c(ethnicity.prim, ethnicity.hes))

### Save the cohort
saveRDS(cohort.var, "data/extraction/cohort_baseline/cohort_var.rds")
print(paste("IMAGE SAVED", Sys.time()))

### Calculate incidence rates by year
print(paste("Incidence rates by year", Sys.time()))
source("code/p2_cohort_baseline/p3_incidence_rates_by_year.R")

### Locate pandoc
Sys.setenv(RSTUDIO_PANDOC="/opt/gridware/apps/binapps/rstudio/0.98.1103/bin/pandoc")

### Render the Rmarkdown document
print("render to markdown")
rmarkdown::render("code/p2_cohort_baseline/p3_table1_cohort_baseline.Rmd")

### Finished
print("FINISHED")