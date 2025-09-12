### This program will create the codelists for the outcome (both primary care and secondary care), and the exclusion criteria (secondary care).
### The exclusion criteria is history of CVD or intracerebral stroke.
### The endeanvour health file with all the code lists, has a code list for the combination of CVD and intracerebral stroke, we
### therefore do not have to create a code list for the exclusion criteria in primary care. This is done in p0.1_create_codelists.

### Clear workspace
rm(list=ls())

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
getwd()

###
### Primary care, outcome
###

### Read in code lists
prim_chd <- data.table::fread(file = paste(getwd(),"/codelists/endeavour_health/edh_chd_emis.csv", sep = ""), 
                              sep = ",", header = TRUE, fill = TRUE, colClasses = "character", 
                              check.names = TRUE, na.strings = c(NA_character_, ""))

prim_stroke_tia <- data.table::fread(file = paste(getwd(),"/codelists/endeavour_health/edh_istroke_tia_emis.csv", sep = ""), 
                              sep = ",", header = TRUE, fill = TRUE, colClasses = "character", 
                              check.names = TRUE, na.strings = c(NA_character_, ""))

### Combine
prim_cvd <- rbind(prim_chd, prim_stroke_tia)
str(prim_cvd)
### Remove NA's, deduplicate and reduce to variable sof interst
prim_cvd <- base::subset(prim_cvd, !is.na(code.id)) |>
  dplyr::distinct(code.id, .keep_all = TRUE)

### Remove excess variable names
prim_cvd <- prim_cvd[,c("term", "subset", "code.id")]
colnames(prim_cvd)[colnames(prim_cvd) == "subset"] <- "condition"
colnames(prim_cvd)[colnames(prim_cvd) == "code.id"] <- "medcodeid"

### Merge with aurum.dictionary
## Get dictionary
aurum.med.dict <- data.table::fread(file = paste(getwd(),"/codelists/202106_emismedicaldictionary.txt", sep = ""), 
                                    sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                                    check.names = TRUE, na.strings = c(NA_character_, ""))
## Reduce to just the medcodeid's
aurum.med.dict <- aurum.med.dict[,c(1)]
colnames(aurum.med.dict) <- c("medcodeid")

## Merge
prim_cvd <- merge(prim_cvd, aurum.med.dict, by.x = "medcodeid", by.y = "medcodeid")

### Save to disk
write.csv(prim_cvd, "codelists/analysis/edh_cvd_medcodeid.csv", row.names= FALSE)

###
### Primary care, exclusion
###

###### CVD exclusion code list was created with "Q code group Cardiovascular disease (original)", from program p0.1
 

###
### Secondary care, outcome and exclusion
###

### These code lists are stored in codelists/icd, and codelists/analysis, and do not need any editing.
### See exclusion_icd.csv and outcome_cvd_icd.csv.





