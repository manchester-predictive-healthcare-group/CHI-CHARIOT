###
### Create Snomed codelists for Evergreen
###

### Clear workspace
rm(list=ls())

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
getwd()

### Medical dictionary
aurum.med.dict <- data.table::fread(file = paste(getwd(),"/codelists/202106_emismedicaldictionary.txt", sep = ""), 
                                    sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                                    check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::rename(medcodeid = MedCodeId)


aurum.prod.dict <- data.table::fread(file = paste(getwd(),"/codelists/202106_emisproductdictionary.txt", sep = ""), 
                                    sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                                    check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::rename(prodcodeid = ProdCodeId)


###################
### MEDCODEID's ###
###################

###
### Write a function which will convert the EDH codelists to Snomed
###
convert_edh <- function(outname){
  
  ### Read in an analysis code list
  codelist.medcodeid <- data.table::fread(paste("codelists/analysis/", outname, "_medcodeid.csv", sep = ""), colClasses = "character")
  
  ### Get Snomed code
  codelist.snomed <- merge(codelist.medcodeid, aurum.med.dict, by.x = "medcodeid", by.y = "medcodeid") |>
    dplyr::select(condition, Term, SnomedCTConceptId, SnomedCTDescriptionId)
  
  ### Save snomed codelist
  write.csv(codelist.snomed, paste("codelists/analysis_snomed/", outname, "_snomed.csv", sep = ""), row.names= FALSE)
  
}

### Convert edh codelists
###

### Exclusion Criteria
convert_edh(outname = "edh_exclusion")

### Outcome
convert_edh(outname = "edh_cvd")

### Demographics ###
convert_edh(outname = "edh_bangladeshi")
convert_edh(outname = "edh_black_african")
convert_edh(outname = "edh_caribbean")
convert_edh(outname = "edh_chinese")
convert_edh(outname = "edh_indian")
convert_edh(outname = "edh_irish")
convert_edh(outname = "edh_not_recorded")
convert_edh(outname = "edh_oth_asian")
convert_edh(outname = "edh_oth_black")
convert_edh(outname = "edh_oth_ethnic")
convert_edh(outname = "edh_oth_mixed")
convert_edh(outname = "edh_oth_white")
convert_edh(outname = "edh_pakistani")
convert_edh(outname = "edh_white_asian")
convert_edh(outname = "edh_white_black_african")
convert_edh(outname = "edh_white_black_caribbean")
convert_edh(outname = "edh_white_british")

### Medical history ###
convert_edh(outname = "edh_hypertension")
convert_edh(outname = "edh_ra")
convert_edh(outname = "edh_af")
convert_edh(outname = "edh_smi")
convert_edh(outname = "edh_fhcvd")
convert_edh(outname = "edh_migraine")
convert_edh(outname = "edh_sle")
convert_edh(outname = "edh_t1dia")
convert_edh(outname = "edh_t2dia")
convert_edh(outname = "edh_ckd")
convert_edh(outname = "edh_impotence")
convert_edh(outname = "edh_smoking_non")
convert_edh(outname = "edh_smoking_ex")
convert_edh(outname = "edh_smoking_light")
convert_edh(outname = "edh_smoking_mod")
convert_edh(outname = "edh_smoking_heavy")

### Test data ###
convert_edh(outname = "edh_sbp")
convert_edh(outname = "edh_bmi")

###
### Write function to convert the AH code lists
###

### These files already contain the Snomed IDs, so its just a case of reformatting so that they are in the same format
### as the edh ones
convert_ah <- function(outname){
  
  ### Read in an analysis code list
  codelist.snomed <- data.table::fread(paste("codelists/analysis/", outname, ".csv", sep = ""), colClasses = "character") |>
    dplyr::rename(condition = disease, Term = descr, SnomedCTConceptId = snomedctconceptid, SnomedCTDescriptionId = snomedctdescriptionid) |>
    dplyr::select(condition, Term, SnomedCTConceptId, SnomedCTDescriptionId)
  
  ### Save snomed codelist
  write.csv(codelist.snomed, paste("codelists/analysis_snomed/", outname, "_snomed.csv", sep = ""), row.names= FALSE)
  
}

### Convert edh codelists
###

convert_ah(outname = "ah_oral_cancer")
convert_ah(outname = "ah_brain_cancer")
convert_ah(outname = "ah_lung_cancer")
convert_ah(outname = "ah_blood_cancer")
convert_ah(outname = "ah_copd")
convert_ah(outname = "ah_downs_syndrome")
convert_ah(outname = "ah_intellectual_disability")

###
### Write function to convert the uom code lists
###

### This is same as the convert_edh function, but condition isn't in there, so we add that manually
convert_uom <- function(outname){
  
  ### Read in an analysis code list
  codelist.medcodeid <- data.table::fread(paste("codelists/analysis/", outname, "_medcodeid.csv", sep = ""), colClasses = "character")
  
  ### Get Snomed code
  codelist.snomed <- merge(codelist.medcodeid, aurum.med.dict, by.x = "medcodeid", by.y = "medcodeid") |>
    dplyr::mutate(condition = outname) |>
    dplyr::select(condition, Term, SnomedCTConceptId, SnomedCTDescriptionId)
  
  ### Save snomed codelist
  write.csv(codelist.snomed, paste("codelists/analysis_snomed/", outname, "_snomed.csv", sep = ""), row.names= FALSE)
  
}


### Convert uom codelists
###
convert_uom(outname = "uom_cholhdl_ratio")
convert_uom(outname = "uom_hdl")
convert_uom(outname = "uom_chol")
convert_uom(outname = "uom_pre_eclampsia")
convert_uom(outname = "uom_postnatal_depression")
convert_uom(outname = "weight")
convert_uom(outname = "height")

####################
### PRODCODEID's ###
####################

###
### Write function to convert the prodcodeid code lists
###

### This is same as the convert_edh function, but condition isn't in there, so we add that manually
convert_drug <- function(outname){

  ### Read in an analysis code list
  codelist.prodcodeid <- data.table::fread(paste("codelists/analysis/", outname, "_prodcodeid.csv", sep = ""), colClasses = "character")
  
  ### Get Snomed code
  codelist.out <- merge(dplyr::select(codelist.prodcodeid, prodcodeid), aurum.prod.dict, by.x = "prodcodeid", by.y = "prodcodeid") |>
    dplyr::mutate(condition = outname)
  
  ### Save snomed codelist
  write.csv(codelist.out, paste("codelists/analysis_snomed/", outname, "_drug.csv", sep = ""), row.names= FALSE)
  
}

### Convert edh codelists
###
convert_drug(outname = "uom_statins")
convert_drug(outname = "uom_antihypertensives")
convert_drug(outname = "uom_antihypertensives_expanded")
convert_drug(outname = "uom_antipsychotics")
convert_drug(outname = "uom_oral_corticosteroids")
convert_drug(outname = "uom_impotence")
convert_drug(outname = "uom_antipsychotics")
convert_drug(outname = "uom_antihypertensives_less_common")
