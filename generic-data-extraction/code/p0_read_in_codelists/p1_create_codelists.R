### This reads in the .csv files from the endeavour health site, and gets into a usable format for extracting CPRD data
### There are three files containing ICD, medcodeid and prodcodeid for all the conditions/variables usd in qrisk3
### Some conditions will not have any codes of a certain type (e.g. medical hstory variables will not have prodcodeid's)

### This program will create codelists for the predictors variables, which are all contained within one file.
### Code lists for the outcome will be created seperately.

### Clear workspace
rm(list=ls())

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
getwd()

## Read in master code list file
qrisk3_all <- data.table::fread(file = paste(getwd(),"/codelists/endeavour_health/edh_qrisk3_all_emis.csv", sep = ""), 
                                 sep = ",", header = TRUE, fill = TRUE, colClasses = "character", 
                                 check.names = TRUE, na.strings = c(NA_character_, ""))

### Get a list of all subset names/variable categories
varnames <- names(table(qrisk3_all$subset))
varnames

### Also read in the Aurum dictionary, some of the EMIS original codes aren't in the Aurum medical dictionary, and we will remove these

### Medical dictionary
## Get dictionary
aurum.med.dict <- data.table::fread(file = paste(getwd(),"/codelists/202106_emismedicaldictionary.txt", sep = ""), 
                               sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                               check.names = TRUE, na.strings = c(NA_character_, ""))
## Reduce to just the medcodeid's
aurum.med.dict <- aurum.med.dict[,c(1)]
colnames(aurum.med.dict) <- c("medcodeid")

### Product dictionary
## Get dictionary
aurum.prod.dict <- data.table::fread(file = paste(getwd(),"/codelists/202106_emisproductdictionary.txt", sep = ""), 
                                   sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                                   check.names = TRUE, na.strings = c(NA_character_, ""))

## Reduce to just the prodcodeid's
aurum.prod.dict <- aurum.prod.dict[,c(1)]
colnames(aurum.prod.dict) <- c("prodcodeid")


### ICD dictionary
## Get dictionary
icd.dict <- data.table::fread(file = paste(getwd(),"/codelists/icd_dict.csv", sep = ""), 
                                     sep = ",", header = TRUE, fill = TRUE, colClasses = "character", 
                                     check.names = TRUE, na.strings = c(NA_character_, ""))
## Reduce to relevant observations
icd.dict <- icd.dict[,c("CODE", "DESCRIPTION")]

### Function to extract a codelist for a specific variable from the parent codelist
create_codelist_edh <- function(outname, 
                                subset_condition, 
                                type, 
                                save.disk = TRUE, 
                                return.output = TRUE){
  
#   outname = "edh_cvd_hist"
#   subset_condition = "Q code group Cardiovascular disease (original)"
#   type = "icd"
#   outname = "edh_statins"
#   subset_condition = "Q code group Statins"
#   type = "prodcodeid"
#   type = "icd"
#   save.disk = FALSE
  
  ### Load in parent code list

  ### Read in master code list files
  if (type %in% c("medcodeid", "prodcodeid")){
    parent_list <- data.table::fread(file = paste(getwd(),"/codelists/endeavour_health/edh_qrisk3_all_emis.csv", sep = ""), 
                                     sep = ",", header = TRUE, fill = TRUE, colClasses = "character", 
                                     check.names = TRUE, na.strings = c(NA_character_, ""))
  } else if (type == "icd"){
    parent_list <- data.table::fread(file = paste(getwd(),"/codelists/endeavour_health/edh_qrisk3_all_icd.csv", sep = ""), 
                                     sep = ",", header = TRUE, fill = TRUE, colClasses = "character", 
                                     check.names = TRUE, na.strings = c(NA_character_, ""))
  }

  ### Check that subset_condition is there
  if(!(subset_condition %in% names(table(parent_list$subset)))){
    stop("subset_condition is not present in the parent code list")
  }
  
  ### If type = "ICD" remove 'code.id' then rename 'legacy.code' to 'code.id', to match that from the edh_qrisk3_all_emis.csv file
  if (type == "icd"){
    parent_list <- dplyr::select(parent_list, -code.id) |>
      dplyr::rename(code.id = legacy.code)
  }
  
  ### Extract codes that have the specified subset_condition and deduplicate
  codelist <- base::subset(parent_list, subset == subset_condition & !is.na(code.id)) |>
    dplyr::distinct(code.id, .keep_all = TRUE)
  
  ### Remove excess variable names
  codelist <- codelist[,c("term", "subset", "code.id")]
  colnames(codelist)[colnames(codelist) == "subset"] <- "condition"
  colnames(codelist)[colnames(codelist) == "code.id"] <- type
  
  ### Merge with the appropriate dictionary, and only keep codes that are in the Aurum medical/product dictionary
  if (type == "medcodeid"){
    codelist <- merge(codelist, aurum.med.dict, by.x = "medcodeid", by.y = "medcodeid")
  } else if (type == "prodcodeid"){
    codelist <- merge(codelist, aurum.prod.dict, by.x = "prodcodeid", by.y = "prodcodeid")
  } else if (type == "icd"){
    codelist <- merge(codelist, icd.dict, by.x = "icd", by.y = "CODE")
    codelist <- dplyr::select(codelist, -term) |> 
      dplyr::rename(term = DESCRIPTION) |> 
      dplyr::select(icd, term, condition)
  }
  
  ### Arrange by term
  codelist <- dplyr::arrange(codelist, term)
  if (type == "icd"){
    codelist <- dplyr::arrange(codelist, icd)
  }
  
  ### Save this code list
  if (save.disk == TRUE){
    write.csv(codelist, paste("codelists/analysis/", outname, "_", type, ".csv", sep = ""), row.names= FALSE)
  }

  ### Return codelist
  if (return.output == TRUE){
    return(codelist)
  }
  
}

########################
### Create codelists ###
########################

### CVD exclusion code lists
create_codelist_edh(outname = "edh_exclusion", subset_condition = "Q code group Cardiovascular disease (original)", type = "medcodeid")

### Demographics ###
create_codelist_edh(outname = "edh_bangladeshi", subset_condition = "Q code group Ethnicity - Bangladeshi", type = "medcodeid")
create_codelist_edh(outname = "edh_black_african", subset_condition = "Q code group Ethnicity - Black African" , type = "medcodeid")
create_codelist_edh(outname = "edh_caribbean", subset_condition = "Q code group Ethnicity - Caribbean" , type = "medcodeid")
create_codelist_edh(outname = "edh_chinese", subset_condition = "Q code group Ethnicity - Chinese" , type = "medcodeid")
create_codelist_edh(outname = "edh_indian", subset_condition = "Q code group Ethnicity - indian" , type = "medcodeid")
create_codelist_edh(outname = "edh_irish", subset_condition = "Q code group Ethnicity - Irish" , type = "medcodeid")
create_codelist_edh(outname = "edh_not_recorded", subset_condition = "Q code group Ethnicity - not recorded" , type = "medcodeid")
create_codelist_edh(outname = "edh_oth_asian", subset_condition = "Q code group Ethnicity - Other Asian" , type = "medcodeid")
create_codelist_edh(outname = "edh_oth_black", subset_condition = "Q code group Ethnicity - other black"   , type = "medcodeid")
create_codelist_edh(outname = "edh_oth_ethnic", subset_condition = "Q code group Ethnicity - Other ethnic group"  , type = "medcodeid")
create_codelist_edh(outname = "edh_oth_mixed", subset_condition = "Q code group Ethnicity - other mixed" , type = "medcodeid")
create_codelist_edh(outname = "edh_oth_white", subset_condition = "Q code group Ethnicity - Other white background" , type = "medcodeid")
create_codelist_edh(outname = "edh_pakistani", subset_condition = "Q code group Ethnicity - Pakistani"  , type = "medcodeid")
create_codelist_edh(outname = "edh_white_asian", subset_condition = "Q code group Ethnicity - White & Asian"    , type = "medcodeid")
create_codelist_edh(outname = "edh_white_black_african", subset_condition = "Q code group Ethnicity - White & Black African"  , type = "medcodeid")
create_codelist_edh(outname = "edh_white_black_caribbean", subset_condition = "Q code group Ethnicity - White & Black Caribbean"  , type = "medcodeid")
create_codelist_edh(outname = "edh_white_british", subset_condition = "Q code group Ethnicity - White British (English, Welsh, Scottish, Northern Irish or British)", type = "medcodeid")

### Medical history ###
create_codelist_edh(outname = "edh_hypertension", subset_condition = "Q code group Hypertension ", type = "medcodeid")
create_codelist_edh(outname = "edh_ra", subset_condition = "Q code group Rheumatoid Arthritis", type = "medcodeid")
create_codelist_edh(outname = "edh_af", subset_condition = "Q code group Atrial Fibrillation or Flutter", type = "medcodeid")
create_codelist_edh(outname = "edh_smi", subset_condition = "Q code group Severe mental illness", type = "medcodeid")
create_codelist_edh(outname = "edh_fhcvd", subset_condition = "Q code group Family history coronary heart disease", type = "medcodeid")
create_codelist_edh(outname = "edh_migraine", subset_condition = "Q code group Migraine", type = "medcodeid")
create_codelist_edh(outname = "edh_sle", subset_condition = "Q code group Systemic lupus erythematosis", type = "medcodeid")
create_codelist_edh(outname = "edh_t1dia", subset_condition = "Q code group Type 1 diabetes ", type = "medcodeid")
create_codelist_edh(outname = "edh_t2dia", subset_condition = "Q code group Type 2 diabetes ", type = "medcodeid")
create_codelist_edh(outname = "edh_ckd", subset_condition = "Q code group Chronic renal disease", type = "medcodeid")
create_codelist_edh(outname = "edh_impotence", subset_condition = "Q code group erectile dysfunction/Impotence", type = "medcodeid")
create_codelist_edh(outname = "edh_smoking_non", subset_condition = "Q code group Non-smoker", type = "medcodeid")
create_codelist_edh(outname = "edh_smoking_ex", subset_condition = "Q code group Ex-smoker", type = "medcodeid")
create_codelist_edh(outname = "edh_smoking_light", subset_condition = "Q code group Light smoker < 10/day or amount not specified", type = "medcodeid")
create_codelist_edh(outname = "edh_smoking_mod", subset_condition = "Q code group Moderate smoker 10-19/day", type = "medcodeid")
create_codelist_edh(outname = "edh_smoking_heavy", subset_condition = "Q code group Heavy smoker 2-+/day", type = "medcodeid")

### Test data ###
create_codelist_edh(outname = "edh_cholhdl_ratio", subset_condition = "Q code group cholesterol/HDL ratio", type = "medcodeid")
create_codelist_edh(outname = "edh_hdl", subset_condition = "Q code group HDL cholesterol ", type = "medcodeid")
create_codelist_edh(outname = "edh_chol", subset_condition = "Q code group Total cholesterol level (value)", type = "medcodeid")
create_codelist_edh(outname = "edh_sbp", subset_condition = "Q code group systolic blood pressure", type = "medcodeid")
create_codelist_edh(outname = "edh_bmi", subset_condition = "Q code group body mass index", type = "medcodeid")

### Prescription data ###
# create_codelist_edh(outname = "edh_statins", subset_condition = "Q code group Statins", type = "prodcodeid")
# create_codelist_edh(outname = "edh_cortico", subset_condition = "Q code group Oral corticosteroids", type = "prodcodeid")
# create_codelist_edh(outname = "edh_antipsy", subset_condition = "Q code group Atypical Antipsychotics", type = "prodcodeid")
# create_codelist_edh(outname = "edh_impotence", subset_condition = "Q code group 7.4.5 Drugs for erectile dysfunction", type = "prodcodeid")
# create_codelist_edh(outname = "edh_thiazides", subset_condition = "Q code group 2.2.1       Thiazides And Related Diuretics", type = "prodcodeid")
# create_codelist_edh(outname = "edh_diuretics", subset_condition = "Q code group 2.2.8       Diuretics With Potassium", type = "prodcodeid")
# create_codelist_edh(outname = "edh_beta_adren", subset_condition = "Q code group 2.4         Beta-Adrenoceptor Blocking Drugs", type = "prodcodeid")
# create_codelist_edh(outname = "edh_angio1", subset_condition = "Q code group 2.5.5.1     Angiotensin-Converting Enzyme Inhibitors", type = "prodcodeid")
# create_codelist_edh(outname = "edh_angio2", subset_condition = "Q code group 2.5.5.2     Angiotensin-II Receptor Antagonists", type = "prodcodeid")
# create_codelist_edh(outname = "edh_calcium_block", subset_condition = "Q code group 2.6.2       Calcium Channel Blockers", type = "prodcodeid")

