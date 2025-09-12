### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load the cohort so far
cohort <- readRDS("data/extraction/cohort_exclu1.rds")
cohort <- data.table::as.data.table(cohort)
str(cohort)

###
### Exclude history of CVD in primary care
###

### Start by identifying which individuals have a history of CVD in primary care
print(paste("ra", Sys.time()))
cvd_hist <- extract_ho(cohort = cohort, 
                       varname = "prior_CVD_prim", 
                       codelist = "edh_exclusion_medcodeid", 
                       indexdt = "fup_start", 
                       db = "aurum", 
                       tab = "obs",
                       out.save.disk = FALSE, 
                       return.output = TRUE)

### Merge with cohort.base
cohort <- merge(cohort, cvd_hist, by.x = "patid", by.y = "patid", all.x = TRUE)
print("AFTER MERGE")
str(cohort)

### Run the exclusion
cohort <- dplyr::filter(cohort, prior_CVD_prim == 0)
print(paste("Cohort after exclusion of CVD in primary care", Sys.time()))
str(cohort)

###
### Exclude history of CVD in secondary care
###

### Load the linkage file
linkage.type1 <- extract_txt("data/unzip/Linkage_type1/Type_1/22_002333_icd_aurum_hesapc.txt", colClasses = "character")
linkage.type1 <- data.table::as.data.table(linkage.type1)
str(linkage.type1)

### Load the ICD10 code list for which exclusion will be applied
cvd_hist_icd_type1_exclusion <- data.table::fread(file = paste(getwd(),"/codelists/analysis/cvd_hist_icd_type1_exclusion.csv", sep = ""), 
                                                  sep = ",", header = TRUE, colClasses = "character")
str(cvd_hist_icd_type1_exclusion)

### Only retain values in linkage.type1, that are in our exclusion code list
linkage.type1$exclude <- as.integer(!is.na(fastmatch::fmatch(linkage.type1$ICD, cvd_hist_icd_type1_exclusion$CODE)))
linkage.type1 <- linkage.type1[exclude == 1]
table(linkage.type1$ICD)

### Format eventdate into a numeric date value
linkage.type1$eventdate <- as.Date(linkage.type1$eventdate, format = "%d/%m/%Y")
str(linkage.type1)
testthat::expect_true(as.Date(as.numeric(linkage.type1$eventdate[1]), origin = "1970-01-01") == "2014-05-30")

### Merge
paste("merge", Sys.time())
patids.exclude <- merge(cohort[,c("patid", "fup_start")], linkage.type1, by.x = "patid", by.y = "patid")|> 
  data.table::as.data.table()

### Subset data to observations that happened prior to fup_start
paste("subset observations that happened prior to fup_start", Sys.time())
patids.exclude <- patids.exclude[eventdate <= fup_start] %>%
  dplyr::pull(patid) %>% 
  base::unique()
length(patids.exclude)

### Add a variable to cohort, depending on if an individual is in patids.exclude
cohort$prior_CVD_hes <- as.integer(!is.na(fastmatch::fmatch(cohort$patid, patids.exclude)))

### Run the exclusion
cohort <- dplyr::filter(cohort, prior_CVD_hes == 0)
print(paste("Cohort after exclusion of CVD in HES", Sys.time()))
str(cohort)


###
### Save final list of patids

### Save the cohort with all variables from the patient file
saveRDS(cohort, "data/extraction/cohort_exclu2.rds")

### Save just the patient ids
saveRDS(cohort$patid, "data/extraction/cohort_exclu2_patid.rds")
write(cohort$patid, "data/extraction/cohort_exclu2_patid.txt")

### Create and save the file required to send to CPRD for type 1 linkage request
## Create
cohort.linkage.request <- cohort[,c("patid", "hes_apc_e", "ons_death_e", "lsoa_e")]
print("cohort for linkage request")
str(cohort.linkage.request)
## Save as tab delimited file with appropriate name
write.table(cohort.linkage.request, "data/extraction/22_002333_university_of_manchester_patientlist.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
print(paste("CREATED TYPE2 PATID FILE", Sys.time()))


### There is a limit of file size 20MB when requesting linkage. Therefore wan to split patients up into groups of size 1 million
cohort.linkage.request <- dplyr::arrange(cohort.linkage.request, patid)
cohort.linkage.request.list <- split(cohort.linkage.request, rep(1:21, each = 1000000)[1:nrow(cohort.linkage.request)])
str(cohort.linkage.request.list)

### Save each list element as a text file
lapply(1:21, function(x) {write.table(cohort.linkage.request.list[[x]], 
                                      paste("data/extraction/22_002333_university_of_manchester_patientlist_set", x, ".txt", sep = ""),
                                      sep = "\t", quote = FALSE, row.names = FALSE)})
print(paste("FINISHED", Sys.time()))
