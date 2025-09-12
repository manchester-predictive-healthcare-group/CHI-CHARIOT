### This program will finalise the cohort
### 1) Apply exclusion criteria for prior CVD using secondary care data, using finalised code list
### 2) Update end of follow up to also be based on end of study period, March 2020, which is the point in time we have finalised study data for.

### Clear workspace
rm(list=ls())

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Read in cohort
cohort <- readRDS("data/extraction/cohort_exclu2.rds")

###
### Exclude individuals with history of CVD in secondary care
###

### Start by identifying which individuals have a history of CVD in primary care
print(paste("ra", Sys.time()))
cvd_hist <- extract_ho(cohort = cohort, 
                       varname = "prior_CVD_hes", 
                       codelist = "uom_exclusion_icd", 
                       indexdt = "fup_start", 
                       db = "aurum_linked", 
                       tab = "hes_primary",
                       out.save.disk = FALSE, 
                       return.output = TRUE)

### Remove prior_CVD_hes from cohort, from when we applied exclusion criteria with type 1 linked data
cohort <- cohort[,-c("prior_CVD_hes")]
str(cohort)

### Merge with cohort.base
cohort <- merge(cohort, cvd_hist, by.x = "patid", by.y = "patid", all.x = TRUE)
print("AFTER MERGE")
str(cohort)

### Run the exclusion
cohort <- dplyr::filter(cohort, prior_CVD_hes == 0)
print(paste("Cohort after exclusion of CVD in secondary care", Sys.time()))
str(cohort)

###
### Extract death dates and compare  ONS/CPRD
###

### Create connection to linked data database
mydb <- RSQLite::dbConnect(RSQLite::SQLite(), "data/sql/aurum_linked.sqlite")

### We now have new death data from ONS, I guess we want to compare this with CPRD data
ons_death <- RSQLite::dbGetQuery(mydb, 'SELECT * FROM death')
ons_death <- ons_death[,c("patid", "dod", "cause")]

### Merge with cohort
cohort <- merge(cohort, ons_death, by.x = "patid", by.y = "patid", all.x = TRUE)

###
### Compare the difference between ddate from CPRD and date from ONS
colnames(cohort)

### 2.2% only have death in ONS
sum(is.na(cohort$ddate_comb) & !is.na(cohort$dod))/nrow(cohort)

### 0.16% only have death in CPRD
sum(!is.na(cohort$ddate_comb) & is.na(cohort$dod))/nrow(cohort)

### When both are recorded, 74% are the same
sum(cohort$ddate_comb[!is.na(cohort$ddate_comb) & !is.na(cohort$dod)] == cohort$dod[!is.na(cohort$ddate_comb) & !is.na(cohort$dod)])/nrow(cohort[!is.na(cohort$ddate_comb) & !is.na(cohort$dod),])

### When both are recorded, 2.9% have CPRD before ONS
sum(cohort$ddate_comb[!is.na(cohort$ddate_comb) & !is.na(cohort$dod)] < cohort$dod[!is.na(cohort$ddate_comb) & !is.na(cohort$dod)])/nrow(cohort[!is.na(cohort$ddate_comb) & !is.na(cohort$dod),])

### When both are recorded, 22.8% have ONS before CPRD
sum(cohort$ddate_comb[!is.na(cohort$ddate_comb) & !is.na(cohort$dod)] > cohort$dod[!is.na(cohort$ddate_comb) & !is.na(cohort$dod)])/nrow(cohort[!is.na(cohort$ddate_comb) & !is.na(cohort$dod),])

### When both are recorded, the absolute difference is bigger than 30 for 1.8%
sum(abs(as.numeric(cohort$ddate_comb[!is.na(cohort$ddate_comb) & !is.na(cohort$dod)]) - cohort$dod[!is.na(cohort$ddate_comb) & !is.na(cohort$dod)]) > 30)/nrow(cohort[!is.na(cohort$ddate_comb) & !is.na(cohort$dod),])


###
### Set follow up end to be min of 01032020 (last follow up in HES), fup_end and dod
### Also define death date to be minimum of ddate_comb (from CPRD AURUM) and death date in ONS
### Also exclude individuals who are older than 85 at start of follow up
###

### Create new start and end of follow dates, and death
cohort <- dplyr::mutate(cohort,
                        ### Define follow up end
                        fup_end = pmin(fup_end, # combined death date
                                       as.Date(paste('01032020', sep = ""), format = "%d%m%Y"), # end of follow up in HES
                                       dod, # death date in ONS
                                       na.rm = T),
                        ### Create a combined death death from primary care and secondary care data
                        ddate_comb = pmin(ddate_comb, # combined death date from CPRD AURUM
                                       dod, # death date in ONS
                                       na.rm = T),
                        do86 = dob + round(365.25*86))

### Subset to those where fup_start < fup_end
cohort <- subset(cohort, fup_start < fup_end)
print(paste("number of patients after requiring >= 1 day valid follow up = ", nrow(cohort), sep = ""))
str(cohort)

### Subset to those where start of follow up is younger than 86
cohort <- subset(cohort, fup_start < do86)
str(cohort)

###
### Save the cohort with all variables from the patient file
saveRDS(cohort, "data/extraction/cohort_exclu3.rds")

print(paste("FINISHED", Sys.time()))