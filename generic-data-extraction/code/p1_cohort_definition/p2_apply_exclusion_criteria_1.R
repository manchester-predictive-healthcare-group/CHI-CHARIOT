### Clear workspace
rm(list=ls())

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

###
### Extract the relevant data

### Extract filepaths for patient and practice
filenames.prac <- list.files("data/unzip/Practice", pattern = ".txt", full.names = TRUE)
filenames.pat <- list.files("data/unzip/Patient", pattern = ".txt", full.names = TRUE)
paste("extract", Sys.time())

### Extract patients from each practice and combine
pat.comb <- do.call("rbind", lapply(filenames.pat, extract_txt_pat, set = TRUE))
# pat.comb <- do.call("rbind", lapply(filenames.pat, extract_txt_pat, nrows = 500))
str(pat.comb)

### Extract practice data
prac.comb <- do.call("rbind", lapply(filenames.prac, extract_txt))

### Note that there are duplicates in the pracid files, i.e.
tail(extract_txt(filenames.prac[1]))
head(extract_txt(filenames.prac[2]))

### Therefore we should deduplicate
prac.comb <- subset(prac.comb, !duplicated(prac.comb$pracid))

### Merge these
paste("merge with practice file", Sys.time())
print(paste(nrow(pat.comb), "rows prior to merge with practice file"))
# cohort <- dplyr::inner_join(pat.comb, prac.comb, by = dplyr::join_by("pracid"))
cohort <- merge(pat.comb, prac.comb, by.x = "pracid", by.y = "pracid")
print(paste(nrow(cohort), "rows after merge with practice file"))

### Convert dates
cohort$regstartdate <- as.Date(cohort$regstartdate, format = "%d/%m/%Y")
cohort$regenddate <- as.Date(cohort$regenddate, format = "%d/%m/%Y")
cohort$emis_ddate <- as.Date(cohort$emis_ddate, format = "%d/%m/%Y")
cohort$cprd_ddate <- as.Date(cohort$cprd_ddate, format = "%d/%m/%Y")
cohort$lcd <- as.Date(cohort$lcd, format = "%d/%m/%Y")

### Create a date turned 18
cohort <- dplyr::mutate(cohort,
                        dob = as.Date(paste('0107',yob, sep = ""), format = "%d%m%Y"),
                        do18 = dob + round(365.25*18))

### Create first day of followup
cohort <- dplyr::mutate(cohort,
                        fup_start = pmax(do18, #date turned 18
                                         as.Date(paste('01012005', sep = ""), format = "%d%m%Y"), # start of study
                                         regstartdate + 365, # 1 year prior follow up
                                         na.rm = T),
                        ddate_comb = pmin(emis_ddate, cprd_ddate, na.rm = T), # create combined death date
                        fup_end = pmin(ddate_comb, # combined death date
                                       lcd, # last collection date of practice
                                       regenddate, # date end of registration with practice
                                       na.rm = T)
)

###
### Apply exclusion criteria for not having >= 1 valid days follow up in the study 

### Remove those without acceptable uts

## Check if any NA's
sum(is.na(cohort$acceptable))
## Remove individuals
cohort <- subset(cohort, acceptable == 1)
print(paste("number of patients with acceptable follow up = ", nrow(cohort), sep = ""))

### Remove individuals where they don't have >= 1 valid days follow up in the study 
### (i.e. exit criteria for cohort occurs prior to entry criteria are all satisfied)
cohort <- subset(cohort, fup_start < fup_end)
print(paste("number of patients after requiring >= 1 day valid follow up = ", nrow(cohort), sep = ""))

###
### Merge with linkage eligibility file and exclude those patients not with HES, ONS and lsoa linkage

### Extract lnikage file
linkage <- read.table("data/unzip/Linkage22/Aurum_enhanced_eligibility_January_2022.txt", sep = "\t", header = TRUE,
                      colClasses = c("character", "integer", "character", "integer", "integer", "integer",
                                     "integer", "integer", "integer","integer", "integer", "integer",
                                     "integer", "integer", "integer","integer"))
linkage <- linkage[,c("patid", "hes_apc_e", "ons_death_e", "lsoa_e")]
str(linkage)

### Apply merge
paste("merge linkage", Sys.time())
print(paste(nrow(cohort), "rows prior to merge with linkage"))
# cohort <- dplyr::inner_join(cohort, linkage, by = dplyr::join_by("patid"))
cohort <- merge(cohort, linkage, by.x = "patid", by.y = "patid")
print(paste(nrow(cohort), "rows after merging with linkage"))
str(cohort)


### Remove individuals not eligibile for linkage
cohort <- subset(cohort, hes_apc_e == 1 & ons_death_e == 1 & lsoa_e == 1)
print(paste(nrow(cohort), "rows after excluding those not eligible for linkage"))


###
### Run some sense checks on the cohort

### See if any people meeting the following criteria are in the dataset (sense checks):

### Died before 1st Jan 2005
sum(cohort$ddate_comb < as.Date(paste('01012005', sep = ""), format = "%d%m%Y"), na.rm = T)

### De-registered with practice before Jan 2005
sum(cohort$regenddate < as.Date(paste('01012005', sep = ""), format = "%d%m%Y"), na.rm = T)

### Practice last upload date before Jan 2005
sum(cohort$lcd < as.Date(paste('01012005', sep = ""), format = "%d%m%Y"), na.rm = T)

### < 18 years old at end of follow up
sum(cohort$do18 > cohort$regenddate, na.rm = T)


###
### Explore differences between cprd_ddate and emis_ddate
n <- nrow(cohort)
print(paste("number of patients = ", n, sep = ""))

## How many individuals have NA for emis_ddate and cprd_ddate
n1 <- sum(is.na(cohort$emis_ddate) & is.na(cohort$cprd_ddate) == TRUE)
print(paste(n1, "(", round(100*n1/n), "%) patients have missing emis_ddate and cprd_ddate", sep = ""))

## How many individuals have NA for emis_ddate and not cprd_ddate
n1 <- sum(is.na(cohort$emis_ddate) & !is.na(cohort$cprd_ddate) == TRUE)
print(paste(n1, "(", round(100*n1/n), "%) patients have missing emis_ddate and recorded cprd_ddate", sep = ""))

## How many individuals have NA for cprd_ddate and not emis_ddate
n1 <- sum(!is.na(cohort$emis_ddate) & is.na(cohort$cprd_ddate) == TRUE)
print(paste(n1, "(", round(100*n1/n), "%) patients have recorded emis_ddate and missing cprd_ddate", sep = ""))

## How many individuals have recorded for cprd_ddate and emis_ddate
n1 <- sum(!is.na(cohort$emis_ddate) & !is.na(cohort$cprd_ddate) == TRUE)
print(paste(n1, "(", round(100*n1/n), "%) patients have recorded emis_ddate and cprd_ddate", sep = ""))

## Of individuals with both recorded, how many are different
n2 <- sum(!is.na(cohort$emis_ddate) & !is.na(cohort$cprd_ddate) & (cohort$emis_ddate != cohort$cprd_ddate) == TRUE)
print(paste(n2, "(", round(100*n2/n1), "%) of individuals with both recorded, have different values recorded", sep = ""))

## This is the distribution of the differences
diff <- (cohort$emis_ddate - cohort$cprd_ddate)[!is.na(cohort$emis_ddate) & !is.na(cohort$cprd_ddate) & (cohort$emis_ddate != cohort$cprd_ddate)]
print("emis_ddate - cprd_ddate")
min(as.numeric(diff))
max(as.numeric(diff))
quantile(as.numeric(diff))
# hist(as.numeric(diff))

## Create a reduced cohort, which is just patids and the variable set
cohort.reduced <- cohort[,c("patid", "set")]

## Save these cohorts
print(Sys.time())
saveRDS(cohort, "data/extraction/cohort_exclu1.rds")
saveRDS(cohort.reduced, "data/extraction/cohort_exclu1_reduced.rds")
print("FINISHED")
