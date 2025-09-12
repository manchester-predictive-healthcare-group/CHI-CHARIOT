### I convert code lists to snomed
### This is done using the medical dictionary from the data extract
### However there is another dictionary in the code browser.
### Want to check I'm not missing anything.

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)


### Read medical dictionary from the data extract
med.dict.data <- data.table::fread(file = paste(getwd(),"/codelists/202106_emismedicaldictionary.txt", sep = ""), 
                                    sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                                    check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::rename(medcodeid = MedCodeId)
colnames(med.dict.data)

### Read medical dictionary from the code browser
med.dict.browser <- data.table::fread(file = paste(getwd(),"/codelists/CPRDAurumMedical.txt", sep = ""), 
                                   sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                                   check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::rename(medcodeid = MedCodeId)
colnames(med.dict.data)

### Get the ones just in the browser
just.cb <- dplyr::anti_join(med.dict.browser, med.dict.data, by = dplyr::join_by(medcodeid))

### Its a very small numner of codes
head(med.dict.browser)

### Lets see if they actually appear in the data
just.cb.data <- db_query(codelist = NULL,
                                     db.open = NULL,
                                     db = NULL,
                                     db.filepath = "/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/data/sql/aurum.sqlite",
                                     tab = "obs",
                                     codelist.vector = just.cb$medcodeid)

### Return the structure of the data
print(str(just.cb.data))
### Check for some in med.dict.data
subset(med.dict.data, SnomedCTConceptId == "816177009")
