### Clear workspace
rm(list=ls())

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Start by creating/joining a database
mydb <- RSQLite::dbConnect(RSQLite::SQLite(), "data/sql/aurum_linked.sqlite")

### Get filename for HES_primary
filename <- "data/unzip/Linkage_type2/22_002333_Type2/Aurum_linked/Final/hes_primary_diag_hosp_22_002333_DM.txt"

### Add to database
add_to_database(filename,
                filetype = "hes_primary", 
                db = mydb,
                overwrite = TRUE)
paste("hes_primary added", Sys.time())

### Get filename for death data
filename <- "data/unzip/Linkage_type2/22_002333_Type2/Aurum_linked/Final/death_patient_22_002333_DM.txt"

### Add to database
add_to_database(filename,
                filetype = "death", 
                db = mydb,
                overwrite = TRUE)

### Get filename for patient data
filename <- "data/unzip/Linkage_type2/22_002333_Type2/Aurum_linked/Final/hes_patient_22_002333_DM.txt"

### Add to database
add_to_database(filename,
                filetype = "hes_patient", 
                db = mydb,
                overwrite = TRUE)

### Get filename for patient data
filename <- "data/unzip/Linkage_type2/22_002333_Type2/Aurum_linked/Final/hes_diagnosis_hosp_22_002333_DM.txt"

### Add to database
add_to_database(filename,
                filetype = "hes_hosp", 
                db = mydb,
                overwrite = TRUE)

### Get filename for patient data
filename <- "data/unzip/Linkage_type2/22_002333_Type2/Aurum_linked/Final/patient_2019_imd_22_002333_DM.txt"

### Add to database
add_to_database(filename,
                filetype = "IMD", 
                db = mydb,
                overwrite = TRUE)

RSQLite::dbListTables(mydb)

### Disconnect
RSQLite::dbDisconnect(mydb)

paste("FINISHED", Sys.time())
