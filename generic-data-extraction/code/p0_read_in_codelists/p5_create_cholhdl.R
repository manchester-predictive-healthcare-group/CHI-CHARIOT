###
### Need a seperate program to create the code list for cholesterol/HDL ratio, as one of the key codes was missing.
### The "term" for the code: Serum cholesterol/high density lipoprotein ratio was there, but the medcodeid was a code
### with only one observation. This "term" was in the codebrowser twice, the other code, which wasn't in the EDH had 40,000,000 observationcausing 
### there to be too much missing data in the variable. Note, both the codes had snomed concept ID 1015681000000109.
###

### The same is true for HDL, although the missing codes are not as impactful, and arguably were not included deliberately,
### whereas the omission of the code for Serum cholesterol/high density lipoprotein ratio was important.

### Clear workspace
rm(list=ls())

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
getwd()

###
### Create new code list for chol/hdl ratio
###

### Read in code list
edh_cholhdl_ratio <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_cholhdl_ratio_medcodeid.csv", sep = ""), 
                                sep = ",", header = TRUE, fill = TRUE, colClasses = "character", 
                                check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::select(medcodeid, term)

### Read in the extra code list I created manually using the code browser, and search *cholesterol*ratio*
uom_cholhdl_ratio <- data.table::fread(file = paste(getwd(),"/codelists/uom/cholhdl_ratio_codebrowser.txt", sep = ""), 
                                       sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                                       check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::rename(medcodeid = MedCodeId, term = Term) |>
  dplyr::select(medcodeid, term)

### Concatenate
cholhdl_ratio <- rbind(edh_cholhdl_ratio, uom_cholhdl_ratio)
### De-duplicate
cholhdl_ratio <- cholhdl_ratio[!duplicated(cholhdl_ratio$medcodeid), ]
### Write to disk
write.csv(cholhdl_ratio, paste("codelists/analysis/uom_cholhdl_ratio_medcodeid.csv", sep = ""), row.names= FALSE)


###
### Create new code list for hdl
###

### Read in code list
edh_hdl <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_hdl_medcodeid.csv", sep = ""), 
                                       sep = ",", header = TRUE, fill = TRUE, colClasses = "character", 
                                       check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::select(medcodeid, term)

### Read in the extra code list I created manually using the code browser, and search *cholesterol*ratio*
uom_hdl <- data.table::fread(file = paste(getwd(),"/codelists/uom/hdl_codebrowser.txt", sep = ""), 
                                       sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                                       check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::rename(medcodeid = MedCodeId, term = Term) |>
  dplyr::select(medcodeid, term)

### Concatenate
hdl <- rbind(edh_hdl, uom_hdl)
### De-duplicate
hdl <- hdl[!duplicated(hdl$medcodeid), ]
### Write to disk
write.csv(hdl, paste("codelists/analysis/uom_hdl_medcodeid.csv", sep = ""), row.names= FALSE)

###
### Create new code list for chol
###

### Read in code list
edh_chol <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_chol_medcodeid.csv", sep = ""), 
                             sep = ",", header = TRUE, fill = TRUE, colClasses = "character", 
                             check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::select(medcodeid, term)

### Read in the extra code list I created manually using the code browser, and search *cholesterol*ratio*
uom_chol <- data.table::fread(file = paste(getwd(),"/codelists/uom/chol_codebrowser.txt", sep = ""), 
                             sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                             check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::rename(medcodeid = MedCodeId, term = Term) |>
  dplyr::select(medcodeid, term)

### Concatenate
chol <- rbind(edh_chol, uom_chol)
### De-duplicate
chol <- chol[!duplicated(chol$medcodeid), ]
### Write to disk
write.csv(chol, paste("codelists/analysis/uom_chol_medcodeid.csv", sep = ""), row.names= FALSE)


###
### Finally going to check for any code lists I may have missed
###

### Read in code lists i am using (which have been derived above)
anal.chol <- read.csv(paste("codelists/analysis/uom_chol_medcodeid.csv", sep = ""), colClasses = c("character", "character"))
anal.hdl <- read.csv(paste("codelists/analysis/uom_hdl_medcodeid.csv", sep = ""), colClasses = c("character", "character"))
anal.cholhdl <- read.csv(paste("codelists/analysis/uom_cholhdl_ratio_medcodeid.csv", sep = ""), colClasses = c("character", "character"))

### Concatenate
anal.all <- rbind(anal.chol, anal.hdl, anal.cholhdl)

### Read in file from codebrowser to see what I have missed
uom.all <- data.table::fread(file = paste(getwd(),"/codelists/uom/all_chol_codebrowser.txt", sep = ""), 
                             sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                             check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::rename(medcodeid = MedCodeId, term = Term) |>
  dplyr::select(medcodeid, term)
str(uom.all)
str(anal.all)

### Lets do an anti_join
just.uom <- dplyr::anti_join(uom.all, anal.all, by = dplyr::join_by(medcodeid))
just.anal <- dplyr::anti_join(anal.all, uom.all, by = dplyr::join_by(medcodeid))

### One partiular code I am interseted in
### 279462018, High density lipoprotein/total cholesterol ratio
### I have checked and they look like normal ratios, so going to add this to a new code list
anal.cholhdl.extended <- rbind(anal.cholhdl, subset(just.uom, medcodeid == "279462018"))

### Write to disk
write.csv(subset(anal.cholhdl.extended, medcodeid == "279462018"), paste("codelists/analysis/uom_hdlchol_ratio_medcodeid.csv", sep = ""), row.names= FALSE)
