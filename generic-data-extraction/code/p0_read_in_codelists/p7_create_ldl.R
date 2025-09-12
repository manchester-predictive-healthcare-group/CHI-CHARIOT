###
### the LDL codelist is creaetd seperately.
### Wasn't intiially needed, but may be needed for prototype 3
###

### Clear workspace
rm(list=ls())

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
getwd()

###
### Create new code list for chol/hdl ratio
###


### Read in the extra code list I created manually using the code browser, and search *cholesterol*ratio*
uom_ldl <- data.table::fread(file = paste(getwd(),"/codelists/uom/ldl_codebrowser.txt", sep = ""), 
                                       sep = "\t", header = TRUE, fill = TRUE, colClasses = "character", 
                                       check.names = TRUE, na.strings = c(NA_character_, "")) |>
  dplyr::rename(medcodeid = MedCodeId, term = Term) |>
  dplyr::select(medcodeid, term)

### Write to disk
write.csv(uom_ldl, paste("codelists/analysis/uom_ldl_medcodeid.csv", sep = ""), row.names= FALSE)


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

### One particular code I am interseted in
### 279462018, High density lipoprotein/total cholesterol ratio
### I have checked and they look like normal ratios, so going to add this to a new code list
anal.cholhdl.extended <- rbind(anal.cholhdl, subset(just.uom, medcodeid == "279462018"))

### Write to disk
write.csv(subset(anal.cholhdl.extended, medcodeid == "279462018"), paste("codelists/analysis/uom_hdlchol_ratio_medcodeid.csv", sep = ""), row.names= FALSE)
