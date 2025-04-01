###
### This file will check extract_ho and the manual approach will give the same results
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project3/")
getwd()

### Read in extract_ho output
rcprd_ho <- readRDS("data/extract/ho_rcprd.rds")
manual_ho <- readRDS("data/extract/ho_manual.rds")

### Test they they are equal
testthat::expect_equal(rcprd_ho, manual_ho)

### Return output
mean(rcprd_ho$ho)*100
mean(manual_ho$ho)*100