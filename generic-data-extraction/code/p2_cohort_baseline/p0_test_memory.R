### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Save as .rds objects
test1 <- readRDS("data/extraction/cohort_baseline/query_sbp.rds")
str(test1)
test2 <- readRDS("data/extraction/cohort_baseline/query_antihypertensives.rds")
str(test2)
print("FINISHED")

