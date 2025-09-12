###
### Explore extraction of smoking
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Smoking
smoking_var1 <- readRDS(paste(getwd(), "/data/extraction/cohort_baseline/var_", "smoking", ".rds", sep = ""))
smoking_var2 <- readRDS(paste(getwd(), "/data/extraction/cohort_baseline_followup/var_", "smoking_t365", ".rds", sep = ""))
table(smoking_var1$smoking)
table(smoking_var2$smoking)

sum(is.na(smoking_var1$smoking))
sum(is.na(smoking_var2$smoking))

# > table(smoking_var1$smoking)

# CHECK
# Non-smoker  Ex-smoker      Light   Moderate      Heavy 
# 9300984    3668009    2363856    1091425     599448 
# > table(smoking_var2$smoking)
# 
# Non-smoker  Ex-smoker      Light   Moderate      Heavy 
# 9458454    3986046    2408563    1086603     591932 
# > sum(is.na(smoking_var1$smoking))
# [1] 2386953
# > sum(is.na(smoking_var2$smoking))
# [1] 1879077

### Does this match smoking from the cohort?
cohort_var <- readRDS("data/extraction/cohort_baseline/cohort_var.rds")
table(cohort_var$smoking)

### Somehow, in cohort_var, there are no problems?
sum(is.na(smoking_var$smoking))
sum(is.na(cohort_var$smoking))

### My new version of smoking still doesn't match that from the initial cohort extraction.
### I just need to figure out what the difference is, so I can decide which is correct.

### Currently, cohort_var has less missingness, so that seems good?
### Or maybe there was a problem I overadjusted for.
### Let's take a look at which individuals don't match

### So the new code seems to have loads of missing data, as well as basically no non-smokers...
### Lets merge to see where the differences are
smoking_var$smoking_new <- smoking_var$smoking
cohort_var$smoking_old <- cohort_var$smoking

### Combine
smoking_comb <- dplyr::full_join(cohort_var[,c("patid", "fup_start", "fup_end", "smoking_old")],
                                 smoking_var[,c("patid", "smoking_new")], 
                                 by = dplyr::join_by(patid))
head(smoking_comb)

### Get the ones where there is difference
smoking_comb_discordant <- subset(smoking_comb, 
                                  (smoking_old != smoking_new) & !is.na(smoking_old) & !is.na(smoking_new) | 
                                    !is.na(smoking_old) & is.na(smoking_new) |
                                    is.na(smoking_old) & !is.na(smoking_new)) |>
  dplyr::mutate(fup_start = as.numeric(fup_start), fup_end = as.numeric(fup_end))

### Read in the smoking db queries that I have saved
smoking_comb_discordant_red <- smoking_comb_discordant[1:100, ]


### Run a database query for type 1 and type 2
# db.qry.non <- db_query(codelist.non,
#                        db = db, 
#                        db.filepath = db.filepath,
#                        tab = "obs")
# 
# db.qry.ex <- db_query(codelist.ex,
#                       db = db, 
#                       db.filepath = db.filepath,
#                       tab = "obs")
# 
# db.qry.light <- db_query(codelist.light,
#                          db = db, 
#                          db.filepath = db.filepath,
#                          tab = "obs")
# 
# db.qry.mod <- db_query(codelist.mod,
#                        db = db, 
#                        db.filepath = db.filepath,
#                        tab = "obs")
# 
# db.qry.heavy <- db_query(codelist.heavy,
#                          db = db, 
#                          db.filepath = db.filepath,
#                          tab = "obs")

### Read in the db_qry of the prescriptions of interest (contains all prescription dates of interest)
db.qry.non <- readRDS(paste("../alex/project2/data/p4/db_qry_", "smoking_non", ".rds", sep = "")) |>
  dplyr::arrange(patid, obsdate)
db.qry.ex <- readRDS(paste("../alex/project2/data/p4/db_qry_", "smoking_ex", ".rds", sep = ""))|>
  dplyr::arrange(patid, obsdate)
db.qry.light <- readRDS(paste("../alex/project2/data/p4/db_qry_", "smoking_light", ".rds", sep = ""))|>
  dplyr::arrange(patid, obsdate)
db.qry.mod <- readRDS(paste("../alex/project2/data/p4/db_qry_", "smoking_mod", ".rds", sep = ""))|>
  dplyr::arrange(patid, obsdate)
db.qry.heavy <- readRDS(paste("../alex/project2/data/p4/db_qry_", "smoking_heavy", ".rds", sep = ""))|>
  dplyr::arrange(patid, obsdate)


### lets look at the records for some of these individuals and see whos right
print_records <- function(patid_in){
  print(subset(smoking_comb_discordant, patid == patid_in))
  print("NON-SMOKER")
  print(subset(db.qry.non, patid == patid_in))
  print("EX-SMOKER")
  print(subset(db.qry.ex, patid == patid_in))
  print("LIGHT-SMOKER")
  print(subset(db.qry.light, patid == patid_in))
  print("MODERATE-SMOKER")
  print(subset(db.qry.mod, patid == patid_in))
  print("HEAVY-SMOKER")
  print(subset(db.qry.heavy, patid == patid_in))
}

head(smoking_comb_discordant)
print_records(1000002320274)

db.qry.non <- readRDS("../alex/project2/data/p4/DEBUG_smoking_db_qry_non.rds") |>
  dplyr::arrange(patid, obsdate)
db.qry.ex <- readRDS("../alex/project2/data/p4/DEBUG_smoking_db_qry_ex.rds")|>
  dplyr::arrange(patid, obsdate)
db.qry.light <- readRDS("../alex/project2/data/p4/DEBUG_smoking_db_qry_light.rds")|>
  dplyr::arrange(patid, obsdate)
db.qry.mod <- readRDS("../alex/project2/data/p4/DEBUG_smoking_db_qry_mod.rds")|>
  dplyr::arrange(patid, obsdate)
db.qry.heavy <- readRDS("../alex/project2/data/p4/DEBUG_smoking_db_qry_heavy.rds")|>
  dplyr::arrange(patid, obsdate)




smoking_comb <- dplyr::mutate(smoking_comb, fup_start = as.numeric(fup_start), fup_end = as.numeric(fup_end))

### Lets look at patient 1000000420274
subset(smoking_comb, patid == 1000000420274) 
subset(smoking.non, patid == 1000000420274) |> dplyr::arrange(obsdate)

### Create cohort with appropriate variables
cohort$indexdt <- cohort$fup_start
cohort <- cohort_var[1:50,]
str(cohort)
### Combine queries with cohort, retaining all smoking records prior to the index date
### We treat this as test data, because smoking status may be identified through number of cigarettes smoked per day
### We specify value.na.rm = FALSE, as we want to keep the NA values, because smoking status can also be identified through
### the medcodeid itself.
smoking.non <- combine_query(cohort,
                             db.qry.non,
                             query.type = "test",
                             time.post = 0,
                             numobs = 100,
                             value.na.rm = FALSE)

smoking.ex <- combine_query(cohort,
                            db.qry.ex,
                            query.type = "test",
                            time.post = 0,
                            numobs = 100,
                            value.na.rm = FALSE)

smoking.light <- combine_query(cohort,
                               db.qry.light,
                               query.type = "test",
                               time.post = 0,
                               numobs = 100,
                               value.na.rm = FALSE)

smoking.mod <- combine_query(cohort,
                             db.qry.mod,
                             query.type = "test",
                             time.post = 0,
                             numobs = 100,
                             value.na.rm = FALSE)

smoking.heavy <- combine_query(cohort,
                               db.qry.heavy,
                               query.type = "test",
                               time.post = 0,
                               numobs = 100,
                               value.na.rm = FALSE)

### Currently heavy and moderate have no number of cigarettes smoked per day data
### Light smoker has lots of data on this, as the codes actually include "light or not stated"
### Ex smoker contains lots of values, but cannot be used as looks like they are either used to state
### how many a day used to be smoked, or the year at which smoking was given up.
### Non-smoker currently has no values bigger than zero
#   sum(!is.na(smoking.heavy$value))
#   sum(!is.na(smoking.mod$value))
#   sum(!is.na(smoking.light$value))
#   sum(!is.na(smoking.non$value) & smoking.non$value > 0)
#   sum(!is.na(smoking.ex$value) & smoking.ex$value > 0)

### Add the smoking variable to each dataset
### Assign a smoking value (0 = non, 1 = ex, 2 = light, 3 = moderate, 4 = heavy) to every observation in each query, 
### defined solely by the code lists and medical codes
smoking.non$smoking <- 0
smoking.ex$smoking <- 1
smoking.light$smoking <- 2
smoking.mod$smoking <- 3
smoking.heavy$smoking <- 4

### Change smoking depending on the test value for smoking.light, smoking.mod and smoking.heavy
### We set to NA if more than 100
smoking.light <- dplyr::mutate(smoking.light,
                               smoking = dplyr::case_when(is.na(value) ~ smoking,
                                                          value == 0 ~ 0,
                                                          value > 0 & value < 10 ~ 2,
                                                          value >= 10 & value < 20 ~ 3,
                                                          value >= 20 & value <= 100 ~ 4,
                                                          value > 100 ~ NA)
)

smoking.mod <- dplyr::mutate(smoking.mod,
                             smoking = dplyr::case_when(is.na(value) ~ smoking,
                                                        value == 0 ~ 0,
                                                        value > 0 & value < 10 ~ 2,
                                                        value >= 10 & value < 20 ~ 3,
                                                        value >= 20 & value <= 100 ~ 4,
                                                        value > 100 ~ NA)
)

smoking.heavy <- dplyr::mutate(smoking.heavy,
                               smoking = dplyr::case_when(is.na(value) ~ smoking,
                                                          value == 0 ~ 0,
                                                          value > 0 & value < 10 ~ 2,
                                                          value >= 10 & value < 20 ~ 3,
                                                          value >= 20 & value <= 100 ~ 4,
                                                          value > 100 ~ NA)
)

### Remove the NA values
smoking.light <- smoking.light[!is.na(smoking)]
smoking.mod <- smoking.mod[!is.na(smoking)]
smoking.heavy <- smoking.heavy[!is.na(smoking)]

### Only retain the most recent observation for each
smoking.non <- smoking.non |> 
  dplyr::group_by(patid) |>
  dplyr::filter(dplyr::row_number(dplyr::desc(obsdate)) == 1) 

smoking.ex <- smoking.ex |> 
  dplyr::group_by(patid) |>
  dplyr::filter(dplyr::row_number(dplyr::desc(obsdate)) == 1) 

smoking.light <- smoking.light |> 
  dplyr::group_by(patid) |>
  dplyr::filter(dplyr::row_number(dplyr::desc(obsdate)) == 1) 

smoking.mod <- smoking.mod |> 
  dplyr::group_by(patid) |>
  dplyr::filter(dplyr::row_number(dplyr::desc(obsdate)) == 1) 

smoking.heavy <- smoking.heavy |> 
  dplyr::group_by(patid) |>
  dplyr::filter(dplyr::row_number(dplyr::desc(obsdate)) == 1) 

### Concatenate
variable.dat <- rbind(smoking.non, smoking.ex, smoking.light, smoking.mod, smoking.heavy)

### Arrange so that the first observation is the most recent
### If there are multiple on the same day, we take the most severe smoking status
variable.dat <-variable.dat |> 
  dplyr::arrange(patid, dplyr::desc(obsdate), dplyr::desc(smoking)) |>
  dplyr::group_by(patid)

### Identify those with a smoking history. Given we have only retained one observations from each category,
### this means this individual with > 1 observation must have some sort of smoking history.
smoking.history <- variable.dat |>
  dplyr::summarise(count = dplyr::n()) |> 
  dplyr::filter(count > 1) |>
  dplyr::select(patid) |>
  dplyr::mutate(smoking.history = 1)

### Reduce variable.dat to most recent observation only
variable.dat <- dplyr::slice(variable.dat, 1)

### If their most recent value is non-smoker and they have a smoking history, we must change non-smoker to ex-smoker.
variable.dat <- merge(variable.dat, smoking.history, all.x = TRUE)
variable.dat <- dplyr::mutate(variable.dat,
                              smoking = dplyr::case_when(smoking == 0 & !is.na(smoking.history) ~ 1,
                                                         TRUE ~ smoking))
### Turn into factor variable
variable.dat$smoking <- factor(variable.dat$smoking, 
                               levels = c(0,1,2,3,4), 
                               labels = c("Non-smoker", "Ex-smoker", "Light", "Moderate", "Heavy"))

### Create dataframe of cohort and the variable of interest
variable.dat <- merge(dplyr::select(cohort, patid), variable.dat, by.x = "patid", by.y = "patid", all.x = TRUE)

### Reduce to variables of interest
variable.dat <- variable.dat[,c("patid", "smoking")]

### Change name of variable to varname
colnames(variable.dat)[colnames(variable.dat) == "smoking"] <- varname