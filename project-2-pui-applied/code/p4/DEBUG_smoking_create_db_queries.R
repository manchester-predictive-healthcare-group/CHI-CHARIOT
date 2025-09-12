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

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
variable.name <- "smoking"
cohort.base.name <- "cohort_exclu3"
cohort.base <- readRDS(paste("data/extraction/", cohort.base.name ,".rds", sep = ""))
cohort.base <- cohort.base[1:5000, ]
cohort.base <- data.table::as.data.table(cohort.base)
indexdt.name <- "fup_start"
out.subdir <- "cohort_baseline"
t <- as.numeric(0)
if (t == 0){t <- NULL}

indexdt = indexdt.name
db = "aurum"
db.filepath = NULL

codelist.non = "edh_smoking_non_medcodeid"
codelist.ex = "edh_smoking_ex_medcodeid"
codelist.light = "edh_smoking_light_medcodeid"
codelist.mod = "edh_smoking_mod_medcodeid"
codelist.heavy = "edh_smoking_heavy_medcodeid"

### Run a database query for type 1 and type 2
db.qry.non <- db_query(codelist.non,
                       db = db,
                       db.filepath = db.filepath,
                       tab = "obs")

db.qry.ex <- db_query(codelist.ex,
                      db = db,
                      db.filepath = db.filepath,
                      tab = "obs")

db.qry.light <- db_query(codelist.light,
                         db = db,
                         db.filepath = db.filepath,
                         tab = "obs")

db.qry.mod <- db_query(codelist.mod,
                       db = db,
                       db.filepath = db.filepath,
                       tab = "obs")

db.qry.heavy <- db_query(codelist.heavy,
                         db = db,
                         db.filepath = db.filepath,
                         tab = "obs")

saveRDS(db.qry.non, paste("../alex/project2/data/p4/DEBUG_smoking_db_qry_non.rds", sep = ""))
saveRDS(db.qry.ex, paste("../alex/project2/data/p4/DEBUG_smoking_db_qry_ex.rds", sep = ""))
saveRDS(db.qry.light, paste("../alex/project2/data/p4/DEBUG_smoking_db_qry_light.rds", sep = ""))
saveRDS(db.qry.mod, paste("../alex/project2/data/p4/DEBUG_smoking_db_qry_mod.rds", sep = ""))
saveRDS(db.qry.heavy, paste("../alex/project2/data/p4/DEBUG_smoking_db_qry_heavy.rds", sep = ""))
print("SAVED")
cohort <- cohort.base

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

saveRDS(smoking.non, paste("../alex/project2/data/p4/DEBUG_smoking_combine_qry_non.rds", sep = ""))
saveRDS(smoking.ex, paste("../alex/project2/data/p4/DEBUG_smoking_combine_qry_ex.rds", sep = ""))
saveRDS(smoking.light, paste("../alex/project2/data/p4/DEBUG_smoking_combine_qry_light.rds", sep = ""))
saveRDS(smoking.mod, paste("../alex/project2/data/p4/DEBUG_smoking_combine_qry_mod.rds", sep = ""))
saveRDS(smoking.heavy, paste("../alex/project2/data/p4/DEBUG_smoking_combine_qry_heavy.rds", sep = ""))

print("FINISHED")