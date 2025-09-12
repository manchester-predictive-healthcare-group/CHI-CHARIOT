### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Define inputs
cohort <- readRDS(paste("data/extraction/cohort_exclu3.rds", sep = ""))
cohort <- data.table::as.data.table(cohort)
indexdt <- "fup_start"
out.subdir <- "cohort_baseline"
t <- as.numeric(0)
if (t == 0){t <- NULL}

### Define other variables
time.prev = 365.25*5
time.post = 0
lower.bound <- 1
upper.bound <- 12
db = "aurum"

### Assign code lists
codelist.ratio = "muzambi_ratio_medcodeid"
codelist.chol = "muzambi_cholesterol_medcodeid"
codelist.hdl = "muzambi_hdl_medcodeid"
codelist.drug = "uom_statins_prodcodeid"

### Preparation
## Add index date variable to cohort and change indexdt based on t
cohort <- prep_cohort(cohort, indexdt, t)

### Need to run three database queries, one for ratio, one for chol and one for hdl
print(paste("run queries ratio", Sys.time()))
db.qry.ratio <- db_query(codelist.ratio,
                         db = db, 
                         db.filepath = db.filepath,
                         tab = "obs")
print(paste("run queries chol", Sys.time()))
db.qry.chol <- db_query(codelist.chol,
                        db = db, 
                        db.filepath = db.filepath,
                        tab = "obs")
print(paste("run queries hdl", Sys.time()))
db.qry.hdl <- db_query(codelist.hdl,
                       db = db, 
                       db.filepath = db.filepath,
                       tab = "obs")

### Query for prescriptions
print(paste("run queries drug", Sys.time()))
db.qry.drug <- db_query(codelist.drug,
                        db = db, 
                        db.filepath = db.filepath,
                        tab = "drug")

### Combine drug query with cohort
print(paste("combine queries drug", Sys.time()))
variable.dat.drug <- combine_query(cohort, 
                                   db.qry.drug, 
                                   query.type = "drug",
                                   time.prev = Inf, 
                                   time.post = time.post, 
                                   numobs = Inf)

### Combine ratio query with cohort
print(paste("combine queries ratio", Sys.time()))
variable.dat.ratio <- combine_query(cohort, 
                                    db.qry.ratio, 
                                    query.type = "test",
                                    time.prev = time.prev, 
                                    time.post = time.post, 
                                    lower.bound = lower.bound, 
                                    upper.bound = upper.bound, 
                                    numobs = Inf)

### Combine chol query with cohort
print(paste("combine queries chol", Sys.time()))
variable.dat.chol <- combine_query(cohort, 
                                   db.qry.chol, 
                                   query.type = "test",
                                   time.prev = time.prev, 
                                   time.post = time.post, 
                                   numobs = Inf)

### Combine hdl query with cohort
print(paste("combine queries hdl", Sys.time()))
variable.dat.hdl <- combine_query(cohort, 
                                  db.qry.hdl, 
                                  query.type = "test",
                                  time.prev = time.prev, 
                                  time.post = time.post, 
                                  numobs = Inf)

### Save as .rds objects
saveRDS(variable.dat.ratio, "data/extraction/cohort_baseline/query_cholhdl.rds")
saveRDS(variable.dat.chol, "data/extraction/cohort_baseline/query_chol.rds")
saveRDS(variable.dat.hdl, "data/extraction/cohort_baseline/query_hdl.rds")
saveRDS(variable.dat.drug, "data/extraction/cohort_baseline/query_statins.rds")

### Save these queries as .rds objects, so that we can run the removed_exposed function on lower memory nodes.
