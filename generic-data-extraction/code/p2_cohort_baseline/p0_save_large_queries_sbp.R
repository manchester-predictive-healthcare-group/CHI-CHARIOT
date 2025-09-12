### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
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
lower.bound <- 70
upper.bound <- 210
db = "aurum"

### Assign code lists
codelist.med = "edh_sbp_medcodeid"
codelist.drug = "uom_antihypertensives_prodcodeid"

### Preparation
## Add index date variable to cohort and change indexdt based on t
cohort <- prep_cohort(cohort, indexdt, t)

### Run database queries
print(paste("running queries med", Sys.time()))
db.qry.med <- db_query(codelist.med,
                       db = db, 
                       db.filepath = db.filepath,
                       tab = "obs")
print(paste("running queries drug", Sys.time()))
db.qry.drug <- db_query(codelist.drug,
                        db = db, 
                        db.filepath = db.filepath,
                        tab = "drug")

### Combine query with cohort
print(paste("combine queries med", Sys.time()))
variable.dat.med <- combine_query(cohort, 
                                  db.qry.med, 
                                  query.type = "test",
                                  time.prev = time.prev, 
                                  time.post = time.post, 
                                  lower.bound = lower.bound, 
                                  upper.bound = upper.bound, 
                                  numobs = Inf)

### Combine query with cohort
print(paste("combine queries drug", Sys.time()))
variable.dat.drug <- combine_query(cohort, 
                                   db.qry.drug, 
                                   query.type = "drug",
                                   time.prev = Inf, 
                                   time.post = time.post, 
                                   numobs = Inf)

### Save as .rds objects
saveRDS(variable.dat.med, "data/extraction/cohort_baseline/query_sbp.rds")
saveRDS(variable.dat.drug, "data/extraction/cohort_baseline/query_antihypertensives.rds")

### Save these queries as .rds objects, so that we can run the removed_exposed function on lower memory nodes.
