### Lets get a couple of things straight
### Unless otherwise specified, filename and varname will follow a default naming convention, specified purely through the varname
### We will assume when people are making multiple cohorts, they will put them in seperate directories within "extraction"
### This means we don't have to include
### There will be a default naming convention, based off the index date, and time from index date at which the varaible is extracted
### var_varname_indexdate_tX
### Write a function which will query the database and look for diagnoses mathcing the codelist
extract_ho <- function(cohort, 
                       varname = NULL, 
                       codelist, 
                       indexdt, 
                       t = NULL, 
                       t.varname = TRUE,
                       time.prev = Inf,
                       time.post = 0,
                       numobs = 1,
                       db, 
                       tab,
                       db.filepath = NULL,
                       out.save.disk = TRUE, 
                       out.subdir = NULL,
                       out.filepath = NULL, 
                       return.output = FALSE){
  
  #       varname = "ho_ra_test10"
  #       codelist = "edh_ra_medcodeid"
  #       cohort = cohortA
  #       indexdt = "fup_start"
  #       t = NULL
  #       db = "aurum_nosubset_randset"
  #       db.filepath = NULL
  #       out.save.disk = TRUE
  #       out.filepath = NULL
  #       out.subdir = NULL
  #       return.output = TRUE
  
  ### Preparation
  ## Add index date variable to cohort and change indexdt based on t
  cohort <- prep_cohort(cohort, indexdt, t)
  ## Assign variable name if unspecified
  if (is.null(varname)){
    varname <- "ho"
  }
  ## Change variable name based off time point specified for extraction
  varname <- prep_varname(varname, t, t.varname)
  ## Create named subdirectory if it doesn't exist
  prep_subdir(out.subdir)
  
  ### Run a database query
  db.qry <- db_query(codelist = codelist,
                     db = db, 
                     db.filepath = db.filepath,
                     tab = tab)
  print(str(db.qry))
  
  ### Identify which individuals have a history of XXX
  cohort[,"ho"] <- combine_query_boolean(cohort, 
                                         db.query = db.qry,
                                         query.type = tab,
                                         time.prev = time.prev, 
                                         time.post = time.post, 
                                         numobs = numobs)
  
  ### Reduce to variables of interest
  variable.dat <- cohort[,c("patid", "ho")]
  
  ### Change name of variable to varname
  colnames(variable.dat)[colnames(variable.dat) == "ho"] <- varname
  
  ### Implement output
  implement_output(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output)
  
}