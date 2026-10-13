### Write a function which will extract the SBP variable
extract_sbp <- function(cohort, 
                        varname = NULL,
                        codelist,
                        indexdt, 
                        t = NULL, 
                        t.varname = TRUE,
                        time.prev = 365.25*5, 
                        time.post = 0, 
                        lower.bound = -Inf, 
                        upper.bound = Inf,
                        db, 
                        db.filepath = NULL,
                        out.save.disk = TRUE, 
                        out.subdir = NULL,
                        out.filepath = NULL, 
                        return.output = FALSE){
  
  #           varname = NULL
  #           codelist = "edh_sbp_medcodeid"
  #           cohort = cohortZ
  #           indexdt = "fup_start"
  #           t = 0
  #           t.varname = TRUE
  #           time.prev = 365.25*5
  #           time.post = 0
  #           lower.bound = -Inf
  #           upper.bound = Inf
  #           db = "aurum_small"
  #           db.filepath = NULL
  #           out.save.disk = FALSE
  #           out.filepath = NULL
  #           out.subdir = NULL
  #           return.output = TRUE
  
  ### Preparation
  ## Add index date variable to cohort and change indexdt based on t
  cohort <- prep_cohort(cohort, indexdt, t)
  ## Assign variable name if unspecified
  if (is.null(varname)){
    varname <- "sbp"
  }
  ## Change variable name based off time point specified for extraction
  varname <- prep_varname(varname, t, t.varname)
  ## Create named subdirectory if it doesn't exist
  prep_subdir(out.subdir)
  
  ### Run a database query
  ### THINKING I MAY HAVE TO SPLIT THE PATIDS UP INTO GROUPS BEFORE QUERYING THE DATABASE FOR SBP CODES AS THERE ARE SO MANY
  db.qry <- db_query(codelist,
                     db = db, 
                     db.filepath = db.filepath,
                     tab = "obs")
  
  ### Get test data for individuals in cohort, within time range and remove outliers
  variable.dat <- combine_query(cohort, 
                                db.qry, 
                                query.type = "test",
                                time.prev = time.prev, 
                                time.post = time.post, 
                                lower.bound = lower.bound, 
                                upper.bound = upper.bound)
  
  print("merge")
  ### Create dataframe of cohort and the variable of interest
  variable.dat <- merge(dplyr::select(cohort, patid), variable.dat, by.x = "patid", by.y = "patid", all.x = TRUE)
  
  ### Reduce to variables of interest
  variable.dat <- variable.dat[,c("patid", "value")]
  
  ### Change name of variable to varname
  colnames(variable.dat)[colnames(variable.dat) == "value"] <- varname
  
  implement_output(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output)
  
}