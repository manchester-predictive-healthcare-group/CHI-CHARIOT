### Write a function which will extract the most recent weight
extract_weight <- function(cohort, 
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

  ### Preparation
  ## Add index date variable to cohort and change indexdt based on t
  cohort <- prep_cohort(cohort, indexdt, t)
  ## Assign variable name if unspecified
  if (is.null(varname)){
    varname <- "value"
  }
  ## Change variable name based off time point specified for extraction
  varname <- prep_varname(varname, t, t.varname)
  ## Create named subdirectory if it doesn't exist
  prep_subdir(out.subdir)
  
  ### Run a database query
  db.qry <- db_query(codelist,
                     db = db, 
                     db.filepath = db.filepath,
                     tab = "obs")
  
  ### Get test data for individuals in cohort, within time range and remove outliers
  variable.dat <- combine_query(cohort, 
                                db.qry, 
                                query.type = "test",
                                time.prev = time.prev, 
                                time.post = time.post)
  
  ### For the weight query, we need to rescale those with numunitid = 1691, 2318, 2997 or 6265 from stones to kg
  variable.dat <- dplyr::mutate(variable.dat, value = dplyr::case_when(numunitid %in% c(1691, 2318, 2997, 6265) ~ 6.35029*value,
                                                                       !(numunitid %in% c(1691, 2318, 2997, 6265)) ~ value))
  
  ### If values are missing, < lower.bound or > upper.bound then delete
  if (!is.null(lower.bound) & !is.null(upper.bound)){
    variable.dat <- variable.dat[value > lower.bound & value < upper.bound]
  } else if (is.null(lower.bound) & !is.null(upper.bound)){
    variable.dat <- variable.dat[value < upper.bound]
  } else if (!is.null(lower.bound) & is.null(upper.bound)){
    variable.dat <- variable.dat[value > lower.bound]
  }
  
  print("merge")
  ### Create dataframe of cohort and the variable of interest
  variable.dat <- merge(dplyr::select(cohort, patid), variable.dat, by.x = "patid", by.y = "patid", all.x = TRUE)
  
  ### Reduce to variables of interest
  variable.dat <- variable.dat[,c("patid", "value")]
  
  ### Change name of variable to varname
  colnames(variable.dat)[colnames(variable.dat) == "value"] <- varname
  
  implement_output(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output)
  
}