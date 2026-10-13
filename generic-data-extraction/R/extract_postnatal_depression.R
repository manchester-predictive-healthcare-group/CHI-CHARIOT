### Program to extract postnatal depression
### Want to either identify this through history of medical code, or EPDS score
extract_postnatal_depression <- function(cohort, 
                                         varname = NULL, 
                                         codelist.med,
                                         codelist.test,
                                         indexdt, 
                                         t = NULL, 
                                         t.varname = TRUE,
                                         db, 
                                         db.filepath = NULL,
                                         out.save.disk = TRUE, 
                                         out.subdir = NULL,
                                         out.filepath = NULL, 
                                         return.output = FALSE){
  
  #   cohortZ <- readRDS("data/extraction/cohortZ.rds")
  #           varname = NULL
  #           codelist.med = "uom_postnatal_depression_medcodeid"
  #           codelist.test = "uom_postnatal_depression_score_medcodeid"
  #           cohort = cohortZ
  #           indexdt = "fup_start"
  #           t = NULL
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
    varname <- "postnatal_depression"
  }
  ## Change variable name based off time point specified for extraction
  varname <- prep_varname(varname, t, t.varname)
  ## Create named subdirectory if it doesn't exist
  prep_subdir(out.subdir)
  
  ### Run a database query for type 1 and type 2
  db.qry.med <- db_query(codelist.med,
                         db = db, 
                         db.filepath = db.filepath,
                         tab = "obs")
  
  db.qry.test <- db_query(codelist.test,
                          db = db, 
                          db.filepath = db.filepath,
                          tab = "obs")
  
  ### Identify which individuals have a history of postnatal depression through the medical codes
  cohort[,"med"] <- combine_query_boolean(cohort, 
                                          db.qry.med,
                                          query.type = "obs")
  
  ### Get test data for individuals in cohort
  ### No time range, as interested in histry of, but only interested in those with a score higher than 14, indicating moderate depression
  postnatal_test <- combine_query(cohort, 
                                  db.qry.test, 
                                  query.type = "test",
                                  lower.bound = 14, 
                                  upper.bound = 31)
  
  ### Only retain the most recent observation for each
  postnatal_test <- postnatal_test |> 
    dplyr::group_by(patid) |>
    dplyr::filter(dplyr::row_number(dplyr::desc(obsdate)) == 1) 
  
  ### Create variable for history of postnatal depression identified through the test data
  cohort[,"test"] <- as.integer(!is.na(fastmatch::fmatch(cohort$patid, postnatal_test$patid)))
  
  ### Create the varaible of interest from these
  ### if an individual has history of type 1 and type 2, they will be classed as type 1
  variable.dat <- dplyr::mutate(cohort, postnatal_depression = pmax(med, test))
  
  ### Reduce to variables of interest
  variable.dat <- variable.dat[,c("patid", "postnatal_depression")]
  
  ### Change name of variable to varname
  colnames(variable.dat)[colnames(variable.dat) == "postnatal_depression"] <- varname
  
  ### Implement output
  implement_output(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output)
  
}