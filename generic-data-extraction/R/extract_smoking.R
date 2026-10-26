### Program to extract smoking status
extract_smoking <- function(cohort, 
                            varname = NULL, 
                            codelist.non,
                            codelist.ex,
                            codelist.light,
                            codelist.mod,
                            codelist.heavy,
                            indexdt, 
                            t = NULL, 
                            t.varname = TRUE,
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
    varname <- "smoking"
  }
  ## Change variable name based off time point specified for extraction
  varname <- prep_varname(varname, t, t.varname)
  ## Create named subdirectory if it doesn't exist
  prep_subdir(out.subdir)
  
  print("PROGRESS1")
  
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
  
  print("PROGRESS1")
  
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
  
  print("PROGRESS2")
  
  ### Currently heavy and moderate have no number of cigarettes smoked per day data
  ### Light smoker has lots of data on this, as the codes actually include "light or not stated"
  ### Ex smoker contains lots of values, but cannot be used as looks like they are either used to state
  ### how many a day used to be smoked, or the year at which smoking was given up.
  ### Non-smoker currently has no values bigger than zero
  
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
                                                            value > 100 ~ NA_real_)
  )
  
  smoking.mod <- dplyr::mutate(smoking.mod,
                               smoking = dplyr::case_when(is.na(value) ~ smoking,
                                                          value == 0 ~ 0,
                                                          value > 0 & value < 10 ~ 2,
                                                          value >= 10 & value < 20 ~ 3,
                                                          value >= 20 & value <= 100 ~ 4,
                                                          value > 100 ~ NA_real_)
  )
  
  smoking.heavy <- dplyr::mutate(smoking.heavy,
                                 smoking = dplyr::case_when(is.na(value) ~ smoking,
                                                            value == 0 ~ 0,
                                                            value > 0 & value < 10 ~ 2,
                                                            value >= 10 & value < 20 ~ 3,
                                                            value >= 20 & value <= 100 ~ 4,
                                                            value > 100 ~ NA_real_)
  )
  
  ### Remove the NA values
  smoking.light <- smoking.light[!is.na(smoking)]
  smoking.mod <- smoking.mod[!is.na(smoking)]
  smoking.heavy <- smoking.heavy[!is.na(smoking)]
  
  ###
  ### Combine smoking records and assign final smoking status
  ###
  
  ### Combine without unnecessary copying
  variable.dat <- data.table::rbindlist(
    list(
      smoking.non,
      smoking.ex,
      smoking.light,
      smoking.mod,
      smoking.heavy
    ),
    use.names = TRUE,
    fill = FALSE
  )
  
  ### Ensure data.table format
  data.table::setDT(variable.dat)
  
  ### Record whether each patient has ever been an ex-smoker or current smoker
  variable.dat[, smoking.history := any(smoking %in% c(1, 2, 3, 4)), by = patid]
  
  ### Order by patient, most recent date, then most severe status
  data.table::setorder(
    variable.dat,
    patid,
    -obsdate,
    -smoking
  )
  
  ### Retain the most recent record per patient
  ### The ordering ensures that the most severe status is used for tied dates
  variable.dat <- variable.dat[!duplicated(patid)]
  
  ### A known previous smoker cannot subsequently be classified as never-smoker
  variable.dat[
    smoking == 0 & smoking.history,
    smoking := 1
  ]
  
  ### Check that the recoding worked
  stopifnot(
    variable.dat[
      smoking == 0 & smoking.history,
      .N
    ] == 0L
  )
  
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
  
  ### Implement output
  implement_output(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output)
  
}