### Extract the time until something occurs in primary care data
extract_time_until <- function(cohort,
                               varname.time = NULL,
                               varname.indicator = NULL,
                               codelist = NULL,
                               codelist.vector = NULL,
                               indexdt,
                               censdt,
                               censdt.lag = 0,
                               t = NULL,
                               t.varname = TRUE,
                               db.open = NULL,
                               db = NULL,
                               db.filepath = NULL,
                               tab = c("obs", "drug", "hes_primary", "death"),
                               out.save.disk = FALSE,
                               out.subdir = NULL,
                               out.filepath = NULL,
                               return.output = FALSE){
  
  #         varname.time = NULL
  #         varname.indicator = NULL
  #         codelist = "edh_cvd_hist_medcodeid"
  #         cohort = cohortZ
  #         indexdt = "fup_start"
  #         t = NULL
  #         db = "aurum_small"
  #         tab = "obs"
  #         db.filepath = NULL
  #         out.save.disk = TRUE
  #         out.filepath = NULL
  #         out.subdir = NULL
  #         return.output = TRUE
  
  # cohort = pat
  # varname.time = NULL
  # varname.indicator = NULL
  # codelist.vector = codelist
  # indexdt = "fup_start"
  # censdt = "fup_end"
  # db.open = aurum_extract3
  # tab = "observation"
  # return.output = TRUE
  #         db.filepath = NULL
  #         out.save.disk = FALSE
  #         out.filepath = NULL
  #         out.subdir = NULL
  
  ### Preparation
  ## Add index date variable to cohort and change indexdt based on t
  cohort <- prep_cohort(cohort, indexdt, t, reduce = FALSE)
  ### Change name of censoring date variable to "censdt" so we can easily refer to it
  colnames(cohort)[colnames(cohort) == censdt] <- "censdt"
  
  ## Reduce cohort to variables of interest
  cohort <- cohort[,c("patid", "indexdt", "censdt")]
  ## Assign variable name if unspecified
  if (is.null(varname.time)){
    varname.time <- "var_time"
  }
  if (is.null(varname.indicator)){
    varname.indicator <- "var_indicator"
  }
  
  ## Change variable name based off time point specified for extraction
  varname.time <- prep_varname(varname.time, t, t.varname)
  varname.indicator <- prep_varname(varname.indicator, t, t.varname)
  
  ## Create named subdirectory if it doesn't exist
  prep_subdir(out.subdir)
  
  ### Run a database query
  db.qry <- db_query(codelist,
                     db.open = db.open,
                     db = db,
                     db.filepath = db.filepath,
                     tab = tab,
                     codelist.vector = codelist.vector)
  
  ### Identify the first CVD event happening after the indexdt
  ## If tab = "obs", this could be a query.type of "med" or "test", choose "med" as not interested in test results themselves
  query.type <- tab
  if (query.type == "obs"){query.type <- "med"}
  
  ## Combine query
  ## reduce.output = FALSE because we want access to censdt and
  variable.dat <- combine_query(cohort,
                                db.query = db.qry,
                                query.type = query.type,
                                time.prev = 0,
                                time.post = Inf,
                                numobs = 1,
                                value.na.rm = FALSE,
                                earliest.values = TRUE,
                                reduce.output = FALSE)
  
  ### Calculate the time until event of interest, set to NA and remove if beyond censdt
  variable.dat <-
    dplyr::mutate(variable.dat,
                  var_time = dplyr::case_when(obsdate <= censdt + censdt.lag ~ pmin(obsdate, censdt) - as.numeric(indexdt),
                                              obsdate > censdt + censdt.lag ~ NA),
                  var_indicator = dplyr::case_when(!is.na(var_time) ~ 1,
                                                   TRUE ~ NA)) |>
    dplyr::filter(!is.na(var_time))
  
  ### Reduce to variables of interst
  variable.dat <- variable.dat[,c("patid", "var_time", "var_indicator")]
  
  ### Merge back with cohort
  variable.dat <- merge(dplyr::select(cohort, patid, indexdt, censdt), variable.dat, by.x = "patid", by.y = "patid", all.x = TRUE)
  
  ### If the event has NA, se the time to censdt, and indicator to 0
  variable.dat <- dplyr::mutate(variable.dat,
                                var_indicator = dplyr::case_when(!is.na(var_time) ~ var_indicator,
                                                                 is.na(var_time) ~ 0),
                                var_time = dplyr::case_when(!is.na(var_time) ~ var_time,
                                                            is.na(var_time) ~ as.numeric(censdt - indexdt))
  )
  
  ### Reduce to variables of interest
  variable.dat <- variable.dat[,c("patid", "var_time", "var_indicator")]
  
  ### Change name of variable to varname
  colnames(variable.dat)[colnames(variable.dat) == "var_time"] <- varname.time
  colnames(variable.dat)[colnames(variable.dat) == "var_indicator"] <- varname.indicator
  
  ### Implement output
  implement_output(variable.dat, varname.time, out.save.disk, out.subdir, out.filepath, return.output)
  
}


# extract_time_until <- function(cohort, 
#                                     varname.time = NULL, 
#                                     varname.indicator = NULL,
#                                     codelist, 
#                                     codelist.vector = NULL,
#                                     indexdt, 
#                                     t = NULL, 
#                                     t.varname = TRUE,
#                                     db, 
#                                     tab,
#                                     db.filepath = NULL,
#                                     out.save.disk = TRUE, 
#                                     out.filename = NULL,
#                                     out.subdir = NULL,
#                                     out.filepath = NULL,
#                                     return.output = FALSE){
#   
# #         varname.time = NULL
# #         varname.indicator = NULL
# #         codelist = "edh_cvd_hist_medcodeid"
# #         cohort = cohortZ
# #         indexdt = "fup_start"
# #         t = NULL
# #         db = "aurum_small"
# #         tab = "obs"
# #         db.filepath = NULL
# #         out.save.disk = TRUE
# #         out.filepath = NULL
# #         out.subdir = NULL
# #         return.output = TRUE
#   
#   ### Preparation
#   ## Add index date variable to cohort and change indexdt based on t
#   cohort <- prep_cohort(cohort, indexdt, t, reduce = FALSE)
#   ## Reduce cohort to variables of interest
#   cohort <- cohort[,c("patid", "indexdt", "fup_end")]
#   ## Assign variable name if unspecified
#   if (is.null(varname.time)){
#     varname.time <- "var_time"
#   }
#   if (is.null(varname.indicator)){
#     varname.indicator <- "var_indicator"
#   }
#   if (is.null(out.filename)){
#     out.filename <- "time_until"
#   }
#   
#   ## Change variable name based off time point specified for extraction
#   varname.time <- prep_varname(varname.time, t, t.varname)
#   varname.indicator <- prep_varname(varname.indicator, t, t.varname)
#   out.filename <- prep_varname(out.filename, t, t.varname)
#   
#   ## Create named subdirectory if it doesn't exist
#   prep_subdir(out.subdir)
#   
#   ### Run a database query
#   db.qry <- db_query(codelist,
#                      db, 
#                      db.filepath,
#                      tab = tab,
#                      codelist.vector = codelist.vector)
#   
#   ### Identify the first CVD event happening after the indexdt
#   ## Define query.type
#   query.type <- tab
#   if (query.type == "observation"){query.type = "med"}
#   
#   ## Combine query
#   ## reduce.output = FALSE because we want access to censdt and
#   variable.dat <- combine_query(cohort,
#                                 db.query = db.qry,
#                                 query.type = tab,
#                                 time.prev = 0,
#                                 time.post = Inf,
#                                 numobs = 1,
#                                 value.na.rm = FALSE,
#                                 earliest.values = TRUE,
#                                 reduce.output = FALSE)
#   
#   ### Calculate the time until event of interest, set to NA and remove if beyond fup_end
#   variable.dat <- 
#     dplyr::mutate(variable.dat, 
#                   var_time = dplyr::case_when(obsdate <= fup_end ~ obsdate - as.numeric(indexdt),
#                                                    obsdate > fup_end ~ NA),
#                   var_indicator = dplyr::case_when(!is.na(var_time) ~ 1,
#                                                    TRUE ~ NA)) |>
#     dplyr::filter(!is.na(var_time))
#   
#   ### Reduce to variables of interst
#   variable.dat <- variable.dat[,c("patid", "var_time", "var_indicator")]
#   
#   ### Merge back with cohort
#   variable.dat <- merge(dplyr::select(cohort, patid, indexdt, fup_end), variable.dat, by.x = "patid", by.y = "patid", all.x = TRUE)
#   
#   ### If the event has NA, se the time to fup_end, and indicator to 0
#   variable.dat <- dplyr::mutate(variable.dat, 
#                                 var_indicator = dplyr::case_when(!is.na(var_time) ~ var_indicator,
#                                                                  is.na(var_time) ~ 0),
#                                 var_time = dplyr::case_when(!is.na(var_time) ~ var_time,
#                                                                  is.na(var_time) ~ as.numeric(fup_end - indexdt))
#                                 )
#                                 
#   ### Reduce to variables of interest
#   variable.dat <- variable.dat[,c("patid", "var_time", "var_indicator")]
#   
#   ### Change name of variable to varname
#   colnames(variable.dat)[colnames(variable.dat) == "var_time"] <- varname.time
#   colnames(variable.dat)[colnames(variable.dat) == "var_indicator"] <- varname.indicator
#   
#   ### Implement output
#   implement_output(variable.dat, out.filename, out.save.disk, out.subdir, out.filepath, return.output)
#   
# }