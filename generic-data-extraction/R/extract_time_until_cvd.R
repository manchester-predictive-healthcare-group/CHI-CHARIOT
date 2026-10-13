### Program to extract CVD

### This function still needs work if wanting to utilise the t and t.varname options, db.filepath, varname.time or varname.indicator options
extract_time_until_cvd <- function(cohort, 
                                   varname.time = NULL, 
                                   varname.indicator = NULL,
                                   codelist.primary,
                                   codelist.hes,
                                   indexdt, 
                                   censdt,
                                   censdt.lag = 30,
                                   t = NULL, 
                                   t.varname = TRUE,
                                   db.primary,
                                   db.hes,
                                   db.filepath = NULL,
                                   out.save.disk = TRUE, 
                                   out.subdir = NULL,
                                   out.filepath = NULL, 
                                   return.output = FALSE){
  
  
  ### Calculate time until CVD in either primary care data, HES or ONS death data
  cvd_prim <- extract_time_until(cohort = cohort,
                                 varname.time = "cvd_time_prim", 
                                 varname.indicator = "cvd_indicator_prim", 
                                 codelist = codelist.primary, 
                                 indexdt = indexdt,  
                                 censdt = censdt,
                                 censdt.lag = 0,
                                 db = db.primary, 
                                 tab = "obs", 
                                 out.save.disk = FALSE, 
                                 return.output = TRUE)
  
  cvd_hes <- extract_time_until(cohort = cohort,
                                varname.time = "cvd_time_hes", 
                                varname.indicator = "cvd_indicator_hes", 
                                codelist = codelist.hes, 
                                indexdt = indexdt,  
                                censdt = censdt,
                                censdt.lag = 0,
                                db = db.hes, 
                                tab = "hes_primary", 
                                out.save.disk = FALSE, 
                                return.output = TRUE)
  
  cvd_death <- extract_time_until(cohort = cohort,
                                  varname.time = "cvd_time_death", 
                                  varname.indicator = "cvd_indicator_death", 
                                  codelist = codelist.hes, 
                                  indexdt = indexdt, 
                                  censdt = censdt,
                                  censdt.lag = censdt.lag,
                                  db = db.hes, 
                                  tab = "death", 
                                  out.save.disk = FALSE, 
                                  return.output = TRUE)
  
  ### Merge
  variable.dat <- Reduce(function(df1, df2) merge(df1, df2, by = "patid", all.x = TRUE), list(cohort, cvd_prim, cvd_hes, cvd_death))
  
  ### Create the varaible of interest from these
  ### if an individual has history of type 1 and type 2, they will be classed as type 1
  variable.dat <- dplyr::mutate(variable.dat, 
                                cvd_time = pmin(cvd_time_prim, cvd_time_hes, cvd_time_death),
                                cvd_indicator = pmax(cvd_indicator_prim, cvd_indicator_hes, cvd_indicator_death))
  
  ### Reduce to variables of interest
  variable.dat <- dplyr::select(variable.dat, patid, 
                                cvd_time, cvd_indicator, 
                                cvd_time_prim, cvd_indicator_prim, 
                                cvd_time_hes, cvd_indicator_hes, 
                                cvd_time_death, cvd_indicator_death)
  
  ### Implement output
  implement_output(variable.dat, "time_until_cvd", out.save.disk, out.subdir, out.filepath, return.output)
  
  #   ### Turn into factor variable
  #   variable.dat$diabetes <- factor(variable.dat$diabetes, 
  #                                   levels = c(0,1,2), 
  #                                   labels = c("Absent", "Type1", "Type2"))
  #   
  #   ### Reduce to variables of interest
  #   variable.dat <- variable.dat[,c("patid", "diabetes")]
  #   
  #   ### Change name of variable to varname
  #   colnames(variable.dat)[colnames(variable.dat) == "diabetes"] <- varname
  #   
  #   ### Implement output
  #   implement_output(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output)
  
}