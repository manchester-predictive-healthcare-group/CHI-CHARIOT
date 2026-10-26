### Write a function which will extract the non-HDL cholesterol variable
extract_nonhdl <- function(cohort, 
                           varname = NULL, 
                           codelist.chol,
                           codelist.hdl,
                           codelist.ldl,
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
    varname <- "nonhdl"
  }
  ## Change variable name based off time point specified for extraction
  varname <- prep_varname(varname, t, t.varname)
  ## Create named subdirectory if it doesn't exist
  prep_subdir(out.subdir)
  
  ### Need to run three database queries, one for ratio, one for chol and one for hdl
  db.qry.chol <- db_query(codelist.chol,
                          db = db, 
                          db.filepath = db.filepath,
                          tab = "obs")
  
  db.qry.hdl <- db_query(codelist.hdl,
                         db = db, 
                         db.filepath = db.filepath,
                         tab = "obs")
  
  db.qry.ldl <- db_query(codelist.ldl,
                         db = db, 
                         db.filepath = db.filepath,
                         tab = "obs")
  
  ### Get latest test result for ratio, chol and hdl from the last five years
  variable.dat.chol <- combine_query(cohort, 
                                     db.qry.chol, 
                                     query.type = "test",
                                     time.prev = time.prev, 
                                     time.post = time.post,
                                     remove.zeros = TRUE)
  
  variable.dat.hdl <- combine_query(cohort, 
                                    db.qry.hdl, 
                                    query.type = "test",
                                    time.prev = time.prev, 
                                    time.post = time.post,
                                    remove.zeros = TRUE)
  
  variable.dat.ldl <- combine_query(cohort, 
                                    db.qry.ldl, 
                                    query.type = "test",
                                    time.prev = time.prev, 
                                    time.post = time.post,
                                    remove.zeros = TRUE)
  
  ### Covert ldl to non-hdl
  variable.dat.ldl <- dplyr::mutate(variable.dat.ldl, nonhdl = 1.24*value, type = 2) |>
    dplyr::select(patid, obsdate, nonhdl, type) |>
    dplyr::filter(nonhdl > lower.bound & nonhdl < upper.bound) |>
    as.data.frame()
  
  ### Merge cholesterol and hdl
  variable.dat.chol.hdl <- merge(dplyr::select(variable.dat.chol, patid, obsdate, value), 
                                 dplyr::select(variable.dat.hdl, patid, obsdate, value), 
                                  by.x = "patid", by.y = "patid", allow.cartesian = TRUE)
  
  ### Estimate nonhdl from cholesterol-hdl
  variable.dat.chol.hdl <- dplyr::mutate(variable.dat.chol.hdl, nonhdl = value.x - value.y, type = 1) |>
    dplyr::filter(nonhdl > lower.bound & nonhdl < upper.bound) |>
    dplyr::arrange(patid, dplyr::desc(obsdate.x), dplyr::desc(obsdate.y)) |>
    dplyr::group_by(patid, obsdate.x) |>
    dplyr::filter(dplyr::row_number() == 1) |>
    dplyr::rename(obsdate = obsdate.x) |>
    dplyr::select(patid, obsdate, nonhdl, type) |>
    as.data.frame()
  
  ### Now combine and take most recent, preferring type == 1
  variable.dat <- rbind(variable.dat.chol.hdl, variable.dat.ldl) |> 
    dplyr::mutate(type = factor(type, levels = c(1,2), labels = c("chol_hdl", "ldl"))) |>
    dplyr::arrange(patid, desc(obsdate), type, desc(nonhdl)) |>
    dplyr::group_by(patid) |>
    dplyr::filter(dplyr::row_number() == 1) |>
    as.data.frame() 
  
  ### Create dataframe of cohort and the variable of interest
  variable.dat <- merge(dplyr::select(cohort, patid), variable.dat, by.x = "patid", by.y = "patid", all.x = TRUE)
  
  ### Reduce to variables of interest
  variable.dat <- variable.dat[,c("patid", "nonhdl")]
  
  ### Change name of variable to varname
  colnames(variable.dat)[colnames(variable.dat) == "nonhdl"] <- varname
  
  ### Implement output
  implement_output(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output)
  
}