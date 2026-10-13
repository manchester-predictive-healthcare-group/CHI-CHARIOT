### Write a function which will extract test data, but remove observations if there is a prescription
### within a specified number of days prior to the test
extract_test_unexposed <- function(cohort, 
                                   varname = NULL,
                                   codelist.med,
                                   codelist.drug,
                                   query.med = NULL,
                                   query.drug = NULL,
                                   task.id = NULL,
                                   task.id.n = NULL,
                                   range.max = 180,
                                   range.min = 0,
                                   range.num = 1,
                                   keep.exposed = FALSE,
                                   exposed.type = c("mr", "np", "md", "all"),
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
  
  # ### Set wd
  # setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/Aurum_Jun2021_extract/")
  # getwd()
  # 
  # ### Source functions
  # ### Source functions
  # R.func.sources = list.files("R", full.names = TRUE)
  # sapply(R.func.sources, source)
  # cohortZ <- readRDS(paste("data/extraction/cohortZ.rds", sep = ""))
  # out.subdir <- "cohortZ"
  # db <- "aurum_small"
  # varname <- "sbp_unexposed"
  # indexdt <- "fup_start"
  # t <- 0
  # cohort <- cohortZ
  # codelist = "edh_sbp_medcodeid"
  # t.varname = TRUE
  # time.prev = 365.25*5
  # time.post = 0
  # lower.bound = -Inf
  # upper.bound = Inf
  # db.filepath = NULL
  # out.save.disk = FALSE
  # out.filepath = NULL
  # out.subdir = NULL
  # return.output = TRUE
  # codelist.med <- "edh_sbp_medcodeid"
  # codelist.drug <- "uom_antihypertensives_prodcodeid"
  # 
  # keep.exposed = TRUE
  # keep.metadata = TRUE
  # query.med = NULL
  # query.drug = NULL
  # task.id = NULL
  # task.id.n = NULL
  # range.min <- 14
  # range.max <- 180
  
  ### Preparation
  ## Match arg for exposed.type
  exposed.type <- match.arg(exposed.type)
  
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
  
  ### If task.id has been specified, we want to reduce the cohort to just these individuals.
  if (!is.null(task.id)){
    # Create groups for splitting
    cohort$groups <- cut(1:nrow(cohort), task.id.n, labels = FALSE)
    print(table(cohort$groups))
    
    # Reduce cohort to group of interest
    cohort <- cohort[!is.na(fastmatch::fmatch(cohort$group, task.id)), ]
    print(str(cohort))
    
  }
  
  ### Run database queries and combine with cohort if they haven't been specified in input
  ### Med
  if (is.null(query.med)){
    # Run query
    db.qry.med <- db_query(codelist.med,
                           db = db, 
                           db.filepath = db.filepath,
                           tab = "obs")
    # Combine query with cohort
    variable.dat.med <- combine_query(cohort, 
                                      db.qry.med, 
                                      query.type = "test",
                                      time.prev = time.prev, 
                                      time.post = time.post, 
                                      lower.bound = lower.bound, 
                                      upper.bound = upper.bound, 
                                      numobs = Inf)
    
  } else {
    # Assign query from input
    variable.dat.med <- query.med
    
    # Reduce db.query to group of interest
    if (!is.null(task.id)){
      variable.dat.med <- variable.dat.med[!is.na(fastmatch::fmatch(variable.dat.med$patid, cohort$patid)), ]
      print(str(variable.dat.med))
    }
  }
  
  ### Drug
  if (is.null(query.drug)){
    # Run query
    db.qry.drug <- db_query(codelist.drug,
                            db = db, 
                            db.filepath = db.filepath,
                            tab = "drug")
    # Combine query with cohort
    variable.dat.drug <- combine_query(cohort, 
                                       db.qry.drug, 
                                       query.type = "drug",
                                       time.prev = Inf, 
                                       time.post = time.post, 
                                       numobs = Inf)
  } else {
    # Assign query from input
    variable.dat.drug <- query.drug
    
    # Reduce db.query to group of interest
    if (!is.null(task.id)){
      variable.dat.drug <- variable.dat.drug[!is.na(fastmatch::fmatch(variable.dat.drug$patid, cohort$patid)), ]
      print(str(variable.dat.drug))
    }
  }
  
  ### Group data
  variable.dat.med <- dplyr::arrange(variable.dat.med, patid) |>
    dplyr::group_by(patid)
  
  ### Run function to remove exposed observations
  variable.dat.med.unexposed <- dplyr::group_modify(.data = variable.dat.med, .f = remove_exposed, 
                                                    range.max = range.max, range.min = range.min, range.num = range.num, data.pres = variable.dat.drug, 
                                                    keep.exposed = keep.exposed, keep.metadata = TRUE)
  
  ### Reduce to the most recent observation, or if keep.exposed == TRUE, options to keep the one with most number of prescriptions in the window,
  ### or the minimum difference from the prescription.
  if (exposed.type == "all"){
    variable.dat.med.unexposed <- variable.dat.med.unexposed |>
      dplyr::arrange(patid, desc(obsdate))
  } else {
    if (keep.exposed == FALSE){
      variable.dat.med.unexposed <- variable.dat.med.unexposed |>
        dplyr::arrange(patid, desc(obsdate)) |>
        dplyr::filter(dplyr::row_number(patid) == 1)
    } else {
      if (exposed.type == "mr"){
        variable.dat.med.unexposed <- variable.dat.med.unexposed |>
          dplyr::arrange(patid, desc(obsdate)) |>
          dplyr::filter(dplyr::row_number(patid) == 1)
      } else if (exposed.type == "np"){
        variable.dat.med.unexposed <- variable.dat.med.unexposed |>
          dplyr::arrange(patid, desc(np), desc(obsdate)) |>
          dplyr::filter(dplyr::row_number(patid) == 1)
      } else if (exposed.type == "md"){
        variable.dat.med.unexposed <- variable.dat.med.unexposed |>
          dplyr::arrange(patid, md, desc(obsdate)) |>
          dplyr::filter(dplyr::row_number(patid) == 1)
      }
    }
  }
  
  rm(variable.dat.med, variable.dat.drug, db.qry.med, db.qry.drug)
  
  ### Create data frame of cohort and the variable of interest
  variable.dat <- merge(dplyr::select(cohort, patid), variable.dat.med.unexposed, by.x = "patid", by.y = "patid", all.x = TRUE)
  
  ### Reduce to variables of interest
  if (keep.exposed == TRUE){
    variable.dat <- variable.dat[,c("patid", "value", "obsdate", "np", "md", "medcodeid")]
  } else if (keep.exposed == FALSE){
    variable.dat <- variable.dat[,c("patid", "value", "obsdate", "medcodeid")]
  }
  
  ### Change name of variable to varname
  colnames(variable.dat)[colnames(variable.dat) == "value"] <- varname
  
  implement_output_taskid(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output, task.id)
  
}