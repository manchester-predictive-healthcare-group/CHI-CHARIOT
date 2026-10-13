### Write a function which will extract cholesterol/hdl values data, but remove observations if there is a prescription
### within a specified number of days prior to the test
extract_cholhdl_ratio_unexposed <- function(cohort, 
                                            varname = NULL,
                                            codelist.ratio,
                                            codelist.chol,
                                            codelist.hdl,
                                            codelist.drug,
                                            query.ratio = NULL,
                                            query.chol = NULL,
                                            query.hdl = NULL,
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
  ### Preparation
  ## Match arg for exposed.type
  exposed.type <- match.arg(exposed.type)
  print(paste("exposed.type = ", exposed.type))
  ## Add index date variable to cohort and change indexdt based on t
  cohort <- prep_cohort(cohort, indexdt, t)
  ## Assign variable name if unspecified
  if (is.null(varname)){
    varname <- "cholhdl_ratio"
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
  
  ### Run queries and combine with cohort if they haven't been inputted as function arguments
  
  ### Ratio
  if (is.null(query.ratio)){
    # Run query
    db.qry.ratio <- db_query(codelist.ratio,
                             db = db, 
                             db.filepath = db.filepath,
                             tab = "obs")
    # Combine query with cohort
    variable.dat.ratio <- combine_query(cohort, 
                                        db.qry.ratio, 
                                        query.type = "test",
                                        time.prev = time.prev, 
                                        time.post = time.post, 
                                        lower.bound = lower.bound, 
                                        upper.bound = upper.bound, 
                                        numobs = Inf)
  } else {
    # Assign query from input
    variable.dat.ratio <- query.ratio
    
    # Reduce db.query to group of interest
    if (!is.null(task.id)){
      variable.dat.ratio <- variable.dat.ratio[!is.na(fastmatch::fmatch(variable.dat.ratio$patid, cohort$patid)), ]
      print(str(variable.dat.ratio))
    }
    
  }
  
  ### cholesterol
  if (is.null(query.chol)){
    # Run query
    db.qry.chol <- db_query(codelist.chol,
                            db = db, 
                            db.filepath = db.filepath,
                            tab = "obs")
    # Combine query with cohort
    variable.dat.chol <- combine_query(cohort, 
                                       db.qry.chol, 
                                       query.type = "test",
                                       time.prev = time.prev, 
                                       time.post = time.post, 
                                       numobs = Inf)
  } else {
    # Assign query from input
    variable.dat.chol <- query.chol
    
    # Reduce db.query to group of interest
    if (!is.null(task.id)){
      variable.dat.chol <- variable.dat.chol[!is.na(fastmatch::fmatch(variable.dat.chol$patid, cohort$patid)), ]
      print(str(variable.dat.chol))
    }
  }
  
  ### HDL
  if (is.null(query.hdl)){
    # Run query
    db.qry.hdl <- db_query(codelist.hdl,
                           db = db, 
                           db.filepath = db.filepath,
                           tab = "obs")
    # Combine query with cohort
    variable.dat.hdl <- combine_query(cohort, 
                                      db.qry.hdl, 
                                      query.type = "test",
                                      time.prev = time.prev, 
                                      time.post = time.post, 
                                      numobs = Inf)
  } else {
    # Assign query from input
    variable.dat.hdl <- query.hdl
    
    # Reduce db.query to group of interest
    if (!is.null(task.id)){
      variable.dat.hdl <- variable.dat.hdl[!is.na(fastmatch::fmatch(variable.dat.hdl$patid, cohort$patid)), ]
      print(str(variable.dat.hdl))
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
  variable.dat.ratio <- dplyr::arrange(variable.dat.ratio, patid) |>
    dplyr::group_by(patid)
  
  variable.dat.chol <- dplyr::arrange(variable.dat.chol, patid) |>
    dplyr::group_by(patid)
  
  variable.dat.hdl <- dplyr::arrange(variable.dat.hdl, patid) |>
    dplyr::group_by(patid)
  
  ### Run function to remove exposed observations for ratio and total cholesterol
  variable.dat.ratio.unexposed <- dplyr::group_modify(.data = variable.dat.ratio, .f = remove_exposed, 
                                                      range.max = range.max, range.min = range.min, range.num = range.num, data.pres = variable.dat.drug, 
                                                      keep.exposed = keep.exposed, keep.metadata = TRUE)
  
  variable.dat.chol.unexposed <- dplyr::group_modify(.data = variable.dat.chol, .f = remove_exposed, 
                                                     range.max = range.max, range.min = range.min, range.num = range.num, data.pres = variable.dat.drug, 
                                                     keep.exposed = keep.exposed, keep.metadata = TRUE)
  
  ### For ratio, and total cholesterol:
  ### Reduce to the most recent observation, or if keep.exposed == TRUE, 
  ### options to keep the one with most number of prescriptions in the window, or the minimum difference from the prescription.
  
  ### NEW
  if (exposed.type == "all"){
    variable.dat.ratio.unexposed <- variable.dat.ratio.unexposed |>
      dplyr::arrange(patid, desc(obsdate))
  } else {
    if (keep.exposed == FALSE){
      variable.dat.ratio.unexposed <- variable.dat.ratio.unexposed |>
        dplyr::arrange(patid, desc(obsdate)) |>
        dplyr::filter(dplyr::row_number(patid) == 1)
    } else {
      if (exposed.type == "mr"){
        variable.dat.ratio.unexposed <- variable.dat.ratio.unexposed |>
          dplyr::arrange(patid, desc(obsdate)) |>
          dplyr::filter(dplyr::row_number(patid) == 1)
      } else if (exposed.type == "np"){
        variable.dat.ratio.unexposed <- variable.dat.ratio.unexposed |>
          dplyr::arrange(patid, desc(np), desc(obsdate)) |>
          dplyr::filter(dplyr::row_number(patid) == 1)
      } else if (exposed.type == "md"){
        variable.dat.ratio.unexposed <- variable.dat.ratio.unexposed |>
          dplyr::arrange(patid, md, desc(obsdate)) |>
          dplyr::filter(dplyr::row_number(patid) == 1)
      }
    }
  }
  
  if (exposed.type == "all"){
    variable.dat.chol.unexposed <- variable.dat.chol.unexposed |>
      dplyr::arrange(patid, desc(obsdate))
  } else {
    if (keep.exposed == FALSE){
      variable.dat.chol.unexposed <- variable.dat.chol.unexposed |>
        dplyr::arrange(patid, desc(obsdate)) |>
        dplyr::filter(dplyr::row_number(patid) == 1)
    } else {
      if (exposed.type == "mr"){
        variable.dat.chol.unexposed <- variable.dat.chol.unexposed |>
          dplyr::arrange(patid, desc(obsdate)) |>
          dplyr::filter(dplyr::row_number(patid) == 1)
      } else if (exposed.type == "np"){
        variable.dat.chol.unexposed <- variable.dat.chol.unexposed |>
          dplyr::arrange(patid, desc(np), desc(obsdate)) |>
          dplyr::filter(dplyr::row_number(patid) == 1)
      } else if (exposed.type == "md"){
        variable.dat.chol.unexposed <- variable.dat.chol.unexposed |>
          dplyr::arrange(patid, md, desc(obsdate)) |>
          dplyr::filter(dplyr::row_number(patid) == 1)
      }
    }
  }
  
  # ### Ratio
  # if (keep.exposed == FALSE){
  #   variable.dat.ratio.unexposed <- variable.dat.ratio.unexposed |>
  #     dplyr::arrange(patid, desc(obsdate)) |>
  #     dplyr::filter(dplyr::row_number(patid) == 1)
  # } else {
  #   if (exposed.type == "mr"){
  #     variable.dat.ratio.unexposed <- variable.dat.ratio.unexposed |>
  #       dplyr::arrange(patid, desc(obsdate)) |>
  #       dplyr::filter(dplyr::row_number(patid) == 1)
  #   } else if (exposed.type == "np"){
  #     variable.dat.ratio.unexposed <- variable.dat.ratio.unexposed |>
  #       dplyr::arrange(patid, desc(np), desc(obsdate)) |>
  #       dplyr::filter(dplyr::row_number(patid) == 1)
  #   } else if (exposed.type == "md"){
  #     variable.dat.ratio.unexposed <- variable.dat.ratio.unexposed |>
  #       dplyr::arrange(patid, md, desc(obsdate)) |>
  #       dplyr::filter(dplyr::row_number(patid) == 1)
  #   }
  # }
  
  # ### Total cholesterol
  # if (keep.exposed == FALSE){
  #   variable.dat.chol.unexposed <- variable.dat.chol.unexposed |>
  #     dplyr::arrange(patid, desc(obsdate)) |>
  #     dplyr::filter(dplyr::row_number(patid) == 1)
  # } else {
  #   if (exposed.type == "mr"){
  #     variable.dat.chol.unexposed <- variable.dat.chol.unexposed |>
  #       dplyr::arrange(patid, desc(obsdate)) |>
  #       dplyr::filter(dplyr::row_number(patid) == 1)
  #   } else if (exposed.type == "np"){
  #     variable.dat.chol.unexposed <- variable.dat.chol.unexposed |>
  #       dplyr::arrange(patid, desc(np), desc(obsdate)) |>
  #       dplyr::filter(dplyr::row_number(patid) == 1)
  #   } else if (exposed.type == "md"){
  #     variable.dat.chol.unexposed <- variable.dat.chol.unexposed |>
  #       dplyr::arrange(patid, md, desc(obsdate)) |>
  #       dplyr::filter(dplyr::row_number(patid) == 1)
  #   }
  # }

  ### For HDL, just keep most recent, because we have not applied any exposed/criteria
  variable.dat.hdl<- variable.dat.hdl |>
    dplyr::arrange(patid, obsdate) |>
    dplyr::filter(dplyr::row_number(dplyr::desc(obsdate)) == 1)
  
  ###
  ### If extracting all, concatenate these datasets and prepare to be outputted
  ###
  if (exposed.type == "all"){
    ### Rbind
    variable.dat <- rbind(variable.dat.hdl, variable.dat.chol.unexposed, variable.dat.ratio.unexposed, fill = TRUE)
    ### Reduce to appropriate variables
    variable.dat[,c("patid", "value", "obsdate", "np", "md", "medcodeid")]
  } else {
    ### Calculate ratio's estimated from chol/hdl
    ## Merge chol and hdl datasets
    variable.dat.manual <- merge(variable.dat.chol.unexposed, variable.dat.hdl, by.x = "patid", by.y = "patid")
    ## Calculate ratio
    variable.dat.manual$value <- variable.dat.manual$value.x/variable.dat.manual$value.y
    ## Take furthest away date
    variable.dat.manual$obsdate <- pmin(variable.dat.manual$obsdate.x, variable.dat.manual$obsdate.y)
    ## Remove values outside of range
    ### If values are missing, < lower.bound or > upper.bound then delete
    if (!is.null(lower.bound) & !is.null(upper.bound)){
      variable.dat.manual <- subset(variable.dat.manual, value > lower.bound & value < upper.bound)
    } else if (is.null(lower.bound) & !is.null(upper.bound)){
      variable.dat.manual <- subset(variable.dat.manual, value < upper.bound)
    } else if (!is.null(lower.bound) & is.null(upper.bound)){
      variable.dat.manual <- subset(variable.dat.manual, value > lower.bound)
    }
    variable.dat.manual <- variable.dat.manual[,c("patid", "value", "obsdate")]
    rm(variable.dat.chol, variable.dat.hdl)
    
    ### Merge the two
    variable.dat <- merge(variable.dat.ratio.unexposed, variable.dat.manual, by.x = "patid", by.y = "patid", all.x = TRUE, all.y = TRUE)
    
    ### Take most recent of the two
    variable.dat <- dplyr::mutate(variable.dat, cholhdl_ratio_unexposed = dplyr::case_when(is.na(value.x) & !is.na(value.y) ~ value.y,
                                                                                           !is.na(value.x) & is.na(value.y) ~ value.x,
                                                                                           !is.na(value.x) & !is.na(value.y) & obsdate.y > obsdate.x ~ value.y,
                                                                                           !is.na(value.x) & !is.na(value.y) & obsdate.y <= obsdate.x ~ value.x))
    
    ### Create dataframe of cohort and the variable of interest
    variable.dat <- merge(dplyr::select(cohort, patid), variable.dat, by.x = "patid", by.y = "patid", all.x = TRUE)
    
    ### Reduce to variables of interest
    variable.dat <- variable.dat[,c("patid", "cholhdl_ratio_unexposed")]
    
    ### Change name of variable to varname
    colnames(variable.dat)[colnames(variable.dat) == "cholhdl_ratio_unexposed"] <- varname
    
  }

  ### Implement output
  implement_output_taskid(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output, task.id)
  
}