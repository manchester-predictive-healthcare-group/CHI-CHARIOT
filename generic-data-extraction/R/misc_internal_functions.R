###
### Write a function which will add the data from one of the CPRD flatfiles, to a sql database, specified via the filename argument.
###
add_to_database <- function(filename, filetype, nrows = -1, select = NULL, subset.patids = NULL, use.set = FALSE, db, ...){
  
  ### Use the extract_txt function for the relevant filetype
  extract_txt_func <- list("obs" = extract_txt_obs,
                           "drug" = extract_txt_drug,
                           "ref" = extract_txt_ref,
                           "prob" = extract_txt_prob,
                           "cons" = extract_txt_cons,
                           "hes_primary" = extract_txt_hes_primary,
                           "hes_patient" = extract_txt_char,
                           "hes_hosp" = extract_txt_char,
                           "IMD" = extract_txt_char,
                           "death" = extract_txt_death)[[filetype]]
  
  ### Read data into R and turn file into data.table for quicker subsetting
  ext.dat <- extract_txt_func(filepath = filename, nrows = nrows, select = select)
  ext.dat <- data.table::as.data.table(ext.dat)
  
  ### Reduce to only patids from cohort of interest
  if (!is.null(subset.patids)){
    ### Turn into data.table
    subset.patids <- data.table::as.data.table(subset.patids)
    
    ### Reduce subset.patids input to just those from the appropriate 'set' if the filenames are compatible
    ### This will make subsetting ext.dat to the patients %in% subset.patids a lot quicker
    if (use.set == TRUE){
      ### Extract the 'set' from the filename
      set.filename <- as.numeric(stringr::str_match(filename, "set\\s*(.*?)\\s*_")[,2])
      ### Apply subsetting
      subset.patids <- subset.patids[set == set.filename, patid]
    }
    ### Apply subsettings
    ext.dat <- ext.dat[patid %in% subset.patids]
  }
  
  print("START WRITE TABLE")
  ### Add to sqlite database
  RSQLite::dbWriteTable(db, filetype, ext.dat, ...)
  print("FINISHED WRITE TABLE")
  
}

###
### A function to extract and create an sqlite database for all files of a particular file type
###
cprd_extract <- function(dbname, filetype, nrows = -1, select = NULL, subset.patids = NULL, use.set = FALSE){
  
  print(paste("nrows", nrows))
  print(paste("select", select))
  print(paste("use.set", use.set))
  
  ### Spit error out if subset.patids not in correct form
  if (!is.null(subset.patids)){
    if (use.set == FALSE){
      if (!is.vector(subset.patids)){
        stop("When use.set = FALSE, subset.patids should be a vector")
      }
    } else if (use.set == TRUE){
      if (!("patid" %in% colnames(subset.patids) & "set" %in% colnames(subset.patids))){
        stop("When use.set = TRUE, subset.patids should contain columns 'patid' and 'set'")
      }
    }
  }
  
  ### Start by creating/joining a database
  mydb <- RSQLite::dbConnect(RSQLite::SQLite(), paste("data/sql/", dbname, ".sqlite", sep = ""))
  
  ### Get filenames
  directory <- c("cons" = "Consultation", "obs" = "Observation", "drug" = "DrugIssue", "ref" = "Referral", "prob" = "Problem")
  filenames <- list.files(paste("data/unzip/", directory[filetype], sep = ""), pattern = ".txt", full.names = TRUE)
  
  ### Apply the add_to_database function to each file
  ## Overwrite for first file
  add_to_database(filenames[1], 
                  filetype = filetype, 
                  nrows = nrows,
                  select = select,
                  subset.patids = subset.patids,
                  use.set = use.set,
                  db = mydb, 
                  overwrite = TRUE)
  ## Append for all subsequent files
  for (filename in filenames[-1]){
    print(paste(filename, Sys.time()))
    add_to_database(filename, 
                    filetype = filetype,
                    nrows = nrows,
                    select = select,
                    subset.patids = subset.patids,
                    use.set = use.set,
                    db = mydb, 
                    append = TRUE)
  }
  
  ### Disconnect
  RSQLite::dbDisconnect(mydb)
  
}

###
### Prep cohort inputted into a variable extraction function
###
prep_cohort <- function(cohort, indexdt, t, reduce = TRUE){
  
  ### Change name of indexdt variable to "indexdt" so we can easily refer to it
  colnames(cohort)[colnames(cohort) == indexdt] <- "indexdt"
  
  ### Reduce cohort to variables of interst
  if (reduce == TRUE){
    cohort <- cohort[,c("patid", "indexdt")]
  }

  ### If t has been specified, round it to nearest integer, add to indexdt to aextract data at appropriate time
  if (!is.null(t)){
    t <- round(t)
    cohort$indexdt <- cohort$indexdt + t
  }
  
  return(cohort)
  
}

###
### Add 't' to varname
###
prep_varname <- function(varname, t, t.varname){
  if (!is.null(t)){
    if (t.varname == TRUE){
      varname <- paste(varname, "_t", t, sep = "")
    }
  }
  return(varname)
}

###
### If out.subdir has been specified, check if it exists. If not, create directory
###
prep_subdir <- function(out.subdir){
  if (!is.null(out.subdir)){
    if (file.exists(paste("data/extraction/", out.subdir, sep = "")) == FALSE){
      dir.create(paste("data/extraction/", out.subdir, sep = ""))
      print(paste("Directory /data/extraction/", out.subdir, " has been created", sep = ""))
    }
  }
}


###
### Function to format output for variable extraction functions
### This has no utility outside of these variable extraction functions
###
implement_output <- function(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output){
  ### We now have the final variable file, we just need to decide how and where to save it
  ### Save data frame to disc
  if (out.save.disk == TRUE){
    ## Use standard naming convention if out.filepath is not specified
    if (is.null(out.filepath)){
      ## Save directly into /extraction if a subdirectory is not specified
      if (is.null(out.subdir)){
        saveRDS(variable.dat, paste("data/extraction/var_", varname, ".rds", sep = ""))
        ## Save into the specified subdirectory if specified
      } else if (!is.null(out.subdir)){
        saveRDS(variable.dat, paste("data/extraction/", out.subdir, "/var_", varname, ".rds", sep = ""))
      }  
      ## Use the full user specified filepath if specified
    } else if (!is.null(out.filepath)){
      saveRDS(variable.dat, out.filepath)
    }
  }
  
  ### Return as output of function if requested
  if (return.output == TRUE){
    return(variable.dat)
  }
  
}


###
### Function to format output for variable extraction functions
### This has no utility outside of these variable extraction functions
### This function also allows specification of task.id, for when parallelising
###
implement_output_taskid <- function(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output, task.id = NULL){
  ### We now have the final variable file, we just need to decide how and where to save it
  ### Save data frame to disc
  if (out.save.disk == TRUE){
    ## Use standard naming convention if out.filepath is not specified
    if (is.null(out.filepath)){
      ## Save directly into /extraction if a subdirectory is not specified
      if (is.null(out.subdir)){
        if (!is.null(task.id)){
          saveRDS(variable.dat, paste("data/extraction/var_", varname, "_id", task.id, ".rds", sep = ""))
        } else {
          saveRDS(variable.dat, paste("data/extraction/var_", varname, ".rds", sep = ""))
        }
        ## Save into the specified subdirectory if specified
      } else if (!is.null(out.subdir)){
        if (!is.null(task.id)){
          saveRDS(variable.dat, paste("data/extraction/", out.subdir, "/var_", varname, "_id", task.id, ".rds", sep = ""))
        } else {
          saveRDS(variable.dat, paste("data/extraction/", out.subdir, "/var_", varname, "_id", task.id, ".rds", sep = ""))
        }
      }  
      ## Use the full user specified filepath if specified
    } else if (!is.null(out.filepath)){
      saveRDS(variable.dat, out.filepath)
    }
  }
  
  ### Return as output of function if requested
  if (return.output == TRUE){
    return(variable.dat)
  }
  
}


###
### Write a function for querying an SQLite database
###
db_query <- function(codelist = NULL,
                     db.open = NULL,
                     db = NULL,
                     db.filepath = NULL,
                     tab = c("obs", "drug", "hes_primary", "death"),
                     codelist.vector = NULL){
  
  ### Extract codelist
  if (is.null(codelist.vector)){
    if (!is.null(codelist)){
      codelist <- data.table::fread(file = paste(getwd(),"/codelists/analysis/", codelist, ".csv", sep = ""),
                                    sep = ",", header = TRUE, colClasses = "character")
      if (tab == "obs"){
        codelist <- codelist$medcodeid
      } else if (tab == "drug"){
        codelist <- codelist$prodcodeid
      } else if (tab %in% c("hes_primary", "death")){
        codelist <- codelist$ICD10
      }
    }
  } else {
    codelist <- codelist.vector
  }
  
  ### Connect to SQLite database
  if (is.null(db.open)){
    if (!is.null(db)){
      mydb <- RSQLite::dbConnect(RSQLite::SQLite(), paste("data/sql/", db, ".sqlite", sep = ""))
    } else if (!is.null(db.filepath)){
      mydb <- RSQLite::dbConnect(RSQLite::SQLite(), db.filepath)
    }
  } else {
    mydb <- db.open
  }
  
  ### Create the query
  if (!is.null(codelist)){
    if (tab == "obs"){
      where_clause <- paste0("`medcodeid` IN (", paste("'", codelist, "'", sep = "", collapse = ","), ")")
      qry <- paste("SELECT * FROM", tab, "WHERE", where_clause)
    } else if (tab == "drug"){
      where_clause <- paste0("`prodcodeid` IN (", paste("'", codelist, "'", sep = "", collapse = ","), ")")
      qry <- paste("SELECT * FROM", tab, "WHERE", where_clause)
    } else if (tab == "hes_primary"){
      where_clause <- paste0("`ICD_PRIMARY` IN (", paste("'", codelist, "'", sep = "", collapse = ","), ")")
      qry <- paste("SELECT * FROM", tab, "WHERE", where_clause)
    } else if (tab == "death"){
      where_clause <- paste0("`cause` IN (", paste("'", codelist, "'", sep = "", collapse = ","), ")")
      qry <- paste("SELECT * FROM", tab, "WHERE", where_clause)
    }
  } else if (is.null(codelist)){
    qry <- paste("SELECT * FROM", tab)
  }
 
  ### Run the query and turn into data.table
  db.query <- RSQLite::dbGetQuery(mydb, qry)
  db.query <- data.table::as.data.table(db.query)
  
  ### Disconnect
  if (is.null(db.open)){
    RSQLite::dbDisconnect(mydb)
  }
  
  ### Return query
  return(db.query)
  
}


###
### ### Write a function for combining a database query with a cohort, and returing a 0/1 vector, for whether there is a medical or drug record
### prior to the index date for individuals in cohort.
###
combine_query_boolean <- function(cohort, 
                                  db.query,
                                  query.type = c("obs", "drug", "hes_primary"),
                                  time.prev = Inf, 
                                  time.post = 0, 
                                  numobs = 1){
  
  ### Merge cohort with the database query keeping observations that are in both 
  cohort.qry <- merge(cohort, db.query, by.x = "patid", by.y = "patid")
  
  ### Reduce to variables of interest
  if (query.type == c("obs")){
    cohort.qry <- cohort.qry[,c("patid", "indexdt", "obsdate")]
  } else if (query.type == "drug"){
    cohort.qry <- cohort.qry[,c("patid", "indexdt", "issuedate")]
    ## rename issuedate to obsdate so we can use same code for medical or drug queries
    colnames(cohort.qry)[colnames(cohort.qry) == "issuedate"] <- "obsdate"
  } else if (query.type == "hes_primary"){
    cohort.qry <- cohort.qry[,c("patid", "indexdt", "admidate")]
    ## rename issuedate to obsdate so we can use same code for medical or drug queries
    colnames(cohort.qry)[colnames(cohort.qry) == "admidate"] <- "obsdate"
  }
  
  ### Remove values outside of specified time range
  cohort.qry <- cohort.qry[obsdate <= indexdt + time.post & obsdate > indexdt - time.prev]
  print("COHORT QUERY 1")
  print(str(cohort.qry))
  
  ### Identify which patients have had the required number of events in the specified time period
  cohort.qry <- cohort.qry |>
    dplyr::group_by(patid) |>
    dplyr::summarise(total = dplyr::n()) |>
    dplyr::filter(total >= numobs) |>
    dplyr::pull(patid)
  print("COHORT QUERY 2")
  print(str(cohort.qry))
  
  ### Create a 0/1 vector for whether patients are in this list and add to cohort
  print(paste("create variable", Sys.time()))
  cohort.qry <- as.integer(!is.na(fastmatch::fmatch(cohort$patid, cohort.qry)))
  print("COHORT QUERY 3")
  print(sum(cohort.qry))
  
  ### Return this vector
  return(cohort.qry)
  
}



###
### Write a function for combining a database query with a cohort, and only retain observations within an
### appropriate time frame, and lower/upper bounds for test data
###
combine_query <- function(cohort,
                          db.query,
                          query.type = c("med", "test", "drug", "hes_primary", "death"),
                          time.prev = Inf,
                          time.post = Inf,
                          lower.bound = -Inf,
                          upper.bound = Inf,
                          remove.zeros = FALSE,
                          numobs = 1,
                          value.na.rm = TRUE,
                          earliest.values = FALSE,
                          reduce.output = TRUE){
  
  ### Merge cohort with the database query keeping observations that are in both
  cohort.qry <- merge(cohort, db.query, by.x = "patid", by.y = "patid")
  cohort.qry <- data.table::as.data.table(cohort.qry)
  
  ### Rename event date variables to all be "obsdate"
  if (query.type == c("med")){
    
  } else if (query.type == "test"){
    
  } else if (query.type == "drug"){
    ## rename issuedate to obsdate so we can use same naming for all queries
    colnames(cohort.qry)[colnames(cohort.qry) == "issuedate"] <- "obsdate"
  } else if (query.type == "hes_primary"){
    ## rename admidate to obsdate so we can use same naming for all queries
    colnames(cohort.qry)[colnames(cohort.qry) == "admidate"] <- "obsdate"
  }  else if (query.type == "death"){
    ## rename dod to obsdate so we can use same naming for all queries
    colnames(cohort.qry)[colnames(cohort.qry) == "dod"] <- "obsdate"
  }
  
  ### Remove values outside of specified time range
  cohort.qry <- cohort.qry[obsdate <= indexdt + time.post & obsdate > indexdt - time.prev]
  
  ### For test data, remove values outside of value range
  if (query.type == "test"){
    ### If values are missing, < lower.bound or > upper.bound then delete
    if (value.na.rm == TRUE){
      cohort.qry <- cohort.qry[!is.na(value) & value > lower.bound & value < upper.bound]
    }
    ### If specified, remove zeros
    if (remove.zeros == TRUE){
      cohort.qry <- cohort.qry[value != 0]
    }
  }
  
  ### Reduce to variables of interest
  if (reduce.output == TRUE){
    if (query.type == c("med")){
      cohort.qry <- cohort.qry[,c("patid", "medcodeid", "obsdate")]
    } else if (query.type == "test"){
      cohort.qry <- cohort.qry[,c("patid", "medcodeid", "obsdate", "value", "numunitid", "numrangelow", "numrangehigh")]
    } else if (query.type == "drug"){
      cohort.qry <- cohort.qry[,c("patid", "prodcodeid", "obsdate")]
    } else if (query.type == "hes_primary"){
      cohort.qry <- cohort.qry[,c("patid", "obsdate")]
    }  else if (query.type == "death"){
      cohort.qry <- cohort.qry[,c("patid", "obsdate")]
    }
  }
  
  ### Group by patid and obsdate and keep the most recent 'numobs' number of observations
  if (earliest.values == FALSE){
    cohort.qry <- dplyr::group_by(cohort.qry, patid) |>
      dplyr::filter(dplyr::row_number(dplyr::desc(obsdate)) <= numobs) |>
      data.table::as.data.table()
  } else if (earliest.values == TRUE){
    cohort.qry <- dplyr::group_by(cohort.qry, patid) |>
      dplyr::filter(dplyr::row_number(obsdate) <= numobs) |>
      data.table::as.data.table()
  }
  
  ### Return cohort.qry
  return(cohort.qry)
  
}

###
### Write a function to remove test data observations, if there is a prescription prior to this test data observation.
### Prescription will only exclude the observation if recorded within a specified window before the test, where is it
### deemed the prescription may have impacted the test result.
### Writing a function that can be utilised with dplyr::group_modify and dplyr::group_map
###
remove_exposed <- function(data, id, range.max = 180, range.min = 0, range.num = 1, data.pres, keep.exposed = FALSE, keep.metadata = TRUE){
  
  ### Reduce prescription data to just patient of interest
  data.pres <- data.pres[!is.na(fastmatch::fmatch(data.pres$patid, id)), ]
  
  ### Identify which row/medical observations are exposed 
  exposed <- unlist(lapply(data$obsdate, function(x) {
    ### calculate diffs
    diffs <- x - data.pres$obsdate
    ### check if any are within range 14 days to 180 days
    (sum(as.numeric(range.min < diffs & diffs < range.max)) >= range.num)
  }))
  
  ### Reduce data by removing or keeping these observations
  ### If keep.exposed == TRUE, and specified to do so, extract number of prescriptions in the window, 
  ### and the minimum distance from a prescription in the window
  if (keep.exposed == TRUE){
    data <- data[exposed, ]
    if (keep.metadata == TRUE){
      ### Get number of prescriptions within window
      np <- unlist(lapply(data$obsdate, function(x) {
        ### calculate diffs
        diffs <- x - data.pres$obsdate
        ### get number within range 14 days to 180 days
        sum(range.min < diffs & diffs < range.max)
      }))
      data$np <- np
      
      ### Get the min distance between observations and prescriptions
      md <- unlist(lapply(data$obsdate, function(x) {
        ### calculate diffs
        diffs <- x - data.pres$obsdate
        ### only keep diffs in time range
        diffs <- diffs[range.min < diffs & diffs < range.max]
        ### keep the minimum (suppressing warnings because some will be NULL)
        suppressWarnings(min(diffs))
      }))
      data$md <- md
    }
  } else if (keep.exposed == FALSE){
    data <- data[!exposed, ]
  }
  
  return(data)
  
}
