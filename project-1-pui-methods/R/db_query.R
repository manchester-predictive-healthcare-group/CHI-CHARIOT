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