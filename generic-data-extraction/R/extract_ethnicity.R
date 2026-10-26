### This program will extract Ethnicity from primary care records (note this will also be recorded in secondary care and have to be extracted from there)

### No need for codelist.in because the code lists for this variable are fixed/specifc
### No need for indexdt, t and t.varname as we are looking over entire patient follow up
extract_ethnicity <- function(cohort, 
                              varname = NULL, 
                              db, 
                              db.linked,
                              db.filepath = NULL,
                              out.save.disk = TRUE, 
                              out.subdir = NULL,
                              out.filepath = NULL, 
                              return.output = FALSE){
  
  ### Preparation
  ## Reduce cohort to important variable
  cohort <- cohort[, "patid"]
  ## Assign variable name if unspecified
  if (is.null(varname)){
    varname <- "ethnicity"
  }
  ## Create named subdirectory if it doesn't exist
  prep_subdir(out.subdir)
  
  #################################################
  ### EXTRACT ETHNICITY USING PRIMARY CARE DATA ###
  #################################################
  
  ### There are a large number of codelists for this variable, extract them all
  ### Note that the user is not specifying the codelists, so they must have the appropriate codelists with names in the codelists/analysis folder
  ### to be able to run this function
  
  ### bangladeshi ADDED
  edh.bangladeshi <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_bangladeshi_medcodeid.csv", sep = ""), 
                                       sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.bangladeshi <- edh.bangladeshi$medcodeid
  
  ### edh.black.african ADDED
  edh.black.african <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_black_african_medcodeid.csv", sep = ""), 
                                         sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.black.african <- edh.black.african$medcodeid
  
  ### edh.black.caribbean ADDED
  edh.black.caribbean <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_caribbean_medcodeid.csv", sep = ""), 
                                           sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.black.caribbean <- edh.black.caribbean$medcodeid
  
  ### edh.chinese ADDED
  edh.chinese <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_chinese_medcodeid.csv", sep = ""), 
                                   sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.chinese <- edh.chinese$medcodeid
  
  ### edh.indian ADDED
  edh.indian <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_indian_medcodeid.csv", sep = ""), 
                                  sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.indian <- edh.indian$medcodeid
  
  ### edh.irish ADDED
  edh.irish <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_irish_medcodeid.csv", sep = ""), 
                                 sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.irish <- edh.irish$medcodeid
  
  ### edh.not.recorded ADDED
  edh.not.recorded <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_not_recorded_medcodeid.csv", sep = ""), 
                                        sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.not.recorded <- edh.not.recorded$medcodeid
  
  ### edh.oth.asian ADDED
  edh.oth.asian <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_oth_asian_medcodeid.csv", sep = ""), 
                                     sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.oth.asian <- edh.oth.asian$medcodeid
  
  ### edh.oth.black ADDED
  edh.oth.black <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_oth_black_medcodeid.csv", sep = ""), 
                                     sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.oth.black <- edh.oth.black$medcodeid
  
  ### edh.oth.ethnic ADDED
  edh.oth.ethnic <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_oth_ethnic_medcodeid.csv", sep = ""), 
                                      sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.oth.ethnic <- edh.oth.ethnic$medcodeid
  
  ### edh.oth.mixed ADDED
  edh.oth.mixed <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_oth_mixed_medcodeid.csv", sep = ""), 
                                     sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.oth.mixed <- edh.oth.mixed$medcodeid
  
  ### edh.oth.white ADDED
  edh.oth.white <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_oth_white_medcodeid.csv", sep = ""), 
                                     sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.oth.white <- edh.oth.white$medcodeid
  
  ### edh.pakistani ADDED
  edh.pakistani <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_pakistani_medcodeid.csv", sep = ""), 
                                     sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.pakistani <- edh.pakistani$medcodeid
  
  ### edh.white.asian ADDED
  edh.white.asian <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_white_asian_medcodeid.csv", sep = ""), 
                                       sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.white.asian <- edh.white.asian$medcodeid
  
  ### edh.white.black.african ADDED
  edh.white.black.african <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_white_black_african_medcodeid.csv", sep = ""), 
                                               sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.white.black.african <- edh.white.black.african$medcodeid
  
  ### edh.white.black.caribbean ADDED
  edh.white.black.caribbean <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_white_black_caribbean_medcodeid.csv", sep = ""), 
                                                 sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.white.black.caribbean <- edh.white.black.caribbean$medcodeid
  
  ### edh.white.british ADDED
  edh.white.british <- data.table::fread(file = paste(getwd(),"/codelists/analysis/edh_white_british_medcodeid.csv", sep = ""), 
                                         sep = ",", header = TRUE, colClasses = rep("character", 3))
  edh.white.british <- edh.white.british$medcodeid
  
  ### Combine into codelist categories
  codes.bangladeshi <- edh.bangladeshi
  codes.black.african <- edh.black.african
  codes.black.caribbean <- edh.black.caribbean
  codes.chinese <- edh.chinese
  codes.indian <- edh.indian
  codes.oth.asian <- c(edh.oth.asian, edh.white.asian)
  codes.oth.ethnic <- c(edh.oth.ethnic, edh.oth.black, edh.oth.mixed, edh.white.black.african, edh.white.black.caribbean)
  codes.pakistani <- edh.pakistani
  codes.white <- c(edh.white.british, edh.irish, edh.oth.white)
  
  ### And one codelist for extracting data
  codelist <- c(codes.bangladeshi, codes.black.african, codes.black.caribbean, codes.chinese, codes.indian,
                codes.oth.asian, codes.oth.ethnic, codes.pakistani, codes.white)
  
  ### Query database
  db.qry <- db_query(NULL,
                     db = db, 
                     db.filepath = db.filepath,
                     tab = "obs",
                     codelist.vector = codelist)
  str(db.qry)
  ### Merge cohort with the database query
  cohort.qry <- merge(cohort, db.qry, by.x = "patid", by.y = "patid")
  
  ### Reduce to variables of interest
  cohort.qry <- cohort.qry[,c("patid", "medcodeid", "obsdate")]
  
  ### Create ethnicity variable
  cohort.qry <- dplyr::mutate(cohort.qry,
                              ethnicity = dplyr::case_when(medcodeid %in% codes.bangladeshi ~ "bangladeshi",
                                                           medcodeid %in% codes.black.african ~ "black.african",
                                                           medcodeid %in% codes.black.caribbean ~ "black.caribbean",
                                                           medcodeid %in% codes.chinese ~ "chinese",
                                                           medcodeid %in% codes.indian ~ "indian",
                                                           medcodeid %in% codes.oth.asian ~ "other.asian",
                                                           medcodeid %in% codes.oth.ethnic ~ "other.ethnic",
                                                           medcodeid %in% codes.pakistani ~ "pakistani",
                                                           medcodeid %in% codes.white ~ "white"))
  
  ### For now, only keep most recent observation for each individual
  cohort.qry <- dplyr::group_by(cohort.qry, patid) |>
    dplyr::filter(dplyr::row_number(dplyr::desc(obsdate)) == 1) 
  
  ### Merge back with cohort
  variable.dat <- merge(cohort, cohort.qry, by.x = "patid", by.y = "patid", all.x = TRUE)
  
  ### If ethnicity is "not.stated", set to missing
  
  ### Turn into factor variable
  variable.dat$ethnicity.prim <- factor(variable.dat$ethnicity, 
                                        levels = c("bangladeshi", "black.african", "black.caribbean", "chinese", "indian", 
                                                   "other.asian", "other.ethnic", "pakistani", "white"), 
                                        labels = c("bangladeshi", "black african", "black caribbean", "chinese", "indian", 
                                                   "other asian", "other ethnic", "pakistani", "white"))
  
  ### Reduce to variables of interest
  variable.dat <- variable.dat[,c("patid", "ethnicity.prim")]
  
  ###################################################
  ### EXTRACT ETHNICITY USING SECONDARY CARE DATA ###
  ###################################################
  
  ### Read in the patient file 
  db.qry <- db_query(codelist = NULL,
                     db = db.linked, 
                     db.filepath = db.filepath,
                     tab = "hes_patient",
                     codelist.vector = NULL)
  
  ### Merge with variable.dat
  variable.dat <- merge(variable.dat, dplyr::select(db.qry, patid, gen_ethnicity), by.x = "patid", by.y = "patid", all.x = TRUE)
  
  ### Create same factor variable as primary care one
  variable.dat <- dplyr::mutate(variable.dat, ethnicity.hes = dplyr::case_when(gen_ethnicity == "White" ~ "white",
                                                                               gen_ethnicity == "Pakistani" ~ "pakistani",
                                                                               gen_ethnicity == "Other" ~ "other.ethnic",
                                                                               gen_ethnicity == "Oth_Asian" ~ "other.asian",
                                                                               gen_ethnicity == "Mixed" ~ "other.ethnic",
                                                                               gen_ethnicity == "Indian" ~ "indian",
                                                                               gen_ethnicity == "Chinese" ~ "chinese",
                                                                               gen_ethnicity == "Bl_Other" ~ "other.ethnic",
                                                                               gen_ethnicity == "Bl_Carib" ~ "black.caribbean",
                                                                               gen_ethnicity == "Bl_Afric" ~ "black.african",
                                                                               gen_ethnicity == "Bangladesi" ~ "bangladeshi"))
  
  ### Create factor
  variable.dat$ethnicity.hes <- factor(variable.dat$ethnicity.hes, 
                                        levels = c("bangladeshi", "black.african", "black.caribbean", "chinese", "indian", 
                                                   "other.asian", "other.ethnic", "pakistani", "white"), 
                                        labels = c("bangladeshi", "black african", "black caribbean", "chinese", "indian", 
                                                   "other asian", "other ethnic", "pakistani", "white"))
  
  ### Reduce to variables of interest
  variable.dat <- variable.dat[,c("patid", "ethnicity.prim", "ethnicity.hes")]
  
  ############################################
  ### Create a combined ethnicity variable ###
  ############################################
  variable.dat <- dplyr::mutate(variable.dat, ethnicity = dplyr::case_when(!is.na(ethnicity.hes) ~ ethnicity.hes,
                                                                           !is.na(ethnicity.prim) ~ ethnicity.prim))
  
  ### Change name of variable to varname
  colnames(variable.dat)[colnames(variable.dat) == "ethnicity"] <- varname
  
  ### Implement output
  implement_output(variable.dat, varname, out.save.disk, out.subdir, out.filepath, return.output)
  
}


