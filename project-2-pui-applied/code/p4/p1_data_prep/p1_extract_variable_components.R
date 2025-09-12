###
### I want to get level of recording of weight, height, BMI, HDL, LDL, total, total/HDL ratio, triglycerides and non-hdl
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd()
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Create variable list
variable.list <- c(
  ### BMI-related
  "height"#1
  ,"weight"#2
  ,"bmi"#3
  ### Cholesterol-related
  ,"total_chol"#4
  ,"hdl"#5
  ,"total_chol_hdl_ratio"#6
  ,"ldl"#7
  ,"triglycerides"#8
  ,"nonhdl"#9
)

### Extract scenario from command line
args <- commandArgs(trailingOnly = T)
variable.name <- variable.list[as.numeric(args[1])]
cohort.base.name <- args[2]
cohort.base <- readRDS(paste("data/extraction/", cohort.base.name ,".rds", sep = ""))
cohort.base <- data.table::as.data.table(cohort.base)
indexdt.name <- args[3]
out.subdir <- args[4]
t <- as.numeric(args[5])
if (t == 0){t <- NULL}

print(paste("variable = ", variable.name, sep = ""))
print(paste("cohort.base.name = ", cohort.base.name, sep = ""))
print(paste("indexdt.name = ", indexdt.name, sep = ""))
print(paste("out.subdir = ", out.subdir, sep = ""))
print(paste("t = ", t, sep = ""))

### Apply extraction
if (variable.name == "height"){
  print(paste(variable.name, Sys.time()))
  extract_height(cohort = cohort.base, 
                 varname = "component_height",
                 codelist = "height_medcodeid",
                 lower.bound = 100,
                 upper.bound = 250,
                 indexdt = indexdt.name, 
                 db = "aurum", 
                 out.save.disk = TRUE, 
                 out.filepath = NULL, 
                 out.subdir = out.subdir, 
                 return.output = FALSE)
  print("FINISHED")
} else if (variable.name == "weight"){
  print(paste(variable.name, Sys.time()))
  extract_weight(cohort = cohort.base, 
                 varname = "component_weight",
                 codelist = "weight_medcodeid",
                 lower.bound = 30,
                 upper.bound = 300,
                 indexdt = indexdt.name, 
                 db = "aurum", 
                 out.save.disk = TRUE, 
                 out.filepath = NULL, 
                 out.subdir = out.subdir, 
                 return.output = FALSE)
  print("FINISHED")
} else if (variable.name == "bmi"){
  print(paste(variable.name, Sys.time()))
  extract_test(cohort = cohort.base, 
               varname = "component_bmi",
               codelist = "edh_bmi_medcodeid",
               indexdt = indexdt.name, 
               db = "aurum", 
               out.save.disk = TRUE, 
               out.filepath = NULL, 
               out.subdir = out.subdir, 
               return.output = FALSE)
  print("FINISHED")
} else if (variable.name == "total_chol"){
  print(paste(variable.name, Sys.time()))
  extract_test(cohort = cohort.base, 
               varname = "component_total_chol",
               codelist = "muzambi_cholesterol_medcodeid",
               indexdt = indexdt.name, 
               lower.bound = 0.4,
               upper.bound = 20.7,
               db = "aurum", 
               out.save.disk = TRUE, 
               out.filepath = NULL, 
               out.subdir = out.subdir, 
               return.output = FALSE)
  print("FINISHED")
} else if (variable.name == "hdl"){
  print(paste(variable.name, Sys.time()))
  extract_test(cohort = cohort.base, 
               varname = "component_hdl",
               codelist = "muzambi_hdl_medcodeid",
               indexdt = indexdt.name, 
               lower.bound = 0.4,
               upper.bound = 20.7,
               db = "aurum", 
               out.save.disk = TRUE, 
               out.filepath = NULL, 
               out.subdir = out.subdir, 
               return.output = FALSE)
  print("FINISHED")
} else if (variable.name == "total_chol_hdl_ratio"){
  print(paste(variable.name, Sys.time()))
  extract_test(cohort = cohort.base, 
               varname = "component_total_chol_hdl_ratio",
               codelist = "muzambi_ratio_medcodeid",
               indexdt = indexdt.name, 
               lower.bound = 1,
               upper.bound = 12,
               db = "aurum", 
               out.save.disk = TRUE, 
               out.filepath = NULL, 
               out.subdir = out.subdir, 
               return.output = FALSE)
  print("FINISHED")
} else if (variable.name == "ldl"){
  print(paste(variable.name, Sys.time()))
  extract_test(cohort = cohort.base, 
               varname = "component_ldl",
               codelist = "muzambi_ldl_medcodeid",
               indexdt = indexdt.name, 
               lower.bound = 0.4,
               upper.bound = 20.7,
               db = "aurum", 
               out.save.disk = TRUE, 
               out.filepath = NULL, 
               out.subdir = out.subdir, 
               return.output = FALSE)
  print("FINISHED")
} else if (variable.name == "triglycerides"){
  print(paste(variable.name, Sys.time()))
  extract_triglycerides(cohort = cohort.base, 
                        varname = "component_triglycerides",
                        codelist = "uom_triglycerides_medcodeid",
                        indexdt = indexdt.name, 
                        db = "aurum", 
                        out.save.disk = TRUE, 
                        out.filepath = NULL, 
                        out.subdir = out.subdir, 
                        return.output = FALSE)
  print("FINISHED")
} else if (variable.name == "nonhdl"){
  print(paste(variable.name, Sys.time()))
  extract_nonhdl(cohort = cohort.base, 
                 varname = "component_nonhdl",
                 codelist.chol = "muzambi_cholesterol_medcodeid",
                 codelist.hdl = "muzambi_hdl_medcodeid",
                 codelist.ldl = "muzambi_ldl_medcodeid",
                 indexdt = indexdt.name, 
                 lower.bound = 0.4,
                 upper.bound = 20.7,
                 db = "aurum", 
                 out.save.disk = TRUE, 
                 out.filepath = NULL, 
                 out.subdir = out.subdir, 
                 return.output = FALSE)
  print("FINISHED")
} 