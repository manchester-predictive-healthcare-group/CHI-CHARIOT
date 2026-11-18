### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("./")
getwd()

### Extract female models
female_models <- readRDS(paste0("data/p5/small_fits_model_", "female", ".rds"))
length(female_models[[1]])

### Take first one
fm1 <- female_models[[1]]

### Look at names
names(fm1)

### Element 1 is coefficients
str(fm1[[1]])

### Element 2 is covariance matrix
fm1[[2]]
str(fm1[[2]])

### Element 3 is loglikelihood
fm1[[3]]
str(fm1[[3]])

### Element 4 is score
fm1[[4]]
str(fm1[[4]])

### Element 5 is iter
fm1[[5]]
str(fm1[[5]])

### Element 6 is means (numeric vector, one for each coefficient)
fm1[[6]]
str(fm1[[6]])

### Element 7 is first (numeric vector, one for each coefficient)
fm1[[7]]
str(fm1[[7]])

### Element 8 is info (rank, rescale, step halving, convergence)
fm1[[8]]
str(fm1[[8]])

### Element 9 is method for handling ties
fm1[[9]]
str(fm1[[9]])

### Element 10 is n (number of rows of split survival data)
fm1[[10]]
str(fm1[[10]])

### Element 11 is numer of events
fm1[[11]]
str(fm1[[11]])

### Element 12 is terms (contains information about model formula)
str(fm1)
fm1[[12]]
fm1[["terms"]]
# terms contains an environment, verify this isn't the global one
terms_env <- environment(fm1[["terms"]])
identical(terms_env, globalenv())
# The environment only contains knot locations
ls(terms_env, all.names = TRUE)
vapply(ls(terms_env, all.names = TRUE),
       function(x) paste(get(x, envir = terms_env, inherits = FALSE), collapse = ", "),
       character(1))

### Element 13 is assign, linking model matrix to formula
fm1[[13]]
str(fm1[[13]])

### Element 14 is wald tests
fm1[[14]]
str(fm1[[14]])

### Element 15 is concordance
fm1[[15]]
str(fm1[[15]])

### Element 16 is timefix (TRUE/FALSE)
fm1[[16]]
str(fm1[[16]])

### Element 17 is formula
fm1[[17]]
str(fm1[[17]])
# formula contains an environment, verify this isn't the global one
formula_env <- environment(fm1[["formula"]])
identical(formula_env, globalenv())
# The environment only contains knot locations
ls(formula_env, all.names = TRUE)
vapply(ls(formula_env, all.names = TRUE),
       function(x) paste(get(x, envir = formula_env, inherits = FALSE), collapse = ", "),
       character(1))

### Element 18 is levels of categorical predictor variables
fm1[[18]]
str(fm1[[18]])

### Element 19 is contrasts
fm1[[19]]
str(fm1[[19]])

### Element 20 is the function call
fm1[[20]]
str(fm1[[20]])