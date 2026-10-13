module load apps/gcc/R/4.4.2
#Rscript p1_diagnostics.R > p1_diagnostics.out 2>&1
#Rscript p2_assess_mar_and_support.R > p2_assess_mar_and_support.out 2>&1
Rscript p3.1_mean_imputation.R > p3.1_mean_imputation.out 2>&1
Rscript p3.2_mean_imputation_stratified.R > p3.2_mean_imputation_stratified.out 2>&1
#Rscript p4_mnar_sensitivity.R > p4_mnar_sensitivity.out 2>&1