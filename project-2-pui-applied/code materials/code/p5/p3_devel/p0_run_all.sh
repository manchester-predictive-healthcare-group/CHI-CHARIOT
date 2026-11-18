#!/bin/bash --login
module load apps/gcc/R/4.4.2
Rscript p1_fit_model_parallel.R > p1_fit_model_parallel.out 2>&1
Rscript p2_assess_nonlinear.R > p2_assess_nonlinear.out 2>&1
Rscript p2_assess_nonlinear_render.R > p2_assess_nonlinear_render.out 2>&1
Rscript p3_save_models_redacted.R > p3_save_models_redacted.out 2>&1
Rscript p4_tests_to_compare_risk_estimation_functions.R > p4_tests_to_compare_risk_estimation_functions.out 2>&1
Rscript p5_sampsize_calc.R > p5_sampsize_calc.out 2>&1