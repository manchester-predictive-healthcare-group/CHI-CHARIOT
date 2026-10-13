module load apps/gcc/R/4.4.2
Rscript p1_calculate_cf_surv_times_parallel.R > p1_calculate_cf_surv_times_parallel.out 2>&1
Rscript p2_validate_internal_parallel.R > p2_validate_internal_parallel.out 2>&1