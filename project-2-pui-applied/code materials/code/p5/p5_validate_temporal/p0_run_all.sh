module load apps/gcc/R/4.4.2
#Rscript p1_temporalv_create_individual_cohorts.R > p1_temporalv_create_individual_cohorts.out 2>&1
#Rscript p2_temporalv_recalculate_cohort_split_times.R 1 > p2_temporalv_recalculate_cohort_split_times.out 2>&1
#Rscript p3_temporalv_test_cohort_sizes_match.R > p3_temporalv_test_cohort_sizes_match.out  2>&1
#Rscript p4_temporalv_calculcate_cf_surv_times.R > p4_temporalv_calculcate_cf_surv_times.out  2>&1
Rscript p5_temporalv_validate.R > p5_temporalv_validate.out  2>&1
Rscript p6_temporalv_validate_intervention_layer.R > p6_temporalv_validate_intervention_layer.out 2>&1