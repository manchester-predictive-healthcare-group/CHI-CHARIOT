module load apps/gcc/R/4.4.2
#Rscript p1_create_cohort.R > p1_create_cohort.out
#Rscript p2_create_variable_objects.R > p2_create_variable_objects.out
#Rscript p3.1_layer_cohort_split_times_statins_ah.R > p3.1_layer_cohort_split_times_statins_ah.out
#Rscript p3.2_layer_cohort_split_times_statins_ah_comb.R > p3.2_layer_cohort_split_times_statins_ah_comb.out
#Rscript p3_impute_cohort.R 1 1 > p3_impute_cohort_1.R
#Rscript p3_impute_cohort.R 1 2 > p3_impute_cohort_2.R
#Rscript p4_create_split_sample.R > p4_create_split_sample.out
#Rscript p6.1_extract_prescriptions.R > p6.1_extract_prescriptions.out
Rscript p6.2_extract_medication_status.R statins > p6.2_extract_medication_status_statins.out
Rscript p6.2_extract_medication_status.R antihyerptensives > p6.2_extract_medication_status_antihyerptensives.out
Rscript p6.3_layer_medication_status.R > p6.3_layer_medication_status.out
Rscript p7.1_augment_medication_status.R statins > p7.1_augment_medication_status_statins.out
Rscript p7.1_augment_medication_status.R antihypertensives > p7.1_augment_medication_status_antihypertensives.out
Rscript p7.2_layer_augmented_medication_status.R > p7.2_layer_augmented_medication_status.out
Rscript p5_assess_imputations.R > p5_assess_imputations.out
