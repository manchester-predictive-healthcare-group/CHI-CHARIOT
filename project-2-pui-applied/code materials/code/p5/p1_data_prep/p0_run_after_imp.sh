#!/bin/bash --login
module load apps/gcc/R/4.4.2
Rscript p2_table1_cohort_render.R > p2_table1_cohort_render.out 2>&1
Rscript p5_create_imp_comb.R > p5_create_imp_comb.out 2>&1
Rscript p7_create_dfs_for_analysis.R > p7_create_dfs_for_analysis.out 2>&1
Rscript p8_get_knot_locations.R > p8_get_knot_locations.out 2>&1