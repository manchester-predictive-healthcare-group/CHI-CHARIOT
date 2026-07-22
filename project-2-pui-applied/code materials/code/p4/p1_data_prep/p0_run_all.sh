module load apps/gcc/R/4.4.2 libs/gcc/gmp/4.3.2 libs/gcc/mpfr/4.1.0 libs/gcc/mpc/1.2.1 compilers/gcc/14.2.0
Rscript p5_create_imp_comb.R > p5_create_imp_comb.out
Rscript p7_create_split_sample.R > p7_create_split_sample.out
Rscript p8_get_knot_locations.R > p8_get_knot_locations.out
Rscript p6_assess_imputations.R > p6_assess_imputations.out