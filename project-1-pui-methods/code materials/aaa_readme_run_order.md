CHARIOT METHODS STUDY: RUN ORDER

Before running the code, set setwd() inside each R file to the root directory of the project.
--------------------------------------------------------------------------------------------------------------------------------------------
1. Data preparation

Run:

code/p1_data_prep/p0_run_all.sh

--------------------------------------------------------------------------------------------------------------------------------------------
2. Model development and main validation

Run:

code/p2_devel/p0_run_all.sh

It fits models 0 to 7 for both female and male cohorts, calculates counterfactual survival times, assesses calibration and discrimination, and creates the main development outputs.

--------------------------------------------------------------------------------------------------------------------------------------------
3. Treatment-effect misspecification sensitivity analysis

First run:

code/p3_sensitivity_treatment_effects/p0_run_all.sh

It creates the counterfactual survival-time outputs for all treatment-effect scenarios and both cohorts.

After it has completed, submit these Slurm arrays:

code/p3_sensitivity_treatment_effects/p2_run_assess_performance.slurm
code/p3_sensitivity_treatment_effects/p2_run_assess_performance_stratified.slurm

The two arrays assess overall and treatment-stratified performance.

--------------------------------------------------------------------------------------------------------------------------------------------
4. K&VG artificial-censoring sensitivity analysis

Run the main K&VG validation analyses for both cohorts, followed by the truncation analysis:

code/p4_sensitivity_kvg_validation/p1_create_validation_data_and_assess_calibration.R 1
code/p4_sensitivity_kvg_validation/p1_create_validation_data_and_assess_calibration.R 2
code/p4_sensitivity_kvg_validation/p1_create_validation_data_and_assess_calibration_truncate.R

Then submit this Slurm array from the code/p4_sensitivity_kvg_validation directory:

code/p4_sensitivity_kvg_validation/p2.slurm

Only after every task in p2.slurm has completed successfully, submit:

code/p4_sensitivity_kvg_validation/p3.slurm

The p2 array creates the sustained-treatment cohorts and estimates the weights. The p3 array calculates the treatment-effect estimates, effective sample sizes and weight diagnostics.

--------------------------------------------------------------------------------------------------------------------------------------------
5. Missingness sensitivity analysis

Run:

code/p5_sensitivity_imputation/p1_performance_stratified_by_unexposed_sbp_missing.R

This R file covers both female and male cohorts and both unexposed-SBP missingness strata.

--------------------------------------------------------------------------------------------------------------------------------------------
6. Direct-effect-of-statins sensitivity analysis

Run the main analysis once for each cohort:

code/p6_sensitivity_direct_effect_statins/p1_direct_effect_statins.R 1
code/p6_sensitivity_direct_effect_statins/p1_direct_effect_statins.R 2

After both runs have completed, run:

code/p6_sensitivity_direct_effect_statins/p2_worked_example_direct_effect_statins.R

--------------------------------------------------------------------------------------------------------------------------------------------
7. Temporal validation

Run:

code/p7_sensitivity_temporal/p0_run_all.sh

It covers both female and male cohorts and temporal index-date shifts of 1 to 5 years.

--------------------------------------------------------------------------------------------------------------------------------------------
8. Produce figures and supplementary material

Only run this stage after all preceding analyses and Slurm arrays have completed successfully.

Run:

code/p2_devel/p6_create_manuscript_figures.R
code/p2_devel/p5_render.R
code/p9_sensitivity_summaries/sens_0_run_specific_rmds.R

The sensitivity-summary script renders supplementary files 4 to 8. The p2_devel rendering script produces the main results supplement and associated development outputs.
--------------------------------------------------------------------------------------------------------------------------------------------