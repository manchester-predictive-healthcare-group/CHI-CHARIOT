CHARIOT APPLIED STUDY: RUN ORDER

This is a run-order guide, not an executable script.
Paths/folders mentioned in this run order are relative to code/p5 unless stated otherwise.
All .R scripts require setting the working directory to the project root (folder containing code, R, data, figures).
Before running the analyses, edit setwd() in each R script to the absolute path of the project root (the folder containing code, R, data and figures).
Run each stage to completion before moving to the next stage.
Shell (.sh) scripts and Slurm (.slurm) scripts indicate inputs for running R files.
R/functions.R contains functions used across analyses.
Adapt hard-coded paths and software/module settings to the local environment.
Create the required data/p5, figures/p5 and p7_rshiny/data output directories.

PREREQUISITES
Access to the extracted CPRD baseline cohort, follow-up cohorts at years 1-5, database, codelists and region lookup
referenced by the scripts is required. Initial cohort extraction is external
to this upload.

1. BASELINE COHORT AND IMPUTATION
In p1_data_prep:
  Rscript p2_create_cohort.R
  Rscript p3_create_variable_objects.R
  sbatch p4_impute_cohort.slurm
  WAIT for all 20 array tasks to finish successfully.
  bash p0_run_after_imp.sh
  Rscript p9_calculate_DAG.R

2. FOLLOW-UP SMOKING AND PRESCRIPTION HISTORIES
In p2_data_prep_further_mfr, run in this order:
  Rscript p1_extract_modifiable_risk_factors.R
  Rscript p2_extract_prescriptions.R
  Rscript p3.1_extract_smoking_followup_for_model.R
  Run p3.2A_extract_prescriptions_followup_for_model.R for both medications
    (statins and antihypertensives), burnout 180, and all groups 1-100:
    Rscript p3.2A_extract_prescriptions_followup_for_model.R MEDICATION 180 100 GROUP
  WAIT for all medication/group combinations to finish successfully.
  Rscript p3.2B_combine_prescriptions_followup_for_model.R statins 180 100
  Rscript p3.2B_combine_prescriptions_followup_for_model.R antihypertensives 180 100
  Rscript p3.2C_augment_prescriptions_followup_for_model.R statins
  Rscript p3.2C_augment_prescriptions_followup_for_model.R antihypertensives
  Rscript p4.1_layer_cohort_split_times.R
  Optional diagnostic: Rscript p5_test_smoking.R

3. MODEL DEVELOPMENT
In p3_devel:
  bash p0_run_all.sh

4. INTERNAL VALIDATION
In p4_validate:
  bash p0_run_all.sh

5. TEMPORAL VALIDATION
In p5_validate_temporal:
  bash p0_run_all.sh

6. BOOTSTRAP OPTIMISM AND INSTABILITY
In p6_optimism_and_instability:
  sbatch p2_optimism_loop.slurm
  WAIT for all bootstrap iterations selected for the final analysis.
The current array 1-2000 runs iterations 1-1000 for each sex.
For 500 iterations per sex, submit with --array=1-1000 instead.
Use the same final iteration range in p9_collate_results.

8. IMPUTATION DIAGNOSTICS AND SENSITIVITY ANALYSES
In p8_imputation_sensitivities:
  bash p0_run_all.sh

9. COLLATE RESULTS AND PRODUCE REPORTS
In p9_collate_results:
  bash p0_run_all.sh
  Rscript supp_0_run_specific_rmds.R
The report-rendering script also runs the counterfactual-method simulation
embedded in supp_4_simulation_validate_xu_extension.Rmd. It requires substantial
memory and refs.bib beside that report; refs.bib is not included in this text upload.
