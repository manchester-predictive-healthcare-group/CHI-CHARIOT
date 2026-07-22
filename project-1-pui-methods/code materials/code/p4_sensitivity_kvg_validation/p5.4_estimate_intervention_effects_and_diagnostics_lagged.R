###
### Estimate the relative risk ratio for the three cohorts, which have been created in file 5.2:
###
### Full cohort
### Overlap cohort (based on IPTW only)
### Overlap cohort comb (based on weights of IPTW, IPCW, and IPACW)
###
### Also calculate effective sample sizes for each cohort, and run weight diagnostics
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/projectM1/")
getwd()

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Load survival
library(survival)

### Define medication, gender and read in dataset
### Extract gender and medication 
args <- commandArgs(trailingOnly = T)
gender_in <- as.numeric(args[1])
med <- args[2]
t_fup_int <- as.numeric(args[3])
t_fup <- t_fup_int*365.25
model_flex <- args[4]

### Read in the three cohorts
# Entire cohort
df_valid_ipacw <- readRDS(paste("data/df_valid_ipacw_", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
# Overlap cohort
df_overlap <- readRDS(paste("data/df_valid_overlap_cohort", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
# Overlap cohort based on combined weights
df_overlap_comb <- readRDS(paste("data/df_valid_overlap_comb_lag_cohort", gender_in, "_", med, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))

#############################################
### Step 1: Estimate relative risk ratios ###
#############################################

###
### Define functions to estimate risks in pseudo-population
###

### Horvitz-Thompson estimator
estimate_risk_pseudo_ht <- function(df, group){
  
  ### Subset data
  df_subset <- dplyr::filter(df, med_status == group)
  
  if (group == 0){
    
    ### Assign weights and indicator
    w <- df_subset$comb_weight_baseline_lag
    y <- df_subset$status_lag == 1 & df_subset$time_lag <= t_fup
    
  } else if (group == 1){
    
    ### Assign weights and indicator
    w <- df_subset$comb_weight_baseline
    y <- df_subset$status == 1 & df_subset$time <= t_fup
    
  }
  
  ### Get risk
  obs_risk <- sum(w * y)/nrow(df_subset)
  
  ### Return it
  return(obs_risk)
  
}

### Hajek-style estimator
estimate_risk_pseudo_hajek <- function(df, group){
  
  ### Subset data
  df_subset <- dplyr::filter(df, med_status == group)
  
  if (group == 0){
    
    ### Exclude individuals that haven't had an event and are censored before t_fup
    df_subset <- dplyr::filter(df_subset, !(status_lag == 0 & time_lag < t_fup))
    
    ### Assign weights and indicator
    w <- df_subset$comb_weight_baseline_lag
    y <- df_subset$status_lag == 1 & df_subset$time_lag <= t_fup
    
  } else if (group == 1){
    
    ### Exclude individuals that haven't had an event and are censored before t_fup
    df_subset <- dplyr::filter(df_subset, !(status == 0 & time < t_fup))
    
    ### Assign weights and indicator
    w <- df_subset$comb_weight_baseline
    y <- df_subset$status == 1 & df_subset$time <= t_fup
    
  }
  
  ### Get risk
  obs_risk <- sum(w * y)/sum(w)
  
  ### Return it
  return(obs_risk)
  
}

###
### Apply functions to estimate risks in each group
###

### Hajek estimators
risk_pseudo_entire_cohort_hajek <- sapply(0:1, function(x){estimate_risk_pseudo_hajek(df = df_valid_ipacw, group = x)})
risk_pseudo_overlap_cohort_hajek <- sapply(0:1, function(x){estimate_risk_pseudo_hajek(df = df_overlap, group = x)})
risk_pseudo_overlap_comb_cohort_hajek <- sapply(0:1, function(x){estimate_risk_pseudo_hajek(df = df_overlap_comb, group = x)})

### Horvitz-Thompsom estimator
risk_pseudo_entire_cohort_ht <- sapply(0:1, function(x){estimate_risk_pseudo_ht(df = df_valid_ipacw, group = x)})
risk_pseudo_overlap_cohort_ht <- sapply(0:1, function(x){estimate_risk_pseudo_ht(df = df_overlap, group = x)})
risk_pseudo_overlap_comb_cohort_ht <- sapply(0:1, function(x){estimate_risk_pseudo_ht(df = df_overlap_comb, group = x)})

### Save output
saveRDS(risk_pseudo_entire_cohort_hajek, 
        paste("data/rr_lagged_entire_cohort_hajek_", med, "_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
saveRDS(risk_pseudo_overlap_cohort_hajek, 
        paste("data/rr_lagged_overlap_cohort_hajek_", med, "_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
saveRDS(risk_pseudo_overlap_comb_cohort_hajek, 
        paste("data/rr_lagged_overlap_comb_cohort_hajek_", med, "_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))

saveRDS(risk_pseudo_entire_cohort_ht, 
        paste("data/rr_lagged_entire_cohort_ht_", med, "_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
saveRDS(risk_pseudo_overlap_cohort_ht, 
        paste("data/rr_lagged_overlap_cohort_ht_", med, "_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
saveRDS(risk_pseudo_overlap_comb_cohort_ht, 
        paste("data/rr_lagged_overlap_comb_cohort_ht_", med, "_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))

print("FINISHED est risk")

###################################
### Step 2: Run ESS diagnostics ###
###################################

############################################################
### ESS FUNCTIONS
############################################################

### Standard Kish formula for effective sample size
### ESS = (sum of weights)^2 / sum of squared weights
### ESS == N when all weights are equal (no variability)
### ESS << N when weights are highly variable
ess <- function(w) sum(w)^2 / sum(w^2)

### Compute ESS overall and by treatment group for a given weight vector
### Arguments:
###   w               - vector of weights
###   label           - descriptive label for this weight stage
###   med_status_vec  - vector of treatment group indicators (0/1), same length as w
### Returns a named list with overall and by-group ESS
compute_ess <- function(w, label, med_status_vec) {
  n_overall   <- length(w)
  ess_overall <- ess(w)
  by_group <- lapply(0:1, function(g) {
    w_g <- w[med_status_vec == g]
    n_g <- length(w_g)
    list(
      group   = g,
      n       = n_g,
      ess     = ess(w_g),
      ess_pct = 100 * ess(w_g) / n_g
    )
  })
  names(by_group) <- paste0("group_", 0:1)
  list(
    label    = label,
    n        = n_overall,
    ess      = ess_overall,
    ess_pct  = 100 * ess_overall / n_overall,
    by_group = by_group
  )
}

### Print a single ESS result object to console
### Prints overall ESS and ESS within each treatment group,
### expressed both as absolute numbers and as % of unweighted n
print_ess <- function(x) {
  cat("\n---", x$label, "---\n")
  cat("  Overall  ESS:", round(x$ess),
      " (", round(x$ess_pct, 1), "% of N =", x$n, ")\n")
  for (g in x$by_group) {
    cat("  Group", g$group, "   ESS:", round(g$ess),
        " (", round(g$ess_pct, 1), "% of n =", g$n, ")\n")
  }
}

### Run full ESS analysis on a cohort, printing results to console and
### returning them as a list that can be saved
### Arguments:
###   df          - data frame containing weights and med_status
###   label       - label for this cohort (used in console output)
### Returns: named list of ESS results at each weighting stage
run_ess_analysis <- function(df, label) {
  
  cat("\n############################################################\n")
  cat("### ESS ANALYSIS:", label, "\n")
  cat("############################################################\n")
  
  ### Compute ESS at each weighting stage
  ### Unweighted entry records N and n by group as a reference point
  results <- list(
    unweighted = list(
      label    = "Unweighted",
      n        = nrow(df),
      by_group = list(
        group_0 = list(group = 0, n = sum(df$med_status == 0)),
        group_1 = list(group = 1, n = sum(df$med_status == 1))
      )
    ),
    ipcw       = compute_ess(df$ipcw_weight,         "IPCW weights only",                     df$med_status),
    ipacw      = compute_ess(df$ipacw_weight_lag,         "IPACW weights only",                    df$med_status),
    ipcw_ipacw = compute_ess(df$comb_weight_censoring_lag,          "IPCW x IPACW combined",                 df$med_status),
    iptw       = compute_ess(df$iptw_weight,          "IPTW weights only",                     df$med_status),
    combined   = compute_ess(df$comb_weight_baseline_lag, "IPCW x IPACW x IPTW (final combined)", df$med_status)
  )
  
  ### Print unweighted cohort size as reference
  cat("\n### Cohort size (unweighted)\n")
  cat("  Overall N:", results$unweighted$n, "\n")
  cat("  Group 0 n:", results$unweighted$by_group$group_0$n, "\n")
  cat("  Group 1 n:", results$unweighted$by_group$group_1$n, "\n")
  
  ### Print ESS at each weighting stage
  for (nm in c("ipcw", "ipacw", "ipcw_ipacw", "iptw", "combined")) {
    print_ess(results[[nm]])
  }
  
  return(results)
}

### Get ESS 
# Entire cohort
ess_results_entire_cohort <- run_ess_analysis(df_valid_ipacw, "ENTIRE COHORT")
saveRDS(ess_results_entire_cohort, paste("data/ess_lagged_", med, "_entire_cohort_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))

# Overlap cohort
ess_results_overlap_cohort <- run_ess_analysis(df_overlap, "OVERLAP COHORT")
saveRDS(ess_results_overlap_cohort, paste("data/ess_lagged_", med, "_overlap_cohort_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))

# Overlap cohort comb
ess_results_overlap_cohort_comb <- run_ess_analysis(df_overlap_comb, "OVERLAP COHORT")
saveRDS(ess_results_overlap_cohort_comb, paste("data/ess_lagged_", med, "_overlap_cohort_comb_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))


###########################################
### Step 3: COMBINED WEIGHT DIAGNOSTICS ###
###########################################
run_weight_diagnostics <- function(df, label) {
  
  df_treated   <- dplyr::filter(df, med_status == 1)
  df_untreated <- dplyr::filter(df, med_status == 0)
  
  cat("\n############################################################\n")
  cat("###", label, "\n")
  cat("############################################################\n")
  
  cat("\nCohort N:", nrow(df), "  Treated:", nrow(df_treated),
      "  Untreated:", nrow(df_untreated), "\n")
  
  ### 1. Distribution of combined weight by treatment group
  cat("\n--- Combined weight distribution ---\n")
  cat("Untreated:\n")
  q_untreated <- quantile(df_untreated$comb_weight_baseline_lag,
                          p = c(0, 0.5, 0.9, 0.95, 0.99, 0.999, 1))
  print(q_untreated)
  cat("Treated:\n")
  q_treated <- quantile(df_treated$comb_weight_baseline_lag,
                        p = c(0, 0.5, 0.9, 0.95, 0.99, 0.999, 1))
  print(q_treated)
  
  ### 2. Weight concentration in treated group
  top1_threshold <- quantile(df_treated$comb_weight_baseline_lag, 0.99)
  top1_share <- sum(df_treated$comb_weight_baseline_lag[df_treated$comb_weight_baseline_lag >= top1_threshold]) /
    sum(df_treated$comb_weight_baseline_lag)
  cat("\n--- Weight concentration in treated group ---\n")
  cat("  Top 1% individuals hold", round(100 * top1_share, 1), "% of total weight\n")
  
  ### 3. Correlations between weight components
  cat("\n--- Spearman correlations between weight components (treated group) ---\n")
  weight_mat <- df_treated[, c("ipcw_weight", "ipacw_weight_lag", "iptw_weight")]
  cor_mat <- round(cor(weight_mat, method = "spearman"), 3)
  print(cor_mat)
  
  ### 4. Profile extreme vs normal combined-weight individuals
  df_extreme_comb <- dplyr::filter(df_treated, comb_weight_baseline_lag >= top1_threshold)
  df_normal_comb  <- dplyr::filter(df_treated, comb_weight_baseline_lag <  top1_threshold)
  
  weight_drivers <- list(
    ipcw  = list(extreme = round(mean(df_extreme_comb$ipcw_weight), 2),
                 normal  = round(mean(df_normal_comb$ipcw_weight),  2)),
    ipacw = list(extreme = round(mean(df_extreme_comb$ipacw_weight_lag), 2),
                 normal  = round(mean(df_normal_comb$ipacw_weight_lag),  2)),
    iptw  = list(extreme = round(mean(df_extreme_comb$iptw_weight), 2),
                 normal  = round(mean(df_normal_comb$iptw_weight),  2))
  )
  
  cat("\n--- What is driving extreme combined weights? (treated, top 1%) ---\n")
  cat("  Mean IPCW weight:  extreme =", weight_drivers$ipcw$extreme,
      "  normal =", weight_drivers$ipcw$normal, "\n")
  cat("  Mean IPACW weight: extreme =", weight_drivers$ipacw$extreme,
      "  normal =", weight_drivers$ipacw$normal, "\n")
  cat("  Mean IPTW weight:  extreme =", weight_drivers$iptw$extreme,
      "  normal =", weight_drivers$iptw$normal, "\n")
  
  ### 5. Covariate profiles
  cat("\n--- Covariate profiles: extreme vs normal combined weights (treated) ---\n")
  covariate_profiles <- lapply(c("age", "sbp", "bmi", "nonhdl", "IMD"), function(var) {
    cat(sprintf("  %-10s  extreme: %.1f (%.1f)   normal: %.1f (%.1f)\n",
                var,
                mean(df_extreme_comb[[var]], na.rm = TRUE),
                sd(df_extreme_comb[[var]],   na.rm = TRUE),
                mean(df_normal_comb[[var]],  na.rm = TRUE),
                sd(df_normal_comb[[var]],    na.rm = TRUE)))
    list(
      extreme_mean = mean(df_extreme_comb[[var]], na.rm = TRUE),
      extreme_sd   = sd(df_extreme_comb[[var]],   na.rm = TRUE),
      normal_mean  = mean(df_normal_comb[[var]],  na.rm = TRUE),
      normal_sd    = sd(df_normal_comb[[var]],    na.rm = TRUE)
    )
  })
  names(covariate_profiles) <- c("age", "sbp", "bmi", "nonhdl", "IMD")
  
  ### 6. Overlap of extreme weights across components
  df_treated$flag_ipcw  <- df_treated$ipcw_weight  >= quantile(df_treated$ipcw_weight,  0.95)
  df_treated$flag_ipacw <- df_treated$ipacw_weight_lag >= quantile(df_treated$ipacw_weight_lag, 0.95)
  df_treated$flag_iptw  <- df_treated$iptw_weight  >= quantile(df_treated$iptw_weight,  0.95)
  df_treated$n_flags    <- df_treated$flag_ipcw + df_treated$flag_ipacw + df_treated$flag_iptw
  
  flag_table <- table(n_components_extreme = df_treated$n_flags)
  n_all_three <- sum(df_treated$n_flags == 3)
  pct_all_three <- round(100 * mean(df_treated$n_flags == 3), 2)
  
  cat("\n--- Overlap of extreme weights across components (treated group) ---\n")
  print(flag_table)
  cat("  In top 5% of all three simultaneously:",
      n_all_three, "(", pct_all_three, "%)\n")
  
  ### Return results as a list
  invisible(list(
    label       = label,
    n           = list(overall = nrow(df), treated = nrow(df_treated),
                       untreated = nrow(df_untreated)),
    quantiles   = list(treated = q_treated, untreated = q_untreated),
    top1_share  = top1_share,
    cor_mat     = cor_mat,
    weight_drivers     = weight_drivers,
    covariate_profiles = covariate_profiles,
    flag_table         = flag_table,
    n_all_three        = n_all_three,
    pct_all_three      = pct_all_three
  ))
}

### Get diagnostics
# Entire cohort
weight_diagnostics_entire_cohort <- run_weight_diagnostics(df_valid_ipacw, "ENTIRE COHORT")
saveRDS(weight_diagnostics_entire_cohort, 
        paste("data/weight_diagnostics_lagged_", med, "_entire_cohort_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))

# Overlap cohort
weight_diagnostics_overlap_cohort <- run_weight_diagnostics(df_overlap,     "OVERLAP COHORT")
saveRDS(weight_diagnostics_overlap_cohort, 
        paste("data/weight_diagnostics_lagged_", med, "_overlap_cohort_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))

# Overlap cohort comb
weight_diagnostics_overlap_cohort_comb <- run_weight_diagnostics(df_overlap_comb,     "OVERLAP COHORT COMB")
saveRDS(weight_diagnostics_overlap_cohort_comb, 
        paste("data/weight_diagnostics_lagged_", med, "_overlap_cohort_comb_", gender_in, "_tfup", t_fup_int, "_", model_flex, ".rds", sep = ""))
