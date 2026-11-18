module load apps/gcc/R/4.4.2
Rscript p1_summarise_optimism_metrics.R > p1_summarise_optimism_metrics.out 2>&1
Rscript p2_summarise_temporal_validation_metrics.R > p2_summarise_temporal_validation_metrics.out 2>&1
Rscript p3_create_validation_plots.R > p3_create_validation_plots.out 2>&1
Rscript p4_summarise_variance.R > p4_summarise_variance.out 2>&1
Rscript p5_discrimination_forestplot.R > p5_discrimination_forestplot.out 2>&1
Rscript p6_instability_plots_synthetic.R > p6_instability_plots_synthetic.out 2>&1
Rscript p7_create_manuscript_figures.R > p7_create_manuscript_figures.out 2>&1
Rscript p8_create_incidence_by_year_tables.R > p8_create_incidence_by_year_tables.out 2>&1