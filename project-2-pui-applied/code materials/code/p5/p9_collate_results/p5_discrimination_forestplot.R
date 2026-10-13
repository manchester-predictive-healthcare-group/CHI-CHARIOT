###
### Create forest plots for optimism-adjusted discrimination
###

### Clear workspace
rm(list=ls())
Sys.time()

### Set wd
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Define filepath to file directory system containing extracted data.
common.data.dir <- file.path("..", "..")

### Source functions
R.func.sources = list.files("R", full.names = TRUE)
sapply(R.func.sources, source)

### Read in the optimism adjusted performance metrics
optimism_df <- readRDS("data/p5/optimism_adjusted_performance_metrics.rds")[[1]] |>
  dplyr::filter(Metric == "Discrimination")

### Create a forest plot for one gender
create_forestplot <- function(gender) {

  ### Set gender labels
  gender_char <- c("male", "female")[gender]
  print(paste("gender =", gender_char))

  ### Create table
  discrimination_combined <- optimism_df |>
    dplyr::filter(Gender == c("Male", "Female")[gender]) |>
    dplyr::mutate(Analysis = factor(Analysis, levels = c("Overall", "Ethnicity", "Age Cat", "Region"))) |>
    dplyr::arrange(Analysis, Subgroup) |>
    dplyr::mutate(
      grouping = dplyr::case_when(
        Analysis == "Overall" ~ "Entire cohort",
        Analysis == "Age Cat" ~ paste0("Aged ", Subgroup),
        TRUE ~ Subgroup
      ),
      grouping = paste0(
        toupper(substr(grouping, 1, 1)),
        substring(grouping, 2)
      ),
      mean = `Optimism-adjusted performance`,
      lower = mean,
      upper = mean,
      Cstat = round(mean, 3)
    ) |>
    dplyr::select(grouping, N, Cstat, mean, lower, upper)

  ### Identical lower and upper values draw only the estimate markers
  ### Create forest plot with the previous formatting
  plot_out <- forestplot::forestplot(
    discrimination_combined,
    labeltext = c(grouping, N, Cstat),
    txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 1)),
    xticks = c(0.55, 0.65, 0.75, 0.85, 0.95),
    boxsize = 0.2,
    zero = discrimination_combined$mean[1],
    graphwidth = grid::unit(2.5, "cm")
  ) |>
    forestplot::fp_add_header(
      "grouping" = "Grouping",
      "N" = "N",
      "Cstat" = "C-index"
    )

  ### Save plot without overwriting the previous forest plot
  ragg::agg_png(
    paste0("figures/p5/discrimination_forestplot_optimism_adjusted_", gender, ".png"),
    width = 6,
    height = 8,
    unit = "in",
    res = 600
  )
  plot(plot_out)
  grDevices::dev.off()
}

### Create forest plots
create_forestplot(gender = 1)
create_forestplot(gender = 2)
