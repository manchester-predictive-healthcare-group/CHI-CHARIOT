###
### Create Figures for manuscript
###

### Prelim
rm(list = ls())
Sys.time()

setwd()
getwd()

results_dir <- "data/p5"

### Load package used for arranging plots
library(patchwork)

###
### Figure 4
###

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
    paste0("figures/p5/Figure4.png"),
    width = 6,
    height = 8,
    unit = "in",
    res = 600
  )
  plot(plot_out)
  grDevices::dev.off()
}

### Create forest plots
create_forestplot(gender = 2)


###
### Figure 5 (temporal validation, baseline layer, entire cohort, 1/3/5 year follow-up visits)
###

### Set gender
gender <- 2

### Extract PH calibration curve and grouped KM data from saved objects
extract_ph <- function(x) x$calib_data_plot
extract_km <- function(x) x$plot.data

### Extract finite patient-level predicted risks
get_risks <- function(calib_object) {
  ### Read risks from the calibration object
  risks <- calib_object$calib_data$pred
  if (is.null(risks)) stop("No patient-level predicted risks in calib_data$pred")
  ### Retain finite risks within the expected range
  risks <- risks[is.finite(risks)]
  if (!length(risks) || any(risks < 0 | risks > 1)) {
    stop("No valid predicted-risk distribution")
  }
  risks
}

### Use at least 0.4 or the 99th percentile as the axis maximum
axis_limit <- function(risks, default_lim = 0.4) {
    limit <- max(default_lim, as.numeric(stats::quantile(risks, 0.99, names = FALSE)))
    if (!is.finite(limit) || limit <= 0) stop("Invalid 99th percentile of predicted risks")
    limit
}

### Create histogram of predicted risks
make_histogram <- function(risks, title, limit) {
  ### Sample at most 200,000 risks for plotting
  if (length(risks) > 200000) {
    risks <- sample(risks, 200000)
  }
  ### Count displayed risks in 50 equal-width bins
  breaks <- seq(0, limit, length.out = 51)
  shown <- risks[risks >= 0 & risks <= limit]
  counts <- hist(shown, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts
  ### Position each bar at its bin midpoint
  bars <- data.frame(pred = head(breaks, -1) + diff(breaks) / 2, n = counts)
  ### Draw the histogram without axes or grid lines
  ggplot2::ggplot(bars, ggplot2::aes(x = pred, y = n)) +
    ggplot2::geom_col(width = limit / 50, fill = "grey65", colour = "white", linewidth = 0.15) +
    ggplot2::coord_cartesian(xlim = c(0, limit), expand = FALSE) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 8, hjust = 0.5),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 0, 5.5)
    )
}

### Create calibration plot
make_calibration_plot <- function(curve, km, limit) {
  ### Label the curve and grouped observed-risk points
  curve$Result <- "Calibration curve"
  km$Result <- "Grouped observed vs expected"
  ### Draw perfect calibration, the PH curve, and grouped KM points
  ggplot2::ggplot() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line(data = curve, ggplot2::aes(pred, pred.obs, colour = Result), linewidth = 0.9) +
    ggplot2::geom_point(data = km, ggplot2::aes(pred, obs, colour = Result), size = 0.8) +
    ### Use the same colours for both figures
    ggplot2::scale_colour_manual(values = c(
      "Calibration curve" = "#0072B2",
      "Grouped observed vs expected" = "#D55E00"
    )) +
    ### Match predicted- and observed-risk axis ranges
    ggplot2::coord_cartesian(xlim = c(0, limit), ylim = c(0, limit), expand = FALSE) +
    ggplot2::labs(x = "Predicted risk", y = "Observed risk", colour = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom",
                   axis.text = ggplot2::element_text(size = 8),
                   axis.title = ggplot2::element_text(size = 8))
}

### Stack histogram above calibration plot at a 1:5 height ratio
with_histogram <- function(histogram, calibration) {
  patchwork::wrap_plots(histogram, calibration, ncol = 1, heights = c(1, 5))  &
    ggplot2::theme(legend.position = "bottom")
}


### Figure 5: baseline layer, female cohort, follow-up visits 1, 3 and 5

### Select the three follow-up visits
follow_up_visits <- c(1, 3, 5)
t_eval_int <- 5
### Construct PH and KM filenames in visit order
ph_files <- file.path(results_dir, paste0("temporalv_calib_ph_", gender, "_t", follow_up_visits, "_teval", t_eval_int, ".rds"))
km_files <- file.path(results_dir, paste0("temporalv_calib_km_group_", gender, "_t", follow_up_visits, "_teval", t_eval_int, ".rds"))

### Build one histogram and calibration panel per follow-up visit
panels <- lapply(seq_along(follow_up_visits), function(i) {
  ### Match the visit to its corresponding file position
  t <- follow_up_visits[i]
  ### Read the PH and grouped KM objects
  ph <- readRDS(ph_files[i])
  km_object <- readRDS(km_files[i])
  ### Extract and order the calibration curve
  curve <- extract_ph(ph)
  curve <- curve[order(curve$pred), ]
  ### Extract grouped KM points and calculate the plot limit
  km <- extract_km(km_object)
  risks <- get_risks(ph)
  limit <- axis_limit(risks, default_lim = 0.31)
  ### Stack the risk histogram over the calibration plot
  with_histogram(
    make_histogram(risks, paste0("Follow-up visit ", t, ", t = 5 years"), limit),
    make_calibration_plot(curve, km, limit)
  )
})

### Arrange Figure 5 in one row
plot_out <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
### Save Figure 5 as a 600 dpi PNG
ragg::agg_png("figures/p5/Figure5.png", width = 9, height = 3.6, units = "in", res = 600)
### Render and close the graphics device
print(plot_out)
grDevices::dev.off()

###
### Figure 6: interventional layer, female cohort, follow-up visit 1, four MFRs
###

### Select visit 1 and the four modifiable risk factors
t <- 1
t_eval_int <- 5
mfrs <- c("sbp", "nonhdl", "bmi", "smoking")

### Build one histogram and calibration panel per risk factor
panels <- lapply(mfrs, function(mfr) {
  ### Read the matching interventional PH and grouped KM objects
  ph <- readRDS(file.path(results_dir, paste0(
    "temporalv_interventional_calib_ph_", gender, "_t", t, "_", mfr, "_teval", t_eval_int, ".rds"
  )))
  km_object <- readRDS(file.path(results_dir, paste0(
    "temporalv_interventional_calib_km_group_", gender, "_t", t, "_", mfr, "_teval", t_eval_int, ".rds"
  )))
  ### Extract and order the calibration curve
  curve <- extract_ph(ph)
  curve <- curve[order(curve$pred), ]
  ### Extract grouped KM points and calculate the plot limit
  km <- extract_km(km_object)
  risks <- get_risks(ph)
  limit <- axis_limit(risks, default_lim = 0.31)
  ### Stack the risk histogram over the calibration plot
  with_histogram(
    make_histogram(risks, paste0("MFR = ", mfr, ", follow-up visit 1, t = 5 years"), limit),
    make_calibration_plot(curve, km, limit)
  )
})

### Arrange Figure 6 in a 2 x 2 grid
plot_out <- patchwork::wrap_plots(panels, ncol = 2, guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
### Save Figure 6 as a 600 dpi PNG
ragg::agg_png("figures/p5/Figure6.png", width = 6, height = 7.2, units = "in", res = 600)
### Render and close the graphics device
print(plot_out)
grDevices::dev.off()