###
### Create internal, bootstrap, instability and temporal validation calibration plots
###
rm(list = ls())
Sys.time()
setwd("/mnt/bmh01-rds/Sperrin_CHARIOT_CPRD/alex/project2/")
getwd()

### Load package used for arranging the plots
library(patchwork)

### Settings
alpha <- 0.15
results_dir <- "data/p5"

### Choose bootstrap iterations to plot 
i <- 1
j <- 500

### Set seed for random sampling for histograms (randomly select 200,000 to speed up computation)
set.seed(1)

### Stop if any requested files are missing
check_files <- function(files, description) {
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop(paste0("Missing ", description, " file(s):\n", paste(missing_files, collapse = "\n")))
  }
}

### Extract data saved by the calibration functions
extract_ph <- function(x) x$calib_data_plot
extract_km <- function(x) x$plot.data

### Interpolate curves onto a common grid and calculate their pointwise median
### This is so we can calculate median calibration curve across bootstrap iterations
median_curve <- function(curves, n = 250) {
  lower <- max(vapply(curves, function(x) min(x$pred, na.rm = TRUE), numeric(1)))
  upper <- min(vapply(curves, function(x) max(x$pred, na.rm = TRUE), numeric(1)))
  pred <- seq(lower, upper, length.out = n)
  pred_obs <- vapply(curves, function(x) {
    stats::approx(x$pred, x$pred.obs, xout = pred, ties = "ordered")$y
  }, numeric(n))
  data.frame(pred = pred, pred.obs = apply(pred_obs, 1, median, na.rm = TRUE))
}

### Calculate medians of grouped KM points across imputations
median_km <- function(km_list) {
  data.frame(
    pred = apply(do.call(cbind, lapply(km_list, function(x) x$pred)), 1, median, na.rm = TRUE),
    obs = apply(do.call(cbind, lapply(km_list, function(x) x$obs)), 1, median, na.rm = TRUE)
  )
}

### Function to extract the full distribution of predictions risks, which is used
### for making the histogram. The plots will be based on a subset of equally spaced points,
### meaning we can't create a histogram from this
get_risks <- function(calib_object) {
  risks <- calib_object$calib_data$pred
  if (is.null(risks)) stop("No patient-level predicted risks in calib_data$pred")
  risks <- risks[is.finite(risks)]
  if (!length(risks) || any(risks < 0 | risks > 1)) {
    stop("No valid predicted-risk distribution (expected risks between 0 and 1)")
  }
  risks
}

### Maximum of 0.4 or 99th percentile for axes
axis_limit <- function(risks) {
  limit <- max(0.5, as.numeric(stats::quantile(risks, 0.99, names = FALSE)))
  if (!is.finite(limit) || limit <= 0) stop("Invalid 99th percentile of predicted risks")
  limit
}

### Create calibration plot with grouped KM results overlaid
make_calibration_plot <- function(curve, km, limit, median_results = TRUE) {
  curve$Result <- if (median_results) "Median calibration curve" else "Calibration curve"
  km$Result <- if (median_results) "Median grouped observed vs expected" else "Grouped observed vs expected"
  ggplot2::ggplot() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line(data = curve, ggplot2::aes(pred, pred.obs, colour = Result), linewidth = 0.9) +
    ggplot2::geom_point(data = km, ggplot2::aes(pred, obs, colour = Result), size = 0.8) +
    ggplot2::scale_colour_manual(values = c(
      "Median calibration curve" = "#0072B2",
      "Median grouped observed vs expected" = "#D55E00",
      "Calibration curve" = "#0072B2",
      "Grouped observed vs expected" = "#D55E00"
    )) +
    ggplot2::coord_cartesian(xlim = c(0, limit), ylim = c(0, limit), expand = FALSE) +
    ggplot2::labs(x = "Predicted risk", y = "Observed risk", colour = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 8),
                   axis.title = ggplot2::element_text(size = 8),
                   legend.position = "bottom", 
                   plot.margin = ggplot2::margin(0, 5.5, 5.5, 5.5))
}

### Bin *all* patients, with fixed breaks spanning the displayed x range.
### Avoids a random sample obscuring a sparse right tail.
make_histogram <- function(risks, title, limit) {
  if (length(risks) > 200000){
    risks <- sample(risks, 200000)
  }
  breaks <- seq(0, limit, length.out = 51)
  shown <- risks[risks >= 0 & risks <= limit]
  counts <- hist(shown, breaks = breaks, plot = FALSE, include.lowest = TRUE)$counts
  stopifnot(sum(counts) == length(shown))
  bars <- data.frame(pred = head(breaks, -1) + diff(breaks) / 2, n = counts)
  ggplot2::ggplot(bars, ggplot2::aes(x = pred, y = n)) +
    ggplot2::geom_col(width = limit / 50, fill = "grey65", colour = "white", linewidth = 0.15) +
    ggplot2::coord_cartesian(xlim = c(0, limit), expand = FALSE) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9, face = "bold", hjust = 0.5),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 0, 5.5)
    )
}

### Stack aligned histogram and calibration panels at a 1:5 height ratio
with_histogram <- function(histogram, calibration) {
  patchwork::wrap_plots(histogram, calibration, ncol = 1,
                        heights = c(1, 5), guides = "auto") &
    ggplot2::theme(legend.position = "bottom")
}

### Save all figures using the ragg PNG device
save_plot <- function(filename, plot, width, height) {
  ggplot2::ggsave(filename, plot, device = ragg::agg_png,
                  width = width, height = height, units = "in", dpi = 150)
}

############################################
### Overall internal-validation calibration
############################################
overall_plots <- vector("list", 2)
for (gender in c(2,1)) {
  # Create filepaths for results (curves and grouped km)
  ph_files <- file.path(results_dir, paste0("calib_ph_", gender, "_", 1:10, ".rds"))
  km_files <- file.path(results_dir, paste0("calib_km_group_", gender, "_", 1:10, ".rds"))
  # Check they exist
  check_files(ph_files, "internal PH calibration")
  check_files(km_files, "internal grouped KM calibration")
  # Read in results
  ph <- lapply(ph_files, readRDS)
  km <- lapply(km_files, readRDS)
  # Extract relevant data
  curve <- median_curve(lapply(ph, extract_ph))
  km_median <- median_km(lapply(km, extract_km))
  # Extracts risks for limits and histogram
  risks <- get_risks(ph[[1]])
  limit <- axis_limit(risks)
  # Create plotes
  overall_plots[[gender]] <- with_histogram(
    make_histogram(risks, c("Male", "Female")[gender], limit),
    make_calibration_plot(curve, km_median, limit)
  )
}
# Wrap plots
overall_plot <- patchwork::wrap_plots(overall_plots[c(2,1)], ncol = 2, guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
# Save
save_plot("figures/p5/internal_validation_calibration.png", overall_plot, 9, 5)

### Save Figure 3 as a 600 dpi PNG
ragg::agg_png("figures/p5/Figure3.png", width = 9, height = 5, units = "in", res = 600)
### Render and close the graphics device
print(overall_plot)
grDevices::dev.off()


#######################################
### Internal validation by subgroup
#######################################
for (gender in 1:2) {
  for (byvar in c("ethnicity", "region", "age_cat")) {
    # Create filepaths for results (curves and grouped km)
    ph_files <- file.path(results_dir, paste0("calib_ph_", gender, "_", 1:10, "_", byvar, ".rds"))
    km_files <- file.path(results_dir, paste0("calib_km_group_", gender, "_", 1:10, "_", byvar, ".rds"))
    # Check they exist
    check_files(ph_files, paste("subgroup PH calibration", byvar))
    check_files(km_files, paste("subgroup grouped KM calibration", byvar))
    # Read in results
    ph <- lapply(ph_files, readRDS)
    km <- lapply(km_files, readRDS)
    # Get levels for 'by' variable (subgroup names)
    level_file <- file.path(results_dir, paste0("discrim_", gender, "_1_", byvar, ".rds"))
    check_files(level_file, paste("subgroup labels", byvar))
    subgroup_names <- as.character(readRDS(level_file)$levels)
    # Create panels
    panels <- lapply(seq_along(subgroup_names), function(k) {
      # Extract relevant data
      curve <- median_curve(lapply(ph, function(x) extract_ph(x[[k]])))
      km_median <- median_km(lapply(km, function(x) extract_km(x[[k]])))
      # Extract risks for histogram and axes
      risks <- get_risks(ph[[1]][[k]])
      limit <- axis_limit(risks)
      # Create plot
      with_histogram(
        make_histogram(risks, subgroup_names[k], limit),
        make_calibration_plot(curve, km_median, limit)
      )
    })
    ### No more than three panels per row
    ncol <- min(3L, length(panels))
    nrow <- ceiling(length(panels) / ncol)
    # Combine plots
    plot_out <- patchwork::wrap_plots(panels, ncol = ncol, guides = "collect") &
      ggplot2::theme(legend.position = "bottom")
    # Save
    save_plot(paste0("figures/p5/internal_validation_calibration_",
                     tolower(c("Male", "Female")[gender]), "_by_", byvar, ".png"),
              plot_out, 3.8 * ncol, 4.1 * nrow)
  }
}

######################################
### Bootstrap test calibration curves
######################################
bootstrap_plots <- vector("list", 2)
iterations <- i:j
for (gender in 1:2) {
  # Get filepaths
  files <- file.path(results_dir, paste0("optimism_boot_calib_ph_data_range_", gender,
                                         "_iter", iterations, "_dforig.rds"))
  # Check they exist
  check_files(files, "bootstrap test calibration")
  # Read in curve data
  curves <- do.call(rbind, lapply(seq_along(files), function(k) {
    x <- readRDS(files[k])
    x$iteration <- iterations[k]
    x
  }))
  ### Bootstrap test curves have no patient-level risk distribution; use the
  ### first internal-validation imputation for creating histograms.
  reference_file <- file.path(results_dir, paste0("calib_ph_", gender, "_1.rds"))
  check_files(reference_file, "bootstrap reference histogram")
  reference_risks <- get_risks(readRDS(reference_file))
  limit <- axis_limit(reference_risks)
  # Plot curves
  calibration <- ggplot2::ggplot(curves, ggplot2::aes(pred, pred.obs, group = iteration)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_line(alpha = alpha, colour = "#0072B2", linewidth = 0.4) +
    ggplot2::coord_cartesian(xlim = c(0, limit), ylim = c(0, limit), expand = FALSE) +
    ggplot2::labs(x = "Predicted risk", y = "Observed risk") +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 5.5, 5.5, 5.5))
  # Combine with histogram
  bootstrap_plots[[gender]] <- with_histogram(
    make_histogram(reference_risks, c("Male", "Female")[gender], limit),
    calibration
  )
}
# Combine plots and save
save_plot("figures/p5/bootstrap_test_calibration_curves.png",
          patchwork::wrap_plots(bootstrap_plots, ncol = 2), 9, 5)

#######################################
### Temporal validation, baseline layer
#######################################
for (gender in 1:2) {
  for (t_eval_int in c(5,10)){
  # Create filepaths
  ph_files <- file.path(results_dir, paste0("temporalv_calib_ph_", gender, "_t", 1:5, "_teval", t_eval_int, ".rds"))
  km_files <- file.path(results_dir, paste0("temporalv_calib_km_group_", gender, "_t", 1:5, "_teval", t_eval_int, ".rds"))
  # Check they exist
  check_files(ph_files, "temporal PH calibration")
  check_files(km_files, "temporal grouped KM calibration")
  # Create panels for t = 1:5
  panels <- lapply(1:5, function(t) {
    # Read in curve data
    ph <- readRDS(ph_files[t])
    # Extract relevant bits
    curve <- extract_ph(ph)
    curve <- curve[order(curve$pred), ]
    # Read in km data
    km_object <- readRDS(km_files[t])
    # Extract relevant bits
    km <- extract_km(km_object)
    # Get risk distribution for histograms and limits
    risks <- get_risks(ph)
    limit <- axis_limit(risks)
    # Create plot
    with_histogram(
      make_histogram(risks, paste0("Follow-up visit ", t, ", t = ", t_eval_int, " years"), limit),
      make_calibration_plot(curve, km, limit, median_results = FALSE)
    )
  })
  # Combine plots
  plot_out <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  # Save
  save_plot(paste0("figures/p5/temporal_validation_calibration_",
                   tolower(c("Male", "Female")[gender]), "_teval", t_eval_int, ".png"), plot_out, 11.4, 8.2)
  }
}


#######################################
### Temporal validation, interventional
#######################################
for (gender in 1:2) {
  for (t_eval_int in c(5,10)){
    for (mfr in c("sbp", "nonhdl", "bmi", "smoking")){
      # Create file paths
      ph_files <- file.path(results_dir, paste0("temporalv_interventional_calib_ph_", gender, "_t", 1:5, "_", mfr, "_teval", t_eval_int, ".rds"))
      km_files <- file.path(results_dir, paste0("temporalv_interventional_calib_km_group_", gender, "_t", 1:5, "_", mfr, "_teval", t_eval_int, ".rds"))
      # Check files exist
      check_files(ph_files, "temporal PH calibration")
      check_files(km_files, "temporal grouped KM calibration")
      # Create panels for t = 1:5
      panels <- lapply(1:5, function(t) {
        # Read in curve data
        ph <- readRDS(ph_files[t])
        # Extract relevant bits
        curve <- extract_ph(ph)
        curve <- curve[order(curve$pred), ]
        # Read in km data
        km_object <- readRDS(km_files[t])
        # Extract relevant bits
        km <- extract_km(km_object)
        # Get risk distribution for histograms and limits
        risks <- get_risks(ph)
        limit <- axis_limit(risks)
        # Create plot
        with_histogram(
          make_histogram(risks, paste0("Follow-up visit ", t, ", t = ", t_eval_int, " years"), limit),
          make_calibration_plot(curve, km, limit, median_results = FALSE)
        )
      })
      # Combie plots
      plot_out <- patchwork::wrap_plots(panels, ncol = 3, guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
      # Save
      save_plot(paste0("figures/p5/temporal_validation_calibration_interventional_",
                       tolower(c("Male", "Female")[gender]), "_", mfr, "_teval", t_eval_int, ".png"), plot_out, 11.4, 8.2)
    }
  }
}
print(paste("FINISHED", Sys.time()))
