# =============================================================================
# CODE 1 — STANDARDIZATION MODEL BUILDING
# =============================================================================
# PURPOSE:
#   Build GAM models to standardize canopy height model (CHM) heights against
#   field-measured dominant height from  forest inventory (FI) plots.
#   Two models are built:
#     - With temporal variables    (for diagnostic / comparison)
#     - Without temporal variables (for operational application)
#
# INPUTS:
#   - FI plot locations           (GeoPackage)
#   - FI tree-level measurements  (GeoPackage)
#   - Annual CHM rasters           (GeoTIFF, one per year)
#   - Flight tile metadata         (GeoPackage with acquisition dates)
#   - Leaf-off masks               (GeoPackage, one per year, optional)
#   - Ancillary rasters: composition, broadleaf ratio, topography, slope, DEM
#
# OUTPUTS:
#   - model_temp_inc.rds  (model with temporal variables)
#   - model_temp_exc.rds  (model without temporal variables — used in Codes 2-3)
#
# NOTE:
#   This script is designed to be run section by section.
#   Adjust file paths in Section 1.2 to match your folder structure.
# =============================================================================


# =============================================================================
# 1. SETUP
# =============================================================================

# --- 1.1 Libraries ---
library(tidyverse)
library(tidyterra)
library(terra)
library(sf)
library(lubridate)
library(mgcv)
library(blockCV)

# --- 1.2 File paths ---
# >> Adjust these to your local folder structure <<

# FI data
path_field_plots    <- "data/field_plots.gpkg"           # FI plot locations
path_field_trees    <- "data/trees.gpkg"                 # tree-level measurements

# CHM data
path_chm_folder     <- "data/chm/"                       # folder with mnh{year}.tif
path_flight_metadata <- "data/flight_tiles.gpkg"         # tile-level acquisition dates

# Leaf-off masks (set to NULL if not available)
path_leafoff_folder <- "data/leafoff_masks/"             # folder with leafoff_mask_{year}.gpkg

# Ancillary rasters (set to NULL if unavailable)
path_composition    <- "data/ancillary/composition.tif"
path_broadleaf_ratio <- "data/ancillary/broadleaf_ratio.tif"
path_topography     <- "data/ancillary/topography.tif"
path_slope          <- "data/ancillary/slope.tif"
path_dem            <- "data/ancillary/dem.tif"

# Output folder
path_output         <- "outputs/data/"

# --- 1.3 Parameters ---

# Coordinate reference system (Belgian Lambert 72)
crs_target <- "EPSG:31370"

# NFI plot radius (m) — used for buffering point plots before extraction
field_plot_radius <- 18

# Maximum allowed time difference (days) between field visit and flight
max_date_diff <- 365

# Years where the source is airborne lidar (vs. DAP for the rest)
lidar_years <- c(2014, 2021)

# Spatial cross-validation
cv_folds  <- 10
cv_repeats <- 5
spatial_block_size <- 16636   # block size in meters (approx. range of spatial autocorrelation)

# Forward-backward stepwise selection
improvement_threshold <- 0.01  # minimum RMSE improvement to keep a variable

# Candidate predictor variables
candidate_vars_temp_inc <- c(
  "CHM_height", "CHM_hsd", "CHM_hmean", "CHM_gap_fraction",
  "compo", "broadleaf_ratio", "topo", "slope", "alt",
  "julian_day", "month_flight", "season_flight", "year_factor"
)

candidate_vars_temp_exc <- c(
  "CHM_height", "CHM_hsd", "CHM_hmean", "CHM_gap_fraction",
  "compo", "broadleaf_ratio", "topo", "slope", "alt", "sensor"
)

# Variables that receive a smooth term s() in the GAM
smooth_vars <- c(
  "CHM_height", "CHM_hsd", "CHM_hmean", "CHM_gap_fraction",
  "broadleaf_ratio", "slope", "alt", "julian_day", "month_flight"
)


# =============================================================================
# 2. FUNCTIONS
# =============================================================================

# --- 2.1 CHM metric extraction ---

# Calculate canopy gap fraction (% of pixels below 50% of the 99th percentile)
calculate_canopy_openness <- function(vals, thr = 0.5) {
  h99 <- quantile(vals, 0.99, na.rm = TRUE)
  if (is.na(h99) || h99 == 0) return(NA_real_)
  threshold <- thr * h99
  sum(vals < threshold, na.rm = TRUE) / sum(!is.na(vals)) * 100
}

# Compute CHM summary statistics within a buffered plot
chm_summary_stats <- function(vals) {
  c(
    CHM_height       = as.numeric(quantile(vals, 0.95, na.rm = TRUE)),
    CHM_hsd          = sd(vals, na.rm = TRUE),
    CHM_hmean        = mean(vals, na.rm = TRUE),
    CHM_gap_fraction = calculate_canopy_openness(vals)
  )
}

# --- 2.2 Helper functions ---

# Extract the four-digit year from a filename
extract_year_from_filename <- function(filename) {
  as.integer(regmatches(basename(filename), regexpr("\\d{4}", basename(filename))))
}

# Map numeric composition codes to species names
map_composition_to_species <- function(code, broadleaf_ratio) {
  case_when(
    code ==1 ~ "Other_broadleaf",
    code == 2  ~ "Birch",
    code == 3  ~ "Oak",
    code == 4  ~ "Douglas_fir",
    code == 5  ~ "Spruce",
    code == 6  ~ "Beech",
    code == 7  ~ "Larch",
    code == 8  ~ "Other_broadleaf", # Because poplars are underrepresented
    code == 9  ~ "Pine",
    code == 10 ~ "Other_coniferous",
    TRUE ~ NA_character_
  )
}

# Safely extract values from a raster for buffered plots
# Returns NULL if the raster file does not exist.
safe_extract <- function(plots, raster_path, fun, col_name, crs) {
  if (is.null(raster_path) || !file.exists(raster_path)) return(NULL)
  rast_obj <- terra::rast(raster_path)
  crs(rast_obj) <- crs
  result <- terra::extract(rast_obj, plots, fun = fun, na.rm = TRUE, bind = TRUE)
  result_df <- as.data.frame(result)
  extracted_col <- names(result_df)[ncol(result_df)]
  result_df %>%
    rename(!!col_name := all_of(extracted_col)) %>%
    select(plot_ID, all_of(col_name))
}

# Modal value (most frequent class) for categorical raster extraction
modal_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# --- 2.3 Model performance ---

# Calculate RMSE, relative RMSE, regression slope and R² between actual and predicted
calculate_performance_metrics <- function(actual, predicted) {
  valid <- !is.na(actual) & !is.na(predicted)
  actual <- actual[valid]
  predicted <- predicted[valid]
  if (length(actual) == 0) return(list(rmse = NA, rel_rmse = NA, slope = NA, r2 = NA))
  rmse     <- sqrt(mean((actual - predicted)^2))
  rel_rmse <- (rmse / mean(actual)) * 100
  lm_fit   <- lm(actual ~ predicted)
  list(
    rmse     = rmse,
    rel_rmse = rel_rmse,
    slope    = coef(lm_fit)[2],
    r2       = summary(lm_fit)$r.squared
  )
}

# --- 2.4 Spatial cross-validation ---

# Create repeated spatial block folds
create_spatial_folds <- function(spatial_data, k, n, block_size) {
  all_folds <- list()
  for (rep in 1:n) {
    set.seed(123 + rep - 1)
    sb <- cv_spatial(
      x = spatial_data, size = block_size,
      k = k, selection = "random", iteration = 100
    )
    all_folds[[rep]] <- sb
  }
  all_folds
}

# Evaluate a GAM formula via spatial cross-validation
# Returns the mean RMSE across all folds and repetitions.
evaluate_formula <- function(formula_str, data, spatial_folds, k, n, bias_var) {
  all_rmse <- numeric()
  for (rep in 1:n) {
    for (fold in 1:k) {
      train_idx <- spatial_folds[[rep]]$folds_list[[fold]][[1]]
      test_idx  <- spatial_folds[[rep]]$folds_list[[fold]][[2]]
      model <- try(
        gam(as.formula(formula_str), data = data[train_idx, ], method = "REML"),
        silent = TRUE
      )
      if (!inherits(model, "try-error")) {
        pred <- predict(model, newdata = data[test_idx, ])
        fold_rmse <- sqrt(mean((data[test_idx, ][[bias_var]] - pred)^2, na.rm = TRUE))
        all_rmse <- c(all_rmse, fold_rmse)
      } else {
        all_rmse <- c(all_rmse, Inf)
      }
    }
  }
  mean(all_rmse, na.rm = TRUE)
}

# --- 2.5 Stepwise variable selection ---

# Forward-backward stepwise variable selection using spatial CV
run_stepwise_selection <- function(data, candidate_vars, smooth_vars, spatial_folds,
                                   k, n, bias_var, improvement_threshold, max_vars) {
  current_vars    <- character(0)
  current_formula <- paste0(bias_var, " ~ 1")
  current_rmse    <- evaluate_formula(current_formula, data, spatial_folds, k, n, bias_var)

  cat("  Baseline RMSE:", round(current_rmse, 4), "\n")

  repeat {
    improved <- FALSE

    # --- Forward step: try adding each remaining variable ---
    remaining_vars <- setdiff(candidate_vars, current_vars)
    if (length(remaining_vars) > 0 && length(current_vars) < max_vars) {
      best_fwd_rmse <- Inf
      best_fwd_var  <- NULL
      for (var in remaining_vars) {
        test_vars <- c(current_vars, var)
        terms <- sapply(test_vars, function(v) {
          if (v %in% smooth_vars) paste0("s(", v, ")") else v
        })
        test_formula <- paste0(bias_var, " ~ ", paste(terms, collapse = " + "))
        test_rmse <- evaluate_formula(test_formula, data, spatial_folds, k, n, bias_var)
        cat("    + ", var, " -> RMSE:", round(test_rmse, 4), "\n")
        if (test_rmse < best_fwd_rmse) {
          best_fwd_rmse <- test_rmse
          best_fwd_var  <- var
        }
      }
      if (best_fwd_rmse < current_rmse - improvement_threshold) {
        current_vars <- c(current_vars, best_fwd_var)
        current_rmse <- best_fwd_rmse
        improved <- TRUE
        cat("    >> Added:", best_fwd_var, "\n")
      }
    }

    # --- Backward step: try removing each current variable ---
    if (length(current_vars) > 0) {
      best_bwd_rmse <- Inf
      best_bwd_var  <- NULL
      for (var in current_vars) {
        test_vars <- setdiff(current_vars, var)
        if (length(test_vars) == 0) {
          test_formula <- paste0(bias_var, " ~ 1")
        } else {
          terms <- sapply(test_vars, function(v) {
            if (v %in% smooth_vars) paste0("s(", v, ")") else v
          })
          test_formula <- paste0(bias_var, " ~ ", paste(terms, collapse = " + "))
        }
        test_rmse <- evaluate_formula(test_formula, data, spatial_folds, k, n, bias_var)
        cat("    - ", var, " -> RMSE:", round(test_rmse, 4), "\n")
        if (test_rmse < best_bwd_rmse) {
          best_bwd_rmse <- test_rmse
          best_bwd_var  <- var
        }
      }
      if (best_bwd_rmse < current_rmse - improvement_threshold) {
        current_vars <- setdiff(current_vars, best_bwd_var)
        current_rmse <- best_bwd_rmse
        improved <- TRUE
        cat("    >> Removed:", best_bwd_var, "\n")
      }
    }

    if (!improved) break
  }

  # Build final formula
  if (length(current_vars) == 0) {
    best_formula <- paste0(bias_var, " ~ 1")
  } else {
    terms <- sapply(current_vars, function(v) {
      if (v %in% smooth_vars) paste0("s(", v, ")") else v
    })
    best_formula <- paste0(bias_var, " ~ ", paste(terms, collapse = " + "))
  }

  list(
    selected_vars = current_vars,
    best_formula  = best_formula,
    final_rmse    = current_rmse
  )
}

# --- 2.6 Final model fitting & evaluation ---

# Fit the final GAM model, compute corrected heights and performance metrics
fit_and_evaluate_model <- function(data, formula_str, predictor_var, response_var, bias_var) {
  final_model <- gam(as.formula(formula_str), data = data, method = "REML")
  data$predicted_bias   <- predict(final_model, newdata = data)
  data$corrected_height <- data[[predictor_var]] + data$predicted_bias

  uncorrected <- calculate_performance_metrics(data[[response_var]], data[[predictor_var]])
  corrected   <- calculate_performance_metrics(data[[response_var]], data$corrected_height)

  list(
    final_model = final_model,
    corrected_data = data,
    evaluation = list(uncorrected = uncorrected, corrected = corrected)
  )
}

# Collect holdout predictions from a final CV pass (for diagnostic tables)
collect_cv_predictions <- function(data, formula_str, spatial_folds, k, n,
                                   predictor_var, response_var, bias_var) {
  cv_predictions <- list()

  for (rep in 1:n) {
    for (fold in 1:k) {
      train_idx <- spatial_folds[[rep]]$folds_list[[fold]][[1]]
      test_idx  <- spatial_folds[[rep]]$folds_list[[fold]][[2]]

      model_cv <- try(
        gam(as.formula(formula_str), data = data[train_idx, ], method = "REML"),
        silent = TRUE
      )

      if (!inherits(model_cv, "try-error")) {
        pred_bias <- predict(model_cv, newdata = data[test_idx, ])
        pred_df <- data[test_idx, ] %>%
          st_drop_geometry() %>%
          mutate(
            predicted_bias   = pred_bias,
            corrected_height = .data[[predictor_var]] + pred_bias,
            residual_bias    = .data[[bias_var]] - pred_bias,
            fold = fold,
            rep  = rep
          )
        cv_predictions[[length(cv_predictions) + 1]] <- pred_df
      }
    }
  }

  bind_rows(cv_predictions)
}

# Compute per-group CV metrics (used for species and year breakdowns)
compute_group_cv_metrics <- function(cv_preds, group_var, predictor_var, response_var) {
  cv_preds %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n_obs            = n(),
      mean_bias_before = mean(.data[[predictor_var]] - .data[[response_var]], na.rm = TRUE),
      rmse_before      = sqrt(mean((.data[[response_var]] - .data[[predictor_var]])^2, na.rm = TRUE)),
      r2_before        = cor(.data[[response_var]], .data[[predictor_var]])^2,
      mean_bias_after  = mean(corrected_height - .data[[response_var]], na.rm = TRUE),
      rmse_after       = sqrt(mean((.data[[response_var]] - corrected_height)^2, na.rm = TRUE)),
      r2_after         = cor(.data[[response_var]], corrected_height)^2,
      .groups = "drop"
    )
}


# =============================================================================
# 3. DATA ASSEMBLY
# =============================================================================
cat("=== 3. DATA ASSEMBLY ===\n")

# --- 3.1 Load NFI field data ---
cat("Loading NFI field data...\n")
plots_raw <- vect(path_field_plots) %>%
  as.data.frame() %>%
  mutate(plot_ID = paste0(IGN, "_", NPL))

trees_raw <- vect(path_field_trees) %>%
  as.data.frame() %>%
  mutate(plot_ID = paste0(IGN, "_", NPL))

# --- 3.2 Compute dominant height ---
# Dominant height = mean height of the 9 tallest trees per plot (in meters)
cat("Computing dominant height...\n")
height_by_plot <- trees_raw %>%
  group_by(plot_ID) %>%
  arrange(desc(HTot)) %>%
  slice_head(n = 9) %>%
  summarise(dominant_height = mean(HTot, na.rm = TRUE) / 100, .groups = "drop") %>%
  filter(!is.nan(dominant_height))

cat("  ", nrow(height_by_plot), "plots with valid dominant height\n")

# --- 3.3 Match plots to closest flight dates ---
cat("Matching plots to flight acquisition dates...\n")
field_plots_vect <- terra::vect(path_field_plots) %>% 
  terra::project(crs_target) %>%
  mutate(plot_ID = paste0(IGN, "_", NPL))
flight_meta_vect <- terra::vect(path_flight_metadata)

plot_flight_df <- terra::intersect(field_plots_vect, flight_meta_vect) %>%
  as.data.frame() %>%
  mutate(
    date_field = as.Date(date),
    date_flight = as.Date(centredate),
    time_diff  = as.numeric(date_field - date_flight)
  ) %>% 
  distinct()

# Keep the closest flight per plot × year, within the allowed time window
plot_flight_matched <- plot_flight_df %>%
  group_by(plot_ID, year) %>%
  filter(abs(time_diff) == min(abs(time_diff))) %>%
  ungroup() %>%
  filter(abs(time_diff) <= max_date_diff)

# Add temporal variables
plot_flight_matched <- plot_flight_matched %>%
  mutate(
    julian_day    = yday(date_flight),
    month_flight  = month(date_flight),
    season_flight = case_when(
      month_flight %in% c(12, 1, 2)  ~ "winter",
      month_flight %in% c(3, 4, 5)   ~ "spring",
      month_flight %in% c(6, 7, 8)   ~ "summer",
      month_flight %in% c(9, 10, 11) ~ "fall"
    ),
    year_factor = as.factor(year),
    sensor      = as.factor(ifelse(year %in% lidar_years, "lidar", "dap"))
  )

# Join dominant height
plot_flight_matched <- plot_flight_matched %>%
  left_join(height_by_plot, by = "plot_ID") %>%
  filter(!is.na(dominant_height))

cat("  ", n_distinct(plot_flight_matched$plot_ID), "plots matched to flights\n")

# --- 3.4 Flag leaf-off plots ---
cat("Flagging leaf-off plots...\n")
leafoff_files <- if (!is.null(path_leafoff_folder) && dir.exists(path_leafoff_folder)) {
  list.files(path_leafoff_folder, pattern = "leafoff_mask_\\d{4}\\.gpkg$", full.names = TRUE)
} else {
  character(0)
}

plot_flight_matched$leaf_off <- FALSE

if (length(leafoff_files) > 0) {
  matched_plots_vect <- field_plots_vect[field_plots_vect$plot_ID %in% plot_flight_matched$plot_ID]

  for (mask_file in leafoff_files) {
    mask_year <- as.integer(regmatches(basename(mask_file), regexpr("\\d{4}", basename(mask_file))))
    leafoff_poly <- vect(mask_file)

    plots_this_year <- plot_flight_matched %>%
      filter(year == mask_year) %>%
      pull(plot_ID) %>%
      unique()

    if (length(plots_this_year) == 0) next

    plots_subset <- matched_plots_vect[matched_plots_vect$plot_ID %in% plots_this_year]
    plots_cropped <- try(crop(plots_subset, ext(leafoff_poly)), silent = TRUE)
    if (inherits(plots_cropped, "try-error") || nrow(plots_cropped) == 0) next

    intersects <- is.related(plots_cropped, leafoff_poly, relation = "intersects")
    leafoff_ids <- plots_cropped[intersects, ]$plot_ID

    plot_flight_matched$leaf_off[
      plot_flight_matched$plot_ID %in% leafoff_ids & plot_flight_matched$year == mask_year
    ] <- TRUE

    cat("  ", mask_year, ":", sum(intersects), "plots flagged as leaf-off\n")
  }
}

# Filter out leaf-off observations
n_before <- nrow(plot_flight_matched)
plot_flight_matched <- plot_flight_matched %>% filter(!leaf_off)
cat("  Removed", n_before - nrow(plot_flight_matched), "leaf-off observations\n")

# --- 3.5 Extract CHM metrics ---
cat("Extracting CHM metrics...\n")
chm_files <- list.files(path_chm_folder, pattern = "\\d{4}\\.tif$", full.names = TRUE)
chm_years <- sapply(chm_files, extract_year_from_filename)

matched_plot_ids <- unique(plot_flight_matched$plot_ID)
matched_plots_vect <- field_plots_vect[field_plots_vect$plot_ID %in% matched_plot_ids]
plots_buffered <- terra::buffer(matched_plots_vect, field_plot_radius)

chm_metrics_all <- data.frame()
for (i in seq_along(chm_files)) {
  yr <- chm_years[i]
  plots_for_year <- unique(plot_flight_matched$plot_ID[plot_flight_matched$year == yr])
  plots_year <- plots_buffered[plots_buffered$plot_ID %in% plots_for_year]
  if (nrow(plots_year) == 0) next
  
  chm <- terra::rast(chm_files[i])
  crs(chm) <- crs_target
  
  # Extract WITHOUT bind — returns a clean data.frame with ID + metric columns
  # The metric columns inherit the raster layer name (e.g. "mnh2006", "mnh2006.1", ...)
  # so we rename them explicitly based on the order in chm_summary_stats().
  metrics_raw <- terra::extract(chm, plots_year, fun = chm_summary_stats) %>%
    as.data.frame()
  
  # Rename: column 1 = ID (row index), columns 2-5 = the 4 metrics in order
  metric_names <- c("CHM_height", "CHM_hsd", "CHM_hmean", "CHM_gap_fraction")
  names(metrics_raw) <- c("ID", metric_names)
  
  # Map the row index back to plot_ID
  metrics_raw$plot_ID <- plots_year$plot_ID[metrics_raw$ID]
  metrics_raw$year    <- yr
  metrics_raw$ID      <- NULL
  
  chm_metrics_all <- bind_rows(chm_metrics_all, metrics_raw)
  cat("  ", yr, ":", nrow(metrics_raw), "plots\n")
}

# --- 3.6 Extract ancillary raster data ---
cat("Extracting ancillary raster data...\n")
ancillary_data <- data.frame(plot_ID = plots_buffered$plot_ID)

# Topographic position (categorical -> modal value)
if (!is.null(path_topography) && file.exists(path_topography)) {
  rast_topo <- terra::rast(path_topography)
  crs(rast_topo) <- crs_target
  topo_data <- terra::extract(rast_topo, plots_buffered, fun = modal_value, bind = TRUE) %>%
    as.data.frame()
  topo_col <- names(topo_data)[ncol(topo_data)]
  topo_data <- topo_data %>%
    rename(topo = all_of(topo_col)) %>%
    mutate(topo = round(topo)) %>%
    select(plot_ID, topo)
  ancillary_data <- left_join(ancillary_data, topo_data, by = "plot_ID")
  cat("  - topography\n")
}

# Slope (continuous -> mean)
slope_data <- safe_extract(plots_buffered, path_slope, mean, "slope", crs_target)
if (!is.null(slope_data)) {
  ancillary_data <- left_join(ancillary_data, slope_data, by = "plot_ID")
  cat("  - slope\n")
}

# Elevation (continuous -> mean)
alt_data <- safe_extract(plots_buffered, path_dem, mean, "alt", crs_target)
if (!is.null(alt_data)) {
  ancillary_data <- left_join(ancillary_data, alt_data, by = "plot_ID")
  cat("  - elevation\n")
}

# Tree species composition (categorical -> modal value)
if (!is.null(path_composition) && file.exists(path_composition)) {
  rast_compo <- terra::rast(path_composition)
  crs(rast_compo) <- crs_target
  compo_data <- terra::extract(rast_compo, plots_buffered, fun = modal_value, bind = TRUE) %>%
    as.data.frame()
  compo_col <- names(compo_data)[ncol(compo_data)]
  compo_data <- compo_data %>%
    rename(compo = all_of(compo_col)) %>%
    mutate(compo = round(compo)) %>%
    select(plot_ID, compo)
  ancillary_data <- left_join(ancillary_data, compo_data, by = "plot_ID")
  cat("  - composition\n")
}

# Broadleaf ratio (continuous -> mean)
br_data <- safe_extract(plots_buffered, path_broadleaf_ratio, mean, "broadleaf_ratio", crs_target)
if (!is.null(br_data)) {
  ancillary_data <- left_join(ancillary_data, br_data, by = "plot_ID")
  cat("  - broadleaf ratio\n")
}

# --- 3.7 Combine all data ---
cat("Combining all data...\n")
coords_df <- terra::crds(matched_plots_vect) %>%
  as.data.frame() %>%
  mutate(plot_ID = matched_plots_vect$plot_ID) %>%
  rename(x = 1, y = 2)

combined_data <- plot_flight_matched %>%
  left_join(chm_metrics_all, by = c("plot_ID", "year")) %>%
  left_join(ancillary_data, by = "plot_ID") %>%
  left_join(coords_df, by = "plot_ID")

# Compute bias (field - CHM)
combined_data$bias <- combined_data$dominant_height - combined_data$CHM_height

# Encode categorical variables
if ("compo" %in% names(combined_data) && "broadleaf_ratio" %in% names(combined_data)) {
  combined_data$compo <- as.factor(
    map_composition_to_species(combined_data$compo, combined_data$broadleaf_ratio)
  )
}
if ("topo" %in% names(combined_data)) {
  combined_data$topo <- as.factor(combined_data$topo)
}

# Remove rows with missing essential variables
combined_data <- combined_data %>%
  filter(!is.na(CHM_height), !is.na(compo), !is.na(broadleaf_ratio))

cat("  ", nrow(combined_data), "complete observations\n")

# --- 3.8 Flag outliers (Tukey's fences on the bias) ---
cat("Flagging outliers...\n")
q1  <- quantile(combined_data$bias, 0.25, na.rm = TRUE)
q3  <- quantile(combined_data$bias, 0.75, na.rm = TRUE)
iqr <- q3 - q1
combined_data$outlier_tukey <- combined_data$bias < (q1 - 1.5 * iqr) |
                               combined_data$bias > (q3 + 1.5 * iqr)
cat("  ", sum(combined_data$outlier_tukey, na.rm = TRUE), "outliers flagged\n")

# --- 3.9 Convert to spatial (sf) for cross-validation ---
valid_coords <- !is.na(combined_data$x) & !is.na(combined_data$y)
prepared_sf <- combined_data[valid_coords, ] %>%
  st_as_sf(coords = c("x", "y"), crs = crs_target)

cat("  ", nrow(prepared_sf), "spatially valid observations ready for modelling\n")


# =============================================================================
# 4. MODEL BUILDING
# =============================================================================
cat("\n=== 4. MODEL BUILDING ===\n")

# Working dataset: outliers excluded
model_data <- prepared_sf %>% filter(!outlier_tukey)
cat("Using", nrow(model_data), "observations (outliers excluded)\n")

# --- 4.1 Create spatial cross-validation folds ---
cat("Creating spatial folds (this may take a few minutes)...\n")
spatial_folds <- create_spatial_folds(model_data, k = cv_folds, n = cv_repeats,
                                      block_size = spatial_block_size)

# --- 4.2 Model WITH temporal variables (diagnostic) ---
cat("\n--- Model WITH temporal variables ---\n")
selection_temp_inc <- run_stepwise_selection(
  data         = model_data,
  candidate_vars = candidate_vars_temp_inc,
  smooth_vars  = smooth_vars,
  spatial_folds = spatial_folds,
  k = cv_folds, n = cv_repeats,
  bias_var     = "bias",
  improvement_threshold = improvement_threshold,
  max_vars     = length(candidate_vars_temp_inc)
)

 result_temp_inc <- fit_and_evaluate_model(
  model_data, selection_temp_inc$best_formula,
  "CHM_height", "dominant_height", "bias"
)

# Collect CV holdout predictions for diagnostic tables
cv_preds_temp_inc <- collect_cv_predictions(
  model_data, selection_temp_inc$best_formula, spatial_folds,
  cv_folds, cv_repeats, "CHM_height", "dominant_height", "bias"
)
cv_by_species_inc <- compute_group_cv_metrics(cv_preds_temp_inc, "compo", "CHM_height", "dominant_height")
cv_by_year_inc    <- compute_group_cv_metrics(cv_preds_temp_inc, "year_factor", "CHM_height", "dominant_height")

# Bundle all results
result_temp_inc$selection          <- selection_temp_inc
result_temp_inc$best_formula       <- selection_temp_inc$best_formula
result_temp_inc$best_variables     <- selection_temp_inc$selected_vars
result_temp_inc$cv_predictions     <- cv_preds_temp_inc
result_temp_inc$cv_metrics_by_species <- cv_by_species_inc
result_temp_inc$cv_metrics_by_year    <- cv_by_year_inc

cat("  Selected:", paste(selection_temp_inc$selected_vars, collapse = ", "), "\n")
cat("  RMSE: ", round(result_temp_inc$evaluation$uncorrected$rmse, 2), " -> ",
    round(result_temp_inc$evaluation$corrected$rmse, 2), "\n")

# --- 4.3 Model WITHOUT temporal variables (operational) ---
cat("\n--- Model WITHOUT temporal variables ---\n")
selection_temp_exc <- run_stepwise_selection(
  data         = model_data,
  candidate_vars = candidate_vars_temp_exc,
  smooth_vars  = smooth_vars,
  spatial_folds = spatial_folds,
  k = cv_folds, n = cv_repeats,
  bias_var     = "bias",
  improvement_threshold = improvement_threshold,
  max_vars     = length(candidate_vars_temp_exc)
)

result_temp_exc <- fit_and_evaluate_model(
  model_data, selection_temp_exc$best_formula,
  "CHM_height", "dominant_height", "bias"
)

# Collect CV holdout predictions
cv_preds_temp_exc <- collect_cv_predictions(
  model_data, selection_temp_exc$best_formula, spatial_folds,
  cv_folds, cv_repeats, "CHM_height", "dominant_height", "bias"
)
cv_by_species_exc <- compute_group_cv_metrics(cv_preds_temp_exc, "compo", "CHM_height", "dominant_height")
cv_by_year_exc    <- compute_group_cv_metrics(cv_preds_temp_exc, "year_factor", "CHM_height", "dominant_height")

result_temp_exc$selection          <- selection_temp_exc
result_temp_exc$best_formula       <- selection_temp_exc$best_formula
result_temp_exc$best_variables     <- selection_temp_exc$selected_vars
result_temp_exc$cv_predictions     <- cv_preds_temp_exc
result_temp_exc$cv_metrics_by_species <- cv_by_species_exc
result_temp_exc$cv_metrics_by_year    <- cv_by_year_exc

cat("  Selected:", paste(selection_temp_exc$selected_vars, collapse = ", "), "\n")
cat("  RMSE: ", round(result_temp_exc$evaluation$uncorrected$rmse, 2), " -> ",
    round(result_temp_exc$evaluation$corrected$rmse, 2), "\n")


# =============================================================================
# 5. SAVE OUTPUTS
# =============================================================================
cat("\n=== 5. SAVING ===\n")

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
saveRDS(result_temp_inc, file.path(path_output, "model_temp_inc.rds"))
saveRDS(result_temp_exc, file.path(path_output, "model_temp_exc.rds"))

cat("Saved:\n")
cat("  ", file.path(path_output, "model_temp_inc.rds"), "(with temporal vars)\n")
cat("  ", file.path(path_output, "model_temp_exc.rds"), "(without temporal vars — use for application)\n")
cat("\n=== COMPLETE ===\n")
