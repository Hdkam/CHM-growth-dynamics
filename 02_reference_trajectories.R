set.seed(42)
# =============================================================================
# CODE 2 — REFERENCE TRAJECTORIES ON RANDOM SUBSAMPLE
# =============================================================================
# PURPOSE:
#   Generate species-specific reference height trajectories by:
#   1. Creating virtual circular plots randomly sampled across the study area
#   2. Extracting CHM metrics for each plot across all available years
#   3. Applying the standardization model (from Code 1) to correct heights
#   4. Cleaning the resulting time series (spike/disturbance detection)
#
# INPUTS:
#   - model_temp_exc.rds    (standardization model without temporal vars)
#   - CHM rasters           (one per year)
#   - Ancillary rasters     (composition, broadleaf ratio, topography, slope, DEM)
#   - Leaf-off masks        (optional)
#   - Flight tile metadata  (for acquisition dates)
#
# OUTPUTS:
#   - reference_trajectories.rds  (cleaned time series for all virtual plots)
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

# --- 1.2 File paths ---
# >> Adjust these to your local folder structure <<

# Model from Code 1
path_model <- "outputs/data/model_temp_exc.rds"

# CHM data
path_chm_folder      <- "data/chm/"
path_flight_metadata <- "data/flight_tiles.gpkg"      

# Leaf-off masks (set to NULL if not available)
path_leafoff_folder <- "data/leafoff_masks/"        

# Ancillary rasters (set to NULL if unavailable)
path_composition     <- "data/ancillary/composition.tif"
path_broadleaf_ratio <- "data/ancillary/broadleaf_ratio.tif"
path_topography      <- "data/ancillary/topography.tif"
path_slope           <- "data/ancillary/slope.tif"
path_dem             <- "data/ancillary/dem.tif"

# Output
path_output <- "outputs/data/"

# --- 1.3 Parameters ---
crs_target         <- "EPSG:31370"
field_plot_radius  <- 18          # same as NFI plot radius (m)
lidar_years        <- c(2014, 2021)

# Subsampling
n_oversample       <- 100000      # initial oversampling (to ensure enough valid points)
n_per_species      <- 1500        # final number of plots per species group
max_na_fraction    <- 0.20        # reject plots with >20% NA pixels in the buffer

# Memory management
terraOptions(memfrac = 0.5)
terraOptions(tempdir = tempdir())


# =============================================================================
# 2. FUNCTIONS
# =============================================================================

# --- 2.1 CHM metrics ---

calculate_canopy_openness <- function(vals, thr = 0.5) {
  h99 <- quantile(vals, 0.99, na.rm = TRUE)
  if (is.na(h99) || h99 == 0) return(NA_real_)
  threshold <- thr * h99
  sum(vals < threshold, na.rm = TRUE) / sum(!is.na(vals)) * 100
}

chm_summary_stats <- function(vals) {
  c(
    CHM_height       = as.numeric(quantile(vals, 0.95, na.rm = TRUE)),
    CHM_hsd          = sd(vals, na.rm = TRUE),
    CHM_gap_fraction = calculate_canopy_openness(vals),
    CHM_hmean        = mean(vals, na.rm = TRUE)
  )
}

# --- 2.2 Composition mapping ---

map_composition_to_species <- function(code, broadleaf_ratio) {
  case_when(
    code %in% c(0, 1) & broadleaf_ratio >= 51 ~ "Other_deciduous",
    code %in% c(0, 1) & broadleaf_ratio <  51 ~ "Other_evergreen",
    code == 2  ~ "Birch",
    code == 3  ~ "Oak",
    code == 4  ~ "Douglas_fir",
    code == 5  ~ "Spruce",
    code == 6  ~ "Beech",
    code == 7  ~ "Larch",
    code == 8  ~ "Other_deciduous",
    code == 9  ~ "Pine",
    code == 10 ~ "Other_evergreen",
    TRUE ~ NA_character_
  )
}

# --- 2.3 Modal value for categorical rasters ---

modal_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# --- 2.4 Time series cleaning functions ---

#' Add temporal context (lags/leads) to a height time series
add_temporal_context <- function(df) {
  df %>%
    arrange(ID, as.numeric(as.character(year_factor))) %>%
    group_by(ID) %>%
    mutate(
      year_num    = as.numeric(as.character(year_factor)),
      prev_height = lag(height, 1),
      next_height = lead(height, 1)
    ) %>%
    ungroup()
}

#' Detect disturbances (large drops not followed by recovery)
detect_disturbances <- function(df, drop_threshold = 0.5,
                                recovery_threshold = 0.66,
                                min_absolute_drop = 10) {
  df %>%
    group_by(ID) %>%
    mutate(
      major_drop = (prev_height - height) > pmax(height * drop_threshold, min_absolute_drop),
      recovered_1yr = (lead(height, 1) > prev_height * recovery_threshold),
      recovered_2yr = (lead(height, 2) > prev_height * recovery_threshold),
      disturbance_trigger = major_drop & (
        (is.na(recovered_1yr) & is.na(recovered_2yr)) |
          (!recovered_1yr & (is.na(recovered_2yr) | !recovered_2yr))
      ),
      disturbance_year = if_else(disturbance_trigger, year_num, NA_real_)
    ) %>%
    summarise(
      has_been_disturbed = any(disturbance_trigger, na.rm = TRUE),
      disturbance_year   = min(disturbance_year, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(disturbance_year = if_else(is.infinite(disturbance_year), NA_real_, disturbance_year))
}

#' Detect spikes (positive and negative outliers)
detect_spikes <- function(df, local_threshold = 50, broad_threshold = 30) {
  df %>%
    arrange(ID, as.numeric(as.character(year_factor))) %>%
    group_by(ID) %>%
    mutate(
      prev_height  = lag(height),
      next_height  = lead(height),
      prev2_height = lag(height, 2),
      next2_height = lead(height, 2)
    ) %>%
    rowwise() %>%
    mutate(
      local_median = median(c(prev_height, height, next_height), na.rm = TRUE),
      broad_mean   = mean(c(prev2_height, prev_height, height, next_height, next2_height), na.rm = TRUE),
      relative_diff_from_neighbors         = 100 * (height - local_median) / local_median,
      relative_diff_from_further_neighbors = 100 * (height - broad_mean) / broad_mean,
      negative_spike = !is.na(prev_height) & !is.na(next_height) &
        relative_diff_from_neighbors <= -local_threshold &
        relative_diff_from_further_neighbors <= -broad_threshold &
        abs(next_height - height) > 2,
      positive_spike = !is.na(prev_height) & !is.na(next_height) &
        relative_diff_from_neighbors >= local_threshold &
        relative_diff_from_further_neighbors >= broad_threshold &
        abs(next_height - height) > 2
    ) %>%
    ungroup() %>%
    mutate(
      negative_spike = replace_na(negative_spike, FALSE),
      positive_spike = replace_na(positive_spike, FALSE)
    )
}

#' Classify detected spikes into true vs. false positives
classify_spikes <- function(df) {
  df %>%
    group_by(ID) %>%
    mutate(
      spike_type = case_when(
        negative_spike ~ "negative",
        positive_spike ~ "positive",
        TRUE ~ "none"
      ),
      true_spike = case_when(
        spike_type == "none" ~ "none",
        spike_type == "negative" &
          (lag(spike_type, 2) == spike_type | lead(spike_type, 2) == spike_type) ~ "none",
        spike_type != "none" &
          !is.na(lag(spike_type)) & !is.na(lead(spike_type)) &
          lag(spike_type) == lead(spike_type) &
          lag(spike_type) != spike_type ~ spike_type,
        spike_type == "positive" & !is.na(lag(spike_type)) & lag(spike_type) == "negative" &
          lead(height) >= 0.7 * height ~ "none",
        spike_type == "positive" &
          !is.na(lead(spike_type)) & lead(spike_type) == "negative" &
          100 * (height - mean(c(lag(height), lead(height, 2), lead(height, 3)), na.rm = TRUE)) /
          mean(c(lag(height), lead(height, 2), lead(height, 3)), na.rm = TRUE) >= 30 &
          abs(lead(height, 2) - height) > 2 ~ spike_type,
        TRUE ~ spike_type
      )
    ) %>%
    ungroup() %>%
    mutate(
      is_negative_spike = true_spike == "negative",
      is_positive_spike = true_spike == "positive",
      spike = is_negative_spike | is_positive_spike
    ) %>%
    select(ID, year_factor, height, spike, true_spike, is_negative_spike, is_positive_spike)
}

#' Reconcile spikes and disturbances (a negative spike at the disturbance year is not a spike)
reconcile_spike_disturbance <- function(spike_df, disturbance_df) {
  spike_df %>%
    left_join(disturbance_df, by = "ID", suffix = c("", "_fp")) %>%
    group_by(ID) %>%
    mutate(
      spike = case_when(
        is_negative_spike & year_factor == disturbance_year ~ FALSE,
        TRUE ~ spike
      ),
      invalid_disturbance = any(
        !is.na(disturbance_year) & year_factor == disturbance_year &
          !is.na(lag(is_positive_spike)) & lag(is_positive_spike) == TRUE,
        na.rm = TRUE
      ),
      disturbance_year = ifelse(invalid_disturbance, NA_real_, disturbance_year)
    ) %>%
    select(-invalid_disturbance) %>%
    ungroup()
}

#' Main time series cleaning pipeline
clean_forest_timeseries <- function(input_data) {
  # Work on corrected heights only
  df_cor <- input_data %>%
    filter(status == "corrected") %>%
    select(-status)

  # Pass 1: initial disturbance detection
  disturbance_fp <- df_cor %>%
    add_temporal_context() %>%
    detect_disturbances(recovery_threshold = 0.66)

  # Pass 2: spike detection
  spike_info <- df_cor %>%
    detect_spikes(local_threshold = 50, broad_threshold = 30) %>%
    classify_spikes()

  # Reconcile
  df_reconciled <- reconcile_spike_disturbance(spike_info, disturbance_fp)

  # Pass 3: re-detect disturbances after removing spikes
  disturbance_sp <- df_reconciled %>%
    filter(spike == FALSE) %>%
    add_temporal_context() %>%
    detect_disturbances(recovery_threshold = 0.66)

  # Final merge
  df_reconciled %>%
    select(-disturbance_year, -has_been_disturbed, -height) %>%
    left_join(disturbance_sp, by = "ID") %>%
    right_join(input_data, by = c("ID", "year_factor"))
}


# =============================================================================
# 3. LOAD MODEL
# =============================================================================
cat("=== 3. LOADING MODEL ===\n")

model_result <- readRDS(path_model)

# Rebuild the GAM from the stored formula and training data
gam_model <- mgcv::gam(
  formula = as.formula(model_result$final_model$formula),
  data    = model_result$final_model$model,
  method  = "REML"
)
cat("  Model formula:", deparse(formula(gam_model)), "\n")


# =============================================================================
# 4. CREATE VIRTUAL PLOTS
# =============================================================================
cat("\n=== 4. CREATING VIRTUAL PLOTS ===\n")

# --- 4.1 Random sampling within the composition raster extent ---
cat("Sampling", n_oversample, "candidate points...\n")
compo_rast <- rast(path_composition)
crs(compo_rast) <- crs_target

points_raw <- terra::spatSample(
  compo_rast, size = n_oversample,
  method = "random", na.rm = TRUE, as.points = TRUE
) %>%
  rename(compo_raw = 1)

# --- 4.2 Reject points too close to raster edges (>20% NA in buffer) ---
cat("Filtering edge plots...\n")
na_prop <- terra::extract(compo_rast, buffer(points_raw, field_plot_radius),
                          fun = function(x) mean(is.na(x)))
points_valid <- points_raw[na_prop[[2]] < max_na_fraction, ]
cat("  ", nrow(points_valid), "points retained after edge filtering\n")

# --- 4.3 Extract modal composition within the buffered plot ---
cat("Extracting plot-level composition...\n")
rast_compo_factor <- compo_rast %>% terra::as.factor()
plots_with_compo <- terra::extract(
  x = rast_compo_factor,
  y = buffer(points_valid, field_plot_radius),
  fun = modal_value,
  bind = TRUE
) %>%
  rename(compo_point = 1, compo_modal = 2)

# Map numeric codes to species names
plots_with_compo$compo <- case_when(
  plots_with_compo$compo_modal == 1  ~ "Other_broadleaf",
  plots_with_compo$compo_modal == 2  ~ "Birch",
  plots_with_compo$compo_modal == 3  ~ "Oak",
  plots_with_compo$compo_modal == 4  ~ "Douglas_fir",
  plots_with_compo$compo_modal == 5  ~ "Spruce",
  plots_with_compo$compo_modal == 6  ~ "Beech",
  plots_with_compo$compo_modal == 7  ~ "Larch",
  plots_with_compo$compo_modal == 8  ~ "Other_broadleaf",
  plots_with_compo$compo_modal == 9  ~ "Pine",
  plots_with_compo$compo_modal == 10 ~ "Other_coniferous",
  TRUE ~ NA_character_
)

# --- 4.4 Stratified subsample: n_per_species plots per species group ---
cat("Stratified subsampling:", n_per_species, "per species...\n")
plots_subsample_vect <- plots_with_compo %>%
  filter(!is.na(compo)) %>%
  group_by(compo) %>%
  slice_sample(n = n_per_species) %>%
  ungroup() %>%
  mutate(plot_ID = row_number())

cat("  ", nrow(plots_subsample_vect), "virtual plots created across",
    n_distinct(plots_subsample_vect$compo), "species groups\n")

# --- 4.5 Extract ancillary data for virtual plots ---
cat("Extracting ancillary data for virtual plots...\n")

# Elevation
rast_alt <- rast(path_dem)
crs(rast_alt) <- crs_target
alt_df <- terra::extract(rast_alt, plots_subsample_vect, fun = mean, bind = TRUE, na.rm = TRUE) %>%
  as.data.frame()
alt_col <- names(alt_df)[ncol(alt_df)]
alt_df <- alt_df %>% rename(alt = all_of(alt_col)) %>% select(plot_ID, alt)

# Broadleaf ratio
rast_br <- rast(path_broadleaf_ratio)
crs(rast_br) <- crs_target
br_df <- terra::extract(rast_br, plots_subsample_vect, fun = mean, bind = TRUE, na.rm = TRUE) %>%
  as.data.frame()
br_col <- names(br_df)[ncol(br_df)]
br_df <- br_df %>% rename(broadleaf_ratio = all_of(br_col)) %>% select(plot_ID, broadleaf_ratio)

# Slope
rast_slope <- rast(path_slope)
crs(rast_slope) <- crs_target
slope_df <- terra::extract(rast_slope, plots_subsample_vect, fun = mean, bind = TRUE, na.rm = TRUE) %>%
  as.data.frame()
slope_col <- names(slope_df)[ncol(slope_df)]
slope_df <- slope_df %>% rename(slope = all_of(slope_col)) %>% select(plot_ID, slope)

# Topographic position
rast_topo <- rast(path_topography)
crs(rast_topo) <- crs_target
topo_df <- terra::extract(rast_topo, plots_subsample_vect, fun = modal_value, bind = TRUE) %>%
  as.data.frame()
topo_col <- names(topo_df)[ncol(topo_df)]
topo_df <- topo_df %>%
  rename(topo = all_of(topo_col)) %>%
  mutate(topo = round(topo)) %>%
  select(plot_ID, topo)

cat("  Done.\n")


# =============================================================================
# 5. EXTRACT CHM TIME SERIES
# =============================================================================
cat("\n=== 5. EXTRACTING CHM TIME SERIES ===\n")

chm_files <- list.files(path_chm_folder, pattern = "\\d{4}\\.tif$", full.names = TRUE)
leafoff_files <- if (!is.null(path_leafoff_folder) && dir.exists(path_leafoff_folder)) {
  list.files(path_leafoff_folder, pattern = "leafoff_mask_\\d{4}\\.gpkg$", full.names = TRUE)
} else {
  character(0)
}

plots_height_df <- data.frame()

for (i in seq_along(chm_files)) {
  year <- as.character(
    as.integer(regmatches(basename(chm_files[i]), regexpr("\\d{4}", basename(chm_files[i]))))
  )
  cat("  Processing year", year, "...\n")

  chm <- rast(chm_files[i])
  crs(chm) <- crs_target

  # Start with all virtual plots
  plots_to_extract <- plots_subsample_vect
  plots_to_extract$leaf_off <- FALSE

  # Flag leaf-off plots for this year
  mask_file <- leafoff_files[grepl(year, leafoff_files)]
  if (length(mask_file) == 1) {
    leafoff_poly <- vect(mask_file)

    # Only beech and oak can be leaf-off in this dataset
    deciduous_ids <- plots_subsample_vect$plot_ID[plots_subsample_vect$compo %in% c("Beech", "Oak")]
    plots_deciduous <- plots_subsample_vect[plots_subsample_vect$plot_ID %in% deciduous_ids]
    plots_cropped <- try(crop(plots_deciduous, ext(leafoff_poly)), silent = TRUE)

    if (!inherits(plots_cropped, "try-error") && nrow(plots_cropped) > 0) {
      intersects <- is.related(plots_cropped, leafoff_poly, relation = "intersects")
      leafoff_ids <- plots_cropped[intersects, ]$plot_ID
      plots_to_extract$leaf_off[plots_to_extract$plot_ID %in% leafoff_ids] <- TRUE
      cat("    Flagged", sum(intersects), "leaf-off plots\n")
    }
  }

  # Extract CHM metrics (without bind — rename columns explicitly)
  metrics_raw <- terra::extract(chm, plots_to_extract, fun = chm_summary_stats) %>%
    as.data.frame()
  metric_names <- c("CHM_height", "CHM_hsd", "CHM_gap_fraction", "CHM_hmean")
  names(metrics_raw) <- c("ID_row", metric_names)

  # Map row index back to plot attributes
  metrics <- metrics_raw %>%
    mutate(
      plot_ID    = plots_to_extract$plot_ID[ID_row],
      compo      = plots_to_extract$compo[ID_row],
      year_factor = as.factor(year),
      sensor     = as.factor(ifelse(year %in% as.character(lidar_years), "lidar", "dap")),
      leaf_off   = plots_to_extract$leaf_off[ID_row]
    ) %>%
    select(-ID_row)

  plots_height_df <- bind_rows(plots_height_df, metrics)
}

# --- 5.1 Combine with ancillary data ---
cat("Combining with ancillary data...\n")
plots_ts_raw <- plots_height_df %>%
  left_join(alt_df, by = "plot_ID") %>%
  left_join(br_df, by = "plot_ID") %>%
  left_join(slope_df, by = "plot_ID") %>%
  left_join(topo_df, by = "plot_ID") %>%
  mutate(
    compo = as.factor(compo),
    topo  = as.factor(topo),
    # Ensure leaf_off is clean boolean (kept as column for downstream inspection)
    leaf_off = !is.na(leaf_off) & leaf_off
  )

cat("  ", nrow(plots_ts_raw), "observations total,",
    sum(plots_ts_raw$leaf_off), "flagged as leaf-off\n")


# =============================================================================
# 6. APPLY STANDARDIZATION MODEL
# =============================================================================
cat("\n=== 6. APPLYING STANDARDIZATION MODEL ===\n")

plots_ts_pred <- plots_ts_raw %>%
  rowwise() %>%
  mutate(
    predicted = list({
      pred <- predict(gam_model, newdata = pick(everything()), se.fit = TRUE)
      list(
        bias_fit     = pred$fit,
        bias_se      = pred$se.fit,
        lower_bound  = pred$fit - 1.96 * pred$se.fit,
        upper_bound  = pred$fit + 1.96 * pred$se.fit
      )
    })
  ) %>%
  ungroup() %>%
  unnest_wider(predicted) %>%
  mutate(
    height_raw       = CHM_height,
    height_corrected = CHM_height + bias_fit
  )

cat("  Model applied to", nrow(plots_ts_pred), "observations\n")

# --- 6.1 Reshape to long format for cleaning pipeline ---
plots_ts_long <- plots_ts_pred %>%
  pivot_longer(
    cols      = c(height_raw, height_corrected),
    names_to  = "status",
    names_prefix = "height_",
    values_to = "height"
  ) %>%
  mutate(
    lower_bound = ifelse(status == "raw", NA_real_, lower_bound),
    upper_bound = ifelse(status == "raw", NA_real_, upper_bound)
  ) %>%
  distinct()


# =============================================================================
# 7. TIME SERIES CLEANING
# =============================================================================
cat("\n=== 7. CLEANING TIME SERIES ===\n")

plots_ts_clean <- clean_forest_timeseries(
  plots_ts_long %>% rename(ID = plot_ID)
) %>%
  rename(plot_ID = ID)

cat("  Cleaning complete.\n")
cat("  Disturbances detected in",
    sum(plots_ts_clean$has_been_disturbed, na.rm = TRUE) /
      n_distinct(plots_ts_clean$year_factor),
    "plots\n")
cat("  Spikes flagged:",
    sum(plots_ts_clean$spike, na.rm = TRUE), "observations\n")

# =============================================================================
# 8. SAVE OUTPUT
# =============================================================================
cat("\n=== 8. SAVING ===\n")

dir.create(path_output, recursive = TRUE, showWarnings = FALSE)
saveRDS(plots_ts_clean, file.path(path_output, "reference_trajectories.rds"))

cat("Saved:", file.path(path_output, "reference_trajectories.rds"), "\n")
cat("\n=== COMPLETE ===\n")

