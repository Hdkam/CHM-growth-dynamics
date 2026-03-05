# =============================================================================
# CODE 3 — APPLICATION TO CASE STUDY ZONES
# =============================================================================
# PURPOSE:
#   Apply the standardization model to specific study zones (tiles) using
#   a hexagonal grid that mimics NFI plot size. For each hexagon:
#   1. Extract CHM metrics across all available years
#   2. Apply the standardization model to correct heights
#   3. Clean the resulting time series (spike/disturbance detection)
#
# INPUTS:
#   - model_temp_exc.rds    (standardization model without temporal vars)
#   - CHM rasters           (one per year)
#   - Ancillary rasters     (composition, broadleaf ratio, topography, slope, DEM, forest mask)
#   - Leaf-off masks        (optional)
#   - Flight tile metadata  (for acquisition dates)
#   - Study zone boundaries (GeoPackage with tile polygons)
#
# OUTPUTS:
#   - tile_{id}_timeseries.rds   (cleaned time series per tile)
#   - tile_{id}_hexgrid.rds      (hexagonal grid with stand delineation)
#
# NOTE:
#   This script processes tiles one at a time inside a loop.
#   For large areas, it can be parallelized or run tile-by-tile.
# =============================================================================


# =============================================================================
# 1. SETUP
# =============================================================================

# --- 1.1 Libraries ---
library(tidyverse)
library(terra)
library(tidyterra)
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

# Leaf-off masks
path_leafoff_folder  <- "data/leafoff_masks/"

# Ancillary rasters
path_composition     <- "data/ancillary/composition.tif"
path_broadleaf_ratio <- "data/ancillary/broadleaf_ratio.tif"
path_topography      <- "data/ancillary/topography.tif"
path_slope           <- "data/ancillary/slope.tif"
path_dem             <- "data/ancillary/dem.tif"
path_forest_mask     <- "data/ancillary/forest_mask.tif"

# Study zone
path_study_tiles     <- "data/study_tiles.gpkg"
tile_ids_to_process  <- c(1, 2, 3)   # IDs of the tiles to process

# --- 1.3 Parameters ---
crs_target    <- "EPSG:31370"
lidar_years   <- c(2014, 2021)

# Hexagon side length (chosen so that hex area ≈ NFI circular plot area = π × 18² ≈ 1018 m²)
hex_side <- sqrt(1018 / (3 * sqrt(3) / 2))

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
    code == 1  ~ "Other_broadleaf",
    code == 2  ~ "Birch",
    code == 3  ~ "Oak",
    code == 4  ~ "Douglas_fir",
    code == 5  ~ "Spruce",
    code == 6  ~ "Beech",
    code == 7  ~ "Larch",
    code == 8  ~ "Other_broadleaf",
    code == 9  ~ "Pine",
    code == 10 ~ "Other_coniferous",
    TRUE ~ NA_character_
  )
}

# --- 2.3 Modal value ---

modal_value <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# --- 2.4 Time series cleaning (self-contained) ---

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

clean_forest_timeseries <- function(input_data) {
  df_cor <- input_data %>%
    filter(status == "corrected") %>%
    select(-status)

  disturbance_fp <- df_cor %>%
    add_temporal_context() %>%
    detect_disturbances(recovery_threshold = 0.66)

  spike_info <- df_cor %>%
    detect_spikes(local_threshold = 50, broad_threshold = 30) %>%
    classify_spikes()

  df_reconciled <- reconcile_spike_disturbance(spike_info, disturbance_fp)

  disturbance_sp <- df_reconciled %>%
    filter(spike == FALSE) %>%
    add_temporal_context() %>%
    detect_disturbances(recovery_threshold = 0.66)

  df_reconciled %>%
    select(-disturbance_year, -has_been_disturbed, -height) %>%
    left_join(disturbance_sp, by = "ID") %>%
    right_join(input_data, by = c("ID", "year_factor"))
}


# =============================================================================
# 3. LOAD MODEL & STUDY ZONES
# =============================================================================
cat("=== 3. LOADING MODEL & STUDY ZONES ===\n")

model_result <- readRDS(path_model)
gam_model <- mgcv::gam(
  formula = as.formula(model_result$final_model$formula),
  data    = model_result$final_model$model,
  method  = "REML"
)
cat("  Model loaded.\n")

tiles_all <- vect(path_study_tiles) 
tiles_sel <- tiles_all[tiles_all$PXID %in% tile_ids_to_process] %>%
  terra::project(crs_target)
cat("  Processing", nrow(tiles_sel), "tile(s):", paste(tile_ids_to_process, collapse = ", "), "\n")

# Pre-load file lists
chm_files <- list.files(path_chm_folder, pattern = "\\d{4}\\.tif$", full.names = TRUE)
leafoff_files <- if (!is.null(path_leafoff_folder) && dir.exists(path_leafoff_folder)) {
  list.files(path_leafoff_folder, pattern = "leafoff_mask_\\d{4}\\.gpkg$", full.names = TRUE)
} else {
  character(0)
}
maillage <- vect(path_flight_metadata)


# =============================================================================
# 4. TILE-BY-TILE PROCESSING LOOP
# =============================================================================

for (tile_idx in seq_len(nrow(tiles_sel))) {

  tile <- tiles_sel[tile_idx, ]
  tile_id <- tile$PXID
  cat("\n==========================================================\n")
  cat("  PROCESSING TILE", tile_id, "\n")
  cat("==========================================================\n")

  # -------------------------------------------------------------------------
  # 4.1 Crop ancillary rasters to tile extent
  # -------------------------------------------------------------------------
  cat("  Cropping ancillary rasters...\n")

  chm_ref <- rast(chm_files[length(chm_files)]) %>% terra::crop(tile)
  crs(chm_ref) <- crs_target

  forest_mask <- rast(path_forest_mask) %>%
    terra::crop(tile) %>%
    terra::resample(chm_ref, method = "bilinear")

  compo_tile <- rast(path_composition) %>%
    terra::crop(tile) %>%
    terra::resample(chm_ref, method = "bilinear") %>%
    terra::mask(forest_mask)

  br_tile <- rast(path_broadleaf_ratio) %>%
    terra::crop(tile) %>%
    terra::resample(chm_ref, method = "bilinear") %>%
    terra::mask(forest_mask)

  topo_tile <- rast(path_topography) %>%
    terra::crop(tile) %>%
    terra::resample(chm_ref, method = "bilinear") %>%
    terra::mask(forest_mask)

  slope_tile <- rast(path_slope) %>%
    terra::crop(tile) %>%
    terra::resample(chm_ref, method = "bilinear") %>%
    terra::mask(forest_mask)

  dem_tile <- rast(path_dem) %>%
    terra::crop(tile) %>%
    terra::resample(chm_ref, method = "bilinear") %>%
    terra::mask(forest_mask)

  # -------------------------------------------------------------------------
  # 4.2 Build hexagonal grid
  # -------------------------------------------------------------------------
  cat("  Building hexagonal grid...\n")

  ext_poly <- st_as_sfc(st_bbox(ext(chm_ref), crs = terra::crs(chm_ref)))
  hex_grid <- st_make_grid(
    ext_poly,
    cellsize = hex_side * sqrt(3),
    square   = FALSE
  ) %>%
    st_sf() %>%
    mutate(ID = row_number()) %>%
    vect()

  cat("    ", nrow(hex_grid), "hexagons created\n")

  # -------------------------------------------------------------------------
  # 4.3 Extract ancillary values per hexagon
  # -------------------------------------------------------------------------
  cat("  Extracting ancillary values per hexagon...\n")

  compo_hex  <- terra::extract(compo_tile %>% as.factor(), hex_grid, fun = modal_value, bind = TRUE)
  br_hex     <- terra::extract(br_tile, hex_grid, fun = mean, bind = TRUE, na.rm = TRUE)
  topo_hex   <- terra::extract(topo_tile %>% as.factor(), hex_grid, fun = modal_value, bind = TRUE)
  slope_hex  <- terra::extract(slope_tile, hex_grid, fun = mean, bind = TRUE, na.rm = TRUE)
  alt_hex    <- terra::extract(dem_tile, hex_grid, fun = mean, bind = TRUE, na.rm = TRUE)

  # Assemble into a single dataframe
  hex_df <- left_join(compo_hex %>% as.data.frame(), br_hex %>% as.data.frame()) %>%
    left_join(topo_hex %>% as.data.frame()) %>%
    left_join(slope_hex %>% as.data.frame()) %>%
    left_join(alt_hex %>% as.data.frame())

  # Rename columns (position-based — adjust if your raster layer names differ)
  col_names <- names(hex_df)
  hex_df <- hex_df %>%
    rename(
      compo_raw       = 2,
      broadleaf_ratio = 3,
      topo            = 4,
      slope           = 5,
      alt             = 6
    ) %>%
    mutate(
      compo = map_composition_to_species(round(compo_raw), round(broadleaf_ratio))
    )

  values(hex_grid) <- hex_df

  # -------------------------------------------------------------------------
  # 4.4 Stack and crop CHM rasters
  # -------------------------------------------------------------------------
  cat("  Stacking CHM rasters...\n")

  base_rast <- rast(ext(chm_ref), resolution = 0.5, crs = crs_target)
  raster_list <- list()

  for (i in seq_along(chm_files)) {
    year <- as.character(
      as.integer(regmatches(basename(chm_files[i]), regexpr("\\d{4}", basename(chm_files[i]))))
    )
    cat("    Year", year, "...\n")

    chm_yr <- rast(chm_files[i])
    crs(chm_yr) <- crs_target
    chm_yr <- chm_yr %>%
      terra::crop(tile) %>%
      terra::resample(base_rast) %>%
      terra::mask(forest_mask)

    raster_list[[i]] <- chm_yr
  }

  chm_stacked <- do.call(c, raster_list)
  gc()

  # -------------------------------------------------------------------------
  # 4.5 Extract CHM metrics per hexagon across all years
  # -------------------------------------------------------------------------
  cat("  Extracting CHM metrics per hexagon...\n")
  
  metrics <- terra::extract(chm_stacked, hex_grid, fun = chm_summary_stats, bind = TRUE)
  
  # Strip the raster layer name prefix (e.g. "mnh2006" → "2006", "mnh2006.1" → "2006.1")
  # The prefix comes from the raster filenames; we only need the year + metric suffix.
  metric_cols <- grep("^[a-zA-Z]+\\d{4}", names(metrics), value = TRUE)
  new_names <- sub("^[a-zA-Z]+(?=\\d{4})", "", metric_cols, perl = TRUE)
  names(metrics)[match(metric_cols, names(metrics))] <- new_names

  # --- Flag leaf-off hexagons ---
  for (i in seq_along(chm_files)) {
    year <- as.character(
      as.integer(regmatches(basename(chm_files[i]), regexpr("\\d{4}", basename(chm_files[i]))))
    )
    mask_file <- leafoff_files[grepl(year, leafoff_files)]

    if (length(mask_file) == 1) {
      leafoff_poly <- try(vect(mask_file) %>% crop(ext(hex_grid)), silent = TRUE)
      if (inherits(leafoff_poly, "try-error") || nrow(leafoff_poly) == 0) next

      hex_deciduous <- hex_grid[hex_grid$compo %in% c("Beech", "Oak", "Birch", "Other_deciduous"), ]
      if (nrow(hex_deciduous) == 0) next

      hex_intersect <- is.related(hex_deciduous, leafoff_poly, "intersects")
      leafoff_ids <- hex_deciduous[hex_intersect, ]$ID

      metrics[[paste0(year, ".leafoff")]] <- metrics$ID %in% leafoff_ids
      cat("    ", year, ":", sum(hex_intersect), "hexagons flagged as leaf-off\n")
    }
  }

  gc()

  # -------------------------------------------------------------------------
  # 4.6 Reshape metrics into a long time-series format
  # -------------------------------------------------------------------------
  cat("  Reshaping to time series format...\n")

  # Intersect with flight metadata (to get acquisition dates)
  maillage_tile <- maillage %>%
    terra::crop(ext(hex_grid)) %>%
    mutate(date_flight = as.Date(centredate))

  tile_ts_raw <- terra::intersect(maillage_tile, metrics) %>%
    as.data.frame()
  
  
  names(tile_ts_raw) <- sub("^[a-zA-Z]+(?=\\d{4})", "", names(tile_ts_raw), perl = TRUE)

  # Separate leaf-off columns
  leafoff_cols <- grep("\\.leafoff$", names(tile_ts_raw), value = TRUE)
  

  if (length(leafoff_cols) > 0) {
    leafoff_df <- tile_ts_raw %>%
      select(ID, all_of(leafoff_cols)) %>%
      distinct() %>%
      pivot_longer(
        cols = all_of(leafoff_cols),
        names_to    = "year",
        names_pattern = "(\\d{4})\\.leafoff$",
        values_to   = "leaf_off"
      ) %>%
      mutate(leaf_off = replace_na(leaf_off, FALSE))
  } else {
    leafoff_df <- NULL
  }

  # Pivot CHM metrics from wide to long
  tile_ts_df <- tile_ts_raw %>%
    select(-any_of(leafoff_cols)) %>%
    rename(year_flight = year) %>%
    pivot_longer(
      cols = matches("^\\d{4}"),
      names_to  = c("year", "stat"),
      names_pattern = "(\\d{4})(.*)"
    ) %>%
    mutate(
      stat = case_when(
        stat == ""   ~ "CHM_height",
        stat == ".1" ~ "CHM_hsd",
        stat == ".2" ~ "CHM_gap_fraction",
        stat == ".3" ~ "CHM_hmean",
        TRUE ~ stat
      )
    ) %>%
    select(-centredate, -date_flight) %>%
    distinct() %>%
    pivot_wider(names_from = stat, values_from = value, values_fn = ~ mean(.x, na.rm = TRUE)) %>%
    filter(year == year_flight) %>%
    mutate(
      year_factor = as.factor(year),
      sensor      = as.factor(ifelse(year %in% as.character(lidar_years), "lidar", "dap")),
      compo       = as.factor(compo),
      topo        = as.factor(round(topo))
    )

  # Join leaf-off flags
  if (!is.null(leafoff_df)) {
    tile_ts_df <- tile_ts_df %>% 
      left_join(leafoff_df, by = c("ID", "year")) %>%
      # Check that leaf-off is flagged accordingly
      mutate(leaf_off = coalesce(leaf_off, FALSE))
  } else {
    tile_ts_df$leaf_off <- FALSE
  }

  # Remove leaf-off observations and incomplete rows
  tile_ts_df <- tile_ts_df %>%
    filter(!leaf_off, !is.na(compo), !is.na(topo))

  cat("    ", nrow(tile_ts_df), "observations ready for model application\n")

  # -------------------------------------------------------------------------
  # 4.7 Apply standardization model
  # -------------------------------------------------------------------------
  cat("  Applying standardization model...\n")

  tile_ts_pred <- tile_ts_df %>%
    rowwise() %>%
    mutate(
      predicted = list({
        pred <- predict(gam_model, newdata = pick(everything()), se.fit = TRUE)
        list(
          bias_fit    = pred$fit,
          bias_se     = pred$se.fit,
          lower_bound = pred$fit - 1.96 * pred$se.fit,
          upper_bound = pred$fit + 1.96 * pred$se.fit
        )
      })
    ) %>%
    ungroup() %>%
    unnest_wider(predicted) %>%
    mutate(
      height_raw       = CHM_height,
      height_corrected = CHM_height + bias_fit
    )

  # Reshape to long format
  tile_ts_long <- tile_ts_pred %>%
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

  # -------------------------------------------------------------------------
  # 4.8 Clean time series
  # -------------------------------------------------------------------------
  cat("  Cleaning time series...\n")

  tile_ts_clean <- clean_forest_timeseries(tile_ts_long)

  cat("    Done.\n")

  # -------------------------------------------------------------------------
  # 4.9 Save tile outputs
  # -------------------------------------------------------------------------
  dir.create(path_output, recursive = TRUE, showWarnings = FALSE)

  saveRDS(tile_ts_clean, file.path(path_output, paste0("tile_", tile_id, "_timeseries.rds")))
  saveRDS(hex_grid,      file.path(path_output, paste0("tile_", tile_id, "_hexgrid.rds")))

  cat("  Saved: tile_", tile_id, "_timeseries.rds + tile_", tile_id, "_hexgrid.rds\n")

  # Clean up memory
  rm(chm_stacked, metrics, tile_ts_raw, tile_ts_df, tile_ts_pred, tile_ts_long, tile_ts_clean)
  gc()
}

cat("\n=== ALL TILES COMPLETE ===\n")
