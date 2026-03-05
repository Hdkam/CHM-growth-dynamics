# =============================================================================
# CODE 4 — RESULTS AND FIGURES
# =============================================================================
# PURPOSE:
#   Produce all tables and figures for the article, organized by manuscript
#   numbering. Each section can be run independently once data are loaded (§2).
#
#   Figures 1-2 and Tables 1, 3 are not code-generated (workflow diagram,
#   study area map, survey table, candidate predictors table).
#
#   This script produces:
#     - Table 2   : Plot distribution by year × species
#     - Figure 3  : Biplots before/after harmonization (both models)
#     - Table S1  : Per-year harmonization metrics
#     - Table S2  : Per-species harmonization metrics
#     - Figure 4  : Reference growth trajectories per species
#     - Figure 5  : Vertical growth and deviation maps (Site A)
#     - Figure 6  : Beech deviation across sites
#     - Table 4   : Moran's I spatial autocorrelation
#
# INPUTS:
#   - model_temp_inc.rds         (model with temporal vars)
#   - model_temp_exc.rds         (model without temporal vars)
#   - reference_trajectories.rds (cleaned subsample time series)
#   - tile_{id}_timeseries.rds   (cleaned tile time series)
#   - tile_{id}_hexgrid.rds      (hexagonal grids)
# =============================================================================


# =============================================================================
# 1. SETUP
# =============================================================================

# --- 1.1 Libraries ---
library(tidyverse)
library(terra)
library(tidyterra)
library(sf)
library(patchwork)
library(ggridges)
library(ggpattern)
library(spdep)

# --- 1.2 File paths ---
path_data <- "outputs/data/"


# --- 1.3 Color palettes ---
species_colors <- c(
  "Other_coniferous"  = "#4c2772",
  "Other_broadleaf"  = "#a6cee3",
  "Birch"            = "#1f78b4",
  "Oak"              = "#b2df8a",
  "Douglas_fir"      = "#33a02c",
  "Spruce"           = "#fb9a99",
  "Beech"            = "#e31a1c",
  "Larch"            = "#fdbf6f",
  "Pine"             = "#cab2d6"
)

status_colors <- c(
  "raw"       = "hotpink2",
  "corrected" = "darkgreen"
)

status_labels <- c(
  "raw"       = "Original",
  "corrected" = "Harmonized\nw/out temporal var."
)

# --- 1.4 Shared theme ---
theme_article <- theme_bw() +
  theme(
    text = element_text(size = 10),
    legend.position = "right"
  )

# --- 1.5 Tile configuration ---
# Map tile IDs to manuscript site labels (adjust to match your data)
tile_config <- tibble(
  tile_id    = c(1, 2, 3),   # <-- adjust these IDs
  site_label = c("Site A", "Site B", "Site C")
)


# =============================================================================
# 2. LOAD DATA
# =============================================================================
cat("=== 2. LOADING DATA ===\n")

# Models
mod_temp_inc <- readRDS(file.path(path_data, "model_temp_inc.rds"))
mod_temp_exc <- readRDS(file.path(path_data, "model_temp_exc.rds"))

# Reference trajectories
ref_traj <- readRDS(file.path(path_data, "reference_trajectories.rds"))

# Tile outputs
tile_ts   <- list()
tile_hex  <- list()
for (i in seq_len(nrow(tile_config))) {
  tid <- as.character(tile_config$tile_id[i])
  lbl <- tile_config$site_label[i]
  ts_file  <- file.path(path_data, paste0("tile_", tid, "_timeseries.rds"))
  hex_file <- file.path(path_data, paste0("tile_", tid, "_hexgrid.rds"))
  if (file.exists(ts_file) & file.exists(hex_file)) {
    tile_ts[[lbl]]  <- readRDS(ts_file)
    tile_hex[[lbl]] <- readRDS(hex_file)
    cat("  Loaded", lbl, "(tile", tid, ")\n")
  } else {
    cat("  SKIPPED", lbl, "(tile", tid, ") — files not found\n")
  }
}

cat("  All available data loaded.\n")


# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

#' Compute CV scores from a model object
compute_cv_scores <- function(model_obj) {
  cv_preds <- model_obj$cv_predictions %>%
    group_by(plot_ID) %>%
    summarise(
      dominant_height  = first(dominant_height),
      CHM_height       = first(CHM_height),
      corrected_height = mean(corrected_height, na.rm = TRUE),
      .groups = "drop"
    )

  cv_preds %>%
    summarise(
      n_plots          = n(),
      mean_bias_before = mean(dominant_height - CHM_height, na.rm = TRUE),
      rmse_before      = sqrt(mean((dominant_height - CHM_height)^2, na.rm = TRUE)),
      r2_before        = 1 - sum((dominant_height - CHM_height)^2) /
                              sum((dominant_height - mean(dominant_height))^2),
      mean_bias_after  = mean(dominant_height - corrected_height, na.rm = TRUE),
      rmse_after       = sqrt(mean((dominant_height - corrected_height)^2, na.rm = TRUE)),
      r2_after         = 1 - sum((dominant_height - corrected_height)^2) /
                              sum((dominant_height - mean(dominant_height))^2)
    )
}

#' Compute growth slope and initial height for each unit × status
compute_growth_metrics <- function(ts_data, id_col = "ID", ref_year = 2006) {
  ts_data %>%
    mutate(
      year_flight = as.numeric(as.character(year_factor)),
      leaf_off    = !is.na(leaf_off) & leaf_off
    ) %>%
    group_by(.data[[id_col]], status) %>%
    arrange(year_flight) %>%
    # Extract h0 BEFORE filtering (ref_year may be flagged as spike or leaf_off)
    mutate(h0 = first(height[year_flight == ref_year], default = NA_real_)) %>%
    # Now filter out spikes and leaf_off for growth slope computation
    filter(
      spike == FALSE | is.na(spike),
      !leaf_off,
      !is.na(h0),
      any(!is.na(height)),
      n() >= 3   # need at least 3 points for a meaningful regression
    ) %>%
    mutate(
      hslop = lm(height ~ year_flight)$coefficients[2]
    ) %>%
    ungroup() %>%
    distinct(.data[[id_col]], status, compo, h0, hslop, has_been_disturbed) %>%
    filter(!is.na(h0), h0 > 2, !is.na(hslop))
}

#' Prepare a tile for mapping: compute growth metrics + deviations
prepare_tile_for_mapping <- function(ts_data, hex_grid, ref_coefs, id_col = "ID") {
  growth <- compute_growth_metrics(ts_data, id_col = id_col)

  growth_dev <- growth %>%
    left_join(ref_coefs, by = c("compo", "status")) %>%
    mutate(
      predicted      = intercept_lm + slope_lm * h0,
      vertical_resid = hslop - predicted
    ) %>%
    select(all_of(id_col), compo, status, hslop, h0, vertical_resid, has_been_disturbed) %>%
    distinct()

  hex_with_data <- hex_grid %>%
    as.data.frame() %>%
    select(-any_of("compo")) %>%
    left_join(growth_dev, by = id_col) %>%
    filter(!is.na(compo))

  hex_with_data
}



# =============================================================================
# TABLE 2 — Plot distribution by reference acquisition year and species
# =============================================================================
cat("\n=== TABLE 2: PLOT DISTRIBUTION ===\n")

table2 <- mod_temp_exc$cv_predictions %>%
  distinct(plot_ID, compo, year_factor) %>%
  count(year_factor, compo) %>%
  pivot_wider(names_from = compo, values_from = n, values_fill = 0) %>%
  mutate(Total = rowSums(across(-year_factor))) %>%
  arrange(year_factor)

cat("  Table 2 — Plots per year × species:\n")
print(table2)


# =============================================================================
# FIGURE 3 — Biplots before/after harmonization (both models)
# =============================================================================
cat("\n=== FIGURE 3: HARMONIZATION BIPLOTS ===\n")

# Aggregate CV holdout predictions per plot
aggregate_cv <- function(model_obj) {
  model_obj$cv_predictions %>%
    group_by(plot_ID) %>%
    summarise(
      dominant_height  = first(dominant_height),
      CHM_height       = first(CHM_height),
      corrected_height = mean(corrected_height, na.rm = TRUE),
      .groups = "drop"
    )
}

cv_inc <- aggregate_cv(mod_temp_inc)
cv_exc <- aggregate_cv(mod_temp_exc)

scores_inc <- compute_cv_scores(mod_temp_inc)
scores_exc <- compute_cv_scores(mod_temp_exc)

# Panel (a): model WITH temporal variables
fig3a <- ggplot(cv_inc) +
  geom_point(aes(x = CHM_height, y = dominant_height, color = "Original"),
             alpha = 0.08) +
  geom_smooth(aes(x = CHM_height, y = dominant_height, color = "Original",
                  fill = "Original"), method = "lm") +
  geom_point(aes(x = corrected_height, y = dominant_height, color = "Harmonized\nw/ temporal var."),
             alpha = 0.08) +
  geom_smooth(aes(x = corrected_height, y = dominant_height, color = "Harmonized\nw/ temporal var.",
                  fill = "Harmonized\nw/ temporal var."), method = "lm") +
  geom_abline(slope = 1, linetype = 2) +
  scale_color_manual(name = "",
    values = c("Original" = "hotpink2", "Harmonized\nw/ temporal var." = "darkgreen")) +
  scale_fill_manual(name = "",
    values = c("Original" = "hotpink2", "Harmonized\nw/ temporal var." = "darkgreen")) +
  annotate("text", x = 5, y = 50, hjust = 0, vjust = 1, size = 3,
           label = paste0("r² = ", round(scores_inc$r2_before, 2),
                          "\nRMSE = ", round(scores_inc$rmse_before, 2),
                          "\nErr (m) = ", round(scores_inc$mean_bias_before, 2))) +
  annotate("text", x = 25, y = 10, hjust = 0, vjust = 0, size = 3,
           label = paste0("r² = ", round(scores_inc$r2_after, 2),
                          "\nRMSE = ", round(scores_inc$rmse_after, 2),
                          "\nErr (m) = ", round(scores_inc$mean_bias_after, 2))) +
  theme_article +
  labs(x = "CHM height (m)", y = "RFI (field) height (m)", title = "(a)") +
  coord_cartesian(xlim = c(5, 50), ylim = c(5, 52))

# Panel (b): model WITHOUT temporal variables
fig3b <- ggplot(cv_exc) +
  geom_point(aes(x = CHM_height, y = dominant_height, color = "Original"),
             alpha = 0.08) +
  geom_smooth(aes(x = CHM_height, y = dominant_height, color = "Original",
                  fill = "Original"), method = "lm") +
  geom_point(aes(x = corrected_height, y = dominant_height, color = "Harmonized\nw/out temporal var."),
             alpha = 0.08) +
  geom_smooth(aes(x = corrected_height, y = dominant_height, color = "Harmonized\nw/out temporal var.",
                  fill = "Harmonized\nw/out temporal var."), method = "lm") +
  geom_abline(slope = 1, linetype = 2) +
  scale_color_manual(name = "",
    values = c("Original" = "hotpink2", "Harmonized\nw/out temporal var." = "#7570b3")) +
  scale_fill_manual(name = "",
    values = c("Original" = "hotpink2", "Harmonized\nw/out temporal var." = "#7570b3")) +
  annotate("text", x = 5, y = 50, hjust = 0, vjust = 1, size = 3,
           label = paste0("r² = ", round(scores_exc$r2_before, 2),
                          "\nRMSE = ", round(scores_exc$rmse_before, 2),
                          "\nErr (m) = ", round(scores_exc$mean_bias_before, 2))) +
  annotate("text", x = 25, y = 10, hjust = 0, vjust = 0, size = 3,
           label = paste0("r² = ", round(scores_exc$r2_after, 2),
                          "\nRMSE = ", round(scores_exc$rmse_after, 2),
                          "\nErr (m) = ", round(scores_exc$mean_bias_after, 2))) +
  theme_article +
  labs(x = "CHM height (m)", y = "RFI (field) height (m)", title = "(b)") +
  coord_cartesian(xlim = c(5, 50), ylim = c(5, 52))

fig3 <- fig3a + fig3b + plot_layout(guides = "collect")
print(fig3)

# Regression summaries
cat("\n  (a) With temporal — original:\n")
print(summary(lm(dominant_height ~ CHM_height, data = cv_inc)))
cat("  (a) With temporal — standardized:\n")
print(summary(lm(dominant_height ~ corrected_height, data = cv_inc)))
cat("  (b) Without temporal — original:\n")
print(summary(lm(dominant_height ~ CHM_height, data = cv_exc)))
cat("  (b) Without temporal — standardized:\n")
print(summary(lm(dominant_height ~ corrected_height, data = cv_exc)))


# =============================================================================
# TABLE S1 — Per-year harmonization metrics
# =============================================================================
cat("\n=== TABLE S1: PER-YEAR HARMONIZATION METRICS ===\n")

table_s1 <- mod_temp_exc$cv_metrics_by_year %>%
  mutate(n_obs = n_obs / 5)   # each plot appears in 5 CV repetitions
cat("  Table S1:\n")
print(table_s1)


# =============================================================================
# TABLE S2 — Per-species harmonization metrics
# =============================================================================
cat("\n=== TABLE S2: PER-SPECIES HARMONIZATION METRICS ===\n")

table_s2 <- mod_temp_exc$cv_metrics_by_species %>%
  mutate(n_obs = n_obs / 5)
cat("  Table S2:\n")
print(table_s2)


# =============================================================================
# FIGURE 4 — Reference growth trajectories per species
# =============================================================================
cat("\n=== FIGURE 4: REFERENCE GROWTH TRAJECTORIES ===\n")

# Compute growth metrics for the subsample
subsample_growth <- compute_growth_metrics(ref_traj, id_col = "plot_ID")

# Keep only plots present in both statuses
both_status <- subsample_growth %>%
  filter(has_been_disturbed == F) %>%
  group_by(plot_ID, compo) %>%
  filter(all(c("raw", "corrected") %in% status)) %>%
  ungroup()

# Subsample for plotting (500 per species, reproducible)
set.seed(999)
plot_sample <- both_status %>%
  distinct(plot_ID, compo) %>%
  slice_sample(n = 500, by = "compo") %>%
  inner_join(subsample_growth, by = c("plot_ID", "compo"))

# Figure 4
fig4 <- plot_sample %>%
  ggplot(aes(x = h0, y = hslop, color = status, group = status)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.1) +
  geom_smooth(aes(fill = status), method = "lm") +
  facet_wrap(~compo, ncol = 3) +
  scale_color_manual(values = status_colors, labels = status_labels, name = "status") +
  scale_fill_manual(values = status_colors, labels = status_labels, name = "status") +
  theme_article +
  labs(
    x = "TCH in 2006 (m)",
    y = "Vertical growth trend over 15 yrs (m)"
  )

print(fig4)

# Growth regression table per species × status
lm_table <- plot_sample %>%
  group_by(compo, status) %>%
  summarise(
    n      = n(),
    lm_fit = list(lm(hslop ~ h0)),
    .groups = "drop"
  ) %>%
  mutate(
    intercept   = sapply(lm_fit, function(m) coef(m)[1]),
    slope       = sapply(lm_fit, function(m) coef(m)[2]),
    r2          = sapply(lm_fit, function(m) summary(m)$r.squared),
    p_intercept = sapply(lm_fit, function(m) summary(m)$coefficients[1, 4]),
    p_slope     = sapply(lm_fit, function(m) summary(m)$coefficients[2, 4])
  ) %>%
  select(-lm_fit) %>%
  mutate(across(where(is.numeric) & !n, ~ round(.x, 4))) %>%
  arrange(compo, status)

cat("\n  Growth regression table:\n")
print(lm_table)

# ANOVA: does harmonization shift the intercept / slope?
cat("\n  ANOVA — does status (raw/corrected) affect growth relationship?\n")
for (sp in unique(plot_sample$compo)) {
  sub <- filter(plot_sample, compo == sp)
  mod0 <- lm(hslop ~ h0, data = sub)
  mod1 <- lm(hslop ~ h0 + status, data = sub)
  mod2 <- lm(hslop ~ h0 * status, data = sub)
  p_int <- anova(mod0, mod1)$`Pr(>F)`[2]
  p_slp <- anova(mod1, mod2)$`Pr(>F)`[2]
  cat("  ", sp, ": intercept p =", formatC(p_int, format = "e", digits = 2),
      ", slope p =", formatC(p_slp, format = "e", digits = 2), "\n")
}

# Paired Wilcoxon test: raw vs. corrected growth slopes
wilcox_results <- plot_sample %>%
  distinct(plot_ID, compo, status, hslop) %>%
  pivot_wider(names_from = status, values_from = hslop) %>%
  group_by(compo) %>%
  summarise(
    wilcox_p = wilcox.test(raw, corrected, paired = TRUE)$p.value,
    .groups = "drop"
  )
cat("\n  Paired Wilcoxon test (raw vs. corrected growth slopes):\n")
print(wilcox_results)

# # Bootstrap stability check
# cat("\n  Bootstrap stability (1000 iterations)...\n")
# stability_check <- plot_sample %>%
#   filter(status == "corrected") %>%
#   group_by(compo) %>%
#   do({
#     boot_slopes <- replicate(1000, {
#       samp <- slice_sample(., n = nrow(.), replace = TRUE)
#       lm(hslop ~ h0, data = samp)$coefficients[2]
#     })
#     tibble(
#       slope_mean     = mean(boot_slopes),
#       slope_se       = sd(boot_slopes),
#       slope_ci_width = diff(quantile(boot_slopes, c(0.025, 0.975)))
#     )
#   })
# print(stability_check)

# Reference growth coefficients (used for Figures 5–6 and Table 4)
ref_coefs_both <- plot_sample %>%
  group_by(compo, status) %>%
  summarise(
    intercept_lm = lm(hslop ~ h0)$coefficients[1],
    slope_lm     = lm(hslop ~ h0)$coefficients[2],
    .groups = "drop"
  )

cat("  Reference growth coefficients:\n")
print(ref_coefs_both)

# =============================================================================
# FIGURE 5 — Vertical growth and deviation maps (Site A)
# =============================================================================
cat("\n=== FIGURE 5: SITE A MAPS ===\n")

# Reusable map function for any site
create_site_map <- function(site_label, ref_coefs,
                            zoom_half_x = 2500, zoom_half_y = 500) {

  ts_data  <- tile_ts[[site_label]]
  hex_grid <- tile_hex[[site_label]]
  if (is.null(ts_data) || is.null(hex_grid)) {
    cat("  SKIPPED", site_label, "— data not loaded\n")
    return(NULL)
  }

  hex_data <- prepare_tile_for_mapping(ts_data, hex_grid, ref_coefs)
  hex_corr <- hex_data %>% filter(status == "corrected")

  # Zoom window centered on tile extent
  bbox <- st_bbox(st_as_sf(hex_grid))
  x_center <- (bbox["xmin"] + bbox["xmax"]) / 2
  y_center <- (bbox["ymin"] + bbox["ymax"]) / 2
  xlim <- c(x_center - zoom_half_x, x_center + zoom_half_x)
  ylim <- c(y_center - zoom_half_y, y_center + zoom_half_y)

  coord_panel <- coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
                          datum = pull_crs(hex_grid))
  theme_panel <- theme_dark() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

  # Helper: merge data onto hex grid for spatial plotting
  hex_merge <- function(cols) hex_grid %>% merge(hex_corr %>% select(ID, all_of(cols)))
  
  # Setup : merge disturbed plots onto hex grid 
  hex_disturbed_agg <- hex_merge("has_been_disturbed") %>%
    terra::aggregate(by = 'has_been_disturbed') %>%
    tidyterra::filter(has_been_disturbed == T)
  
  hex_stands <- hex_merge("compo")%>%
    terra::aggregate(by = 'compo')

  # Panel (b): Dominant species
  p_species <- ggplot() +
    geom_spatvector(data = hex_merge("compo"), aes(fill = compo), color = NA) +
    scale_fill_manual(values = species_colors, na.value = "transparent",
                      name = "Dominant species") +
    coord_panel + 
    theme_panel + labs(title = "(b)")

  # Panel (c): Initial height (2006)
  p_height <- ggplot() +
    geom_spatvector(data = hex_merge("h0"), aes(fill = h0), color = NA) +
    geom_spatvector(data = hex_stands, color = 'grey30', fill = NA) +
    scale_fill_viridis_c(name = "Initial height\n(2006 ; m)") +
    coord_panel + 
    theme_panel + labs(title = "(c)")

  # Panel (d): Growth trend — with disturbance hatching
  hex_growth_v <- hex_merge(c("hslop", "has_been_disturbed"))
  p_growth <- ggplot() +
    geom_spatvector(data = hex_growth_v, aes(fill = hslop), color = NA) +
    scale_fill_gradient2(
      low = "#d7191c", mid = "#ffffff", high = "#1a9641", midpoint = 0,
      name = "Growth\ntrend (m/yr)", limits = c(-1, 1), oob = scales::squish
    ) +
    geom_spatvector(data = hex_stands, color = 'grey30', fill = NA) +
    geom_sf_pattern(
      data = hex_disturbed_agg,
      color = 'black',
      pattern_spacing = 0.04,
      pattern_density = 0.2,
      pattern_fill = 'black',
      fill = NA,
      linewidth = 0.5
    ) +
    coord_panel + 
    theme_panel + labs(title = "(d)")

  # Panel (e): Deviation from reference trajectory — with disturbance hatching
  hex_dev_v <- hex_merge(c("vertical_resid", "has_been_disturbed"))
  p_deviation <- ggplot() +
    geom_spatvector(data = hex_dev_v, aes(fill = vertical_resid), color = NA) +
    scale_fill_gradient2(
      low = "#e66101", mid = "#ffffff", high = "#5e3c99", midpoint = 0,
      name = "Deviation from\nreference\ntrajectory (m/yr)",
      limits = c(-0.5, 0.5), oob = scales::squish
    ) +
    geom_spatvector(data = hex_stands, color = 'grey30', fill = NA) +
    geom_sf_pattern(
      data = hex_disturbed_agg,
      color = 'black',
      pattern_spacing = 0.04,
      pattern_density = 0.2,
      pattern_fill = 'black',
      fill = NA,
      linewidth = 0.5
    ) +
    coord_panel + 
    theme_panel + labs(title = "(e)")

  # Panel (f): Ridgeline distribution of deviation per species
  # Panel (f): Ridgeline distribution of deviation per species
  zoom_ext <- ext(xlim[1], xlim[2], ylim[1], ylim[2])
  hex_corr_vect <- hex_merge(c("vertical_resid", "has_been_disturbed", "compo"))
  hex_cropped <- terra::crop(hex_corr_vect, zoom_ext)
  
  df_ridge <- as.data.frame(hex_cropped) %>%
    filter(!has_been_disturbed, !is.na(vertical_resid))
  
  p_ridge <- ggplot(df_ridge, aes(x = vertical_resid, y = compo, fill = compo)) +
    geom_density_ridges(alpha = 0.75) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_fill_manual(values = species_colors) +
    theme_bw() +
    labs(x = "Deviation from the reference trajectory (m/yr)", y = "",
         title = "(f)") +
    coord_cartesian(xlim = c(-0.6, 0.6)) +
    theme(legend.position = 'none', 
          aspect.ratio = 0.8/4,
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank()) 

  # Combine all panels
  combined <- (p_species / p_height / p_growth / p_deviation / p_ridge) +
    plot_layout(ncol = 1, heights = c(1, 1, 1, 1, 1, 1))

  list(
    panels   = list(species = p_species, height = p_height, growth = p_growth,
                    deviation = p_deviation, ridge = p_ridge),
    combined = combined,
    hex_data = hex_data
  )
}

# Generate maps for all sites
site_maps <- list()
for (lbl in tile_config$site_label) {
  cat("  Processing", lbl, "...\n")
  site_maps[[lbl]] <- create_site_map(lbl, ref_coefs = ref_coefs_both)
}

# Print Figure 5 (Site A)
if (!is.null(site_maps[["Site A"]])) {
  print(site_maps[["Site A"]]$combined)
}


# =============================================================================
# FIGURE 6 — Beech deviation across sites
# =============================================================================
cat("\n=== FIGURE 6: BEECH DEVIATION ACROSS SITES ===\n")

focus_species <- "Beech"

# Panel (a): Per-site maps with beech highlighted, others greyed out
beech_map_panels <- list()
for (lbl in names(site_maps)) {
  if (is.null(site_maps[[lbl]])) next

  hex_data <- site_maps[[lbl]]$hex_data %>% filter(status == "corrected")
  hex_grid <- tile_hex[[lbl]]

  hex_dev <- hex_grid %>%
    merge(hex_data %>% select(ID, compo, vertical_resid, has_been_disturbed))

  # Zoom window
  bbox <- st_bbox(st_as_sf(hex_grid))
  x_center <- (bbox["xmin"] + bbox["xmax"]) / 2
  y_center <- (bbox["ymin"] + bbox["ymax"]) / 2
  xlim <- c(x_center - 2500, x_center + 2500)
  if(lbl == 'Site C'){
  ylim <- c(y_center + 0, y_center + 1000)
  } else {
    ylim <- c(y_center - 500, y_center + 500)
  }

  # Grey out non-focus species
  hex_focus    <- hex_dev[hex_dev$compo == focus_species, ]
  hex_nonfocus <- hex_dev[hex_dev$compo != focus_species, ]

  p <- ggplot() +
    geom_spatvector(data = hex_nonfocus, fill = "grey80", color = NA, alpha = 0.6) +
    geom_spatvector(data = hex_focus, aes(fill = vertical_resid), color = NA) +
    scale_fill_gradient2(
      low = "#e66101", mid = "#ffffff", high = "#5e3c99", midpoint = 0,
      name = "Deviation from\nreference\ntrajectory (m/yr)",
      limits = c(-0.5, 0.5), oob = scales::squish
    ) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
             datum = pull_crs(hex_grid)) +
    theme_dark() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
    labs(title = lbl) + theme(legend.position = 'none')

  beech_map_panels[[lbl]] <- p
}

# Panel (b): Ridgeline distribution of beech deviation per site
df_beech_all <- map_dfr(names(site_maps), function(lbl) {
  if (is.null(site_maps[[lbl]])) return(NULL)
  hex_grid <- tile_hex[[lbl]]
  bbox <- st_bbox(st_as_sf(hex_grid))
  x_center <- (bbox["xmin"] + bbox["xmax"]) / 2
  y_center <- (bbox["ymin"] + bbox["ymax"]) / 2
  if (lbl == "Site C") {
    zoom_ext <- ext(x_center - 2500, x_center + 2500, y_center + 0, y_center + 1000)
  } else {
    zoom_ext <- ext(x_center - 2500, x_center + 2500, y_center - 500, y_center + 500)
  }
  
  hex_data_vect <- hex_grid %>%
    merge(site_maps[[lbl]]$hex_data %>% 
            filter(status == "corrected") %>%
            select(ID, compo, vertical_resid, has_been_disturbed))
  hex_cropped <- terra::crop(hex_data_vect, zoom_ext) %>%
    as.data.frame()
  
  hex_cropped %>%
    filter(
      !has_been_disturbed,
      compo == focus_species,
      !is.na(vertical_resid)
    ) %>%
    mutate(site = lbl)
})

fig6b <- ggplot(
  df_beech_all,
  aes(x = vertical_resid, y = factor(site, levels = rev(unique(site))),
      fill = after_stat(x))
) +
  geom_density_ridges_gradient(
    scale = 1.25, rel_min_height = 0.001, alpha = 0.75
  ) +
  scale_fill_gradient2(
    low = "#e66101", mid = "#ffffff", high = "#5e3c99", midpoint = 0,
    name = "Deviation\n(m/yr)", limits = c(-0.5, 0.5), oob = scales::squish
  ) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  theme_bw() +
  labs(x = "Deviation from\nreference trajectory (m/yr)", y = "") +
  coord_cartesian(xlim = c(-0.7, 0.7)) +
  theme(legend.position = "none", aspect.ratio = 3 / 1.8)

# Combine Figure 6
if (length(beech_map_panels) > 0) {
  fig6a <- wrap_plots(beech_map_panels, ncol = 1) + plot_layout(guides = "collect")
  fig6 <- fig6a | fig6b
  print(fig6)
}


# =============================================================================
# TABLE 4 — Moran's I spatial autocorrelation
# =============================================================================
cat("\n=== TABLE 4: MORAN'S I (ZOOMED AREA ONLY) ===\n")

moran_results <- list()

for (lbl in names(site_maps)) {
  if (is.null(site_maps[[lbl]])) next
  
  # --- Retrieve data ---
  hex_data_full <- site_maps[[lbl]]$hex_data %>%
    filter(status == "corrected",
           !has_been_disturbed,
           !is.na(vertical_resid))
  
  hex_grid <- tile_hex[[lbl]]
  hex_sf   <- st_as_sf(hex_grid)
  
  # --- Recreate zoom window (same as Figure 5) ---
  bbox <- st_bbox(hex_sf)
  
  x_center <- mean(as.numeric(bbox[c("xmin", "xmax")]))
  y_center <- mean(as.numeric(bbox[c("ymin", "ymax")]))
  
  zoom_bbox <- st_bbox(c(
    xmin = x_center - 2500,
    ymin = y_center - 500,
    xmax = x_center + 2500,
    ymax = y_center + 500
  ), crs = st_crs(hex_sf))
  
  zoom_window <- st_as_sfc(zoom_bbox)
  
  # --- Crop to zoom area ---
  hex_zoom <- st_intersection(hex_sf, zoom_window)
  
  # Join metrics
  hex_joined <- hex_zoom %>%
    left_join(hex_data_full %>%
                select(ID, vertical_resid, hslop),
              by = "ID") %>%
    filter(!is.na(vertical_resid))
  
  if (nrow(hex_joined) < 10) {
    cat(" ", lbl, ": too few observations in zoom area\n")
    next
  }
  
  # --- Build neighbors on zoomed centroids ---
  coords <- st_coordinates(st_centroid(hex_joined))
  knn <- knearneigh(coords, k = 8)
  nb  <- knn2nb(knn)
  lw  <- nb2listw(nb, style = "W")
  
  # --- Moran tests ---
  mi_growth <- moran.test(hex_joined$hslop, lw)
  mi_dev    <- moran.test(hex_joined$vertical_resid, lw)
  
  moran_results[[lbl]] <- tibble(
    Site      = lbl,
    Growth    = paste0(round(mi_growth$estimate["Moran I statistic"], 2),
                       " (p ",
                       ifelse(mi_growth$p.value < 0.001, "< 0.001",
                              paste0("= ", round(mi_growth$p.value, 3))),
                       ")"),
    Deviation = paste0(round(mi_dev$estimate["Moran I statistic"], 2),
                       " (p ",
                       ifelse(mi_dev$p.value < 0.001, "< 0.001",
                              paste0("= ", round(mi_dev$p.value, 3))),
                       ")")
  )
  
  cat(" ", lbl,
      "— Growth:", moran_results[[lbl]]$Growth,
      "| Deviation:", moran_results[[lbl]]$Deviation, "\n")
}

table4 <- bind_rows(moran_results)

cat("\nTable 4 (Zoomed area only):\n")
print(table4)

cat("\n=== ALL RESULTS COMPLETE ===\n")
