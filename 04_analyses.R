# =============================================================================
# CODE 4 — RESULTS AND FIGURES
# =============================================================================
# PURPOSE:
#   Produce all code-generated tables and figures for the article, organized by
#   the REVISED manuscript / SI numbering. Each section can be run independently
#   once data are loaded (§2).
#
#   Figures 1-2 and Tables 1, 3 are not code-generated (workflow diagram,
#   study area map, survey table, candidate predictors table). Appendix S1
#   (Figure S1 / Table S1, percentile sensitivity) is produced elsewhere.
#
#   This script produces:
#     - Table 2    : Plot distribution by reference year x species
#     - Figure 3   : Distribution of dominant-height bias before/after harmon.
#     - Table S2   : GAM coefficients (parametric + smooth + fit), both models
#     - Figure S2  : GAM partial effects of retained smooth terms
#     - Table S3   : Per-year harmonization metrics
#     - Table S4   : Per-species harmonization metrics
#     - Figure 4   : Species-specific reference growth trajectories (exp. decay)
#     - Figure 5   : Vertical growth and deviation maps (Site A)
#     - Figure 6   : Beech deviation across sites
#     - Table 4    : Moran's I spatial autocorrelation
#
# INPUTS:
#   - model_temp_inc.rds         (GAM with temporal vars; $final_model + $cv_*)
#   - model_temp_exc.rds         (GAM without temporal vars; $final_model + $cv_*)
#   - reference_trajectories.rds (cleaned subsample time series)
#   - tile_{id}_timeseries.rds   (cleaned tile time series)
#   - tile_{id}_hexgrid.rds      (hexagonal grids)
#   - tile_{id}_ortho.tif        (RGB orthoimage, optional; Site A only)
#
# NOTE ON THE GROWTH MODEL:
#   Reference trajectories are now fitted as an EXPONENTIAL DECAY of the growth
#   trend on initial height:  hslop = Ka * Kb ^ h0   (Eq. 4), via nls with a
#   log-linear fallback. Cell-level deviation is the residual to that curve
#   (Eq. 5):  vertical_resid = hslop - Ka * Kb ^ h0.
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
library(gratia)    # smooth_estimates() for GAM partial effects (Figure S2)
library(ggtext)    # element_markdown() for HTML legend titles (Figure 5)
library(ggpubr)    # stat_compare_means() for the growth diagnostic
library(knitr)     # kable() for Table S2

# --- 1.2 File paths ---
path_data <- "outputs/data/"

# --- 1.3 Color palettes and labels ---
species_colors <- c(
  "Non_Ligneous"     = "#ffffff",
  "Other_coniferous" = "#4c2772",
  "Other_broadleaf"  = "#a6cee3",
  "Birch"            = "#1f78b4",
  "Oak"              = "#b2df8a",
  "Douglas_fir"      = "#33a02c",
  "Spruce"           = "#fb9a99",
  "Beech"            = "#e31a1c",
  "Larch"            = "#fdbf6f",
  "Poplar"           = "#ff7f00",
  "Pine"             = "#cab2d6"
)

species_labels <- c(
  "Beech"            = "Beech",
  "Birch"            = "Birch",
  "Douglas_fir"      = "Douglas-fir",
  "Larch"            = "Larch",
  "Oak"              = "Oak",
  "Other_broadleaf"  = "Other broadleaves",
  "Other_coniferous" = "Other conifers",
  "Pine"             = "Pine",
  "Spruce"           = "Spruce"
)

status_colors <- c("raw" = "hotpink2", "corrected" = "darkgreen")
status_labels <- c("raw" = "Original\n(biased)",
                   "corrected" = "Without temporal\nvariables")

# Harmonization palette (Figures 3 and S2): three conditions
harmon_colors <- c(
  "Original\n(biased)"          = "hotpink2",
  "Without temporal\nvariables" = "darkgreen",
  "With temporal\nvariables"    = "#7570b3"
)

# Diverging scales reused across maps
scale_dev  <- function() scale_fill_gradient2(
  low = "#e66101", mid = "#ffffff", high = "#5e3c99", midpoint = 0,
  name = "Deviation from the<br>reference trajectory<br>(m/yr)",
  limits = c(-0.5, 0.5), oob = scales::squish)

# --- 1.4 Shared theme ---
theme_article <- theme_bw() +
  theme(
    text            = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.text     = element_text(size = 7),
    legend.position = "right"
  )

# --- 1.5 Tile configuration ---
# Map tile IDs to manuscript site labels. The zoom window is x_center +/- 2500 m
# with a site-specific vertical offset (Site C is shifted up). `ortho` names the
# RGB basemap used for panel (a) of Figure 5; NA where absent.
tile_config <- tibble(
  tile_id    = c(153, 45, 484),      # <-- adjust these IDs to your data
  site_label = c("Site A", "Site B", "Site C"),
  y_off_lo   = c(-500, -500,    0),  # vertical zoom offsets (m) from y_center
  y_off_hi   = c( 500,  500, 1000),
  ortho      = c("tile_153_ortho.tif", NA, NA)
)

zoom_half_x <- 2500  # horizontal half-window (m), constant across sites


# =============================================================================
# 2. LOAD DATA
# =============================================================================
cat("=== 2. LOADING DATA ===\n")

# Models (each object carries $final_model and the $cv_* fields)
mod_temp_inc <- readRDS(file.path(path_data, "model_temp_inc.rds"))
mod_temp_exc <- readRDS(file.path(path_data, "model_temp_exc.rds"))

# Reference-trajectory subsample time series
ref_traj <- readRDS(file.path(path_data, "reference_trajectories.rds"))

# Tile outputs
tile_ts  <- list()
tile_hex <- list()
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
    cat("  SKIPPED", lbl, "(tile", tid, ") - files not found\n")
  }
}
cat("  All available data loaded.\n")


# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

#' Aggregate CV holdout predictions to one row per plot (each plot appears in
#' 5 CV repetitions).
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

#' Compute harmonization scores (bias / RMSE / R2, before and after) from a
#' model object's CV predictions.
compute_cv_scores <- function(model_obj) {
  aggregate_cv(model_obj) %>%
    summarise(
      n_plots          = n(),
      mean_bias_before = mean(dominant_height - CHM_height, na.rm = TRUE),
      sd_bias_before   = sd(dominant_height - CHM_height, na.rm = TRUE),
      rmse_before      = sqrt(mean((dominant_height - CHM_height)^2, na.rm = TRUE)),
      r2_before        = 1 - sum((dominant_height - CHM_height)^2) /
        sum((dominant_height - mean(dominant_height))^2),
      mean_bias_after  = mean(dominant_height - corrected_height, na.rm = TRUE),
      sd_bias_after    = sd(dominant_height - corrected_height, na.rm = TRUE),
      rmse_after       = sqrt(mean((dominant_height - corrected_height)^2, na.rm = TRUE)),
      r2_after         = 1 - sum((dominant_height - corrected_height)^2) /
        sum((dominant_height - mean(dominant_height))^2)
    )
}

#' Extract parametric + smooth term tables and fit stats from a gam object.
extract_gam_table <- function(model, model_name) {
  s <- summary(model)
  
  param <- as.data.frame(s$p.table) %>%
    tibble::rownames_to_column("term") %>%
    rename(estimate = Estimate, se = `Std. Error`,
           stat = `t value`, p = `Pr(>|t|)`) %>%
    mutate(type = "parametric", model = model_name)
  
  smooth <- as.data.frame(s$s.table) %>%
    tibble::rownames_to_column("term") %>%
    rename(edf = edf, ref.df = Ref.df, stat = F, p = `p-value`) %>%
    mutate(type = "smooth", model = model_name)
  
  list(param = param, smooth = smooth,
       r2 = s$r.sq, dev_expl = s$dev.expl, n = nrow(model$model))
}

#' Fit the exponential decay hslop = Ka * Kb^h0 via bounded nls, with a
#' log-linear fallback on positive growth only. Returns Ka, Kb, R2, n, method.
fit_exponential <- function(data) {
  fit <- tryCatch(
    nls(hslop ~ Ka * Kb^h0,
        data      = data,
        start     = list(Ka = 0.8, Kb = 0.95),
        lower     = list(Ka = 0.01, Kb = 0.5),
        upper     = list(Ka = 3.0,  Kb = 1.0),
        algorithm = "port",
        control   = nls.control(maxiter = 200)),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    r2 <- 1 - sum(residuals(fit)^2) / sum((data$hslop - mean(data$hslop))^2)
    tibble(Ka = coef(fit)[["Ka"]], Kb = coef(fit)[["Kb"]],
           r2_exp = r2, n_exp = nrow(data), method = "nls")
  } else {
    d_pos <- data %>% filter(hslop > 0)
    if (nrow(d_pos) < 10) {
      return(tibble(Ka = NA_real_, Kb = NA_real_, r2_exp = NA_real_,
                    n_exp = nrow(data), method = "failed"))
    }
    m <- lm(log(hslop) ~ h0, data = d_pos)
    tibble(Ka = exp(coef(m)[1]), Kb = exp(coef(m)[2]),
           r2_exp = summary(m)$r.squared, n_exp = nrow(d_pos),
           method = "log-linear fallback")
  }
}

#' Compute initial height (h0 at ref_year) and growth slope (hslop) for each
#' unit x status. h0 is read BEFORE the spike/leaf-off filtering, since the
#' reference year may itself be flagged.
compute_growth_metrics <- function(ts_data, id_col = "ID", ref_year = 2006) {
  ts_data %>%
    mutate(
      year_flight = as.numeric(as.character(year_factor)),
      leaf_off    = !is.na(leaf_off) & leaf_off
    ) %>%
    group_by(.data[[id_col]], status) %>%
    arrange(year_flight) %>%
    mutate(h0 = first(height[year_flight == ref_year], default = NA_real_)) %>%
    filter(spike == FALSE | is.na(spike), !leaf_off) %>%
    filter(any(year_flight == ref_year), any(!is.na(height))) %>%
    mutate(hslop = lm(height ~ year_flight)$coefficients[2]) %>%
    ungroup() %>%
    group_by(.data[[id_col]], status) %>%
    summarise(
      compo              = as.character(first(compo)),
      h0                 = first(h0),
      hslop              = first(hslop),
      has_been_disturbed = any(has_been_disturbed, na.rm = TRUE),
      spike              = any(spike, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(h0), !is.na(hslop))
}

#' Prepare a tile for mapping: growth metrics + exponential-decay deviations
#' joined onto the hex grid geometry.
prepare_tile_for_mapping <- function(ts_data, hex_grid, ref_coefs, id_col = "ID") {
  growth_dev <- compute_growth_metrics(ts_data, id_col = id_col) %>%
    left_join(ref_coefs, by = c("compo", "status")) %>%
    mutate(
      predicted      = Ka * Kb^h0,
      vertical_resid = hslop - predicted
    ) %>%
    filter(!is.na(Ka)) %>%
    select(all_of(id_col), compo, status, hslop, h0,
           vertical_resid, has_been_disturbed, spike) %>%
    distinct()
  
  hex_grid %>%
    select(-any_of("compo")) %>%
    left_join(growth_dev, by = id_col) %>%
    filter(!is.na(compo))
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

print(table2)


# =============================================================================
# FIGURE 3 — Distribution of dominant-height bias before/after harmonization
# =============================================================================
cat("\n=== FIGURE 3: BIAS DISTRIBUTION ===\n")

cv_inc <- aggregate_cv(mod_temp_inc)
cv_exc <- aggregate_cv(mod_temp_exc)

# In-text harmonization scores (bias / RMSE / R2)
scores_inc <- compute_cv_scores(mod_temp_inc)
scores_exc <- compute_cv_scores(mod_temp_exc)
cat("  With temporal:\n");    print(t(scores_inc))
cat("  Without temporal:\n"); print(t(scores_exc))

# Long form: bias = dominant_height - height, across the three conditions
bias_long <- bind_rows(
  cv_inc %>% transmute(dominant_height, height = CHM_height,
                       group3 = "Original\n(biased)"),
  cv_inc %>% transmute(dominant_height, height = corrected_height,
                       group3 = "With temporal\nvariables"),
  cv_exc %>% transmute(dominant_height, height = corrected_height,
                       group3 = "Without temporal\nvariables")
)

fig3 <- ggplot(bias_long, aes(x = dominant_height - height,
                              color = group3, fill = group3)) +
  geom_density(alpha = 0.2) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
  scale_color_manual(name = "Harmonization", values = harmon_colors) +
  scale_fill_manual(name = "Harmonization", values = harmon_colors) +
  labs(x = expression("Bias = "*H[dom*","*RFI] - H[dom*","*CHM]*" (m)"),
       y = "Density") +
  theme_article

print(fig3)
ggsave(file.path(path_data, "Figure3_bias_density.png"),
       fig3, width = 90, height = 60, units = "mm", dpi = 500)


# =============================================================================
# TABLE S2 + FIGURE S2 — GAM coefficients and partial effects
# =============================================================================
cat("\n=== TABLE S2 / FIGURE S2: GAM PERFORMANCE ===\n")

res_inc <- extract_gam_table(mod_temp_inc$final_model, "with_time")
res_exc <- extract_gam_table(mod_temp_exc$final_model, "no_time")

# --- Table S2 ---
param_combined  <- bind_rows(res_inc$param,  res_exc$param)
smooth_combined <- bind_rows(res_inc$smooth, res_exc$smooth)
fit_summary <- tibble(
  model        = c("with_time", "no_time"),
  r2_adj       = c(res_inc$r2, res_exc$r2),
  dev_expl_pct = c(res_inc$dev_expl, res_exc$dev_expl) * 100,
  n            = c(res_inc$n, res_exc$n)
)

print(kable(param_combined  %>% select(model, term, estimate, se, stat, p), digits = 3))
print(kable(smooth_combined %>% select(model, term, edf, ref.df, stat, p), digits = 3))
print(kable(fit_summary, digits = 3))

# --- Figure S2: partial effects of the three retained smooth terms ---
# Predictor names come from Code 1: broadleaf-to-conifer ratio -> broadleaf_ratio,
# canopy heterogeneity -> CHM_hsd, gap fraction -> CHM_gap_fraction.
smooth_terms <- c("s(broadleaf_ratio)", "s(CHM_hsd)", "s(CHM_gap_fraction)")

sm_combined <- bind_rows(
  smooth_estimates(mod_temp_inc$final_model) %>%
    filter(.smooth %in% smooth_terms) %>% mutate(model = "With temporal\nvariables"),
  smooth_estimates(mod_temp_exc$final_model) %>%
    filter(.smooth %in% smooth_terms) %>% mutate(model = "Without temporal\nvariables")
) %>%
  mutate(x_value = coalesce(broadleaf_ratio, CHM_hsd, CHM_gap_fraction))

make_panel <- function(smooth_name, xlab, show_y = TRUE) {
  ggplot(filter(sm_combined, .smooth == smooth_name),
         aes(x = x_value, y = .estimate, color = model, fill = model)) +
    geom_ribbon(aes(ymin = .estimate - .se, ymax = .estimate + .se),
                alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.8) +
    coord_cartesian(ylim = c(-2, 2)) +
    scale_color_manual(name = "Harmonization", values = harmon_colors) +
    scale_fill_manual(name = "Harmonization", values = harmon_colors) +
    labs(x = xlab, y = if (show_y) "Partial effect on bias (m)" else NULL) +
    theme_bw()
}

figS2 <- (make_panel("s(broadleaf_ratio)",   "Broadleaf-to-conifer ratio (%)") +
            make_panel("s(CHM_hsd)",           "Canopy heterogeneity (m)", show_y = FALSE) +
            make_panel("s(CHM_gap_fraction)",  "Gap fraction (%)",         show_y = FALSE)) +
  plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") &
  theme_article

print(figS2)
ggsave(file.path(path_data, "FigureS2_GAM_partial_effects.png"),
       figS2, width = 190, height = 60, units = "mm", dpi = 500)


# =============================================================================
# TABLE S3 — Per-year harmonization metrics
# =============================================================================
cat("\n=== TABLE S3: PER-YEAR METRICS ===\n")

table_s3 <- mod_temp_exc$cv_metrics_by_year %>%
  mutate(n_obs = n_obs / 5)   # each plot appears in 5 CV repetitions
print(table_s3)


# =============================================================================
# TABLE S4 — Per-species harmonization metrics
# =============================================================================
cat("\n=== TABLE S4: PER-SPECIES METRICS ===\n")

table_s4 <- mod_temp_exc$cv_metrics_by_species %>%
  mutate(n_obs = n_obs / 5)
print(table_s4)


# =============================================================================
# FIGURE 4 — Species-specific reference growth trajectories (exp. decay)
# =============================================================================
cat("\n=== FIGURE 4: REFERENCE TRAJECTORIES ===\n")

set.seed(999)  # reproducible 500-plot-per-species subsample

# Growth metrics for the reference subsample; keep plots present in BOTH
# statuses whose raw initial height exceeds 2 m.
subsample_growth <- compute_growth_metrics(ref_traj, id_col = "plot_ID") %>%
  filter(!has_been_disturbed) %>%
  group_by(plot_ID) %>%
  filter(any(h0[status == "raw"] > 2)) %>%
  ungroup()

both_status <- subsample_growth %>%
  group_by(plot_ID, compo) %>%
  filter(all(c("raw", "corrected") %in% status)) %>%
  ungroup()

plot_sample <- both_status %>%
  distinct(plot_ID, compo) %>%
  slice_sample(n = 500, by = "compo") %>%
  inner_join(subsample_growth, by = c("plot_ID", "compo"))

# Exponential-decay fits per species x status -> reference coefficients (Ka, Kb)
ref_coefs <- plot_sample %>%
  group_by(compo, status) %>%
  group_modify(~ fit_exponential(.x)) %>%
  ungroup()

cat("  Reference coefficients (Ka, Kb):\n")
print(ref_coefs %>% select(compo, status, Ka, Kb, r2_exp, method))

# Prediction curves and equation labels
pred_lines <- plot_sample %>%
  group_by(compo, status) %>%
  summarise(h0 = list(seq(min(h0), max(h0), length.out = 100)), .groups = "drop") %>%
  unnest(h0) %>%
  left_join(ref_coefs, by = c("compo", "status")) %>%
  mutate(hslop_pred = Ka * Kb^h0)

fits_labeled <- ref_coefs %>%
  mutate(
    label  = sprintf('bold(y == "%s" %%*%% "%s"^%s)',
                     formatC(Ka, format = "f", digits = 2),
                     formatC(Kb, format = "f", digits = 2), "x"),
    status = factor(status, levels = c("corrected", "raw"))
  ) %>%
  arrange(compo, status) %>%
  group_by(compo) %>%
  mutate(vjust_val = c(1.2, 2.3)[row_number()]) %>%
  ungroup()

fig4 <- ggplot(plot_sample, aes(x = h0, y = hslop, color = status, group = status)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(alpha = 0.1) +
  geom_line(data = pred_lines, aes(y = hslop_pred),
            linewidth = 1, alpha = 0.75, inherit.aes = TRUE) +
  geom_text(data = fits_labeled,
            aes(x = Inf, y = Inf, label = label, color = status, vjust = vjust_val),
            hjust = 1.05, size = 3, parse = TRUE, show.legend = FALSE,
            inherit.aes = FALSE) +
  facet_wrap(vars(compo), ncol = 3, labeller = labeller(compo = species_labels)) +
  scale_color_manual(values = status_colors, labels = status_labels,
                     name = "Harmonization") +
  labs(x = expression("Initial top-of-canopy height (TCH"["2006"]*"; m)"),
       y = "Vertical growth trend over 2006-2021 (m/yr)") +
  theme_article

print(fig4)
ggsave(file.path(path_data, "Figure4_reference_trajectories.png"),
       fig4, width = 190, height = 150, units = "mm", dpi = 500)

# --- In-text growth statistics (§3.2) ---
plot_sample %>%
  group_by(status) %>%
  summarise(mean_gwth = mean(hslop, na.rm = TRUE),
            sd_gwth   = sd(hslop, na.rm = TRUE)) %>%
  print()

# Paired Wilcoxon: raw vs corrected growth slopes, per species and overall
wilcox_by_species <- plot_sample %>%
  distinct(plot_ID, compo, status, hslop) %>%
  pivot_wider(names_from = status, values_from = hslop) %>%
  group_by(compo) %>%
  summarise(wilcox_p = wilcox.test(raw, corrected, paired = TRUE)$p.value,
            .groups = "drop")
cat("  Paired Wilcoxon (raw vs corrected), by species:\n")
print(wilcox_by_species)

wilcox_overall <- plot_sample %>%
  distinct(plot_ID, status, hslop) %>%
  pivot_wider(names_from = status, values_from = hslop)
print(wilcox.test(wilcox_overall$raw, wilcox_overall$corrected, paired = TRUE))

# Optional diagnostic: nested ANOVA of growth by species / status
print(summary(aov(hslop ~ compo / status, data = plot_sample)))


# =============================================================================
# FIGURE 5 — Vertical growth and deviation maps (Site A)
# =============================================================================
cat("\n=== FIGURE 5: SITE MAPS ===\n")

# Reusable striped-hatch layer marking disturbed cells
disturbance_hatch <- function(hex_vect) {
  geom_sf_pattern(
    data = sf::st_as_sf(
      terra::aggregate(hex_vect, by = "has_been_disturbed") %>%
        tidyterra::filter(has_been_disturbed == TRUE)
    ),
    pattern = "stripe", fill = NA, color = "black",
    pattern_spacing = 0.02, pattern_density = 0.2,
    pattern_fill = "black", linewidth = 0.5
  )
}

#' Build the six-panel site figure: (a) ortho, (b) species, (c) initial height,
#' (d) growth trend, (e) deviation, (f) ridge distribution of deviation.
create_site_map <- function(site_label, ref_coefs) {
  ts_data  <- tile_ts[[site_label]]
  hex_grid <- tile_hex[[site_label]]
  if (is.null(ts_data) || is.null(hex_grid)) {
    cat("  SKIPPED", site_label, "- data not loaded\n"); return(NULL)
  }
  
  cfg      <- tile_config %>% filter(site_label == !!site_label)
  hex_data <- prepare_tile_for_mapping(ts_data, hex_grid, ref_coefs)
  hex_corr <- hex_data %>% filter(status == "corrected")
  hex_stands <- terra::aggregate(hex_corr, by = "compo")
  
  # Zoom window
  bbox <- st_bbox(st_as_sf(hex_grid))
  x_c  <- (bbox["xmin"] + bbox["xmax"]) / 2
  y_c  <- (bbox["ymin"] + bbox["ymax"]) / 2
  xlim <- c(x_c - zoom_half_x, x_c + zoom_half_x)
  ylim <- c(y_c + cfg$y_off_lo, y_c + cfg$y_off_hi)
  
  coord_panel <- coord_sf(xlim = xlim, ylim = ylim, expand = FALSE,
                          datum = pull_crs(hex_grid))
  theme_panel <- theme_dark() +
    theme(panel.border   = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "right",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_markdown())
  
  # (a) Orthoimage basemap (Site A only; skipped if the file is absent)
  p_a <- NULL
  if (!is.na(cfg$ortho) && file.exists(file.path(path_data, cfg$ortho))) {
    baserast <- rast(file.path(path_data, cfg$ortho))
    p_a <- ggplot() +
      geom_spatraster_rgb(data = baserast, stretch = "lin") +
      coord_panel + theme_panel + labs(title = "(a)")
  }
  
  # (b) Dominant species
  p_b <- ggplot() +
    geom_spatvector(data = hex_corr, aes(fill = compo), color = NA) +
    scale_fill_manual(values = species_colors, na.value = "transparent",
                      name = "Dominant species") +
    disturbance_hatch(hex_corr) +
    coord_panel + theme_panel + labs(title = "(b)")
  
  # (c) Initial height (2006)
  p_c <- ggplot() +
    geom_spatvector(data = hex_corr, aes(fill = h0), color = NA) +
    geom_spatvector(data = hex_stands, fill = NA, color = "grey30") +
    scale_fill_viridis_c(name = "Initial top-of-canopy<br>height (TCH<sub>2006</sub>; m)") +
    disturbance_hatch(hex_corr) +
    coord_panel + theme_panel + labs(title = "(c)")
  
  # (d) Growth trend
  p_d <- ggplot() +
    geom_spatvector(data = hex_corr, aes(fill = hslop), color = NA) +
    geom_spatvector(data = hex_stands, fill = NA, color = "grey30") +
    scale_fill_gradient2(low = "#d7191c", mid = "#ffffff", high = "#1a9641",
                         midpoint = 0, name = "Growth trend<br>(m/yr)",
                         limits = c(-1, 1), oob = scales::squish) +
    disturbance_hatch(hex_corr) +
    coord_panel + theme_panel + labs(title = "(d)")
  
  # (e) Deviation from reference trajectory
  p_e <- ggplot() +
    geom_spatvector(data = hex_corr, aes(fill = vertical_resid), color = NA) +
    geom_spatvector(data = hex_stands, fill = NA, color = "grey30") +
    scale_dev() +
    disturbance_hatch(hex_corr) +
    coord_panel + theme_panel + labs(title = "(e)")
  
  # (f) Ridgeline distribution of deviation per species (colored as in (b))
  zoom_ext <- ext(xlim[1], xlim[2], ylim[1], ylim[2])
  ridge_df <- terra::crop(hex_corr, zoom_ext) %>%
    as.data.frame() %>%
    filter(!has_been_disturbed, spike == FALSE) %>%
    group_by(compo) %>% mutate(nplots = n()) %>% ungroup() %>%
    filter(nplots > 25)
  
  p_f <- ggplot(ridge_df, aes(x = vertical_resid, y = compo, fill = compo)) +
    geom_density_ridges(alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_fill_manual(values = species_colors, name = "Composition") +
    labs(x = "Deviation from the reference trajectory (m/yr)", y = "", title = "(f)") +
    coord_cartesian(xlim = c(-0.6, 0.6)) +
    theme_bw() +
    theme(legend.position = "none", aspect.ratio = 0.8 / 4,
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank())
  
  panels  <- Filter(Negate(is.null), list(p_a, p_b, p_c, p_d, p_e, p_f))
  combined <- wrap_plots(panels, ncol = 1) &
    theme(text = element_text(size = 8),
          legend.text = element_text(size = 4.5),
          legend.title = element_markdown(size = 5.5, lineheight = 1),
          legend.key.size = unit(2, "mm"),
          legend.spacing.y = unit(0.5, "mm"))
  
  list(panels = panels, combined = combined, hex_data = hex_data,
       xlim = xlim, ylim = ylim)
}

# Build maps for all sites (reused by Figure 6 and Table 4)
site_maps <- list()
for (lbl in tile_config$site_label) {
  cat("  Processing", lbl, "...\n")
  site_maps[[lbl]] <- create_site_map(lbl, ref_coefs = ref_coefs)
}

# Print / save Figure 5 (Site A)
if (!is.null(site_maps[["Site A"]])) {
  print(site_maps[["Site A"]]$combined)
  ggsave(file.path(path_data, "Figure5_site_maps.png"),
         site_maps[["Site A"]]$combined,
         width = 190, height = 170, units = "mm", dpi = 500)
}


# =============================================================================
# FIGURE 6 — Beech deviation across sites
# =============================================================================
cat("\n=== FIGURE 6: BEECH DEVIATION ACROSS SITES ===\n")

focus_species <- "Beech"

# Panel (a): per-site deviation maps, beech highlighted, other species greyed
beech_panels <- list()
for (lbl in tile_config$site_label) {
  sm <- site_maps[[lbl]]
  if (is.null(sm)) next
  
  hex_corr     <- sm$hex_data %>% filter(status == "corrected")
  hex_focus    <- hex_corr %>% filter(compo == focus_species)
  hex_nonfocus <- hex_corr %>% filter(compo != focus_species)
  hex_stands   <- terra::aggregate(hex_nonfocus, by = "compo")
  
  beech_panels[[lbl]] <- ggplot() +
    geom_spatvector(data = hex_corr, aes(fill = vertical_resid), color = NA) +
    scale_dev() +
    disturbance_hatch(hex_corr) +
    geom_spatvector(data = hex_stands, fill = NA, color = "grey30", linewidth = 0.5) +
    geom_spatvector(data = hex_nonfocus, fill = "grey60", alpha = 0.6, color = NA) +
    coord_sf(xlim = sm$xlim, ylim = sm$ylim, expand = FALSE,
             datum = pull_crs(tile_hex[[lbl]])) +
    scale_y_continuous(position = "right") +
    theme_dark() +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          legend.position = "none",
          axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(title = lbl)
}

# Panel (b): ridgeline of beech deviation per site
df_beech_all <- map_dfr(tile_config$site_label, function(lbl) {
  sm <- site_maps[[lbl]]
  if (is.null(sm)) return(NULL)
  zoom_ext <- ext(sm$xlim[1], sm$xlim[2], sm$ylim[1], sm$ylim[2])
  sm$hex_data %>%
    filter(status == "corrected", compo == focus_species) %>%
    terra::crop(zoom_ext) %>%
    as.data.frame() %>%
    filter(!has_been_disturbed, spike == FALSE, !is.na(vertical_resid)) %>%
    mutate(site = lbl)
})

fig6b <- ggplot(df_beech_all,
                aes(x = vertical_resid,
                    y = factor(site, levels = rev(tile_config$site_label)),
                    fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 1.25, rel_min_height = 0.001, alpha = 0.75) +
  scale_fill_gradient2(low = "#e66101", mid = "#ffffff", high = "#5e3c99",
                       midpoint = 0, name = "Deviation\n(m/yr)",
                       limits = c(-0.5, 0.5), oob = scales::squish) +
  geom_vline(xintercept = 0, linetype = 2, color = "black") +
  scale_y_discrete(position = "right") +
  labs(x = "Deviation from\nreference trajectory (m/yr)", y = "") +
  coord_cartesian(xlim = c(-0.7, 0.7), expand = FALSE) +
  theme_bw() +
  theme(legend.position = "none", aspect.ratio = 3 / 1.8)

if (length(beech_panels) > 0) {
  fig6 <- (wrap_plots(beech_panels, ncol = 1) | fig6b) & theme_article
  print(fig6)
  ggsave(file.path(path_data, "Figure6_beech_deviation.png"),
         fig6, width = 190, height = 95, units = "mm", dpi = 500)
}


# =============================================================================
# TABLE 4 — Moran's I spatial autocorrelation (zoomed windows)
# =============================================================================
cat("\n=== TABLE 4: MORAN'S I ===\n")

fmt_moran <- function(mt) paste0(
  round(mt$estimate["Moran I statistic"], 2), " (p ",
  ifelse(mt$p.value < 0.001, "< 0.001", paste0("= ", round(mt$p.value, 3))), ")")

moran_results <- list()
for (lbl in tile_config$site_label) {
  sm <- site_maps[[lbl]]
  if (is.null(sm)) next
  
  zoom_ext <- ext(sm$xlim[1], sm$xlim[2], sm$ylim[1], sm$ylim[2])
  hex_sf <- sm$hex_data %>%
    filter(status == "corrected", !has_been_disturbed,
           spike == FALSE, !is.na(vertical_resid)) %>%
    terra::crop(zoom_ext) %>%
    st_as_sf()
  
  if (nrow(hex_sf) < 10) { cat(" ", lbl, ": too few observations\n"); next }
  
  coords <- st_coordinates(st_centroid(hex_sf))
  lw <- nb2listw(knn2nb(knearneigh(coords, k = 8)), style = "W")
  
  moran_results[[lbl]] <- tibble(
    Site      = lbl,
    Growth    = fmt_moran(moran.test(hex_sf$hslop, lw)),
    Deviation = fmt_moran(moran.test(hex_sf$vertical_resid, lw))
  )
  cat(" ", lbl, "- Growth:", moran_results[[lbl]]$Growth,
      "| Deviation:", moran_results[[lbl]]$Deviation, "\n")
}

table4 <- bind_rows(moran_results)
cat("\nTable 4:\n"); print(table4)

cat("\n=== ALL RESULTS COMPLETE ===\n")