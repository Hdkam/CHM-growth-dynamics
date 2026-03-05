# CHM-timeseries-harmonization

R code accompanying the paper:

> **de Lame H.**, Favaro A., Lejeune P., Michez A., Plumacker A., Pollet T., Ponette Q., Vander Linden A., Tasseroul M.-P., Messier C. & Bastin J.-F. *Characterizing forest growth dynamics from multi-source time series of canopy height, illustrated in Belgian temperate forests*. Submitted for peer review.

## Overview

This repository contains the four R scripts used to:

1. Build a GAM-based harmonization model that reduces systematic height biases between multi-source canopy height models (CHM) derived from airborne laser scanning (ALS) and digital aerial photogrammetry (DAP);
2. Derive the data required for generating species-specific reference trajectories of vertical growth;
3. Apply the harmonization workflow exhaustively to case study sites using hexagonal grids;
4. Produce all figures and tables presented in the paper.

## Scripts

| Script | Description | Key outputs |
|---|---|---|
| `01_harmonization.R` | Builds GAM standardization models (with/without temporal variables) using forest inventory plots paired with CHM acquisitions | `model_temp_inc.rds`, `model_temp_exc.rds` |
| `02_reference_trajectories.R` | Applies the harmonization model to a random subsample of virtual plots across the study area, to generate species-specific reference trajectories of vertical growth in code 04 | `reference_trajectories.rds` |
| `03_case_study.R` | Applies the harmonization model to case study zones on a hexagonal grid (~1000 m² cells) | `tile_{id}_timeseries.rds`, `tile_{id}_hexgrid.rds` |
| `04_analyses.R` | Produces all tables and figures for the manuscript | Tables 2, 4, S1, S2; Figures 3–6 |

Scripts are designed to be run **sequentially** (01 → 02 → 03 → 04). Each script contains detailed header comments describing inputs, outputs, and usage notes.

## Requirements

- **R** ≥ 4.x
- R packages:

```r
install.packages(c(
  "tidyverse", "sf", "terra", "tidyterra",
  "mgcv", "blockCV", "spdep", "lubridate",
  "patchwork", "ggridges", "ggpattern"
))
```

## Data

- **CHM rasters and ancillary maps** are publicly accessible at 2 m GSD on [Forestimator](https://forestimator.gembloux.ulg.ac.be/).
- **Forest inventory plot data** is accessible upon request from the [IPRFW](http://iprfw.spw.wallonie.be/).

File paths must be adjusted in each script's setup section to match your local folder structure.

## Citation

If you use this code, please cite:

```
de Lame, H., Favaro, A., Lejeune, P., Michez, A., Plumacker, A., Pollet, T., Ponette, Q., Vander Linden, A., Tasseroul, M.-P., Messier, C., & Bastin, J.-F. Characterizing forest growth dynamics from multi-source time series of canopy height, illustrated in Belgian temperate forests. Submitted for peer review.
```

## License

This project is licensed under the MIT License – see [LICENSE](LICENSE) for details.
