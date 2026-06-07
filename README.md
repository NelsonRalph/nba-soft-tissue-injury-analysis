# 🏀 The Soft Tissue Crisis in the NBA
### A Statistical Analysis of Injury Trends, Play Style Factors, and Structural Interventions

> **Full research paper** | R statistical analysis | 2006–2023 | 30,316 IL records | 647,721 player-game observations

---

## 📋 Table of Contents
- [Overview](#overview)
- [Key Findings](#key-findings)
- [Repository Structure](#repository-structure)
- [What to Upload to GitHub](#what-to-upload-to-github)
- [Data Sources](#data-sources)
- [Requirements](#requirements)
- [How to Run](#how-to-run)
- [Figures](#figures)
- [Methods Summary](#methods-summary)
- [Citation](#citation)
- [License](#license)

---

## Overview

This repository contains the full code, figures, and research paper for a statistical investigation into the dramatic rise of soft tissue injuries in the NBA from 2006 to 2023.

**Soft tissue injuries** (hamstring strains, calf injuries, Achilles tendon pathology, quad strains, groin/hip injuries, ankle sprains, knee ligament injuries) now account for **56% of all NBA injured list placements** — up from 27% in 2006. This project asks: *why*, and *what can be done about it*?

I test five hypotheses:

| # | Hypothesis | Finding |
|---|------------|---------|
| H1 | Soft tissue IL placements have increased significantly | ✅ **+94%** from 2006–11 to 2018–23 (p < 0.001) |
| H2 | Modern play style (3PT rate, pace) drives injury | ✅ **3PT rate r = 0.884** (p < 0.001); pace not significant |
| H3 | Shorter season (60/70 games) reduces injuries | ⚠️ **Confounded** by COVID/lockout schedule compression |
| H4 | 40-minute games (FIBA format) reduce injuries | ✅ **−8.3%** projected via exposure reduction |
| H5 | Capping player minutes (33–35 MPG) reduces injuries | ✅ **−32% to −38%** projected via NB model IRR |

---

## Key Findings

- **94% increase** in soft tissue IL placements per 82-game equivalent (2006–2023)
- **Three-point attempt rate**, not pace, is the dominant play style predictor (r = 0.884)
- **Prior soft tissue injury** is the strongest individual predictor: IRR = **2.31** (doubles re-injury risk)
- Each additional **minute per game** carries a **2.4% increase** in expected injury events (IRR = 1.024)
- The **30–35 MPG tier** has the *highest* injury event rate (49.6%) — higher than players averaging 35+ MPG (survivor-selection effect)
- Negative binomial regression strongly preferred over Poisson: **ΔAIC = 636**, overdispersion θ = 1.232
- All VIF values < 1.2 — **no multicollinearity**

---

## Repository Structure

```
nba-soft-tissue-injuries/
│
├── README.md                          ← You are here
├── LICENSE                            ← MIT License
├── .gitignore                         ← Excludes large/private files
│
├── nba_injury_real_data.R             ← MAIN analysis script (run this)
│
├── data/
│   └── NBA_Player_Injury_Stats_1951_-_2023_.csv   ← Injury dataset (you supply)
│
├── figures/
│   ├── fig1_real_injury_trend.png     ← Stacked area: injury types over time
│   ├── fig2_real_playstyle_injury.png ← 3PT rate & pace vs injury load
│   ├── fig3_real_mpg_injury_rate.png  ← MPG vs injury event rate (GAM)
│   ├── fig4_real_irr_forest.png       ← IRR forest plot (NB model)
│   ├── fig5_real_policy_scenarios.png ← Policy scenario projections
│   ├── diag1_resid_vs_fitted.png      ← Residuals vs fitted
│   ├── diag2_qq.png                   ← Q-Q plot of deviance residuals
│   ├── diag3_scale_location.png       ← Scale-location plot
│   └── diag4_resid_by_season.png      ← Temporal residual check
│
├── paper/
│   └── NBA_Soft_Tissue_Research_Paper_v2.docx  ← Full ACM-format paper
│
└── output/
    └── (model output tables saved here when script runs)
```

---


## Data Sources

This project combines three datasets. The R script downloads #2 and #3 automatically on first run.

| # | Dataset | Source | How to get it |
|---|---------|--------|---------------|
| 1 | **NBA Injured List History 1951–2023** | Kaggle (Logan Launton) | [Download here](https://www.kaggle.com/datasets/loganlauton/nba-injury-stats-1951-2023) — place as `data/NBA_Player_Injury_Stats_1951_-_2023_.csv` |
| 2 | **Player box scores 2005–2023** | hoopR-data (sportsdataverse) | Auto-downloaded by script from [GitHub](https://github.com/sportsdataverse/hoopR-data) |
| 3 | **Team box scores 2005–2023** | hoopR-data (sportsdataverse) | Auto-downloaded by script from [GitHub](https://github.com/sportsdataverse/hoopR-data) |

> **Note on the injury CSV:** The Kaggle dataset is free to download with a Kaggle account. Due to GitHub's 100MB file size limit, it is not included in this repository. The script will error with a clear message if it cannot find the file and tell you exactly where to place it.

---

## Requirements

### R packages
The script installs missing packages automatically on first run. For manual installation:

```r
install.packages(c(
  "ggplot2",   # plotting
  "tidyr",     # data reshaping
  "dplyr",     # data manipulation
  "MASS",      # negative binomial regression (glm.nb)
  "lme4",      # mixed-effects models (glmer)
  "mgcv",      # GAM smoother for Figure 3
  "broom",     # tidy model output
  "httr"       # HTTP requests for data download
))
```

### R version
Developed and tested on **R 4.3.3**. Should work on R 4.0+.

### Internet connection
Required on first run only — to download the hoopR-data RDS files (~18 MB total across 38 files). Subsequent runs use the cached `nba_data/` folder.

---

## How to Run

### Step 1 — Set up your folder
```
my-project-folder/
├── nba_injury_real_data.R
└── data/
    └── NBA_Player_Injury_Stats_1951_-_2023_.csv
```

### Step 2 — Set working directory in RStudio
```
Session → Set Working Directory → To Source File Location
```
Or in the R console:
```r
setwd("path/to/your/folder")
```

### Step 3 — Source the script
```r
source("nba_injury_real_data.R")
```
Or press **Ctrl+Shift+S** (Windows/Linux) / **Cmd+Shift+S** (Mac) in RStudio.

### What happens on first run
1. Script checks for the injury CSV (stops with a clear error if missing)
2. Creates `nba_data/player_box/` and `nba_data/team_box/` folders
3. Downloads 19 seasons × 2 files = 38 RDS files from GitHub (~18 MB, ~30 seconds)
4. Runs all analysis sections (1–14), printing results to console
5. Saves 5 figures to `nba_figures/`

### Subsequent runs
Skips download entirely (files already cached). Full analysis runs in ~2–3 minutes.

---

## Figures

| Figure | Description |
|--------|-------------|
| ![Fig 1](figures/fig1_real_injury_trend.png) | **Figure 1** — Stacked area of soft tissue IL placements by injury type, 2006–2023 |
| ![Fig 2](figures/fig2_real_playstyle_injury.png) | **Figure 2** — 3PT attempt rate & pace vs. injury load (season-level scatter) |
| ![Fig 3](figures/fig3_real_mpg_injury_rate.png) | **Figure 3** — MPG vs. soft tissue IL event rate with GAM smoother |
| ![Fig 4](figures/fig4_real_irr_forest.png) | **Figure 4** — IRR forest plot from negative binomial regression |
| ![Fig 5](figures/fig5_real_policy_scenarios.png) | **Figure 5** — Projected % reduction by policy scenario |

---

## Methods Summary

| Component | Detail |
|-----------|--------|
| **Study period** | NBA seasons 2006–2023 (18 seasons) |
| **Injury records** | 30,316 IL placements (post-2005 reporting standard) |
| **Player-game obs** | 647,721 (hoopR-data via ESPN) |
| **Panel size** | 8,301 player-seasons (≥10 GP, ≥5 MPG); 1,702 unique players |
| **Injury classification** | 9-category keyword taxonomy on IL Notes field |
| **Primary model** | Negative binomial GLM (MASS::glm.nb); AIC = 16,512; θ = 1.232 |
| **Model comparison** | NB vs. Poisson: ΔAIC = 636; overdispersion ratio = 1.245 |
| **Multicollinearity** | All VIF < 1.2 (no concern) |
| **Mixed effects** | Poisson GLMM with player random intercept (lme4::glmer) |
| **Natural experiments** | 2012 lockout (66 games), 2020–21 COVID (72 games) for fatigue exponent |
| **Season normalization** | Per-82-game equivalent using actual games played |

---

## Citation

If you use this code or paper in your own research, please cite:

```bibtex
@article{nba_softtissue_2026,
  title   = {The Soft Tissue Crisis in the NBA: A Statistical Analysis of
             Injury Trends, Play Style Factors, and Structural Interventions},
  author  = {[Authors]},
  year    = {2026},
  note    = {GitHub repository: https://github.com/NelsonRalph/nba-soft-tissue-injury-analysis}
}
```

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

The injury dataset is sourced from Kaggle and subject to its own terms of use. The hoopR-data repository data originates from ESPN and is subject to ESPN's terms of service.

---

## Acknowledgements

- **Logan Launton** — NBA Injured List dataset (Kaggle)
- **Saiem Gilani & Billy Hutchinson** — hoopR-data (sportsdataverse)
- **Jeff Stotts** — InStreetClothes injury tracking methodology
- **ESPN Statistics and Information Group** — underlying box score data

---

*Last updated: June 2026*
