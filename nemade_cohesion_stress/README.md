# Neighborhood problems, social cohesion and perceived stress in Accra, Ghana

Analysis code for a multilevel study of perceived neighborhood problems, neighborhood social cohesion and
perceived stress among adults in the Greater Accra Metropolitan Area (GAMA), Ghana.

## Survey

The data come from the CityPulse Accra survey, a cross-sectional, population-based survey of environmental
and psychosocial stressors and health in GAMA conducted between August 2024 and May 2025. The sampling frame
was based on the 2021 Ghana Population and Housing Census. Using a two-stage stratified design, 300
enumeration areas were selected across the districts of GAMA with probability proportional to population
size, and 15 households were to be randomly selected in each (4,500 households in total). One adult aged
18 years or older was interviewed face to face in each household, and the location of each home was
recorded. Completed questionnaires were obtained from 3,157 adults; the analytic sample comprises 3,150
adults living in 22 districts.

Measures used here:

- **Perceived stress**: 10-item Perceived Stress Scale (PSS-10), standardized total score.
- **Perceived neighborhood problems**: mean of 14 items (0-4), with crime and disorder, environmental and
  socioeconomic subscales.
- **Social cohesion**: mean of four trust and collective-action items (1-5); two negatively worded trust
  items form a separate distrust measure.
- **Covariates**: age, sex, education, household wealth quintile, marital status, employment status, years
  at the current dwelling, ethnicity and district.

Enumeration-area identifiers are not available, so sampling clusters are reconstructed from home
coordinates (complete-linkage clustering within district). Models are linear mixed models with a
reconstructed-cluster random intercept and district fixed effects. Reported inference uses CR2
cluster-robust standard errors clustered by reconstructed cluster (Satterthwaite degrees of freedom), with
standard errors clustered by district as a conservative bound.

## Files

| File | Purpose |
|---|---|
| `config.R` | Paths: `data_file` (analytic data set) and `out_dir` (outputs), relative to the project root |
| `01_main_analysis.R` | Scales and reliability, reconstructed clusters, leave-self-out neighbor means, models A-G (model-based inference), model-based sensitivity analyses, descriptives |
| `02_diagnostics.R` | Residual diagnostics for models C-E (residual vs fitted, Q-Q plots, summary statistics) |
| `03_robust_inference.R` | CR2 cluster-robust inference for all reported terms, simple slopes, group-specific slopes, slopes by wealth quintile, predicted values and sensitivity analyses |
| `04_tables_figures.R` | Table 1, Tables 2-3, S1-S2 Tables (.docx) and Figs 1-2 (.tif, with PNG previews) |
| `sessionInfo.txt` | R session information from the run that produced the published results |

## Running the code

1. Place the analytic data set at `data/analytic_data.rds` (or edit `config.R`).
2. From the project root, run the scripts in order:

```r
source("01_main_analysis.R")
source("02_diagnostics.R")
source("03_robust_inference.R")
source("04_tables_figures.R")
```

`02_diagnostics.R` and `03_robust_inference.R` reuse the data preparation and model fits of
`01_main_analysis.R` (they read it with `readLines()`), and both check that their fixed effects match
`outputs/tables/models.csv`, so `01_main_analysis.R` must be run first. Outputs are written to
`outputs/tables/` and `outputs/figures/`.

Required R packages: lme4, performance, psych, lavaan, ordinal, sandwich, lmtest, spdep, sf, ggplot2,
clubSandwich, ragg, flextable and officer. The results were produced with R 4.6.2; package versions are in
`sessionInfo.txt`.

The analytic data set is a data frame with one row per participant, including `Respondent_ID`, `DISTRICT`,
`LONG`, `LAT`, the PSS-10 items (`PSS_1`-`PSS_10`, reverse-coded `revPSS_4`, `revPSS_5`, `revPSS_7`,
`revPSS_8`) and `PSS_score`, the 14 neighborhood problem items (`PN_Noise`, `PN_Air`, `PN_Unemp`,
`PN_Litter`, `PN_Robbery`, `PN_Drugs`, `PN_Pov`, `PN_SexViol`, `PN_Vandalism`, `PN_Housing`, `PN_Disaster`,
`PN_Homeless`, `PN_Prost`, `PN_Safety`), the trust and collective-action items (`PN_TS1`-`PN_TS4`,
`PN_CAC1`, `PN_CAC2`) and the covariates `age`, `female`, `educ_f`, `wealth_quintile`, `Marital`, `Emp`,
`YrsInHome` and `Ethn`.

## Data access

The participant-level data are not included in this repository. They contain the locations of
participants' homes and cannot be shared publicly without risking the identification of participants.
Researchers may request access from the CityPulse Accra study team (contact the corresponding author of
the associated article); access is subject to a data-use agreement and approval by the Noguchi Memorial
Institute for Medical Research Institutional Review Board, which approved the study. The code writes only
aggregate tables and figures, except for a participant ID to cluster lookup in `outputs/derived/`, which
should not be shared.

## License

MIT; see `LICENSE`.
