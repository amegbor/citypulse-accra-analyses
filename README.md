# CityPulse Accra: analysis code

Analysis code for studies using the **CityPulse Accra survey**, a population-based, cross-sectional survey of urban environments, neighborhood stressors, place experiences, and health among adults in the Greater Accra Metropolitan Area (GAMA), Ghana.

This repository contains the code for the statistical analyses, model diagnostics, tables, and figures. Each analysis has its own folder.

> **Data and data-preparation code are not included.** The scripts start from prepared analytic datasets that contain geocoded residential and mapped locations and health information, and these cannot be shared publicly. See [Data access](#data-access).

---

## The survey

| | |
|---|---|
| Setting | Greater Accra Metropolitan Area, Ghana (22 districts) |
| Design | Two-stage stratified sample drawn from the 2021 Ghana Population and Housing Census frame: 301 enumeration areas selected with probability proportional to size, ~15 households per enumeration area (4,500 households targeted); one adult (≥18 years) per household |
| Fieldwork | August 2024 – May 2025 (95% of interviews in 2024) |
| Sample | 3,157 completed questionnaires; analytic sample of 3,150 adults |
| Mode | Face-to-face, tablet-based questionnaire with a public participation GIS (PPGIS) component: participants recorded their home location and could map favourite and stressful places, with open-text descriptions |
| Measures | Perceived Stress Scale (PSS-10); Patient Health Questionnaire-9 (PHQ-9); anxiety symptoms; blood pressure and hypertension; perceived neighborhood problems and stressors; neighborhood social cohesion; sociodemographic and household characteristics; linked satellite-derived and OpenStreetMap environmental exposures |
| Ethics | Noguchi Memorial Institute for Medical Research Institutional Review Board (NMIMR-IRB CPN 007/24-25); written informed consent obtained electronically |

---

## Repository structure

```
citypulse-accra-analyses/
├── p1_restorative_stressful_places/
├── p2_restorative_resources_ambient_stressors/
├── p3_hypertension_depression_geographies/
├── nemade_cohesion_stress/
├── depression_crime_environment/
├── LICENSE
└── README.md
```

Each folder contains `analysis/` (models, diagnostics, sensitivity analyses) and `tables_figures/` (tables and visualizations).

---

## Analyses

### `p1_restorative_stressful_places/`
**Restorative and stressful places in everyday life**

Compares each participant's mapped favourite (restorative) and stressful places to identify which measured environmental features distinguish them, net of who the person is.
- **Text coding:** open-text place descriptions coded into themes with a rule-based dictionary. Agreement with two independent human coders, and between the coders, assessed with Cohen's κ, using bootstrapped 95% confidence intervals.
- **Within-person design:** conditional logistic regression with each participant as a stratum, comparing environmental features (greenness, tree canopy, land surface temperature, roads and traffic proxies, flood exposure, built form, land use) at favourite versus stressful places.
- **Further analyses:**
  - co-location of favourite and stressful places;
  - distance from home by place type;
  - correspondence between reported stressors (e.g. traffic, dust, flooding) and measured environmental conditions.
- **Robustness:** selection weights for mapping, effect modification by sex and poverty, and sensitivity analyses.

*Scripts to be added.*

### `p2_restorative_resources_ambient_stressors/`
**Restorative resources, ambient stressors and mental health**

Examines how restorative resources and ambient stressors in the residential environment relate to mental health outcomes.
- **Exposures:**
  - favourite and stressful places mapped by *other* residents within 500 m of home (leave-self-out, avoiding same-source bias), adjusted for local sampling density;
  - residential greenness, heat and air pollution.
- **Models:** Bayesian hierarchical models with an SPDE spatial random field and district random effects (R-INLA).
- **Further analyses:** effect modification of stressor associations by restorative resources and individual characteristics, with sensitivity analyses on distance bands and exposure definitions.

*Scripts to be added.*

### `p3_hypertension_depression_geographies/`
**Shared geographies of hypertension and depression**

Examines whether hypertension and depressive symptoms share a common spatial pattern across GAMA.
- **Joint model:** a joint Bayesian model (R-INLA) with a shared SPDE spatial field and outcome-specific fields. Alternative structures are compared by WAIC, DIC and LCPO, and the share of spatial variation that is common is estimated.
- **Further analyses:** whether psychological distress relates to a hypertension diagnosis or to measured blood pressure (normotensive, undiagnosed and diagnosed hypertension).

*Scripts to be added.*

### `nemade_cohesion_stress/`
**Neighborhood problems, social cohesion and perceived stress**

Examines whether perceived neighborhood problems are associated with perceived stress, whether neighborhood social cohesion is associated with lower stress or buffers the association with neighborhood problems, and whether these associations reflect residents' own perceptions or conditions reported by their neighbors.
- **Clusters:** sampling clusters reconstructed from home locations (complete-linkage hierarchical clustering within districts).
- **Neighbors' reports:** leave-self-out means of other residents within 500 m, separating each participant's own perception (within-neighborhood component) from neighbors' reports (between-neighborhood component).
- **Models:** multilevel linear models with cluster random intercepts and district fixed effects (Models A–G), including problems × cohesion interactions and effect modification by household wealth and education.
- **Inference:** cluster-robust (CR2) standard errors with Satterthwaite degrees of freedom, clustered by reconstructed cluster; district clustering as a conservative bound.
- **Diagnostics and sensitivity analyses:** residual and Q–Q plots; ordinal and logistic models, problem subscales, alternative cluster definitions, distrust, and residual spatial autocorrelation (Moran's I).

| Script | Purpose |
|---|---|
| `analysis/01_scs_analysis.R` | Cluster reconstruction, neighbor means, Models A–G, sensitivity analyses |
| `analysis/02_scs_diagnostics.R` | Residual diagnostics for Models C–E |
| `analysis/03_scs_robust_se.R` | Initial cluster-robust standard-error check |
| `analysis/04_scs_analysis_CR2.R` | CR2 cluster-robust inference for all reported estimates |
| `tables_figures/01_scs_figures.R` | Figure functions |
| `tables_figures/02_scs_tables_figures.R` | Tables and figures |

### `depression_crime_environment/`
**Neighborhood crime concern, urban environment and depressive symptoms**

Examines associations of perceived neighborhood stressors and satellite-derived environmental conditions (built form, population density, heat stress, PM2.5, greenness) with depressive symptoms (PHQ-9 ≥ 5; PHQ-9 ≥ 10 as sensitivity).
- **Models:** Bayesian hierarchical logistic models (R-INLA), progressing from single-exposure to domain and fully adjusted specifications, with and without an SPDE spatial random field.
- **Residual spatial field:** mapped with posterior evidence classification.

*Scripts to be added.*

---

## Requirements

- R ≥ 4.6
- Main packages: `lme4`, `ordinal`, `sandwich`, `clubSandwich`, `survival`, `INLA`, `fmesher`, `sf`, `terra`, `spdep`, `ggplot2`, `dplyr`, `tidyr`

## Running the code

1. Set the path to the analytic dataset at the top of the first script in each folder.
2. Run the scripts in `analysis/` in numerical order, then those in `tables_figures/`.

Outputs (tables, figures, model objects) are written to local folders excluded from version control.

---

## Data access

The CityPulse Accra data contain geocoded residential and mapped locations linked to health information. Public release is not permitted under the conditions of ethical approval. De-identified data, with locations removed, are available to qualified researchers on request to the Noguchi Memorial Institute for Medical Research Institutional Review Board, subject to approval and a data-use agreement.

## Citation

If you use this code, please cite the associated article(s) and this repository (Zenodo DOI to be added).

## License

Code is released under the MIT License. See `LICENSE`.

## AI assistance

Generative AI (ChatGPT, OpenAI) was used to adapt the authors’ original analysis code for GitHub, including code organization, relative file paths, formatting, and documentation. The authors reviewed and tested the adapted scripts to verify that they preserved the intended analyses and reproduced the original results. The authors take full responsibility for the code and its outputs.


## Contact

Prince M. Amegbor, Department of Global and Environmental Health, School of Global Public Health, New York University.
