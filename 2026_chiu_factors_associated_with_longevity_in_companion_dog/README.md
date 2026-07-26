ABOUT THIS PUBLICATION

Year: 2026

Authors: Weihsueh A. Chiu, Matthew D. Dunbar, Gang Han, Benjamin R Harrison, Amber J. Keyser, Paul E. Litwin, Dog Aging Project Consortium, Daniel Promislow, Kate E. Creevy

Title: Factors Associated with Longevity in Companion Dogs: Initial Findings from the Dog Aging Project

Journal: Geoscience

Citation: 

DOI link: 

ABOUT THIS REPO: This folder contains datasets and/or code used for the analyses reported in the article cited above in fulfillment of the requirements of the publisher and for the sole purpose of replicating the analyses contained in the article. These data represent currently embargoed Dog Aging Project (DAP) data or data otherwise unavailable in a DAP Curated Data Release.

ABOUT DAP DATA: Dog Aging Project data consist of an extensive set of anonymized variables collected from tens of thousands of dogs. Data types include survey data, environmental data, and biospecimen lab results, among other data types. Curated Data are released annually via the Terra platform to researchers worldwide for use in scientific research, curricular, and other nonprofit uses.

**Explore DAP data at** [https://data.dogagingproject.org/](https://data.dogagingproject.org/)

**Apply for DAP data access** [https://dogagingproject.org/data-access](https://dogagingproject.org/data-access)

**IMPORTANT**: By accessing the datasets in this repo, you are agreeing to abide by the terms of the DAP Data Use policy, including:

* The Dog Aging Project retains all rights of ownership to DAP data.
* DAP data are for the sole purposes of research, education, and investigative journalism. 
* DAP data will not be used to train generative AI. 
* DAP data will not be misused or misrepresented to make false or unvalidated claims or to develop commercial products.
* DAP data must not be shared with third parties.
* DAP data must not be used to try to identify or contact any research participants.
* Use of DAP data does not imply collaboration with the Dog Aging Project.
* Users must acknowledge DAP as the source of the data and cite DAP Curated Data in the appropriate format. 

# 2026_chiu_factors_associated_with_longevity_in_companion_dog

The full pipeline is orchestrated by **`00.RunAll.R`** in `code/`, 
which sources the scripts below in order. Each stage runs from a clean workspace (`rm(list=ls())` between calls). Inputs live in `data/`, 
derived outputs in `code/results/`, and plots in `code/figures/`.

## Requirements

- R 4.5.1 (as reported in the manuscript)
- CRAN packages: `survival` (v3.8-3), `flexsurv`, `survminer`, `lubridate`,
  `dplyr`, `tidyverse`, `ggpubr`, `viridisLite`, `stringr`, `haven`,
  `usdata`, `wCorr`, `emmeans`

## Input data (in `data/`)

- `SurvivalAnalysisCuratedDogs_thru_2024_ran_2025_12_16.csv` — DAP curated
  survival release with mortality follow-up through 2024
- `DAP_2024_DogOverview_v1.0.csv` — dog-level demographic overview
- `DAP_2024_CODEBOOK_v1.0.csv` — HLES codebook used to drive the SWAS
- HLES data files loaded by the SWAS scripts
- State-level human mortality / life expectancy files (HDPulse) used by
  `4.2.Geography_DAP_vs_Human.R`

## Pipeline (as sourced by `00.RunAll.R` in `code\`)

| # | Script | Purpose | Manuscript element |
|---|---|---|---|
| 1 | `0.process_DAP_datafiles.R` | Load and clean the curated release; parse dates; compute `first.age`, `last.age`, `event`; factorize Size × Breed_Class × Sex; save `SurvivalData.RData`. | Cohort construction (N = 41,047) |
| 2 | `1.Cohort_Descriptive_Stats.R` | Crude mortality rates and entry/follow-up/death-age quantiles by Breed_Class × Size × Sex strata; prop tests and stratum ANOVA. | **Table 1**; supplemental AOV tables |
| 3 | `2.1.Survival_DAP_demographics.R` | Primary Cox PH: `Surv(first.age, last.age, event) ~ Size + Breed_Class + Sex`; K-M by 20 strata; median/IQR lifespans; stratified size/breed/sex effect p-values. | **Figure 1** (K-M + forest); **Table 2** |
| 4 | `2.2.Survival_DAP_demographics.Alt.R` | Sensitivity analysis — alternative parameterization (10 kg weight bins). | Supplemental sensitivity results |
| 5 | `2.3.Survival_DAP_demographics.Alt2.R` | Sensitivity analysis — alternative reference strata / model form. | Supplemental sensitivity results |
| 6 | `3.Survival_CommonBreeds.R` | Repeats the demographic survival analysis on the 16 most common single AKC breeds (n > 250). | **Table 3**; breed-specific K-M plots |
| 7 | `4.1.Geography_DAP.R` | Cox stratified by Size × Breed × Sex, adjusting for rural/suburban/urban and state; `emmeans` with `method = "eff"` to obtain each state's deviation from the count-weighted grand mean (no reference state pinned to zero). Saves `Geoeffect-Cox-results.Rdata`. | Geographic results text |
| 8 | `4.2.Geography_DAP_vs_Human.R` | Joins state-level dog HRs to human age-adjusted mortality and life expectancy (HDPulse); inverse-variance-weighted regression plus `wCorr` weighted Pearson/Spearman. | **Figure 2** (r = 0.43, ρ = 0.63, p = 0.0017) |
| 9 | `5.1.HLES_SWAS_analysis.R` | Survey-wide association study over the 804 HLES codebook variables (grouped by module: `dd, oc, pa, de, db, df, dt, mp, hs`, …). Univariate Cox stratified by Size × Breed × Sex; modal level as reference for categorical variables; Benjamini–Yekutieli FDR control. | Models underlying Figs 3–4 |
| 10 | `5.2.0.HLES_SWAS_figures.R` | Manhattan-style −log10(q) plot by module and forest plots of selected hits. | **Figure 3**, **Figure 4**, **Table 4** |
| 11 | `5.2.1.HLES_SWAS_figures_SensAn_MatureAdult.R` | SWAS sensitivity analysis restricted to the "Mature Adult" lifestage (n = 20,505). | Supplemental SWAS sensitivity |
| 12 | `5.2.2.HLES_SWAS_figures_SensAn_CommonBreeds.R` | SWAS sensitivity analysis restricted to the 16 most common single breeds (reduced genetic heterogeneity). | Supplemental SWAS sensitivity |

### Helper

- `ggforest2.R` — custom `survminer::ggforest` variant, `source()`-ed by
  `4.1.Geography_DAP.R`, `5.1.HLES_SWAS_analysis.R`, and the SWAS figure
  scripts. Required by the pipeline; not sourced directly from `00.RunAll.R`.

## Outputs

- `code/results/` — every CSV/RData artifact named after the script that produced
  it, so any manuscript number can be traced back to a specific script
  (e.g. `Survival_DAP_demographics.SumStats.csv`,
  `Geoeffect-Cox-Human.MR.csv`, `Supp.HLES_Cox_signif_results.csv`).
- `code/figures/` — PDF/JPG panels used in the main and supplemental figures.

## To reproduce

```r
setwd("<repo root>/code")
source("00.RunAll.R")
```

