# TPHP

This repository contains the core analysis code for the TPHP project, which supports large-scale, spatially resolved human proteome profiling. The dataset quantifies >13,000 proteins across 2856 samples spanning 58 major tissue types (251 tissue subtypes) and 25 cancer types using DIA-MS.

The study is described in the preprint:
**bioRxiv (2025)** — https://www.biorxiv.org/content/10.1101/2025.02.14.638212v1

The code performs per-cancer, per-protein tumor versus paired non-tumor comparisons using linear mixed-effects regression, enabling systematic analysis of oncogenic proteome changes across tissues. Results are exported as RDS objects for downstream statistical analysis, visualization and integration with the public proteome database.

## Repository layout

- `tumor_dysregulation_analysis.R` — main analysis script
- `data/tumor_compare_data.parquet` — input table (paired tumor/non-tumor)
- `output/compare_report_output.rds` — main result
- `output/package_versions.txt` — R and package versions used


# 1. System requirements

## 1.1 Operating systems

The workflow is expected to run on every common desktop/server platforms that support R. The test environment is on Windows 11 x64 system.

## 1.2 Software dependencies

### Required

- R (recommended R 4.0+)
- CRAN packages
  - `tidyverse`
  - `arrow`
  - `lme4`
  - `lmerTest`
  - `plyr`

Exact tested version recorded in `output/package_versions.txt`

## 1.3 Hardware

- CPU-only execution (no GPU required)
- Recommended RAM: ≥8 GB (larger datasets may require more)
- No required non-standard hardware


# 2. Installation guide

## 2.1 Get the code

```bash
git clone <REPO_URL>
cd <REPO_DIR>
```

## 2.2 Install R

Install R for the operating system from [CRAN](https://cran.r-project.org/bin/windows/base/).

Optionally, install packages manually:

```r
install.packages(c("tidyverse","arrow","lmerTest","lme4"))
```

This project is script-based, and no software installation step is required beyond installing dependencies.

# 3. Demo

## Quick start

From the repository root:

```bash
Rscript tumor_dysregulation_analysis.R
```

## 3.1 Input files

Default input path:

- `data/tumor_compare_data.parquet`

## 3.2 Expected output files

After a successful run, the script writes:

- `output/compare_report_output.rds`

  An R list `DEA` with:

  - `DEA$Diff.report`: full per-(cancer, protein) model results
  - `DEA$Diff.report.filter`: filtered results using `Hedges'g >= 0.5` and `p_adj_BH < 0.05`

- `output/package_versions.txt` 
  R version and package versions used.

## 3.3 Expected runtime

Runtime scales primarily with the number of model fits, i.e.,
 *cancer types × proteins × samples*.

As a rough guide, on a standard desktop CPU (8–16 threads, 16–32 GB RAM), throughput is ~25 model fits per second, corresponding to **< 10 minutes per cancer type** under typical settings.


# 4. Instructions

## 4.1 Input data format

The script expects a Parquet table with metadata columns:

- `patient_ID` (string)
- `sample_type` (must include `NT` for non-tumor and `T` for tumor)
- `cancer_abbr` (string)
- `cancer_subtype` (string; may be constant within a subset)
- `Gender` (categorical)
- `Age` (integer-like)
- `Dataset` (categorical)

All remaining columns are treated as numeric features (proteins).


## 4.2 What the script does

For each `cancer_abbr` and each `protein`, the script fits a mixed-effects model:

- `value ~ 1 + sample_type + (optional covariates) + (1 | patient_ID)`

Optional covariates among `{cancer_subtype, Gender, Age_c, Dataset}` are included only when estimable within the cancer/protein subset.

The reported tumor effect is the fixed-effect coefficient for `sample_typeT`.

## 4.3 Output columns

Each row in `DEA$Diff.report` corresponds to one (cancer, protein) fit and includes:

- `effect`, `se`, `t`
- `df`(degree of freedom)
- `p`, `p_adj_BH` (p value)
- `sigma` (residual SD)
- `es_adj = effect / sigma` (standardized effect)
- `g_adj = J(df) * es_adj` (small-sample adjusted standardized effect)
- `formula` (the exact model formula used)
- `is_singular`, `re_var_patient` (fit diagnostics)



