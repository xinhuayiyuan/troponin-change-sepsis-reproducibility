# Troponin change in sepsis — reproducibility package

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Analysis code and SQL for the study:

> **Prognostic Value of Troponin Change Beyond Baseline Troponin and Clinical State in Sepsis:
> A Cohort Study with Multicentre External Validation**

This repository contains the analysis pipeline (R and Python) and the database extraction
SQL required to reproduce every estimate reported in the manuscript and its Supplementary
Information. It contains **no patient-level data**.

> **Version.** The code archive is released as `v1.0.0` (2026-10-13) — the frozen snapshot
> corresponding to the submitted manuscript. Please cite that version rather than the moving
> `main` branch. See [Citation](#7-citation).


---

## 1. Data access

Both source databases are distributed by PhysioNet under a **credentialed data use agreement**
and therefore cannot be redistributed here. Users must obtain access independently:

| Database | Version | Access |
|---|---|---|
| MIMIC-IV | v3.1 | https://physionet.org/content/mimiciv/ (credentialed access + CITI training) |
| eICU Collaborative Research Database | v2.0 | https://physionet.org/content/eicu-crd/ (credentialed access + CITI training) |
| MIMIC-IV-ECHO (for LVEF) | — | https://physionet.org/content/mimic-iv-echo/ |

After access is granted, load both databases into a local PostgreSQL instance:

* MIMIC-IV  -> database `mimic`, schemas `mimiciv_hosp` / `mimiciv_icu`
* eICU-CRD  -> database `eicu`,  schema `eicu_crd`

Database connection parameters are read from environment variables, so no credentials are
stored in the code:

| Variable | Default | Purpose |
|---|---|---|
| `MIMIC_HOST` | `localhost` | PostgreSQL host for MIMIC-IV |
| `MIMIC_PORT` | `5432` | PostgreSQL port |
| `MIMIC_DB` | `mimic` | MIMIC-IV database name |
| `MIMIC_USER` | `postgres` | database user |
| `MIMIC_PASSWORD` | *(empty)* | database password |

Example:

```bash
export MIMIC_HOST=localhost MIMIC_PORT=5432 MIMIC_DB=mimic
export MIMIC_USER=<user> MIMIC_PASSWORD=<password>
export PROJECT_ROOT=/path/to/reproducibility_package
```

The eICU extraction SQL is run with any PostgreSQL client against the `eicu` database
(see the header of `sql/eicu_extract_01_cohort.sql`).

---

## 2. Repository layout

```
reproducibility_package/
├── README.md
├── LICENSE                          # MIT (code only)
├── CITATION.cff                     # how to cite this repository
├── AI_DISCLOSURE.md                 # disclosure of AI-tool use
├── requirements.txt                 # Python dependencies
├── .gitignore                       # keeps data, outputs and figures out of the repo
├── .gitattributes                   # LF line endings
├── meta/
│   ├── variable_definitions.md      # derived-variable definitions used throughout
│   └── r_environment.md             # R version, packages, seeds
├── code/
│   ├── analysis_01..30              # extraction, cohort construction, modelling
│   ├── audit_*.R / *.py             # independent numerical audits
│   ├── gen_*                        # figure and table generation
│   └── build_*.R / *.py             # result-locking helpers
├── sql/
│   ├── eicu_extract_01_cohort.sql   # eICU cohort + site/severity extraction
│   └── eicu_extract_02_export.sql   # eICU export helper
│
├── data/                            # NOT tracked: user-supplied input tables
├── output/                          # created at run time: intermediate results
│   └── lock/                        # created by build_01: frozen result tables
└── figures/                         # created at run time: figures at 300 dpi
```

`data/`, `output/` and `figures/` do not exist in the repository. The first two must be
populated by the pipeline (steps 0a-0c below); `output/` and `figures/` are created automatically.
All three are listed in `.gitignore` so that no patient-level data is ever committed.


All scripts resolve the project root through the environment variable `PROJECT_ROOT`
(default: the working directory). Set it once before running, e.g.

```bash
export PROJECT_ROOT=/path/to/reproducibility_package      # Linux / macOS
setx  PROJECT_ROOT  <path>\to\reproducibility_package  # Windows (new shell)
```

---

## 3. Software environment

**R 4.6.0** with:

| Package | Use |
|---|---|
| survival 3.8-6 | Cox models, Kaplan-Meier, restricted cubic splines, Schoenfeld diagnostics |
| cmprsk 2.2-12 | Fine-Gray competing-risks models |
| mice | multiple imputation (m = 20, chained equations, Rubin pooling) |
| riskRegression, pec | time-dependent (IPCW) AUC |
| survRM2 | restricted mean survival time |
| boot | bootstrap internal validation of the C-index |
| EValue | E-value / tipping-point analysis |
| rms, Hmisc | spline terms and validation utilities |

**Python 3** with: `pandas`, `numpy`, `psycopg2`, `matplotlib`, `Pillow`, `python-docx`.

## 4. Pipeline

Scripts are numbered in execution order and are run from the project root.

> **Entry point.** The pipeline starts from the derived analysis tables
> `v4_cohort_2plus_clean.csv` (MIMIC-IV cohort with dynamic-group assignment) and, for
> eICU, the CSVs produced by `sql/eicu_extract_01_cohort.sql`.
> The definitions needed to rebuild `v4_cohort_2plus_clean.csv` from the raw MIMIC-IV
> tables (Sepsis-3 identification, first ICU admission, age and length-of-stay limits,
> troponin windowing and group assignment) are given in `meta/variable_definitions.md`.
> The first executable script is `analysis_09_event_time_fix.py`, which derives the
> survival-time variables from the recorded dates of death.

### 4.0 Data preparation

| Step | Item | Produces |
|---|---|---|
| 0a | Rebuild the base cohort tables from MIMIC-IV following `meta/variable_definitions.md` | `temp_v4_cohort`, `temp_v4_troponin`, `temp_v4_trop_dynamic`, `temp_v4_clinical`, `temp_v4_comorb` |
| 0b | Export the base cohort to CSV | `data/v4_cohort_2plus_clean.csv` |
| 0c | Run `sql/eicu_extract_01_cohort.sql` (parts 1-5) against the eICU database | `data/eicu_extract_*.csv` |

### 4.1 Cohort construction (MIMIC-IV)

| Step | Script | Purpose |
|---|---|---|
| 1 | `analysis_09_event_time_fix.py` | derive `time_28d_true` / `time_90d_true` from the recorded dates of death -> `data/v4_cohort_2plus_true_time.csv` |
| 2 | `analysis_14_mimic_window_extract.py` | per-stay troponin and covariate window extraction |
| 3 | `analysis_16_trajectory_extract.py` | per-stay, per-window median cTnT for the trajectory figure |
| 4 | `analysis_17_landmark_extract.py` | 72-hour landmark variables |
| 5 | `analysis_18_landmark_dataset.py` | assembly of the landmark analysis dataset |
| 6 | `analysis_27_alt_cohort_extract.py` | alternative (single-measurement-inclusive) cohort variables |
| 7 | `build_01_lock_results.R` | frozen analysis cohort; Tables 1 and 2 and the cohort-flow table |

### 4.2 Primary and sensitivity models

| Step | Script | Purpose |
|---|---|---|
| 8 | `analysis_19_landmark_models.R` | primary 72-hour landmark continuous-rate model (M3a); attenuation ladder; splines; PH diagnostics; RMST |
| 9 | `analysis_20_incremental_value.R` | change in C-index and likelihood-ratio tests |
| 10 | `analysis_21_evalue_bias.R` | E-value, inverse-probability weighting, regression dilution |
| 11 | `analysis_22_legacy_recheck.R` | categorical models (Models 1-3) |
| 12 | `analysis_15_window_sensitivity.R` | measurement-window sensitivity analyses |
| 13 | `analysis_28_sensitivity_altcohort.R` | widened-cohort and negative-control analyses |
| 14 | `analysis_29_sensitivity_core.R` | competing risks, IPCW AUC, multiple imputation, exact landmark |
| 15 | `analysis_30_validation_and_bias.R` | quantitative bias audits |

### 4.3 External validation (eICU-CRD)

| Step | Script | Purpose |
|---|---|---|
| 16 | `analysis_26_eicu_site_gradient.R` | eICU models, measurement-intensity tertiles, hospital-clustered standard errors |

### 4.4 Independent audits

`audit_power_mde.R`, `audit_table1_recompute.py`, `audit_landmark_logrank.py`,
`audit_m3a_baseline_definition.R`, `audit_ST5_excluded.py` recompute reported quantities
independently and verify that every number in the manuscript matches the model output.

### 4.5 Figures and tables

`gen_09_figure_flow.py`, `gen_10_figure_km.R`, `gen_11_figure_forest.py`,
`gen_12_figure_composite.py`, `gen_13_figures_supp.R`, `gen_14_figure_trajectories.py`,
`gen_15_figure_landmark_forest.py`, `gen_18_fig2_km_m3.R`, `gen_19_figure_forest8.py`,
`gen_20_compose_figures.py` and `gen_22_figure_eicu_forest.py` regenerate all main and
supplementary figures at 300 dpi.

Intermediate outputs are written to `output/`, frozen result tables to `output/lock/`, and
figures to `figures/` (all created automatically).


## 5. Key analytical definitions

* **Exposure window**: the first 72 hours after ICU admission.
* **Baseline troponin**: minimum cTnT within the first 6-hour window.
* **Dynamic groups** (assigned in the order listed):
  * *Rising*: peak >= 1.5x baseline **and** peak > 0.03 ng/mL
  * *Falling*: last <= 0.67x baseline **and** baseline > 0.03 ng/mL
  * *Stable-low*: baseline <= 0.03 ng/mL **and** peak <= 0.03 ng/mL
  * *Stable-high*: baseline > 0.03 ng/mL, peak <= 1.5x baseline, last >= 0.67x baseline
* **Rate of change**: r = ln(last / first) x (72 / hours between first and last measurement),
  winsorised at the 1st and 99th percentiles (r_w).
* **Outcomes**: 28-day all-cause mortality (primary); 90-day mortality (secondary).
* **Primary model (M3a)**: fitted in the 72-hour landmark cohort (alive and in hospital at
  72 hours); Cox model adjusted for admission covariates, 72-hour clinical state,
  log(first troponin + 0.001) and measurement count.
* **Prespecified minimum clinically important effect**: hazard ratio 1.20 per SD of r_w.
* **eICU validation**: cardiac troponin I; ICD-coded sepsis; in-hospital death within
  28 days of the landmark; standard errors clustered by hospital.

Further detail is given in `meta/variable_definitions.md`.

---

## 6. Reproducibility notes

* Scripts are deterministic given the same database snapshot and package versions.
* Survival times are derived from the recorded date of death. Calendar dates in MIMIC-IV
  are shifted per patient as part of de-identification, which does not affect intervals.
* No script contains an absolute path; the project root is taken from `PROJECT_ROOT`.
* Reported estimates were frozen in the analysis locks of the submission package; the
  audit scripts re-derive each number from the raw database and the fitted models.

---

## 7. Citation

If you use this code, please cite the manuscript above together with the source databases:

* Johnson AEW, Bulgarelli L, Shen L, et al. MIMIC-IV, a freely accessible electronic
  health record dataset. *Sci Data*. 2023;10(1):1.
* Pollard TJ, Johnson AEW, Raffa JD, et al. The eICU Collaborative Research Database,
  a freely available multi-center database for critical care research. *Sci Data*. 2018;5:180178.

### Citing this repository

Cite the **versioned release**, not the branch. The metadata in `CITATION.cff`
(`version: 1.0.0`, `date-released: 2026-10-13`) matches the Git tag `v1.0.0` and the
archived release:

```bibtex
@software{han_troponin_change_sepsis_2026,
  author  = {Han, Zhongshu},
  title   = {Troponin change in sepsis --- reproducibility package},
  version = {1.0.0},
  year    = {2026},
  url     = {https://github.com/xinhuayiyuan/troponin-change-sepsis-reproducibility}
}
```

To create the release that this citation points to:

```bash
git tag -a v1.0.0 -m "Release v1.0.0 - submitted manuscript"
git push origin v1.0.0      # then publish a GitHub Release from this tag
```

The tag keeps a frozen snapshot of the code exactly as it was when the manuscript was
submitted; cite `v1.0.0` rather than the moving `main` branch.

## 7.1 AI-tool disclosure

`AI_DISCLOSURE.md` records which AI tool was used and, more importantly, which parts of the
work it did and did not touch. It is provided for transparency under ICMJE and Nature
Portfolio policy. The disclosure covers consultation on statistical methodology and
assistance with analysis code only; no AI tool had access to the data or wrote any part of
the manuscript.


---

## 8. Licence

Code is released under the MIT Licence (see `LICENSE`). The source databases remain
subject to the PhysioNet Credentialed Health Data Licence and Data Use Agreement
and are **not** redistributed here.


