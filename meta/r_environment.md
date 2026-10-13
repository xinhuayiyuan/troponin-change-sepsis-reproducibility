# R environment

**R version 4.6.0 (2026-04-24)**, `x86_64-w64-mingw32`, Windows 11 x64.

| Package | Version | Used for |
|---|---|---|
| `survival` | 3.8-6 | Cox proportional-hazards models, Kaplan-Meier curves, restricted cubic splines, Schoenfeld residuals |
| `cmprsk` | 2.2-12 | Fine-Gray competing-risks models |
| `mice` | see note | multiple imputation (m = 20, chained equations, Rubin pooling) |
| `riskRegression`, `pec` | see note | time-dependent (IPCW) AUC |
| `survRM2` | see note | restricted mean survival time |
| `boot` | see note | bootstrap internal validation of the C-index |
| `EValue` | see note | E-value and tipping-point analysis |
| `rms`, `Hmisc` | see note | spline terms and validation utilities |

**Note on versions.** The archived session information of the primary analysis recorded
explicit versions for `survival` and `cmprsk` (given above). The remaining packages are
recorded by name only; to reproduce the exact fit, install the versions that were current
at the time of analysis or the latest release compatible with R 4.6.0.

**Random seeds.** Every resampling procedure (multiple imputation, bootstrap, permutation
tests) sets a fixed seed inside the script, so results are deterministic given the same
input data and package versions.

**Database access.** Connection parameters are read from environment variables
(`MIMIC_HOST`, `MIMIC_PORT`, `MIMIC_DB`, `MIMIC_USER`, `MIMIC_PASSWORD`); see README
section 1.
