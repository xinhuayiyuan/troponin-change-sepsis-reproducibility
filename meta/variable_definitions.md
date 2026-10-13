# Derived variable definitions

All variables are derived from the source tables of MIMIC-IV v3.1 (schemas `mimiciv_hosp`,
`mimiciv_icu`) and the eICU Collaborative Research Database v2.0 (schema `eicu_crd`).
Item identifiers refer to MIMIC-IV unless stated otherwise.

## 1. Cohort

| Variable | Definition |
|---|---|
| `stay_id` | ICU stay identifier (`mimiciv_icu.icustays`) |
| Sepsis-3 | suspected infection plus a SOFA increase >= 2 points, evaluated on the first ICU day |
| First ICU admission | `icustays` record with the lowest `intime` per `subject_id` |
| Age | age at ICU admission; values > 89 years are shifted in MIMIC-IV, hence the 18-90 inclusion range |
| `los_icu` | ICU length of stay in days; inclusion required > 1 day |

## 2. Troponin

| Variable | Definition |
|---|---|
| `ctnt` | cardiac troponin T, `labevents` itemid 51003 (assay generation unspecified; reported as cTnT, not hs-cTnT) |
| Exposure window | first 72 hours from ICU admission |
| `baseline_troponin` | minimum cTnT within the first 6-hour window |
| `peak_troponin` | maximum cTnT within 72 hours |
| `last_troponin` | last cTnT value within 72 hours |
| `delta_troponin` | `peak_troponin` - `baseline_troponin` |
| `n_meas` | number of cTnT measurements within 72 hours |
| `dt_hours` | hours between the first and last measurement within 72 hours |
| `r_w` | rate of change, `ln(last / first) x (72 / dt_hours)`, winsorised at the 1st and 99th percentiles |

### Dynamic groups (precedence order)

| Group | Rule |
|---|---|
| Rising | peak >= 1.5 x baseline AND peak > 0.03 ng/mL |
| Falling | last <= 0.67 x baseline AND baseline > 0.03 ng/mL |
| Stable-low | baseline <= 0.03 ng/mL AND peak <= 0.03 ng/mL |
| Stable-high | baseline > 0.03 ng/mL AND peak <= 1.5 x baseline AND last >= 0.67 x baseline |

Patients satisfying none of the rules are classified as indeterminate and excluded.

## 3. Covariates

| Variable | Definition |
|---|---|
| `age`, `gender` | as recorded at ICU admission |
| `sofa_24h` | SOFA score over the first 24 hours |
| `invasive_vent` | invasive mechanical ventilation, ascertained over the ICU stay |
| `norepinephrine` | norepinephrine infusion, ascertained over the ICU stay |
| `creatinine_max` | maximum serum creatinine in the first 24 hours |
| `lactate_max` | maximum lactate in the first 24 hours |
| `myocardial_infarct` | any coded myocardial infarction (`I21*`, `I22*`, `I252*` in `diagnoses_icd`, `icd_version = 10`) |
| `congestive_heart_failure` | any coded heart failure (`I50*`) |
| `min_lvef_72h` | nadir left ventricular ejection fraction over the index ICU stay (MIMIC-IV-ECHO); values > 100% treated as missing |
| `hadm_id`, `subject_id` | hospital and patient identifiers |

## 4. Outcomes and time variables

| Variable | Definition |
|---|---|
| `mortality_28`, `mortality_90` | all-cause death within 28 / 90 days of ICU admission |
| `time_28d_true`, `time_90d_true` | hours from ICU admission to death, or to the administrative censoring point (672 h / 2,160 h) |
| 72-hour landmark | patient alive and still in hospital 72 hours after ICU admission |
| `hours` (eICU) | ICU length of stay in hours; `unitdischargeoffset` / 60 |

## 5. eICU-CRD specifics

| Variable | Definition |
|---|---|
| `ctni` | cardiac troponin I (`labresult` + `labresultoffset` within 0-72 h) |
| Sepsis | ICD-9/10 codes (038, 995.91, 995.92, 785.52, A40, A41, R65.2) or diagnosis-string keywords |
| Inclusion | age 18-90; ICU stay > 24 h; still in hospital at 72 h (`hospitaldischargeoffset` >= 4,320 min) |
| Outcome | in-hospital death within 28 days of the landmark |
| Site | hospital identifier; sites contributing < 20 eligible stays excluded |
| Covariates | APACHE IV predicted mortality; any vasopressor in the 48-72 h window; 72-hour lactate |
