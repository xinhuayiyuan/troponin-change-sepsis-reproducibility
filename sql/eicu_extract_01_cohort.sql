-- =============================================================================
-- eICU-CRD external validation: data-extraction SQL
-- Created: 2026-10-04
-- Purpose: extraction for the eICU external validation and site-gradient
--          (measurement-intensity) analyses reported in the manuscript.
-- Target:  PostgreSQL database "eicu", schema "eicu_crd" (populated from the
--          PhysioNet eICU Collaborative Research Database v2.0 distribution).
-- Scope:   read-only SELECT statements only; no INSERT/UPDATE/DELETE/CREATE.
-- Time convention: all eICU offsets are minutes relative to ICU admission;
--          labresultoffset may be negative (pre-admission sampling), so the
--          extraction window is [0, 4320] minutes from ICU admission.
-- =============================================================================
-- Run each part separately; each part writes one CSV (see the part headers).
-- Export settings: CSV / UTF-8 / column headers included / comma separated.
-- =============================================================================


-- =============================================================================
-- PART 0 - schema self-check (optional; run first if PART 1 reports a missing column)
-- =============================================================================

-- 0.a column names of the tables used below
SELECT table_name, column_name, data_type
FROM information_schema.columns
WHERE table_schema = 'eicu_crd'
  AND table_name IN ('diagnosis', 'apachepatientresult', 'hospital', 'patient')
ORDER BY table_name, ordinal_position;

-- 0.b coverage of candidate sepsis code sets (to confirm the coding definition)
SELECT icd9code, count(*) AS n
FROM eicu_crd.diagnosis
WHERE icd9code LIKE '038%'  OR icd9code LIKE '99591%' OR icd9code LIKE '99592%'
   OR icd9code LIKE '78552%' OR icd9code LIKE 'A40%'   OR icd9code LIKE 'A41%'
   OR icd9code LIKE 'R652%'  OR icd9code LIKE 'R65.2%'
GROUP BY 1 ORDER BY n DESC LIMIT 50;

-- 0.c textual forms of the cardiac troponin I (cTnI) limit-of-detection (LOD) values
SELECT labresulttext, labmeasurenameinterface, count(*) AS n
FROM eicu_crd.lab
WHERE labname = 'troponin - I'
  AND (labresulttext ILIKE '%<%' OR labresulttext ILIKE '%less%'
       OR labresulttext ILIKE '%detect%')
GROUP BY 1, 2 ORDER BY n DESC LIMIT 50;


-- =============================================================================
-- PART 1 - primary extraction: stay-level cohort
--   :  data/eicu_extract_cohort_stays.csv
--   one row = one eligible ICU stay (including those with n_meas = 1, so that the
-- =============================================================================
WITH base AS (
  SELECT
    p.patientunitstayid, p.uniquepid, p.hospitalid, p.gender, p.age,
    p.unittype, p.unitvisitnumber, p.hospitaladmitsource,
    p.unitdischargeoffset, p.hospitaldischargeoffset,
    p.hospitaldischargestatus, p.hospitaldischargeyear,
    CASE WHEN p.age = '> 89' THEN 90
         WHEN p.age ~ '^[0-9]+$' THEN p.age::numeric
         ELSE NULL END AS age_num
  FROM eicu_crd.patient p
  WHERE p.unitvisitnumber = 1                    -- first ICU stay
),
sep AS (                                          -- sepsis: ICD-code and diagnosis-string pathways
  SELECT DISTINCT d.patientunitstayid
  FROM eicu_crd.diagnosis d
  WHERE d.icd9code LIKE '038%'  OR d.icd9code LIKE '99591%' OR d.icd9code LIKE '99592%'
     OR d.icd9code LIKE '78552%' OR d.icd9code LIKE 'A40%'   OR d.icd9code LIKE 'A41%'
     OR d.icd9code LIKE 'R652%'  OR d.icd9code LIKE 'R65.2%'
     OR lower(d.diagnosisstring) LIKE '%sepsis%'
     OR lower(d.diagnosisstring) LIKE '%septic%'
),
trop AS (                                         -- cTnI within the 0-72 h window plus quality flags
  SELECT
    l.patientunitstayid,
    count(*)                                                    AS n_rows_72h,
    count(*) FILTER (WHERE l.labresult IS NULL)                 AS n_null,
    count(*) FILTER (WHERE l.labresult <= 0)                    AS n_nonpos,
    count(*) FILTER (WHERE lower(coalesce(l.labmeasurenameinterface,''))
                            IN ('mg/ml','mg/dl'))               AS n_unit_implausible,
    count(*) FILTER (WHERE l.labresulttext ILIKE '%<%'
                        OR l.labresulttext ILIKE '%less%'
                        OR l.labresulttext ILIKE '%detect%')    AS n_below_detect,
    count(*) FILTER (WHERE l.labresult > 0
                       AND lower(coalesce(l.labmeasurenameinterface,''))
                           NOT IN ('mg/ml','mg/dl'))            AS n_meas,
    min(l.labresultoffset) FILTER (WHERE l.labresult > 0
                       AND lower(coalesce(l.labmeasurenameinterface,''))
                           NOT IN ('mg/ml','mg/dl'))            AS first_offset_min,
    max(l.labresultoffset) FILTER (WHERE l.labresult > 0
                       AND lower(coalesce(l.labmeasurenameinterface,''))
                           NOT IN ('mg/ml','mg/dl'))            AS last_offset_min,
    (array_agg(l.labresult ORDER BY l.labresultoffset ASC) FILTER (WHERE l.labresult > 0
                       AND lower(coalesce(l.labmeasurenameinterface,''))
                           NOT IN ('mg/ml','mg/dl')))[1]        AS first_val,
    (array_agg(l.labresult ORDER BY l.labresultoffset DESC) FILTER (WHERE l.labresult > 0
                       AND lower(coalesce(l.labmeasurenameinterface,''))
                           NOT IN ('mg/ml','mg/dl')))[1]        AS last_val,
    max(l.labresult) FILTER (WHERE l.labresult > 0
                       AND lower(coalesce(l.labmeasurenameinterface,''))
                           NOT IN ('mg/ml','mg/dl'))            AS max_val
  FROM eicu_crd.lab l
  WHERE l.labname = 'troponin - I'
    AND l.labresultoffset >= 0
    AND l.labresultoffset <= 4320                 -- 0-72 h
  GROUP BY l.patientunitstayid
)
SELECT
  b.patientunitstayid, b.uniquepid, b.hospitalid, b.hospitaldischargeyear,
  b.gender, b.age_num, b.age AS age_raw, b.unittype, b.hospitaladmitsource,
  b.unitdischargeoffset, b.hospitaldischargeoffset, b.hospitaldischargestatus,
  t.n_rows_72h, t.n_null, t.n_nonpos, t.n_unit_implausible, t.n_below_detect,
  t.n_meas, t.first_offset_min, t.last_offset_min, t.first_val, t.last_val, t.max_val,
  (b.hospitaldischargeoffset - 4320) / 1440.0                        AS days_from_landmark_to_discharge,
  CASE WHEN b.hospitaldischargestatus = 'Expired' THEN 1 ELSE 0 END  AS in_hosp_death
FROM base b
JOIN sep  s ON s.patientunitstayid = b.patientunitstayid
JOIN trop t ON t.patientunitstayid = b.patientunitstayid
WHERE b.age_num >= 18
  AND b.age_num <= 90
  AND b.unitdischargeoffset     >= 1440          -- ICU length of stay > 24 h
  AND b.hospitaldischargeoffset >= 4320          -- still in hospital at 72 h (also implies alive at 72 h)
  AND b.hospitaldischargestatus IN ('Alive','Expired')
ORDER BY b.hospitalid, b.patientunitstayid;


-- =============================================================================
-- PART 2 - per-site audit of cTnI units and limits of detection
--   :  data/eicu_extract_site_unit_lod.csv
-- =============================================================================
SELECT
  p.hospitalid,
  count(DISTINCT p.patientunitstayid)                                    AS stays_with_any_ctni,
  count(*)                                                               AS ctni_rows,
  count(DISTINCT lower(coalesce(l.labmeasurenameinterface,'')))          AS n_interface_units,
  string_agg(DISTINCT coalesce(l.labmeasurenameinterface,'(null)'), '|') AS interface_units,
  min(l.labresult)                                                       AS min_result,
  max(l.labresult)                                                       AS max_result,
  count(*) FILTER (WHERE l.labresult <= 0)                               AS n_nonpos,
  count(*) FILTER (WHERE l.labresulttext ILIKE '%<%'
                      OR l.labresulttext ILIKE '%less%'
                      OR l.labresulttext ILIKE '%detect%')               AS n_below_detect,
  min(l.labresultoffset)                                                 AS min_offset,
  max(l.labresultoffset)                                                 AS max_offset
FROM eicu_crd.lab l
JOIN eicu_crd.patient p ON p.patientunitstayid = l.patientunitstayid
WHERE l.labname = 'troponin - I'
  AND l.labresultoffset >= 0 AND l.labresultoffset <= 4320
GROUP BY p.hospitalid
ORDER BY stays_with_any_ctni DESC;


-- =============================================================================
-- PART 3 - site descriptors (hospital size / region / teaching status)
--   :  data/eicu_extract_site_hospital.csv
-- =============================================================================
SELECT * FROM eicu_crd.hospital ORDER BY hospitalid;


-- =============================================================================
-- PART 4 - severity covariates (APACHE IV; optional but recommended)
--   :  data/eicu_extract_apache.csv
--   If this part reports a missing column, the column name can be adjusted accordingly;
--   The analysis remains feasible without it (the adjustment set then reduces to
-- =============================================================================
SELECT
  a.patientunitstayid,
  max(a.apacheversion)               AS apacheversion,
  max(a.acutephysiologyscore)        AS acute_physiology_score,
  max(a.apachescore)                 AS apache_score,
  max(a.predictedhospitalmortality)  AS pred_hosp_mortality,
  max(a.actualhospitalmortality)     AS actual_hosp_mortality
FROM eicu_crd.apachepatientresult a
GROUP BY a.patientunitstayid
ORDER BY a.patientunitstayid;


-- =============================================================================
-- PART 5 - 72-hour clinical-state proxies (to mirror the primary MIMIC model M3a)
--   :  data/eicu_extract_window72.csv
--   Definition: eICU infusiondrug records only an offset (no start/end times), so
--   vaso72 is a conservative snapshot proxy (any vasopressor infusion charted in 48-72 h)
--   and lac72 is the last lactate value within the 48-72 h window (mmol/L).
-- =============================================================================
WITH coh AS (
  SELECT DISTINCT l.patientunitstayid
  FROM eicu_crd.lab l
  WHERE l.labname = 'troponin - I'
    AND l.labresultoffset >= 0 AND l.labresultoffset <= 4320
),
vaso AS (
  SELECT d.patientunitstayid, 1 AS vaso72
  FROM eicu_crd.infusiondrug d
  WHERE d.infusionoffset > 2880 AND d.infusionoffset <= 4320
    AND (lower(d.drugname) LIKE '%norepinephrine%'
      OR lower(d.drugname) LIKE '%epinephrine%'
      OR lower(d.drugname) LIKE '%vasopressin%'
      OR lower(d.drugname) LIKE '%phenylephrine%'
      OR lower(d.drugname) LIKE '%dopamine%'
      OR lower(d.drugname) LIKE '%dobutamine%')
  GROUP BY d.patientunitstayid
),
lac AS (
  SELECT l.patientunitstayid,
         (array_agg(l.labresult ORDER BY l.labresultoffset DESC))[1] AS lac72
  FROM eicu_crd.lab l
  WHERE l.labname = 'lactate'
    AND l.labresult IS NOT NULL
    AND l.labresultoffset > 2880 AND l.labresultoffset <= 4320
  GROUP BY l.patientunitstayid
)
SELECT c.patientunitstayid,
       COALESCE(v.vaso72, 0) AS vaso72,
       lc.lac72              AS lac72
FROM coh c
LEFT JOIN vaso v ON v.patientunitstayid = c.patientunitstayid
LEFT JOIN lac  lc ON lc.patientunitstayid = c.patientunitstayid
ORDER BY c.patientunitstayid;



