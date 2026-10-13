"""Build the landmark analysis dataset.

Exposure:
    r            = (ln(last) - ln(first)) / dt_hours * 72      [log fold-change per 72 h]
    log_baseline = ln(first)
    n_meas       = number of valid cTnT measurements in 0-72 h

Landmark cohort: Sepsis-3 first ICU, >=2 cTnT in 0-72 h, alive at 72 h
AND still in hospital at 72 h.  Outcome: death after the landmark within 28 days
of ICU admission.  Follow-up time for survivors is censored at 28 days.

Every field is either read from the locked cohort file or from the direct
MIMIC-IV extraction (analysis_17); nothing is imputed in this step.

Output: data/landmark_analysis_dataset.csv
"""
import sys
import os
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding='utf-8')
pd.set_option('display.width', 200)

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
coh = pd.read_csv(os.path.join(BASE, 'data', 'v4_cohort_2plus_true_time.csv'))
lm = pd.read_csv(os.path.join(BASE, 'data', 'landmark_vars_mimic.csv'))

d = coh.merge(lm, on='stay_id', how='left', suffixes=('', '_lm'))
assert len(d) == len(coh), 'merge changed row count'

# ---------------------------------------------------------------- timing ----
d['intime'] = pd.to_datetime(d['icu_intime'])
d['first_t'] = pd.to_datetime(d['first_t'])
d['last_t'] = pd.to_datetime(d['last_t'])
d['disch'] = pd.to_datetime(d['dischtime'])
d['hours_to_disch'] = (d['disch'] - d['intime']).dt.total_seconds() / 3600.0
d['dt_hours'] = (d['last_t'] - d['first_t']).dt.total_seconds() / 3600.0

# ------------------------------------------------------------- exposures ----
d['first_tn'] = d['first_val']
d['last_tn'] = d['last_val']
d['log_baseline'] = np.log(d['first_tn'])
d['r'] = (np.log(d['last_tn']) - np.log(d['first_tn'])) / d['dt_hours'] * 72.0
d['n_meas'] = d['n_meas'].astype(int)

# ------------------------------------------------------ landmark cohort -----
d['alive_at_72'] = (d['time_28d_true'] > 72) | (d['mortality_28'] == 0)
d['in_hosp_at_72'] = d['hours_to_disch'] > 72
d['landmark'] = d['alive_at_72'] & d['in_hosp_at_72']

# outcome: death within 28 days of ICU admission, among landmark survivors
d['death28'] = d['mortality_28']
d['time_landmark'] = d['time_28d_true'] - 72.0          # days from landmark to event/censor
d['time_landmark'] = d['time_landmark'].clip(lower=0)

# lactate clearance (only where both available)
d['lactate_clear'] = (d['lactate_bl'] - d['lactate72']) / d['lactate_bl']

# ------------------------------------------------------------- reporting ----
n_all = len(d)
sub = d[d['landmark']].copy()
print('=== cohort ===')
print(f'total cohort                 : {n_all}')
print(f'alive at 72 h                : {int(d.alive_at_72.sum())}')
print(f'still in hospital at 72 h    : {int(d.in_hosp_at_72.sum())}')
print(f'LANDMARK (both)              : {len(sub)}')
print(f'  deaths after landmark (28d): {int(sub.death28.sum())}')
print(f'  survivors                  : {int((sub.death28 == 0).sum())}')
print()
print('=== exposure derivability (landmark cohort) ===')
print(f'dt_hours > 0        : {int((sub.dt_hours > 0).sum())} / {len(sub)}')
print(f'first_tn >= LOD     : {int((sub.first_tn > 0).sum())} / {len(sub)}')
print(f'r finite            : {int(np.isfinite(sub.r).sum())} / {len(sub)}')
print()
print('=== r distribution (landmark cohort) ===')
print(sub['r'].describe(percentiles=[.01, .05, .25, .5, .75, .95, .99]).round(4).to_string())
print()
print('=== covariate availability (landmark cohort) ===')
covs = ['sofa72', 'vaso72', 'lactate72', 'lactate_bl', 'age', 'gender', 'sofa_24h',
        'creatinine_max', 'lactate_max', 'myocardial_infarct',
        'congestive_heart_failure', 'invasive_vent', 'norepinephrine', 'min_lvef_72h']
for c in covs:
    miss = sub[c].isna().sum()
    print(f'  {c:<28} missing {miss:5d} / {len(sub)}  ({miss/len(sub)*100:5.1f}%)')
print()
print('=== definitional note: raw last vs locked window-min last ===')
diff = ~np.isclose(d['last_troponin'], d['last_val'], atol=1e-9, rtol=0)
print(f'  rows differing: {int(diff.sum())} / {len(d)} '
      f'({diff.mean()*100:.1f}%) - legacy rule = min value per 12-h window')

out = os.path.join(BASE, 'data', 'landmark_analysis_dataset.csv')
d.to_csv(out, index=False)
print()
print('saved', out, '| rows:', len(d))
