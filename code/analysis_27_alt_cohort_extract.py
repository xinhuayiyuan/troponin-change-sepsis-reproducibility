"""Sensitivity-cohort extraction.

Builds, from the ORIGINAL cohort pipeline's own intermediate tables in the local
MIMIC-IV PostgreSQL database (public.temp_v4_cohort / temp_v4_troponin /
temp_v4_clinical / temp_v4_comorb), the data needed for analyses the locked
landmark dataset cannot support:

  (a) sensitivity 9 : widen inclusion to >=1 cTnT in 0-72 h (rate missing for
                      single-measurement patients -> imputation sensitivity)
  (b) sensitivity 1 : alternative landmark times (24 h and 48 h)
  (c) negative-control outcome: 30-day readmission with a NON-cardiac principal
                      diagnosis

Anchor check: temp_v4_cohort must contain 24,078 rows, temp_v4_troponin 6,261 and
public.v4_analysis_cohort_2plus 3,204 - exactly the values in the locked
cohort_flow.csv. The script aborts if they differ.

Read-only (SELECT only). Output: data/alt_cohort_vars.csv
"""
import os
import sys

import numpy as np
import pandas as pd
import psycopg2

sys.stdout.reconfigure(encoding='utf-8')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, 'data', 'alt_cohort_vars.csv')
CFG = dict(host=os.environ.get('MIMIC_HOST', 'localhost'),
          port=int(os.environ.get('MIMIC_PORT', 5432)),
          dbname=os.environ.get('MIMIC_DB', 'mimic'),
          user=os.environ.get('MIMIC_USER', 'postgres'),
          password=os.environ.get('MIMIC_PASSWORD', ''))

conn = psycopg2.connect(**CFG)
cur = conn.cursor()

for tbl, expect in [('public.temp_v4_cohort', 24078),
                    ('public.temp_v4_troponin', 6261),
                    ('public.v4_analysis_cohort_2plus', 3204)]:
    cur.execute('select count(*) from %s' % tbl)
    got = cur.fetchone()[0]
    print('%-36s %6d (expected %6d) %s' % (tbl, got, expect, 'OK' if got == expect else 'MISMATCH'))
    if got != expect:
        raise SystemExit('anchor mismatch - refusing to build the sensitivity cohort')

SQL = """
with w as (
  select c.stay_id,
    count(*) filter (where le.charttime < c.icu_intime + interval '24 hours') as n24,
    min(le.charttime) filter (where le.charttime < c.icu_intime + interval '24 hours') as first_t24,
    max(le.charttime) filter (where le.charttime < c.icu_intime + interval '24 hours') as last_t24,
    (array_agg(le.valuenum order by le.charttime asc) filter (
        where le.charttime < c.icu_intime + interval '24 hours'))[1] as first_v24,
    (array_agg(le.valuenum order by le.charttime desc) filter (
        where le.charttime < c.icu_intime + interval '24 hours'))[1] as last_v24,
    count(*) filter (where le.charttime < c.icu_intime + interval '48 hours') as n48,
    min(le.charttime) filter (where le.charttime < c.icu_intime + interval '48 hours') as first_t48,
    max(le.charttime) filter (where le.charttime < c.icu_intime + interval '48 hours') as last_t48,
    (array_agg(le.valuenum order by le.charttime asc) filter (
        where le.charttime < c.icu_intime + interval '48 hours'))[1] as first_v48,
    (array_agg(le.valuenum order by le.charttime desc) filter (
        where le.charttime < c.icu_intime + interval '48 hours'))[1] as last_v48,
    count(*) as n72,
    min(le.charttime) as first_t72,
    max(le.charttime) as last_t72,
    (array_agg(le.valuenum order by le.charttime asc))[1]  as first_v72,
    (array_agg(le.valuenum order by le.charttime desc))[1] as last_v72
  from public.temp_v4_troponin t
  join public.temp_v4_cohort c on c.stay_id = t.stay_id
  join mimiciv_hosp.labevents le
    on le.hadm_id = c.hadm_id and le.itemid = 51003
   and le.valuenum is not null
   and le.charttime >= c.icu_intime
   and le.charttime <  c.icu_intime + interval '72 hours'
  group by c.stay_id
),
nxt as (
  select hadm_id,
         lead(admittime) over (partition by subject_id order by admittime) as next_admittime,
         lead(hadm_id)   over (partition by subject_id order by admittime) as next_hadm_id
  from mimiciv_hosp.admissions
)
select c.stay_id, c.hadm_id, c.subject_id, c.icu_intime, c.icu_72h,
       c.ethnicity, c.gender, c.age, c.hosp_mort, c.mortality_28, c.mortality_90, c.los_icu,
       cl.sofa_24h, cl.creatinine_max, cl.lactate_max, cl.invasive_vent,
       cl.norepinephrine, cl.rrt,
       cm.ckd, cm.congestive_heart_failure, cm.myocardial_infarct,
       a.dischtime, a.deathtime, a.hospital_expire_flag, p.dod,
       w.n24, w.first_t24, w.last_t24, w.first_v24, w.last_v24,
       w.n48, w.first_t48, w.last_t48, w.first_v48, w.last_v48,
       w.n72, w.first_t72, w.last_t72, w.first_v72, w.last_v72,
       n.next_admittime, n.next_hadm_id,
       dx.icd_code as next_principal_icd, dx.icd_version as next_principal_icd_version
from public.temp_v4_troponin t
join public.temp_v4_cohort  c on c.stay_id = t.stay_id
left join public.temp_v4_clinical cl on cl.stay_id = c.stay_id
left join public.temp_v4_comorb   cm on cm.hadm_id = c.hadm_id
left join mimiciv_hosp.admissions a on a.hadm_id = c.hadm_id
left join mimiciv_hosp.patients   p on p.subject_id = c.subject_id
left join w on w.stay_id = c.stay_id
left join nxt n on n.hadm_id = c.hadm_id
left join mimiciv_hosp.diagnoses_icd dx
       on dx.hadm_id = n.next_hadm_id and dx.seq_num = 1
order by c.stay_id
"""

df = pd.read_sql(SQL, conn)
conn.close()
print('rows extracted:', len(df))
print('with >=1 cTnT in 0-72h join:', int(df.n72.notna().sum()))

# ---- derived fields, mirroring analysis_18 exactly --------------------------
hr = lambda a, b: (pd.to_datetime(a) - pd.to_datetime(b)).dt.total_seconds() / 3600.0

df['hours_to_disch'] = hr(df['dischtime'], df['icu_intime'])
dod_h = hr(df['dod'], df['icu_intime'])
df['time_28d_true'] = np.where(df['mortality_28'] == 1, dod_h, 672.0)
df['alive_at_72'] = (df['time_28d_true'] > 72) | (df['mortality_28'] == 0)
df['in_hosp_at_72'] = df['hours_to_disch'] > 72
df['landmark'] = df['alive_at_72'] & df['in_hosp_at_72']
df['death28'] = df['mortality_28']
df['time_landmark'] = (df['time_28d_true'] - 72.0).clip(lower=0)

# in-hospital death as an alternative outcome (prespecified sensitivity analysis)
df['t_hosp'] = (df['hours_to_disch'] - 72.0).clip(lower=0)
df['event_hosp'] = ((df['hospital_expire_flag'] == 1) & (df['t_hosp'] <= 28 * 24)).astype(int)
df['t_hosp'] = df['t_hosp'].clip(upper=28 * 24)

# 30-day readmission with a NON-cardiac principal diagnosis (negative control)
def is_cardiac(code, ver):
    if not isinstance(code, str):
        return False
    c = code.upper().replace('.', '')
    if ver == 10 or (ver is None and c[:1].isalpha()):
        return c.startswith('I')
    return c.startswith(('410', '411', '412', '413', '414', '425', '427', '428'))


df['next_within_30d'] = (
    (pd.to_datetime(df['next_admittime']) - pd.to_datetime(df['dischtime']))
    .dt.total_seconds() / 86400.0).le(30).fillna(False)
df['next_cardiac'] = [is_cardiac(c, v) for c, v in
                      zip(df['next_principal_icd'], df['next_principal_icd_version'])]
df['negctrl_readmit_noncardiac'] = (df['next_within_30d'] & ~df['next_cardiac']).astype(int)

# rate over each window, always standardised to 72 h
for w in (24, 48, 72):
    dt = hr(df['last_t%d' % w], df['first_t%d' % w])
    df['dt_h%d' % w] = dt
    with np.errstate(divide='ignore', invalid='ignore'):
        df['r%d' % w] = (np.log(df['last_v%d' % w]) - np.log(df['first_v%d' % w])) / dt * 72.0

print('landmark cohort (>=1 cTnT):', int(df.landmark.sum()))
print('  of which n72>=2 :', int((df.landmark & (df.n72 >= 2)).sum()))
print('  of which n72==1 :', int((df.landmark & (df.n72 == 1)).sum()))
print('negative-control events (non-cardiac readmission <=30d):',
      int(df.negctrl_readmit_noncardiac.sum()))
df.to_csv(OUT, index=False)
print('saved', OUT)
