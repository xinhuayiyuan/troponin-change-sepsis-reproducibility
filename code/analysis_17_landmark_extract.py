"""Extract landmark (72 h) variables
from the local MIMIC-IV PostgreSQL database for the locked cohort.

Adds, per stay_id:
  first_t / last_t / n_meas  : cTnT (itemid 51003) measurement timing within 0-72 h
  sofa72                     : max 24-h rolling SOFA whose window ends in (48, 72] h
  vaso72                     : any vasoactive agent active across the 48-72 h window
  lactate72 / lactate_bl     : last lactate in (48,72] h / first lactate in [0,24) h
  dischtime / deathtime      : for the "still in hospital at 72 h" landmark rule

No values are imputed or invented here: every column is either read directly
from the database or a documented aggregation of read values.

Output: data/landmark_vars_mimic.csv  (+ copy in the new lock folder)
"""
import sys
import os
import psycopg2
import pandas as pd

sys.stdout.reconfigure(encoding='utf-8')

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COHORT = os.path.join(BASE, 'data', 'v4_cohort_2plus_true_time.csv')
OUT = os.path.join(BASE, 'data', 'landmark_vars_mimic.csv')

CFG = dict(host=os.environ.get('MIMIC_HOST', 'localhost'),
          port=int(os.environ.get('MIMIC_PORT', 5432)),
          dbname=os.environ.get('MIMIC_DB', 'mimic'),
          user=os.environ.get('MIMIC_USER', 'postgres'),
          password=os.environ.get('MIMIC_PASSWORD', ''))

coh = pd.read_csv(COHORT, usecols=['stay_id', 'hadm_id', 'subject_id', 'icu_intime'])
print('cohort rows:', len(coh))

conn = psycopg2.connect(**CFG)
cur = conn.cursor()
cur.execute('drop table if exists tmp_landmark_cohort')
cur.execute('create temp table tmp_landmark_cohort '
            '(stay_id bigint, hadm_id bigint, subject_id bigint, intime timestamp)')
cur.executemany('insert into tmp_landmark_cohort values (%s,%s,%s,%s)',
                [(int(r.stay_id), int(r.hadm_id), int(r.subject_id), r.icu_intime)
                 for r in coh.itertuples()])
print('uploaded to temp table')

SQL = """
with coh as (
  select stay_id, hadm_id, subject_id, intime from tmp_landmark_cohort
),
trop as (
  select c.stay_id,
         min(le.charttime)                                        as first_t,
         max(le.charttime)                                        as last_t,
         count(*)                                                 as n_meas,
         (array_agg(le.valuenum order by le.charttime))[1]         as first_val,
         (array_agg(le.valuenum order by le.charttime desc))[1]    as last_val
  from coh c
  join mimiciv_hosp.labevents le
    on le.hadm_id = c.hadm_id and le.itemid = 51003
   and le.charttime >= c.intime
   and le.charttime <  c.intime + interval '72 hours'
   and le.valuenum is not null
  group by c.stay_id
),
sofa72 as (
  select c.stay_id, max(s.sofa_24hours) as sofa72
  from coh c
  join mimiciv_derived.sofa s on s.stay_id = c.stay_id
  where s.endtime >  c.intime + interval '48 hours'
    and s.endtime <= c.intime + interval '72 hours'
  group by c.stay_id
),
vaso72 as (
  select distinct c.stay_id, 1 as vaso72
  from coh c
  join mimiciv_derived.vasoactive_agent v on v.stay_id = c.stay_id
  where v.starttime < c.intime + interval '72 hours'
    and v.endtime   > c.intime + interval '48 hours'
),
lac72 as (
  select c.stay_id, (array_agg(b.lactate order by b.charttime desc))[1] as lactate72
  from coh c
  join mimiciv_derived.bg b on b.hadm_id = c.hadm_id
  where b.lactate is not null
    and b.charttime >  c.intime + interval '48 hours'
    and b.charttime <= c.intime + interval '72 hours'
  group by c.stay_id
),
lacbl as (
  select c.stay_id, (array_agg(b.lactate order by b.charttime))[1] as lactate_bl
  from coh c
  join mimiciv_derived.bg b on b.hadm_id = c.hadm_id
  where b.lactate is not null
    and b.charttime >= c.intime
    and b.charttime <  c.intime + interval '24 hours'
  group by c.stay_id
)
select c.stay_id,
       t.first_t, t.last_t, t.n_meas, t.first_val, t.last_val,
       s.sofa72,
       coalesce(v.vaso72, 0) as vaso72,
       l72.lactate72, lb.lactate_bl,
       a.dischtime, a.deathtime, a.hospital_expire_flag
from coh c
left join trop   t   on t.stay_id   = c.stay_id
left join sofa72 s   on s.stay_id   = c.stay_id
left join vaso72 v   on v.stay_id   = c.stay_id
left join lac72  l72 on l72.stay_id = c.stay_id
left join lacbl  lb  on lb.stay_id  = c.stay_id
left join mimiciv_hosp.admissions a on a.hadm_id = c.hadm_id
order by c.stay_id
"""

res = pd.read_sql(SQL, conn)
conn.close()

print('fetched rows:', len(res))
print('stay_id coverage:', res['stay_id'].isin(coh['stay_id']).all())
print('n_meas NA:', res['n_meas'].isna().sum())
print('sofa72 NA:', res['sofa72'].isna().sum())
print('vaso72 = 1:', int((res['vaso72'] == 1).sum()))
print('lactate72 NA:', res['lactate72'].isna().sum(),
      '| lactate_bl NA:', res['lactate_bl'].isna().sum())
print('dischtime NA:', res['dischtime'].isna().sum(),
      '| deathtime NA:', res['deathtime'].isna().sum())

res.to_csv(OUT, index=False)
print('saved', OUT)
