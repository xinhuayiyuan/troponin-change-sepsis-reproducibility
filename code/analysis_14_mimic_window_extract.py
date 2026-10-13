"""Extract 72h-window organ-support flags, anchor year, and last-troponin timing
from the local MIMIC-IV PostgreSQL database for the locked cohort (stay_id keyed).
Output: output/mimic_window_vars.csv"""
import sys, os
sys.stdout.reconfigure(encoding='utf-8')
import psycopg2
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, 'output')
os.makedirs(OUT, exist_ok=True)

df = pd.read_csv(os.path.join(BASE, 'data', 'v4_cohort_2plus_true_time.csv'),
                 usecols=['stay_id', 'hadm_id', 'subject_id', 'icu_intime'])
print('cohort:', len(df))

CFG = dict(host=os.environ.get('MIMIC_HOST', 'localhost'),
          port=int(os.environ.get('MIMIC_PORT', 5432)),
          dbname=os.environ.get('MIMIC_DB', 'mimic'),
          user=os.environ.get('MIMIC_USER', 'postgres'),
          password=os.environ.get('MIMIC_PASSWORD', ''))
conn = psycopg2.connect(**CFG)
cur = conn.cursor()

cur.execute('drop table if exists tmp_v4cohort')
cur.execute('create temp table tmp_v4cohort (stay_id bigint, hadm_id bigint, subject_id bigint, intime timestamp)')
rows = [(int(r.stay_id), int(r.hadm_id), int(r.subject_id), r.icu_intime) for r in df.itertuples()]
cur.executemany('insert into tmp_v4cohort values (%s,%s,%s,%s)', rows)
print('uploaded:', len(rows))

SQL = """
with win as (
  select stay_id, hadm_id, subject_id, intime, intime + interval '72 hours' as t72
  from tmp_v4cohort
)
select
  w.stay_id,
  coalesce(v.vent72, 0) as vent72,
  coalesce(n.norepi72, 0) as norepi72,
  coalesce(r.rrt72, 0) as rrt72,
  p.anchor_year,
  p.anchor_year_group,
  t.last_trop_hours
from win w
left join (
  select distinct d.stay_id, 1 as vent72
  from mimiciv_derived.ventilation d
  join win w on d.stay_id = w.stay_id
  where d.ventilation_status = 'InvasiveVent' and d.starttime < w.t72
) v on v.stay_id = w.stay_id
left join (
  select distinct d.stay_id, 1 as norepi72
  from mimiciv_derived.norepinephrine d
  join win w on d.stay_id = w.stay_id
  where d.starttime < w.t72
) n on n.stay_id = w.stay_id
left join (
  select distinct d.stay_id, 1 as rrt72
  from mimiciv_derived.crrt d
  join win w on d.stay_id = w.stay_id
  where d.charttime < w.t72
  union
  select distinct d.stay_id, 1 as rrt72
  from mimiciv_derived.rrt d
  join win w on d.stay_id = w.stay_id
  where d.charttime < w.t72 and d.dialysis_present = 1
) r on r.stay_id = w.stay_id
left join mimiciv_hosp.patients p on p.subject_id = w.subject_id
left join (
  select w.stay_id,
         extract(epoch from max(le.charttime) - w.intime) / 3600.0 as last_trop_hours
  from win w
  join mimiciv_hosp.labevents le
    on le.hadm_id = w.hadm_id and le.itemid = 51003
   and le.charttime >= w.intime and le.charttime < w.t72
   and le.valuenum is not null
  group by w.stay_id, w.intime
) t on t.stay_id = w.stay_id
order by w.stay_id
"""
res = pd.read_sql(SQL, conn)
conn.close()
print('fetched:', len(res), '| stay_id match:', res['stay_id'].isin(df['stay_id']).all())
print('vent72:', res['vent72'].sum(), '| norepi72:', res['norepi72'].sum(),
      '| rrt72:', res['rrt72'].sum())
print('anchor_year NA:', res['anchor_year'].isna().sum(),
      '| last_trop_hours NA:', res['last_trop_hours'].isna().sum())
res.to_csv(os.path.join(OUT, 'mimic_window_vars.csv'), index=False)
print('saved mimic_window_vars.csv')
