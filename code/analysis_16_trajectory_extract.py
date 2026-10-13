"""Extract per-stay per-window median cTnT from MIMIC-IV for trajectory figure
regeneration at 300 dpi. Windows centred at 0,12,...,72 h (+-6 h; first [0,6)).
Output: output/traj_window_values.csv (stay_id, group, win, ctnT)"""
import sys, os
sys.stdout.reconfigure(encoding='utf-8')
import psycopg2
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, 'output')
df = pd.read_csv(os.path.join(BASE, 'data', 'v4_cohort_2plus_true_time.csv'),
                 usecols=['stay_id', 'hadm_id', 'subject_id', 'icu_intime', 'trop_dynamic_group'])

CFG = dict(host=os.environ.get('MIMIC_HOST', 'localhost'),
          port=int(os.environ.get('MIMIC_PORT', 5432)),
          dbname=os.environ.get('MIMIC_DB', 'mimic'),
          user=os.environ.get('MIMIC_USER', 'postgres'),
          password=os.environ.get('MIMIC_PASSWORD', ''))
conn = psycopg2.connect(**CFG)
cur = conn.cursor()
cur.execute('drop table if exists tmp_v4cohort')
cur.execute('create temp table tmp_v4cohort (stay_id bigint, hadm_id bigint, intime timestamp)')
rows = [(int(r.stay_id), int(r.hadm_id), r.icu_intime) for r in df.itertuples()]
cur.executemany('insert into tmp_v4cohort values (%s,%s,%s)', rows)

SQL = """
with le as (
  select w.stay_id,
         floor((extract(epoch from le.charttime - w.intime)/3600.0 + 6) / 12) * 12 as win,
         le.valuenum
  from tmp_v4cohort w
  join mimiciv_hosp.labevents le
    on le.hadm_id = w.hadm_id and le.itemid = 51003
   and le.charttime >= w.intime and le.charttime < w.intime + interval '78 hours'
   and le.valuenum is not null
), per_stay as (
  select stay_id, win, percentile_cont(0.5) within group (order by valuenum) as ctnT
  from le group by stay_id, win
)
select stay_id, win, ctnT from per_stay order by stay_id, win
"""
res = pd.read_sql(SQL, conn)
conn.close()
res = res.merge(df[['stay_id', 'trop_dynamic_group']], on='stay_id', how='left')
print('rows:', len(res), '| windows:', sorted(res['win'].unique()))
res.to_csv(os.path.join(OUT, 'traj_window_values.csv'), index=False)

# group medians/IQR for comparison with the existing figure
g = res.groupby(['trop_dynamic_group', 'win'])['ctnt'].agg(
    med='median', q1=lambda x: x.quantile(0.25), q3=lambda x: x.quantile(0.75), n='count')
for grp in ['Stable-low', 'Falling', 'Stable-high', 'Rising']:
    print('==', grp)
    print(g.loc[grp].round(3).to_string())

