import os
import pandas as pd
import numpy as np
import psycopg2

BASE = os.environ.get("PROJECT_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

DB = dict(host=os.environ.get('MIMIC_HOST', 'localhost'),
          port=int(os.environ.get('MIMIC_PORT', 5432)),
          dbname=os.environ.get('MIMIC_DB', 'mimic'),
          user=os.environ.get('MIMIC_USER', 'postgres'),
          password=os.environ.get('MIMIC_PASSWORD', ''))
conn = psycopg2.connect(**DB)
cur = conn.cursor()

cur.execute("""
SELECT c.stay_id,
    CASE WHEN pa.dod IS NOT NULL
         THEN EXTRACT(EPOCH FROM (pa.dod - c.icu_intime))/3600.0
         ELSE NULL END AS true_hours
FROM temp_v4_cohort c
INNER JOIN mimiciv_hosp.patients pa ON c.subject_id = pa.subject_id;
""")
tt = pd.DataFrame(cur.fetchall(), columns=['stay_id', 'true_hours'])
cur.close()
conn.close()

df = pd.read_csv(os.path.join(BASE, 'data', 'v4_cohort_2plus_clean.csv'))
df = df.merge(tt, on='stay_id', how='left')

# Corrected survival time
df['time_28d_true'] = 672.0
mask = df['true_hours'].notna()
df.loc[mask, 'time_28d_true'] = np.minimum(df.loc[mask, 'true_hours'], 672.0)
df.loc[df['time_28d_true'] <= 0, 'time_28d_true'] = 0.5

# For mortality_90: 90 days = 2160 hours
df['time_90d_true'] = 2160.0
df.loc[mask, 'time_90d_true'] = np.minimum(df.loc[mask, 'true_hours'], 2160.0)
df.loc[df['time_90d_true'] <= 0, 'time_90d_true'] = 0.5

print('Corrected times built.')
print('Deaths 28d by group:')
print(df.groupby('trop_dynamic_group')['mortality_28'].sum())
print()
print('Deaths 90d by group:')
print(df.groupby('trop_dynamic_group')['mortality_90'].sum())

df.to_csv(os.path.join(BASE, 'data', 'v4_cohort_2plus_true_time.csv'), index=False)
print('Saved.')
