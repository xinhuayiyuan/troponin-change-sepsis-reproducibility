import psycopg2
import pandas as pd

DB = dict(host=os.environ.get('MIMIC_HOST', 'localhost'),
          port=int(os.environ.get('MIMIC_PORT', 5432)),
          dbname=os.environ.get('MIMIC_DB', 'mimic'),
          user=os.environ.get('MIMIC_USER', 'postgres'),
          password=os.environ.get('MIMIC_PASSWORD', ''))
conn = psycopg2.connect(**DB)
cur = conn.cursor()

print('=== FLOW RECONCILIATION ===')
cur.execute("SELECT COUNT(*) FROM temp_v4_cohort;")
print('temp_v4_cohort:', cur.fetchone()[0])

cur.execute("SELECT COUNT(*) FROM temp_v4_troponin;")
print('temp_v4_troponin (any measurement):', cur.fetchone()[0])

cur.execute("SELECT COUNT(*) FROM temp_v4_trop_dynamic WHERE trop_valid_count >= 2;")
print('>=2 measurements:', cur.fetchone()[0])

cur.execute("SELECT COUNT(*) FROM temp_v4_trop_dynamic WHERE trop_valid_count >= 2 AND trop_dynamic_group != 'Unknown';")
print('>=2 AND classified:', cur.fetchone()[0])

cur.execute("SELECT COUNT(*) FROM temp_v4_trop_dynamic WHERE trop_valid_count >= 2 AND trop_dynamic_group = 'Unknown';")
print('>=2 AND Unknown:', cur.fetchone()[0])

cur.execute("SELECT COUNT(*) FROM temp_v4_trop_dynamic WHERE trop_valid_count < 2;")
print('<2 measurements (in trop_dynamic):', cur.fetchone()[0])

# Also check v4_analysis_cohort_2plus actual count
cur.execute("SELECT COUNT(*) FROM v4_analysis_cohort_2plus;")
print('v4_analysis_cohort_2plus:', cur.fetchone()[0])

# Where did 4262 come from? Check alternative: >=2 valid in 0-72h strict
cur.execute("SELECT COUNT(*) FROM temp_v4_trop_dynamic WHERE trop_valid_count >= 2 AND baseline_troponin IS NOT NULL;")
print('>=2 AND baseline non-null:', cur.fetchone()[0])

cur.close()
conn.close()
