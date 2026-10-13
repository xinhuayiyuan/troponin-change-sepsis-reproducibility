"""Verify Table ST5 (included vs excluded) from the widened-cohort extract."""
import io, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
alt = pd.read_csv(os.path.join(BASE, 'data', 'alt_cohort_vars.csv'))
inc = pd.read_csv(os.path.join(BASE, 'data', 'v4_cohort_2plus_true_time.csv'),
                  usecols=['stay_id'])

print('widened cohort (>=1 measurement):', len(alt))
inc_ids = set(inc['stay_id'])

for col in ['n24', 'n48', 'n72']:
    ge2 = alt[alt[col].fillna(0) >= 2]
    print(f'  {col} >= 2 -> {len(ge2)}')

sub = alt[alt['n72'].fillna(0) >= 2].copy()
print()
print('>=2 measurements within 72 h (n72>=2):', len(sub))
sub['included'] = sub['stay_id'].isin(inc_ids)
print('  of which included in primary cohort:', int(sub['included'].sum()))
exc = sub[~sub['included']]
print('  excluded (indeterminate):', len(exc))

print()
print('=== EXCLUDED GROUP: recomputed vs Table ST5 ===')
print(f"  age, mean                    : {exc['age'].mean():.1f}   (ST5: 67.8)")
print(f"  SOFA score, mean             : {exc['sofa_24h'].mean():.1f}   (ST5: 7.0)")
print(f"  troponin measurements, mean  : {exc['n72'].mean():.1f}   (ST5: 3.0)")
print(f"  28-day mortality, %          : {100*exc['mortality_28'].mean():.1f}   (ST5: 29.7)")
print()
print('=== INCLUDED GROUP: recomputed vs Table ST5 ===')
i2 = sub[sub['included']]
print(f"  age, mean                    : {i2['age'].mean():.1f}   (ST5: 67.1)")
print(f"  SOFA score, mean             : {i2['sofa_24h'].mean():.1f}   (ST5: 7.4)")
print(f"  troponin measurements, mean  : {i2['n72'].mean():.1f}   (ST5: 3.4)")
print(f"  28-day mortality, %          : {100*i2['mortality_28'].mean():.1f}   (ST5: 32.8)")
