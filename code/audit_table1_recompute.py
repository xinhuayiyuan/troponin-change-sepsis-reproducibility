"""Independent recomputation of manuscript numbers from raw cohort data.
Read-only audit: compares manuscript claims against the raw dataset and lock files.
"""
import io, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
df = pd.read_csv(os.path.join(BASE, 'data', 'v4_cohort_2plus_true_time.csv'))

order = ['Stable-low', 'Falling', 'Stable-high', 'Rising']
df['g'] = pd.Categorical(df['trop_dynamic_group'], categories=order, ordered=True)

print('=== TABLE 1 RECOMPUTATION (mean (SD) / n (%)) ===')
rows = {
    'n': lambda g: f"{len(g)}",
    'age': lambda g: f"{g['age'].mean():.1f} ({g['age'].std():.1f})",
    'male %': lambda g: f"{100*(g['gender']=='M').mean():.1f}",
    'sofa': lambda g: f"{g['sofa_24h'].mean():.1f} ({g['sofa_24h'].std():.1f})",
    'icu_los': lambda g: f"{g['los_icu'].mean():.1f} ({g['los_icu'].std():.1f})",
    'inv_vent %': lambda g: f"{100*g['invasive_vent'].mean():.1f}",
    'norepi %': lambda g: f"{100*g['norepinephrine'].mean():.1f}",
    'rrt %': lambda g: f"{100*g['rrt'].mean():.1f}",
    'mi %': lambda g: f"{100*g['myocardial_infarct'].mean():.1f}",
    'chf %': lambda g: f"{100*g['congestive_heart_failure'].mean():.1f}",
    'ihd %': lambda g: f"{100*g['ischemic_heart_disease'].mean():.1f}",
    'htn %': lambda g: f"{100*g['hypertension'].mean():.1f}",
    'dm %': lambda g: f"{100*g['diabetes'].mean():.1f}",
    'ckd %': lambda g: f"{100*g['ckd'].mean():.1f}",
    'creat_max': lambda g: f"{g['creatinine_max'].mean():.1f} ({g['creatinine_max'].std():.1f})",
    'lactate_max': lambda g: f"{g['lactate_max'].mean():.1f} ({g['lactate_max'].std():.1f})",
    'plt_min': lambda g: f"{g['platelets_min'].mean():.1f} ({g['platelets_min'].std():.1f})",
    'lvef': lambda g: f"{g['min_lvef_72h'].mean():.1f} ({g['min_lvef_72h'].std():.1f})",
    'echo %': lambda g: f"{100*g['has_echo'].mean():.1f}",
    'mort28 %': lambda g: f"{100*g['mortality_28'].mean():.1f}",
}
hdr = f"{'variable':14s}" + ''.join(f"{g:>16s}" for g in order)
print(hdr)
for name, fn in rows.items():
    line = f"{name:14s}"
    for g in order:
        line += f"{fn(df[df['g'] == g]):>16s}"
    print(line)

print()
print('=== TROPONIN MEDIANS [IQR] ===')
for col in ['baseline_troponin', 'peak_troponin', 'delta_troponin']:
    line = f"{col:20s}"
    for g in order:
        v = df.loc[df['g'] == g, col]
        line += f"{v.median():>10.2f} [{v.quantile(.25):.2f}-{v.quantile(.75):.2f}]  "
    print(line)

print()
print('=== MISSINGNESS ===')
tot = len(df)
for col in ['lactate_max', 'min_lvef_72h', 'creatinine_max', 'sofa_24h',
            'myocardial_infarct', 'congestive_heart_failure']:
    print(f"  {col:22s} missing {df[col].isna().sum():5d} ({100*df[col].isna().mean():.1f}%)")
print(f"  echo performed          {int(df['has_echo'].sum())} ({100*df['has_echo'].mean():.1f}%)")
