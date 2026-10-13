"""Audit: landmark cohort size, log-rank statistic, 90-day estimates."""
import io, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

df = pd.read_csv(os.path.join(BASE, 'data', 'v4_cohort_2plus_true_time.csv'))
print('full cohort n =', len(df), ' 28d deaths =', int(df['mortality_28'].sum()))

print()
print('=== LANDMARK COHORT DEFINITIONS ===')
alive_or_late = df[(df['mortality_28'] == 0) | (df['time_28d_true'] >= 72)]
print('  (mortality_28==0 OR time>=72)          n =', len(alive_or_late),
      ' deaths =', int(alive_or_late['mortality_28'].sum()))
died_before = df[(df['mortality_28'] == 1) & (df['time_28d_true'] < 72)]
print('  died before 72h                        n =', len(died_before))
print('  => subsequent deaths =', int(df['mortality_28'].sum()) - len(died_before))

# optional: per-stay landmark dataset produced earlier in the pipeline
lp = os.path.join(BASE, 'data', 'landmark_analysis_dataset.csv')
if os.path.exists(lp):
    ldf = pd.read_csv(lp)
    print()
    print('  locked landmark dataset shape:', ldf.shape)
    for c in ldf.columns:
        if 'mort' in c.lower() or 'death' in c.lower() or 'time' in c.lower():
            print(f'    {c}: sum={ldf[c].sum() if ldf[c].dtype!=object else "n/a"}')

print()
print('=== LOG-RANK (28-day, 4 groups) ===')
try:
    from lifelines.statistics import multivariate_logrank_test
    res = multivariate_logrank_test(df['time_28d_true'], df['trop_dynamic_group'],
                                    df['mortality_28'])
    print(f'  chi2 = {res.test_statistic:.2f}   p = {res.p_value:.3e}')
except ImportError:
    print('  lifelines not installed; computing manually')
    from scipy import stats
    # manual log-rank
    groups = df['trop_dynamic_group'].unique()
    df2 = df.copy()
    times = np.sort(df2['time_28d_true'].unique())
    O = {g: 0.0 for g in groups}
    E = {g: 0.0 for g in groups}
    V = np.zeros((len(groups), len(groups)))
    gi = {g: i for i, g in enumerate(groups)}
    for t in times:
        at_risk = df2['time_28d_true'] >= t
        n = at_risk.sum()
        if n < 2:
            continue
        d = ((df2['time_28d_true'] == t) & (df2['mortality_28'] == 1)).sum()
        if d == 0:
            continue
        for g in groups:
            in_g = df2['trop_dynamic_group'] == g
            n_g = (at_risk & in_g).sum()
            d_g = ((df2['time_28d_true'] == t) & (df2['mortality_28'] == 1) & in_g).sum()
            O[g] += d_g
            E[g] += d * n_g / n
            for g2 in groups:
                n_g2 = (at_risk & (df2['trop_dynamic_group'] == g2)).sum()
                if g == g2:
                    V[gi[g], gi[g2]] += d * (n_g / n) * (1 - n_g / n) * (n - d) / (n - 1) if n > 1 else 0
                else:
                    V[gi[g], gi[g2]] += -d * (n_g / n) * (n_g2 / n) * (n - d) / (n - 1) if n > 1 else 0
    Ovec = np.array([O[g] - E[g] for g in groups[:-1]])
    chi2 = float(Ovec @ np.linalg.pinv(V[:-1, :-1]) @ Ovec)
    print(f'  chi2 = {chi2:.2f}   p = {stats.chi2.sf(chi2, len(groups)-1):.3e}')
