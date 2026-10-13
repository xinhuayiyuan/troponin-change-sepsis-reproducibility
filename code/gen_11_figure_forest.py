"""Regenerate forest plot (Rising vs Stable-low HR across models/sensitivity)
with left-side row labels, dark colours, no overlaps. Reads locked CSVs."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCK = os.path.join(BASE, 'output', 'lock')
OUT = os.path.join(BASE, 'figures', 'fig3b_forest.png')

t2 = pd.read_csv(os.path.join(LOCK, 'table2_cox.csv'))
t4 = pd.read_csv(os.path.join(LOCK, 'table4_sensitivity.csv'))

def grab(df, model, n=None, events=None):
    r = df[(df['model'] == model) & (df['term'] == 'trop_dynamic_groupRising')].iloc[0]
    return (n if n is not None else int(r['n']),
            events if events is not None else int(r['events']),
            r['estimate'], r['ci_low'], r['ci_high'])

# display strings exactly as printed in manuscript Tables 2/4
TXT = {
    'cox_model1_crude': '1.68 (1.30\u20132.17)',
    'cox_model2_main': '1.50 (1.00\u20132.25)',
    'cox_model3_full': '1.00 (0.59\u20131.69)',
    'sens_ge4_measurements': '1.38 (0.18\u201310.49)',
    'sens_excluding_prior_mi': '0.62 (0.32\u20131.19)',
    'sens_landmark_72h': '1.17 (0.62\u20132.22)',
    'sens_ph_stratified': '1.09 (0.63\u20131.87)',
    'sens_mortality_90d': '1.07 (0.66\u20131.74)',
    'sens_landmark_72h_90d': '1.24 (0.70\u20132.18)',
    'sens_window72_organ': '1.03 (0.61\u20131.76)',
    'sens_anchor_year': '1.01 (0.59\u20131.70)',
    'sens_excl_early_falling': '1.01 (0.60\u20131.72)',
}

# numeric values for the three raw-database sensitivity analyses
NEW = {
    'sens_window72_organ': (902, 349, 1.03, 0.61, 1.76),
    'sens_anchor_year': (902, 349, 1.01, 0.59, 1.70),
    'sens_excl_early_falling': (776, 309, 1.01, 0.60, 1.72),
}

rows = [
    ('Model 1 (crude)', 'cox_model1_crude', t2, True),
    ('Model 2 (main adjustment)', 'cox_model2_main', t2, True),
    ('Model 3 (+ baseline troponin, LVEF)', 'cox_model3_full', t2, True),
    ('\u22654 troponin measurements', 'sens_ge4_measurements', t4, False),
    ('Excluding any coded myocardial infarction', 'sens_excluding_prior_mi', t4, False),
    ('72-hour landmark, categorical', 'sens_landmark_72h', t4, False),
    ('Proportional-hazards-stratified', 'sens_ph_stratified', t4, False),
    ('72h-window organ support', 'sens_window72_organ', None, False),
    ('Model 3 + admission year', 'sens_anchor_year', None, False),
    ('Excluding early-last Falling', 'sens_excl_early_falling', None, False),
    ('90-day mortality', 'sens_mortality_90d', t4, False),
    ('90-day mortality, 72-hour landmark', 'sens_landmark_72h_90d', t4, False),
]

DARK = '#111111'
GREY = '#5a5a5a'

from matplotlib.transforms import blended_transform_factory

fig, ax = plt.subplots(figsize=(12.5, 7.8), dpi=300)
fig.subplots_adjust(left=0.30, right=0.66, bottom=0.13, top=0.90)
ys = list(range(len(rows), 0, -1))

# blended transform: x in axes fraction, y in data units (for columns outside axes)
bx = blended_transform_factory(ax.transAxes, ax.transData)

ax.axvline(1.0, ls='--', lw=1.2, color='#333333')

for y, (label, key, df, bold) in zip(ys, rows):
    n, ev, hr, lo, hi = NEW[key] if df is None else grab(df, key)
    c = DARK if bold else GREY
    ax.plot([lo, hi], [y, y], color=c, lw=2.2 if bold else 1.8, solid_capstyle='round')
    ax.plot(hr, y, 's', color=c, markersize=9 if bold else 7.5)
    ax.text(-0.02, y, label, ha='right', va='center', fontsize=12,
            fontweight='bold' if bold else 'normal', color=c, transform=bx, clip_on=False)
    ax.text(1.04, y, f'{n:,} ({ev:,})', ha='left', va='center', fontsize=11,
            color=c, transform=bx, clip_on=False)
    ax.text(1.30, y, TXT[key], ha='left', va='center', fontsize=11,
            color=c, transform=bx, clip_on=False)

ax.text(1.04, len(rows) + 0.9, 'n (deaths)', ha='left', va='center',
        fontsize=11.5, fontweight='bold', color=DARK, transform=bx, clip_on=False)
ax.text(1.30, len(rows) + 0.9, 'HR (95% CI)', ha='left', va='center',
        fontsize=11.5, fontweight='bold', color=DARK, transform=bx, clip_on=False)

ax.set_xscale('log')
ax.set_xlim(0.12, 16)
ax.set_ylim(0.3, len(rows) + 1.5)
ax.set_yticks([])
ax.set_xticks([0.2, 0.5, 1, 2, 4, 8])
ax.set_xticklabels(['0.2', '0.5', '1', '2', '4', '8'], fontsize=12)
ax.set_xlabel('Hazard ratio, Rising vs Stable-low (log scale)', fontsize=13, color=DARK)
ax.set_title('Hazard ratios for Rising trajectory across models and sensitivity analyses',
             fontsize=14.5, fontweight='bold', color=DARK, pad=14)
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('#333333')
ax.tick_params(axis='x', colors=DARK)

fig.text(0.02, 0.025,
         'Bold rows: the three nested Cox models; grey rows: prespecified sensitivity analyses. '
         'n (deaths) are complete-case model samples.',
         fontsize=10.5, color=DARK)

fig.savefig(OUT, facecolor='white')
print('saved', OUT)
