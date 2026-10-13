"""Regenerate the trajectory panel (Figure 3A) at 300 dpi from the original
windowed summary data (data/trajectory_window_summary.csv).
Layout mirrors the previously verified panel; note wording aligned with the
'classification threshold' terminology."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, 'data', 'trajectory_window_summary.csv')
OUT = os.path.join(BASE, 'figures', 'fig3a_trajectories.png')

d = pd.read_csv(DATA)
groups = ['Stable-low', 'Falling', 'Stable-high', 'Rising']
colors = {'Stable-low': '#2ca02c', 'Falling': '#1f1fb8', 'Stable-high': '#e69f00', 'Rising': '#d62728'}
ns = {'Stable-low': 284, 'Falling': 673, 'Stable-high': 1457, 'Rising': 790}

fig, ax = plt.subplots(figsize=(11.5, 8.0), dpi=300)
fig.subplots_adjust(left=0.145, right=0.985, top=0.845, bottom=0.30)

for g in groups:
    sub = d[d['trop_dynamic_group'] == g].sort_values('window_h')
    x, m = sub['window_h'].values, sub['median'].values
    ax.fill_between(x, sub['q1'].values, sub['q3'].values, color=colors[g], alpha=0.18, lw=0)
    ax.plot(x, m, '-o', color=colors[g], lw=2.4, ms=6.5,
            label=f"{g} (n={ns[g]:,})")

ax.axhline(0.03, ls=':', lw=1.4, color='#555555')
ax.text(70.5, 0.0345, '0.03 ng/mL', fontsize=11, color='#555555', ha='right')

ax.set_yscale('log')
ax.set_ylim(0.011, 7)
ax.set_yticks([0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
ax.set_yticklabels(['0.02', '0.05', '0.10', '0.20', '0.50', '1.00', '2.00', '5.00'], fontsize=12.5)
ax.set_xticks([0, 12, 24, 36, 48, 60, 72])
ax.set_xticklabels(['0', '12', '24', '36', '48', '60', '72'], fontsize=12.5)
ax.set_xlim(-3, 75)
ax.set_ylabel('Cardiac troponin T, ng/mL (log scale)', fontsize=14.5)
ax.set_xlabel('Hours from ICU admission', fontsize=15, fontweight='bold', labelpad=10)
ax.set_title('Median cardiac troponin T trajectories by dynamic group',
             fontsize=17, fontweight='bold', pad=26)
ax.text(0.5, 1.045, 'First 72 hours of ICU admission; bands show interquartile range (25th-75th percentile)',
        transform=ax.transAxes, ha='center', fontsize=12, color='#555555')
leg = ax.legend(loc='upper left', fontsize=12, title='Dynamic group', title_fontsize=12.5,
                frameon=False, borderaxespad=0.1)
for t, g in zip(leg.get_texts(), groups):
    t.set_color(colors[g])

note = [
    'Note: cTnT measured by itemid 51003 (MIMIC-IV); assay generation not specified in the database.',
    'The Stable-low trajectory ends at 60 h because no patient in this group had a measurement in the 66-72 h window,',
    'reflecting cessation of serial testing when values remain at or below the classification threshold',
    '(mean 4.7 measurements in Rising versus 2.5 in Stable-low; Supplementary Table ST2).',
]
y = 0.145
for ln in note:
    fig.text(0.115, y, ln, fontsize=10.5, color='#222222')
    y -= 0.032

fig.savefig(OUT, facecolor='white')
print('saved', OUT)
