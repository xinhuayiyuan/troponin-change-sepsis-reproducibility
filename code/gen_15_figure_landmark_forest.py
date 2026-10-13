"""Generate Figure S4: landmark continuous-rate sensitivity forest plot.
Based on Table S14 values."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, 'figures', 'figS4_landmark_forest.png')

# (label, n, deaths, hr, lo, hi)
rows = [
    ('S0  Primary model (M3a)', 2107, 616, 1.00, 0.97, 1.03),
    ('S2  Exact landmark (6h grace)', 2079, 591, 1.00, 0.97, 1.03),
    ('S4  Multiple imputation (m=20)', 2834, 735, 1.00, 0.97, 1.02),
    ('S6  Excluding CKD or RRT', 1412, 381, 1.00, 0.96, 1.03),
    ('S7  In-hospital death', 2107, 521, 0.99, 0.96, 1.02),
    ('S9  Echo subgroup (M3b)', 1499, 454, 1.00, 0.97, 1.03),
    ('S10 LOD/√2 substitution', 2107, 616, 1.00, 0.97, 1.03),
    ('S11 dt≥6h between measurements', 1914, 551, 1.00, 0.97, 1.04),
    ('eICU external validation', 3111, 574, 1.02, 1.00, 1.04),
]

DARK = '#111111'
GREY = '#5a5a5a'
BLUE = '#1f4e79'

from matplotlib.transforms import blended_transform_factory

fig, ax = plt.subplots(figsize=(11, 6.5), dpi=300)
fig.subplots_adjust(left=0.32, right=0.64, bottom=0.14, top=0.88)
ys = list(range(len(rows), 0, -1))
bx = blended_transform_factory(ax.transAxes, ax.transData)

ax.axvline(1.0, ls='--', lw=1.2, color='#333333')

for y, (label, n, ev, hr, lo, hi) in zip(ys, rows):
    c = BLUE if 'Primary' in label else GREY
    lw = 2.6 if 'Primary' in label else 1.8
    ms = 9 if 'Primary' in label else 7.5
    ax.plot([lo, hi], [y, y], color=c, lw=lw, solid_capstyle='round')
    ax.plot(hr, y, 's', color=c, markersize=ms)
    ax.text(-0.02, y, label, ha='right', va='center', fontsize=12,
            fontweight='bold' if 'Primary' in label else 'normal', color=c, transform=bx, clip_on=False)
    ax.text(1.04, y, f'{n:,} ({ev:,})', ha='left', va='center', fontsize=11,
            color=c, transform=bx, clip_on=False)
    ax.text(1.34, y, f'{hr:.2f} ({lo:.2f}\u2013{hi:.2f})', ha='left', va='center', fontsize=11,
            color=c, transform=bx, clip_on=False)

ax.text(1.04, len(rows) + 0.8, 'n (deaths)', ha='left', va='center',
        fontsize=11.5, fontweight='bold', color=DARK, transform=bx, clip_on=False)
ax.text(1.34, len(rows) + 0.8, 'HR (95% CI)', ha='left', va='center',
        fontsize=11.5, fontweight='bold', color=DARK, transform=bx, clip_on=False)

ax.set_xlim(0.92, 1.09)
ax.set_ylim(0.3, len(rows) + 1.4)
ax.set_yticks([])
ax.set_xticks([0.93, 0.95, 0.97, 1.00, 1.03, 1.05, 1.07])
ax.set_xticklabels(['0.93', '0.95', '0.97', '1.00', '1.03', '1.05', '1.07'], fontsize=12)
ax.set_xlabel('Hazard ratio per 1 unit of r (linear scale)', fontsize=13, color=DARK)
ax.set_title('Landmark continuous-rate sensitivity analyses',
             fontsize=14.5, fontweight='bold', color=DARK, pad=14)
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('#333333')
ax.tick_params(axis='x', colors=DARK)

fig.text(0.02, 0.025,
         'Bold row: primary model (S0). All models fitted in the 72-hour landmark cohort unless otherwise specified.',
         fontsize=10.5, color=DARK)

fig.savefig(OUT, facecolor='white')
print('saved', OUT)
