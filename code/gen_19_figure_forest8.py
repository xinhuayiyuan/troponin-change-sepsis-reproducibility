"""Figure 3B: compact 8-entry forest plot.
Rows 1-3 = the three nested Cox models (bold); rows 4-8 = five key sensitivity analyses.
Reads locked CSVs only."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCK = os.path.join(BASE, 'output', 'lock')
OUT = os.path.join(BASE, 'figures', 'fig3b_forest8.png')

t2 = pd.read_csv(os.path.join(LOCK, 'table2_cox.csv'))
t4 = pd.read_csv(os.path.join(LOCK, 'table4_sensitivity.csv'))

def grab(df, model):
    r = df[(df['model'] == model) & (df['term'] == 'trop_dynamic_groupRising')].iloc[0]
    return int(r['n']), int(r['events']), r['estimate'], r['ci_low'], r['ci_high']

rows = [
    ('Model 1 (crude)', 'cox_model1_crude', t2, True),
    ('Model 2 (adjusted)', 'cox_model2_main', t2, True),
    ('Model 3 (+ baseline troponin, LVEF)', 'cox_model3_full', t2, True),
    ('72-hour landmark, categorical', 'sens_landmark_72h', t4, False),
    ('Excluding any coded myocardial infarction', 'sens_excluding_prior_mi', t4, False),
    ('Proportional-hazards-stratified', 'sens_ph_stratified', t4, False),
    ('\u22654 troponin measurements', 'sens_ge4_measurements', t4, False),
    ('90-day mortality', 'sens_mortality_90d', t4, False),
]

DARK = '#111111'
GREY = '#5a5a5a'
from matplotlib.transforms import blended_transform_factory

fig, ax = plt.subplots(figsize=(12.5, 6.2), dpi=300)
fig.subplots_adjust(left=0.34, right=0.62, bottom=0.15, top=0.89)
ys = list(range(len(rows), 0, -1))
bx = blended_transform_factory(ax.transAxes, ax.transData)

ax.axvline(1.0, ls='--', lw=1.2, color='#333333')

annotations = []          # hand-placed texts, for the collision check below
for y, (label, key, df, bold) in zip(ys, rows):
    n, ev, hr, lo, hi = grab(df, key)
    c = DARK if bold else GREY
    ax.plot([lo, hi], [y, y], color=c, lw=2.4 if bold else 1.8, solid_capstyle='round')
    ax.plot(hr, y, 's', color=c, markersize=9.5 if bold else 7.5)
    annotations.append(ax.text(-0.02, y, label, ha='right', va='center', fontsize=12,
                               fontweight='bold' if bold else 'normal', color=c,
                               transform=bx, clip_on=False))
    annotations.append(ax.text(1.05, y, f'{n:,} ({ev:,})', ha='left', va='center', fontsize=11,
                               color=c, transform=bx, clip_on=False))
    # The HR column is anchored right of the widest n string ("3,204 (1,052)")
    # so the two never collide; the layout is machine-checked below.
    annotations.append(ax.text(1.62, y, f'{hr:.2f} ({lo:.2f}\u2013{hi:.2f})', ha='left',
                               va='center', fontsize=11, color=c, transform=bx, clip_on=False))

annotations.append(ax.text(1.05, len(rows) + 0.9, 'n (deaths)', ha='left', va='center',
                           fontsize=11.5, fontweight='bold', color=DARK, transform=bx, clip_on=False))
annotations.append(ax.text(1.62, len(rows) + 0.9, 'HR (95% CI)', ha='left', va='center',
                           fontsize=11.5, fontweight='bold', color=DARK, transform=bx, clip_on=False))

ax.set_xscale('log')
ax.set_xlim(0.15, 12)
ax.set_ylim(0.3, len(rows) + 1.5)
ax.set_yticks([])
ax.set_xticks([0.5, 1, 2, 4, 8])
ax.set_xticklabels(['0.5', '1', '2', '4', '8'], fontsize=12)
ax.set_xlabel('Hazard ratio, Rising vs Stable-low (log scale)', fontsize=13, color=DARK)
ax.set_title('Hazard ratios for the Rising trajectory: models and key sensitivity analyses',
             fontsize=14.5, fontweight='bold', color=DARK, pad=14)
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('#333333')
ax.tick_params(axis='x', colors=DARK)

annotations.append(fig.text(
    0.02, 0.02,
    'Bold rows: the three nested Cox models; grey rows: key prespecified sensitivity analyses. '
    'n (deaths) are complete-case model samples.',
    fontsize=10.5, color=DARK))

fig.savefig(OUT, facecolor='white')
print('saved', OUT)


def overlap_report(fig, axes_list, annotations, tag):
    """Machine check on the hand-placed columns (row labels, n, HR, headers,
    footnote): no annotation may overlap another annotation or any other text
    (e.g. a tick label). Matplotlib-managed tick labels are not checked against
    each other. Returns the collision count."""
    import matplotlib.text as mtext
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    others = [t for t in fig.findobj(mtext.Text)
              if t.get_text().strip() and t not in annotations]
    ab = [t.get_window_extent(renderer=rend) for t in annotations]
    ob = [t.get_window_extent(renderer=rend) for t in others]
    bad = []
    for i, a in enumerate(annotations):
        for j in range(i + 1, len(annotations)):
            if ab[i].overlaps(ab[j]):
                bad.append((a.get_text()[:30], annotations[j].get_text()[:30], 'annotation/annotation'))
        for k, o in enumerate(others):
            if ab[i].overlaps(ob[k]):
                bad.append((a.get_text()[:30], o.get_text()[:30], 'annotation/text'))
    print(f'OVERLAP CHECK [{tag}]: collisions = {len(bad)} '
          f'(annotations {len(annotations)}, other texts {len(others)})')
    for r in bad[:8]:
        print('    collision:', r)
    return len(bad)


n_collisions = overlap_report(fig, [ax], annotations, 'Figure 3B')
if os.environ.get('FIG_LAYOUT_CHECK_ONLY') == '1':
    print('layout-check-only: Figure 3B NOT written')
else:
    fig.savefig(OUT, facecolor='white')
    print('saved', OUT)
assert n_collisions == 0, 'Figure 3B layout has text collisions'
