"""Supplementary Figure SF6: external validation forest plot (MIMIC-IV vs eICU-CRD).

Reads locked outputs only:
  landmark/landmark_attenuation_ladder.csv  (MIMIC-IV M3a, per unit)
  landmark/... -> per-SD from sensitivity/sens_main_table.csv (S0 row)
  eicu/eicu_site_gradient_models.csv             (eICU, cluster-robust)
Output: figures/figS6_eicu_forest.png
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCK = os.path.join(BASE, 'output', 'lock')
OUT = os.path.join(BASE, 'figures', 'figS6_eicu_forest.png')

lad = pd.read_csv(os.path.join(LOCK, 'landmark', 'landmark_attenuation_ladder.csv'))
m3a = lad[lad['model'] == 'M3a +landmark state (main)'].iloc[0]
sens = pd.read_csv(os.path.join(LOCK, 'sensitivity', 'sens_main_table.csv'))
s0 = sens[sens['analysis'].str.contains('S0 primary model')].iloc[0]
eicu = pd.read_csv(os.path.join(LOCK, 'eicu', 'eicu_site_gradient_models.csv'))

rows_unit = [('MIMIC-IV primary (M3a)', m3a['HR'], m3a['lower'], m3a['upper'])]
rows_sd = [('MIMIC-IV primary (M3a)', s0['HR_per_SD'], s0['lo_per_SD'], s0['hi_per_SD'])]
labels = {'eICU crude (cluster by hospital)': 'eICU crude',
          'eICU + admission-level covariates (cluster)': 'eICU + admission covariates',
          'eICU + admission + 72h clinical state (cluster)': 'eICU + 72 h clinical state'}
for _, r in eicu.iterrows():
    if r['model'] in labels:
        rows_unit.append((labels[r['model']], r['HR'], r['lo'], r['hi']))
        rows_sd.append((labels[r['model']], r['HR_per_SD'], r['lo_per_SD'], r['hi_per_SD']))

# ---- layout -----------------------------------------------------------------
# Every text column is anchored at a FIXED figure coordinate, so the annotation
# columns of one panel cannot run into the labels of the next panel.
FIG_W, FIG_H = 16.5, 4.9
LAYOUT = [
    dict(rows=rows_unit, title='A  Per 1 unit of r_w',
         rect=[0.160, 0.175, 0.226, 0.615], label_x=0.150, hr_x=0.400, mcid=False,
         xlim=(0.955, 1.065), ticks=[0.96, 0.98, 1.00, 1.02, 1.04, 1.06]),
    dict(rows=rows_sd, title='B  Per 1 standard deviation of r_w',
         rect=[0.650, 0.175, 0.226, 0.615], label_x=0.640, hr_x=0.895, mcid=True,
         xlim=(0.900, 1.235), ticks=[0.95, 1.00, 1.05, 1.10, 1.15, 1.20]),
]
fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300)
panel_axes = []
annotations = []          # every hand-placed text, for the collision check below
for spec in LAYOUT:
    ax = fig.add_axes(spec['rect'])
    x0, w = spec['rect'][0], spec['rect'][2]
    fx = lambda x, x0=x0, w=w: (x - x0) / w          # figure x -> axes fraction
    rows = spec['rows']
    ys = list(range(len(rows), 0, -1))
    for y, (lab, hr, lo, hi) in zip(ys, rows):
        primary = lab.startswith('MIMIC')
        color = '#1f4e79' if primary else '#5a5a5a'
        ax.plot([lo, hi], [y, y], color=color, lw=2.6 if primary else 2.0)
        ax.plot([hr], [y], 'o', color=color, ms=9 if primary else 7.5)
        annotations.append(ax.text(
            fx(spec['label_x']), y, lab, transform=ax.get_yaxis_transform(), ha='right',
            va='center', fontsize=11.5, color=color,
            fontweight='bold' if primary else 'normal'))
        annotations.append(ax.text(
            fx(spec['hr_x']), y, f'{hr:.2f} ({lo:.2f}-{hi:.2f})',
            transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize=11, color=color))
    ax.axvline(1.0, ls='--', lw=1.2, color='#333333')
    if spec['mcid']:
        ax.axvline(1.20, ls=':', lw=1.6, color='#b2182b')
        annotations.append(ax.text(1.196, len(rows) + 0.55, 'MCID 1.20', ha='right', va='center',
                                   fontsize=10.5, color='#b2182b'))
    ax.set_yticks(ys); ax.set_yticklabels([])
    ax.set_ylim(0.4, len(rows) + 0.85)
    ax.set_xlim(*spec['xlim'])
    ax.set_xticks(spec['ticks']); ax.set_xticklabels([f'{t:.2f}' for t in spec['ticks']])
    ax.set_xlabel('Hazard ratio (95% CI)', fontsize=12)
    ax.set_title(spec['title'], fontsize=13, fontweight='bold', pad=12)
    ax.tick_params(axis='x', labelsize=11)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    panel_axes.append(ax)


def overlap_report(fig, axes_list, annotations, tag):
    """Machine check on the hand-placed annotation columns: no annotation may
    overlap another annotation or any other text (e.g. a tick label), and no
    annotation may sit on a different panel's axes box. Matplotlib-managed
    tick labels are not checked against each other. Returns the collision count."""
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
                bad.append((a.get_text()[:26], annotations[j].get_text()[:26], 'annotation/annotation'))
        for k, o in enumerate(others):
            if ab[i].overlaps(ob[k]):
                bad.append((a.get_text()[:26], o.get_text()[:26], 'annotation/text'))
        for k, ax in enumerate(axes_list):
            if a.axes is ax:
                continue
            if ab[i].overlaps(ax.get_window_extent(renderer=rend)):
                bad.append((a.get_text()[:26], f'panel{k} box', 'annotation/panel'))
    print(f'OVERLAP CHECK [{tag}]: collisions = {len(bad)} '
          f'(annotations {len(annotations)}, other texts {len(others)})')
    for r in bad[:8]:
        print('    collision:', r)
    return len(bad)


n_collisions = overlap_report(fig, panel_axes, annotations, 'SF6')
if os.environ.get('FIG_LAYOUT_CHECK_ONLY') == '1':
    print('layout-check-only: SF6 NOT written')
else:
    fig.savefig(OUT, bbox_inches='tight')
    print('saved', OUT)
print('panel A rows:', [(r[0], round(r[1], 3), round(r[2], 3), round(r[3], 3)) for r in rows_unit])
print('panel B rows:', [(r[0], round(r[1], 3), round(r[2], 3), round(r[3], 3)) for r in rows_sd])
assert n_collisions == 0, 'SF6 layout has text collisions'

