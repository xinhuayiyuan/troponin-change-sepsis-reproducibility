"""Redraw the STROBE cohort flow diagram: larger fonts, darker colours, no overlaps."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, 'figures', 'fig1_flow.png')

NAVY = '#0e2a47'      # main border / heading
DARK = '#111111'      # main text
GREY = '#3d3d3d'      # exclusion border
GREY_T = '#1f1f1f'    # exclusion text

fig, ax = plt.subplots(figsize=(9.2, 10.0), dpi=300)
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis('off')

ax.text(50, 97.5, 'Cohort assembly', ha='center', va='top',
        fontsize=23, fontweight='bold', color=NAVY)
ax.text(50, 92.8, 'STROBE flow diagram', ha='center', va='top',
        fontsize=16, color=DARK)

def box(x, y, w, h, lines, border, tcolor, dashed=False, bold_last=True, fs=14.5):
    b = FancyBboxPatch((x - w/2, y - h/2), w, h,
                       boxstyle='round,pad=0.4,rounding_size=1.0',
                       linewidth=1.8, edgecolor=border, facecolor='white',
                       linestyle='--' if dashed else '-')
    ax.add_patch(b)
    n = len(lines)
    for i, ln in enumerate(lines):
        yy = y + (n - 1) * 1.5 - i * 3.0
        bold = bold_last and (i == n - 1)
        ax.text(x, yy, ln, ha='center', va='center',
                fontsize=fs + (1.5 if bold else 0),
                fontweight='bold' if bold else 'normal', color=tcolor)

def arrow(x1, y1, x2, y2, style='-', color=DARK):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',
                                 mutation_scale=22, linewidth=1.8,
                                 color=color, linestyle=style))

MX, MW = 29, 50   # main column centre / width
EX, EW = 80, 36   # exclusion column centre / width

box(MX, 82.5, MW, 8.5, ['Adults with sepsis-3,', 'first ICU admission', 'n = 24,078'], NAVY, DARK)
box(EX, 71.5, EW, 8.0, ['Excluded: no troponin', 'measurement', '(n = 17,817)'], GREY, GREY_T, dashed=True, bold_last=False, fs=13.5)
box(MX, 60.5, MW, 11.5, ['Any cardiac troponin T', 'measurement within 72 h', 'of ICU admission', 'n = 6,261'], NAVY, DARK)
box(EX, 49.5, EW, 8.0, ['Excluded: fewer than', '2 measurements', '(n = 1,873)'], GREY, GREY_T, dashed=True, bold_last=False, fs=13.5)
box(MX, 38.5, MW, 8.5, ['At least 2 troponin', 'measurements within 72 h', 'n = 4,388'], NAVY, DARK)
box(78, 27.5, 40, 12.5, ['Excluded (n = 1,184; 27.0%):',
                         'no 0-6 h baseline, n = 1,163',
                         'group unclassifiable, n = 21'],
    GREY, GREY_T, dashed=True, bold_last=False, fs=12)
box(MX, 16.5, MW, 8.5, ['Primary analysis cohort', 'n = 3,204'], NAVY, DARK)

arrow(MX, 78.2, MX, 65.6)
arrow(MX, 55.4, MX, 43.1)
arrow(MX, 33.9, MX, 21.1)
for y in (71.5, 49.5, 27.5):
    ax.plot([MX + MW/2 + 1.5, EX - EW/2 - 1.5], [y, y], ls=':', lw=1.6, color=GREY)

note = ('Note. Sepsis-3 defined per the Third International Consensus Definitions. Cardiac troponin T measured by itemid 51003.\n'
        'Dynamic classification used predefined criteria (Rising, Falling, Stable-high, Stable-low). "No 0-6 h baseline" denotes patients\n'
        'without any cTnT measurement in the 0-6 h baseline window, in whom the baseline value is undefined; "group unclassifiable"\n'
        'denotes a defined baseline but serial values satisfying none of the four group definitions.')
ax.text(50, 6.0, note, ha='center', va='top', fontsize=11.5, color=DARK, linespacing=1.6)

fig.savefig(OUT, bbox_inches='tight', facecolor='white')
print('saved', OUT)
