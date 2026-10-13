"""Compose the v3 submission figures from regenerated panels.

Figure 1 = cohort flow diagram (single panel, fig1_flow.png).
Figure 2 = fig2_km_m3.png (two panels + shared number-at-risk strip).
Figure 3 = trajectory panel (A) stacked over the 8-entry forest plot (B).
"""
from PIL import Image, ImageDraw, ImageFont
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(BASE, 'figures')

def font(sz, bold=False):
    name = 'arialbd.ttf' if bold else 'arial.ttf'
    for d in (os.path.join(os.environ.get('WINDIR', ''), 'Fonts'),
              '/usr/share/fonts/truetype/msttcorefonts',
              '/Library/Fonts'):
        p = os.path.join(d, name) if d else ''
        if p and os.path.exists(p):
            return ImageFont.truetype(p, sz)
    return ImageFont.load_default()

# ---- Figure 3: (A) trajectories over (B) forest plot ----
imgA = Image.open(os.path.join(FIG, 'fig3a_trajectories.png')).convert('RGB')
imgB = Image.open(os.path.join(FIG, 'fig3b_forest8.png')).convert('RGB')

W = 3000
imgA = imgA.resize((W, int(imgA.height * W / imgA.width)), Image.LANCZOS)
imgB = imgB.resize((W, int(imgB.height * W / imgB.width)), Image.LANCZOS)

MARGIN_L, GAP = 190, 110
canvas = Image.new('RGB', (W + MARGIN_L, imgA.height + GAP + imgB.height), 'white')
canvas.paste(imgA, (MARGIN_L, 0))
canvas.paste(imgB, (MARGIN_L, imgA.height + GAP))

d = ImageDraw.Draw(canvas)
f_lab = font(110, bold=True)
d.text((20, 20), 'a', fill='#111111', font=f_lab)
d.text((20, imgA.height + GAP + 20), 'b', fill='#111111', font=f_lab)

out = os.path.join(FIG, 'fig3_composite.png')
canvas.save(out, dpi=(300, 300))
print('saved', out, canvas.size)

# ---- tidy: retire the superseded composite drafts ----
for stale in ('fig1_composite.png', 'fig2_composite.png'):
    p = os.path.join(FIG, stale)
    if os.path.exists(p):
        os.remove(p)
        print('removed superseded', stale)

# ---- Figure 1 stays as the flow diagram alone ----
print('Figure 1 =', os.path.join(FIG, 'fig1_flow.png'))
print('Figure 2 =', os.path.join(FIG, 'fig2_km_m3.png'))
