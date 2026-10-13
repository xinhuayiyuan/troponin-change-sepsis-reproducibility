"""Fix trajectory figure bottom overlap (re-add axis title + note cleanly),
then merge trajectories (panel A) + forest plot (panel B) into fig3_composite.png."""
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

# --- 1) repair trajectory figure (skip if legacy source absent; reuse fig3a) ---
_src_path = os.path.join(FIG, 'fig2_trajectories.png')
if not os.path.exists(_src_path):
    fixed = Image.open(os.path.join(FIG, 'fig3a_trajectories.png')).convert('RGB')
    print('reusing existing fig3a_trajectories.png', fixed.size)
    imgA = fixed
    imgB = Image.open(os.path.join(FIG, 'fig3b_forest.png')).convert('RGB')
    W = 3000
    imgA = imgA.resize((W, int(imgA.height * W / imgA.width)), Image.LANCZOS)
    imgB = imgB.resize((W, int(imgB.height * W / imgB.width)), Image.LANCZOS)
    MARGIN_L, GAP = 190, 100
    canvas = Image.new('RGB', (W + MARGIN_L, imgA.height + GAP + imgB.height), 'white')
    canvas.paste(imgA, (MARGIN_L, 0))
    canvas.paste(imgB, (MARGIN_L, imgA.height + GAP))
    d = ImageDraw.Draw(canvas)
    f_lab = font(116, bold=True)
    d.text((20, 20), 'A', fill='#111111', font=f_lab)
    d.text((20, imgA.height + GAP + 20), 'B', fill='#111111', font=f_lab)
    out = os.path.join(FIG, 'fig3_composite.png')
    canvas.save(out, dpi=(300, 300))
    print('saved', out, canvas.size)
    raise SystemExit(0)
src = Image.open(_src_path).convert('RGB')
W0, H0 = src.size                      # 1150 x 820
crop = src.crop((0, 0, W0, 715))       # keep plot + tick labels; drop overlapped bottom
EXTRA = 175
fixed = Image.new('RGB', (W0, 715 + EXTRA), 'white')
fixed.paste(crop, (0, 0))
d = ImageDraw.Draw(fixed)
# axis title (centred, clear of tick labels)
t = 'Hours from ICU admission'
f_t = font(26, bold=True)
tw = d.textlength(t, font=f_t)
d.text(((W0 - tw) / 2, 728), t, fill='#111111', font=f_t)
# note (3 lines, well spaced)
note = [
    'Note: cTnT measured by itemid 51003 (MIMIC-IV); assay generation not specified in the database.',
    'The Stable-low trajectory ends at 60 h because no patient in this group had a measurement in the 66-72 h window,',
    'reflecting cessation of serial testing when values remain at the classification threshold (mean 2.5 vs 4.7 measurements in Rising; ST3).',
]
f_n = font(18)
y = 775
for ln in note:
    d.text((70, y), ln, fill='#222222', font=f_n)
    y += 32
fixed_path = os.path.join(FIG, 'fig3a_trajectories.png')
fixed.save(fixed_path)
print('saved', fixed_path, fixed.size)

# --- 2) composite A (trajectories) + B (forest) ---
imgA = fixed
imgB = Image.open(os.path.join(FIG, 'fig3b_forest.png')).convert('RGB')
W = 1500
imgA = imgA.resize((W, int(imgA.height * W / imgA.width)), Image.LANCZOS)
imgB = imgB.resize((W, int(imgB.height * W / imgB.width)), Image.LANCZOS)
MARGIN_L, GAP = 95, 50
canvas = Image.new('RGB', (W + MARGIN_L, imgA.height + GAP + imgB.height), 'white')
canvas.paste(imgA, (MARGIN_L, 0))
canvas.paste(imgB, (MARGIN_L, imgA.height + GAP))
d = ImageDraw.Draw(canvas)
f_lab = font(58, bold=True)
d.text((10, 10), 'A', fill='#111111', font=f_lab)
d.text((10, imgA.height + GAP + 10), 'B', fill='#111111', font=f_lab)
out = os.path.join(FIG, 'fig3_composite.png')
canvas.save(out, dpi=(300, 300))
print('saved', out, canvas.size)
