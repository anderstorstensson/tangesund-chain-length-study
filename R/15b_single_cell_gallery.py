# =============================================================================
# Single-cell Chaetoceros exemplars (Fig. S3), stage 2: detector check and strip.
# Keeps the first four candidates per depth from 14b with exactly one detector box,
# foreground (> 25 gray levels from the median background) outside the box dilated by
# 3 px below 2% of the box area, and a box margin >= 2 px. Composited at 600 dpi in
# Noto Sans.
#
# Reads:  outputs/tangesund/single_cell_work/{selection.csv, pngs/} (from 14b);
#         data/chains/chains_chaetoceros_yolo11s/weights/best.pt
# Writes: outputs/tangesund/15b_chaeto_single_cells.png, figures/FigS3.png,
#         single_cell_work/15b_candidate_scan.csv, 15b_selection_final.csv
# Deps:   ultralytics, torch, pillow, matplotlib, pandas
# Run:    python3 R/15b_single_cell_gallery.py   (from the repository root)
# =============================================================================
import os, shutil
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from PIL import Image
from ultralytics import YOLO

plt.rcParams["font.family"] = "Noto Sans"
WORK    = "outputs/tangesund/single_cell_work"
WEIGHTS = "data/chains/chains_chaetoceros_yolo11s/weights/best.pt"
OUT     = "outputs/tangesund/15b_chaeto_single_cells.png"
FIG_OUT = "outputs/tangesund/figures/FigS3.png"
MF      = 1/2.77                    # um per pixel
CONF, IOU = 0.25, 0.30              # production Chaetoceros detector settings
K       = 4                         # exemplars per layer
FG_DIFF, DILATE, MAX_CLUTTER, MIN_MARGIN = 25, 3, 0.02, 2   # frame-cleanliness thresholds
DPI     = 600

sel = pd.read_csv(f"{WORK}/selection.csv").sort_values(["layer", "rank"])
def imgpath(s, roi): return f"{WORK}/pngs/{s}/{s}_{roi:05d}.png"
model = YOLO(WEIGHTS)

def frame_check(path, box):
    """Clutter outside the box as a fraction of the box area, and the box's
    smallest distance to the frame edge in pixels."""
    im = np.asarray(Image.open(path).convert("L")).astype(float)
    h, w = im.shape
    x1, y1, x2, y2 = [int(round(v)) for v in box]
    fg = np.abs(im - np.median(im)) > FG_DIFF
    inside = np.zeros_like(fg)
    inside[max(0, y1 - DILATE):y2 + DILATE, max(0, x1 - DILATE):x2 + DILATE] = True
    clutter = (fg & ~inside).sum() / max(1, (x2 - x1) * (y2 - y1))
    margin = min(x1, y1, w - x2, h - y2)
    return clutter, margin

boxes, keep, scan = {}, [], []
for layer in ["surface_3m", "deep_16m"]:
    n = 0
    for _, r in sel[sel.layer == layer].iterrows():
        path = imgpath(r["sample"], r["roi_number"])
        res = model(path, conf=CONF, iou=IOU, verbose=False)[0]
        b = res.boxes.xyxy.cpu().numpy().tolist()
        key = f"{r['sample']}_{r['roi_number']:05d}"
        clutter, margin = frame_check(path, b[0]) if len(b) == 1 else (np.nan, np.nan)
        if len(b) != 1:
            status = f"{len(b)} boxes, rejected"
        elif clutter >= MAX_CLUTTER:
            status = f"other object in frame (clutter {clutter:.3f}), rejected"
        elif margin < MIN_MARGIN:
            status = f"box at frame edge (margin {margin} px), rejected"
        else:
            status = "kept" if n < K else "eligible, surplus"
        print(f"  {layer:11s} rank {int(r['rank']):2d} {key} {r['pc_esd']:.1f} um: {status}")
        scan.append(dict(layer=layer, rank=int(r["rank"]), key=key, n_boxes=len(b),
                         clutter=round(clutter, 4), margin=margin, status=status))
        if status == "kept":
            boxes[key] = b[0]; keep.append(dict(r, key=key, clutter=round(clutter, 4), margin=margin)); n += 1
    if n < K:
        raise SystemExit(f"{layer}: only {n} eligible candidates; raise k_candidates in stage 1")
pd.DataFrame(scan).to_csv(f"{WORK}/15b_candidate_scan.csv", index=False)
fin = pd.DataFrame(keep)
fin.to_csv(f"{WORK}/15b_selection_final.csv", index=False)

SURF, DEEP = "#c0392b", "#2471a3"
SCALE = 0.031; GAP = 14.0; LBL_IN = 1.18; LAB_H = 0.24
PADL, PADR, PADT, PADB = 0.10, 0.72, 0.28, 0.24
def dims_um(s, roi):
    w, h = Image.open(imgpath(s, roi)).size
    return w * MF, h * MF
def place(fig, fw, fh, x_in, ytop_in, wi, hi):
    return fig.add_axes([x_in / fw, 1 - (ytop_in + hi) / fh, wi / fw, hi / fh])
def draw(ax, key, s, roi, boxcol, lw=1.3):
    im = np.asarray(Image.open(imgpath(s, roi)).convert("L")); h, w = im.shape
    ax.imshow(im, cmap="gray", origin="upper", vmin=0, vmax=255, extent=[0, w*MF, h*MF, 0])
    x1, y1, x2, y2 = boxes[key]
    ax.add_patch(Rectangle((x1*MF, y1*MF), (x2-x1)*MF, (y2-y1)*MF, fill=False, ec=boxcol, lw=lw, zorder=3))
    ax.set_xlim(0, w*MF); ax.set_ylim(h*MF, 0); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
def scalebar(fig, fw, fh, x_in, y_in, um=20):
    L = um * SCALE
    fig.add_artist(Line2D([x_in/fw, (x_in+L)/fw], [y_in/fh, y_in/fh], color="k", lw=2.6, solid_capstyle="butt"))
    fig.text((x_in + L/2)/fw, (y_in + 0.05)/fh, f"{um} µm", ha="center", va="bottom", fontsize=8.5)
def band_w(df): return (sum(dims_um(r["sample"], r["roi_number"])[0] for _, r in df.iterrows()) + GAP * (len(df) - 1)) * SCALE
def band_h(df): return max(dims_um(r["sample"], r["roi_number"])[1] for _, r in df.iterrows()) * SCALE

bands = [fin[fin.layer == "surface_3m"].reset_index(drop=True), fin[fin.layer == "deep_16m"].reset_index(drop=True)]
cols = [SURF, DEEP]; labels = ["Surface\n3 m", "Deep\n16 m"]
bh = [band_h(b) for b in bands]; bw = [band_w(b) for b in bands]
GAP_ROW, SB_STRIP = 0.34, 0.36
fw = PADL + LBL_IN + max(bw) + PADR
fh = PADT + (LAB_H + bh[0]) + GAP_ROW + (LAB_H + bh[1]) + SB_STRIP + PADB
fig = plt.figure(figsize=(fw, fh), dpi=DPI)
xcrop = PADL + LBL_IN; y = PADT
for bi, df in enumerate(bands):
    h = bh[bi]; x = xcrop
    for _, r in df.iterrows():
        wi, hi = [d*SCALE for d in dims_um(r["sample"], r["roi_number"])]
        ax = place(fig, fw, fh, x, y + LAB_H + (h - hi)/2, wi, hi)
        draw(ax, r["key"], r["sample"], r["roi_number"], cols[bi])
        fig.text((x + wi/2)/fw, 1 - (y + LAB_H - 0.03)/fh, f"{r['pc_esd']:.1f} µm", ha="center", va="bottom",
                 fontsize=8.5, color=cols[bi])
        x += wi + GAP*SCALE
    fig.text((PADL + LBL_IN*0.46)/fw, 1 - (y + LAB_H + h/2)/fh, labels[bi], ha="center", va="center",
             color=cols[bi], fontsize=12, fontweight="bold", linespacing=1.4)
    y += LAB_H + h + GAP_ROW
scalebar(fig, fw, fh, fw - PADR - 20*SCALE - 0.05, PADB + SB_STRIP*0.3)
fig.savefig(OUT, dpi=DPI, facecolor="white"); plt.close(fig)
os.makedirs(os.path.dirname(FIG_OUT), exist_ok=True); shutil.copy(OUT, FIG_OUT)
print("wrote", OUT, "and", FIG_OUT, f"({fw:.1f}x{fh:.1f} in at {DPI} dpi)")
print("chosen per-cell ESD: surface", [round(v, 1) for v in bands[0].pc_esd], "deep", [round(v, 1) for v in bands[1].pc_esd])
