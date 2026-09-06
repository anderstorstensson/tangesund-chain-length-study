# =============================================================================
# Chaetoceros depth gallery, stage 2: YOLO recount of the ROIs selected by 14 and the
# composite at a common micron scale (Fig. 4d).
#
# Reads:  outputs/tangesund/chaeto_gallery_work/{selection.csv, pngs/} (from 14);
#         data/chains/chains_chaetoceros_yolo11s/weights/best.pt
# Writes: outputs/tangesund/15_chaeto_gallery_inset.png (embedded in Fig. 4 by 16),
#         15_chaeto_single_cells.png (unused; Fig. S3 comes from 15b),
#         chaeto_gallery_work/yolo_boxes.json
# Deps:   ultralytics, torch, pillow, matplotlib, pandas
# Run:    python3 R/15_chaeto_gallery.py   (from the repository root)
# =============================================================================
import os, glob, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
from ultralytics import YOLO

WORK    = "outputs/tangesund/chaeto_gallery_work"
WEIGHTS = "data/chains/chains_chaetoceros_yolo11s/weights/best.pt"
OUT        = "outputs/tangesund/15_chaeto_gallery_inset.png"    # Fig. 4 panel (d), embedded by 16_percell_size.R
OUT_SINGLE = "outputs/tangesund/15_chaeto_single_cells.png"      # unused; Fig. S3 comes from 15b
QUAD_ASPECT = 1.30            # width/height of the Fig. 4d quadrant
MF      = 1/2.77                      # um per pixel (project micron_factor)
CONF, IOU = 0.25, 0.30               # production Chaetoceros detector settings

sel = pd.read_csv(f"{WORK}/selection.csv")
sel["key"] = sel.apply(lambda r: f"{r['sample']}_{r['roi_number']:05d}.png", axis=1)
def imgpath(s, roi): return f"{WORK}/pngs/{s}/{s}_{roi:05d}.png"

# ---- YOLO: count cells + record boxes, then verify vs stored count ----------
model = YOLO(WEIGHTS)
boxes = {}
for _, r in sel.iterrows():
    res = model(imgpath(r["sample"], r["roi_number"]), conf=CONF, iou=IOU, verbose=False)[0]
    boxes[r["key"]] = res.boxes.xyxy.cpu().numpy().tolist()
sel["yolo_n"] = sel["key"].map(lambda k: len(boxes[k]))
# boxes only on crops whose recount matches cell_count_resolved; a mismatch does not abort
BOXED = set(sel.loc[sel.yolo_n == sel.ncell, "key"])
agree = len(BOXED)
print(f"YOLO re-count agreement with stored cell_count_resolved: {agree}/{len(sel)}")
mism = sel[sel.yolo_n != sel.ncell]
if len(mism):
    print("  boxes suppressed on miscounted crops:\n",
          mism[["layer", "key", "ncell", "yolo_n"]].to_string(index=False))

# ---- layout: images placed at a common micron scale, sized to fit (no letterbox)
SURF, DEEP = "#c0392b", "#2471a3"        # red = surface, blue = deep (boxes + labels)
SCALE  = 0.031                            # inch per um (common physical scale, all crops)
GAP    = 14.0                             # horizontal gap between crops (um)
LBL_IN = 1.18                             # left column for row labels (inch)
LAB_H  = 0.24                             # strip above each crop for its text label (inch)
PADL, PADR, PADT, PADB = 0.10, 0.72, 0.28, 0.24   # wide PADR keeps end-crop labels from clipping
from matplotlib.lines import Line2D

def dims_um(s, roi):
    w, h = Image.open(imgpath(s, roi)).size
    return w * MF, h * MF

def place(fig, fw, fh, x_in, ytop_in, wi, hi):
    """add_axes from an inch box measured from the top-left."""
    return fig.add_axes([x_in / fw, 1 - (ytop_in + hi) / fh, wi / fw, hi / fh])

def draw(ax, key, s, roi, boxcol, lw=1.3, show_box=True):
    im = np.asarray(Image.open(imgpath(s, roi)).convert("L")); h, w = im.shape
    ax.imshow(im, cmap="gray", origin="upper", vmin=0, vmax=255, extent=[0, w*MF, h*MF, 0])
    if show_box and key in BOXED:                     # only draw boxes where the recount matched
        for x1, y1, x2, y2 in boxes[key]:
            ax.add_patch(Rectangle((x1*MF, y1*MF), (x2-x1)*MF, (y2-y1)*MF,
                         fill=False, ec=boxcol, lw=lw, zorder=3))
    ax.set_xlim(0, w*MF); ax.set_ylim(h*MF, 0)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

def scalebar(fig, fw, fh, x_in, y_in, um=20):
    L = um * SCALE
    fig.add_artist(Line2D([x_in/fw, (x_in+L)/fw], [y_in/fh, y_in/fh],
                          color="k", lw=2.6, solid_capstyle="butt"))
    fig.text((x_in + L/2)/fw, (y_in + 0.05)/fh, f"{um} µm", ha="center", va="bottom", fontsize=8.5)

def lay(name): return sel[sel.layer == name].reset_index(drop=True)

def pad_to_aspect(path, target, top_align=False):
    """Pad a PNG with white to a width/height aspect; top_align keeps content at the top."""
    im = Image.open(path).convert("RGB"); w, h = im.size
    if w / h > target: W, H = w, round(w / target)
    else:              W, H = round(h * target), h
    canvas = Image.new("RGB", (W, H), "white")
    oy = 0 if top_align else (H - h) // 2
    canvas.paste(im, ((W - w) // 2, oy)); canvas.save(path)

def _band_w(df):
    return (sum(dims_um(df.iloc[i]["sample"], df.iloc[i]["roi_number"])[0] for i in range(len(df)))
            + GAP * (len(df) - 1)) * SCALE
def _band_h(df):
    return max(dims_um(df.iloc[i]["sample"], df.iloc[i]["roi_number"])[1] for i in range(len(df))) * SCALE

# ---- chain gallery: surface stacked 2x2 (wider crops), deep in one row of 4 ----
# the Fig. 4 panel letter is added in R (16)
def build_chains(out, pad_aspect=None, tag=None):
    surf = lay("surface_3m").iloc[:4].reset_index(drop=True)
    deep = lay("deep_16m").iloc[:4].reset_index(drop=True)
    bands = [surf.iloc[0:2].reset_index(drop=True),      # surface row 1
             surf.iloc[2:4].reset_index(drop=True),      # surface row 2
             deep]                                        # deep row (4)
    cols  = [SURF, SURF, DEEP]
    band_h = [_band_h(b) for b in bands]; band_w = [_band_w(b) for b in bands]
    TAG_H    = 0.30 if tag else 0.0                       # small strip for the baked panel letter
    GAP_SURF, GAP_DEEP, SB_STRIP = 0.26, 0.44, 0.34      # gap between surface rows / before deep / scale-bar strip
    gaps = [GAP_SURF, GAP_DEEP, 0.0]
    fw = PADL + LBL_IN + max(band_w) + PADR
    fh = PADT + TAG_H + sum(LAB_H + h for h in band_h) + GAP_SURF + GAP_DEEP + SB_STRIP + PADB
    fig = plt.figure(figsize=(fw, fh), dpi=200)
    xcrop = PADL + LBL_IN                                 # left edge of the crops
    # optional baked letter, ~1.25x size because the inset is downscaled ~0.8x in Fig. 4
    if tag:
        fig.text(xcrop/fw, 1 - 0.05/fh, tag, ha="left", va="top",
                 fontsize=19, fontweight="bold", color="#000")
    y = PADT + TAG_H
    surf_top = y + LAB_H; surf_bot = None; deep_top = deep_h = None
    for bi, df in enumerate(bands):
        h = band_h[bi]; x = xcrop
        for i in range(len(df)):
            r = df.iloc[i]; wi, hi = [d*SCALE for d in dims_um(r["sample"], r["roi_number"])]
            ax = place(fig, fw, fh, x, y + LAB_H + (h - hi)/2, wi, hi)
            draw(ax, r["key"], r["sample"], r["roi_number"], cols[bi])
            fig.text((x + wi/2)/fw, 1 - (y + LAB_H - 0.03)/fh,
                     f"{int(r['ncell'])} cells, {r['pc_esd']:.1f} µm/cell",
                     ha="center", va="bottom", fontsize=8.5, color="#222")
            x += wi + GAP*SCALE
        if bi == 1: surf_bot = y + LAB_H + h
        if bi == 2: deep_top, deep_h = y + LAB_H, h
        y += LAB_H + h + gaps[bi]
    # group labels (Surface spans the two surface rows; Deep on its row)
    fig.text((PADL + LBL_IN*0.46)/fw, 1 - (surf_top + (surf_bot - surf_top)/2)/fh,
             "Surface\n3 m", ha="center", va="center", color=SURF, fontsize=12, fontweight="bold", linespacing=1.4)
    fig.text((PADL + LBL_IN*0.46)/fw, 1 - (deep_top + deep_h/2)/fh,
             "Deep\n16 m", ha="center", va="center", color=DEEP, fontsize=12, fontweight="bold", linespacing=1.4)
    scalebar(fig, fw, fh, fw - PADR - 20*SCALE - 0.05, PADB + SB_STRIP*0.4)
    fig.savefig(out, dpi=200, facecolor="white"); plt.close(fig)
    print(f"  gallery natural size {fw:.2f}x{fh:.2f} in (aspect {fw/fh:.2f})")
    if pad_aspect: pad_to_aspect(out, pad_aspect, top_align=True)
    print("wrote", out, f"({fw:.1f}x{fh:.1f} in)")

# ---- single-cell strip (unused; Fig. S3 comes from 15b) -----------------------
def build_singles(out):
    bands  = [lay("single_surface_3m").reset_index(drop=True),
              lay("single_deep_16m").reset_index(drop=True)]
    cols   = [SURF, DEEP]
    labels = ["Surface\n3 m", "Deep\n16 m"]
    band_h = [_band_h(b) for b in bands]
    band_w = [_band_w(b) for b in bands]
    GAP_ROW, SB_STRIP = 0.34, 0.36
    fw = PADL + LBL_IN + max(band_w) + PADR
    fh = PADT + (LAB_H + band_h[0]) + GAP_ROW + (LAB_H + band_h[1]) + SB_STRIP + PADB
    fig = plt.figure(figsize=(fw, fh), dpi=200)
    xcrop = PADL + LBL_IN
    y = PADT
    for bi, df in enumerate(bands):
        h = band_h[bi]; x = xcrop
        for i in range(len(df)):
            r = df.iloc[i]; wi, hi = [d*SCALE for d in dims_um(r["sample"], r["roi_number"])]
            ax = place(fig, fw, fh, x, y + LAB_H + (h - hi)/2, wi, hi)
            draw(ax, r["key"], r["sample"], r["roi_number"], cols[bi])   # boxes where the recount matched
            fig.text((x + wi/2)/fw, 1 - (y + LAB_H - 0.03)/fh, f"{r['pc_esd']:.1f} µm",
                     ha="center", va="bottom", fontsize=8.5, color=cols[bi])
            x += wi + GAP*SCALE
        fig.text((PADL + LBL_IN*0.46)/fw, 1 - (y + LAB_H + h/2)/fh, labels[bi],
                 ha="center", va="center", color=cols[bi], fontsize=12, fontweight="bold", linespacing=1.4)
        y += LAB_H + h + GAP_ROW
    scalebar(fig, fw, fh, fw - PADR - 20*SCALE - 0.05, PADB + SB_STRIP*0.3)
    fig.savefig(out, dpi=200, facecolor="white"); plt.close(fig)
    print("wrote", out, f"({fw:.1f}x{fh:.1f} in)")

build_chains(OUT, pad_aspect=QUAD_ASPECT)             # surface 2x2, deep 1x4; letter added in R
build_singles(OUT_SINGLE)                             # unused
json.dump(boxes, open(f"{WORK}/yolo_boxes.json", "w"))
