# =============================================================================
# Chaetoceros count error by per-cell size tertile on the 500 held-out test chains of
# Groves et al. (2026), with trend tests vs size and chain length.
#
# Reads:  data/chains_validation/chaetoceros_test/ (images + test_image_counts.csv),
#         data/chains/chains_chaetoceros_yolo11s/weights/best.pt
# Writes: outputs/tangesund/ 21_chaeto_size_bias.csv, 21_chaeto_size_validation_raw.csv
# Deps:   numpy, pandas, pillow, scipy, ultralytics
# Run:    python3 R/21_chaeto_size_validation.py   (from the repo root)
# =============================================================================
import os
import numpy as np
import pandas as pd
from PIL import Image
from ultralytics import YOLO

WEIGHTS = "data/chains/chains_chaetoceros_yolo11s/weights/best.pt"
TEST = "data/chains_validation/chaetoceros_test"
MICRON = 1 / 2.77  # um per pixel
CONF, IOU = 0.25, 0.30
OUT = "outputs/tangesund"

df = pd.read_csv(os.path.join(TEST, "test_image_counts.csv")).dropna(axis=1, how="all")
df = df[df["file_name"].apply(lambda f: os.path.exists(os.path.join(TEST, f)))].reset_index(drop=True)
files = [os.path.join(TEST, f) for f in df["file_name"]]
manual = df["cell_count"].to_numpy().astype(int)

# image geometry -> per-cell size proxies
W = np.array([Image.open(f).size[0] for f in files])
H = np.array([Image.open(f).size[1] for f in files])
area_px = (W * H).astype(float)
percell_area = area_px / manual
percell_esd = 2 * np.sqrt(percell_area / np.pi) * MICRON   # approx per-cell ESD (um)
mindim_um = np.minimum(W, H) * MICRON                       # ROI short axis ~ cell width (um)

# inference
model = YOLO(WEIGHTS)
pred = []
for i in range(0, len(files), 64):
    for r in model(files[i:i+64], conf=CONF, iou=IOU, verbose=False):
        pred.append(len(r.boxes))
pred = np.array(pred)
err = pred - manual

R = pd.DataFrame(dict(manual=manual, pred=pred, err=err, percell_esd=percell_esd,
                      mindim_um=mindim_um, W=W, H=H))
R.to_csv(os.path.join(OUT, "21_chaeto_size_validation_raw.csv"), index=False)

def block(mask=None, name="ALL"):
    e = err if mask is None else err[mask]
    return dict(bin=name, n=len(e), exact_pct=round(100*np.mean(e == 0), 1),
                within1_pct=round(100*np.mean(np.abs(e) <= 1), 1),
                mae=round(np.mean(np.abs(e)), 3), mean_bias=round(np.mean(e), 3))

print(f"CHAETOCEROS count validation | {len(files)} chains | conf={CONF} iou={IOU}")
print("Overall:", block())
print("  (paper Table S6 chaetoceros: 63.5% exact, MAE 0.57 -> sanity check)\n")

# stratify by per-cell size tertiles
q = np.quantile(percell_esd, [1/3, 2/3])
tert = np.where(percell_esd <= q[0], "small", np.where(percell_esd <= q[1], "mid", "large"))
rows = []
for lab in ["small", "mid", "large"]:
    m = tert == lab
    b = block(m, f"percell-{lab}")
    b["esd_um_range"] = f"{percell_esd[m].min():.1f}-{percell_esd[m].max():.1f}"
    b["esd_um_median"] = round(np.median(percell_esd[m]), 1)
    rows.append(b)
tab = pd.DataFrame(rows)
print("BY PER-CELL SIZE TERTILE (ESD proxy):")
print(tab.to_string(index=False))
tab.to_csv(os.path.join(OUT, "21_chaeto_size_bias.csv"), index=False)

# trend tests: does bias depend on size? on chain length?
from scipy.stats import spearmanr, pearsonr
rs_esd, ps_esd = spearmanr(percell_esd, err)
rs_cnt, ps_cnt = spearmanr(manual, err)
rp_esd, _ = pearsonr(percell_esd, err)
print(f"\nTREND  err vs per-cell ESD : Spearman rho={rs_esd:+.3f} (p={ps_esd:.2g}) | Pearson r={rp_esd:+.3f}")
print(f"Trend    err vs chain length : Spearman rho={rs_cnt:+.3f} (p={ps_cnt:.2g})")

# partial Spearman of err vs per-cell ESD given chain length: correlate rank residuals
from scipy.stats import rankdata
def rank_residual(y, x):
    ry, rx = rankdata(y), rankdata(x)
    slope, intercept = np.polyfit(rx, ry, 1)
    return ry - (slope * rx + intercept)
resid_err = rank_residual(err, manual)
resid_esd = rank_residual(percell_esd, manual)
rs_partial, ps_partial = pearsonr(resid_esd, resid_err)  # Pearson r on rank-residuals = partial Spearman
print(f"Partial  err vs per-cell ESD | chain length : partial Spearman rho={rs_partial:+.3f} (p={ps_partial:.2g})")

# bias in size windows around the Tangesund surface (~9.3 um) and deep (~7.4 um) per-cell ESD
for tag, lo, hi in [("deep~7.4um", 6.5, 8.3), ("surface~9.3um", 8.4, 10.5)]:
    m = (percell_esd >= lo) & (percell_esd < hi)
    if m.sum() > 5:
        print(f"  window {tag:14s} (proxy {lo}-{hi}um, n={m.sum()}): mean_bias={err[m].mean():+.3f}, MAE={np.abs(err[m]).mean():.3f}")
print("\nDONE. -> outputs/tangesund/21_chaeto_size_bias.csv (+ _raw.csv)")
