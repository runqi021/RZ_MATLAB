"""
Combine per-ROI calcium videos (row 1..n) with a behavior row (row n+1)
into a single stacked video.

Row 1..n : existing Fmovie_perROI/ROI_NN_movie.mp4  (patch + dFF trace)
Row n+1  : behavior video (left)  +  breath trace with yellow time line (right)

Alignment:
  - calcium F has TossFrames dropped at the start (default 30).
  - behavior AVI + breath trace are on the raw timebase (6000 frames).
  - We drop the same TossFrames from the behavior + breath so all rows share
    the same post-toss frame index.

Originals are read only; nothing is modified in place.
"""
import argparse, json
import cv2
import numpy as np
from pathlib import Path
from scipy.io import loadmat

# ---- config ---------------------------------------------------------------
parser = argparse.ArgumentParser(description="Combine per-ROI + behavior + breath into one video.")
parser.add_argument("--fov",  default=r"D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing\roi5_1400-1230-0_x4.4_15lp_6000f_00001",
                    help="FOV folder containing the AVI, breath_peak_data.mat and Fmovie_perROI/")
parser.add_argument("--rois", default="3,5,7", help="Comma-separated ROI indices, e.g. 4,10,11")
parser.add_argument("--toss", type=int, default=30, help="Frames dropped from calcium at start")
args = parser.parse_args()

fov       = Path(args.fov)
roi_dir   = fov / "Fmovie_perROI"
roi_list  = [int(x) for x in args.rois.split(",") if x.strip()]
toss      = args.toss

# Auto-discover AVI + breath mat (prefer ones with matching DLC breath_peak_data suffix)
avi_hits = sorted(fov.glob("Basler_*.avi"))
assert avi_hits, f"No Basler_*.avi in {fov}"
avi_path = avi_hits[0]

bm_hits = sorted(fov.glob("*breath_peak_data.mat"))
assert bm_hits, f"No *breath_peak_data.mat in {fov}"
breath_mat = bm_hits[0]

out_path = roi_dir / f"combined_ROI_{'_'.join(map(str, roi_list))}_with_behavior.mp4"
print(f"  AVI:    {avi_path.name}")
print(f"  breath: {breath_mat.name}")

# ---- open ROI videos ------------------------------------------------------
caps, W, H, N, fps = [], [], [], [], []
for r in roi_list:
    fn = roi_dir / f"ROI_{r:02d}_movie.mp4"
    assert fn.exists(), f"Missing {fn}"
    c = cv2.VideoCapture(str(fn))
    caps.append(c)
    W.append(int(c.get(cv2.CAP_PROP_FRAME_WIDTH)))
    H.append(int(c.get(cv2.CAP_PROP_FRAME_HEIGHT)))
    N.append(int(c.get(cv2.CAP_PROP_FRAME_COUNT)))
    fps.append(c.get(cv2.CAP_PROP_FPS))
    print(f"  ROI {r:02d}: {W[-1]}x{H[-1]}, {N[-1]} frames @ {fps[-1]:.2f} fps")

rowW, rowH = W[0], H[0]
assert all(w == rowW and h == rowH for w, h in zip(W, H)), "ROI videos differ in size"
assert len(set(round(f, 2) for f in fps)) == 1, "Frame rates differ across ROI videos"
fps0 = fps[0]

# ROI-row dimensions (from generate_perROI_videos.m)
patch_sz   = rowH                        # 300
trace_W    = rowW - patch_sz - 4         # 1574 - 300 - 4 = 1270
trace_DrawW= trace_W - 70                # scaleMargin=70 → 1200 drawable

# ---- open behavior AVI ----------------------------------------------------
cap_b = cv2.VideoCapture(str(avi_path))
assert cap_b.isOpened(), f"Cannot open {avi_path}"
bW = int(cap_b.get(cv2.CAP_PROP_FRAME_WIDTH))
bH = int(cap_b.get(cv2.CAP_PROP_FRAME_HEIGHT))
bN = int(cap_b.get(cv2.CAP_PROP_FRAME_COUNT))
bFPS = cap_b.get(cv2.CAP_PROP_FPS)
print(f"  Behavior AVI: {bW}x{bH}, {bN} frames @ {bFPS:.2f} fps")

# Drop first `toss` frames to align with calcium
for _ in range(toss):
    cap_b.read()

# Optional crop (from pick_behavior_crop.py)
crop_file = fov / "behavior_crop.json"
if crop_file.exists():
    crop = json.loads(crop_file.read_text())
    cx, cy, cw, ch = crop["x"], crop["y"], crop["w"], crop["h"]
    print(f"  Cropping behavior to [x={cx}, y={cy}, w={cw}, h={ch}]")
else:
    cx, cy, cw, ch = 0, 0, bW, bH

# Letterbox cropped behavior into patch_sz x patch_sz
scale_b = min(patch_sz / cw, patch_sz / ch)
bW2, bH2 = int(round(cw * scale_b)), int(round(ch * scale_b))
bx0 = (patch_sz - bW2) // 2
by0 = (patch_sz - bH2) // 2

# ---- load + trim breath trace --------------------------------------------
bm = loadmat(str(breath_mat))
breath_full = np.asarray(bm["breath"], dtype=np.float64).ravel()
t_breath_full = np.asarray(bm["t_breath"], dtype=np.float64).ravel()
# CAP_PROP_FRAME_COUNT can be off by 1; trust whichever is smaller.
breath = breath_full[toss:]
t_breath = t_breath_full[toss:] - t_breath_full[toss]  # re-base to 0
T = min(min(N), breath.size, bN - toss)
print(f"  Using T = {T} frames after toss (calcium={min(N)}, breath={breath.size}, beh={bN-toss})")
breath = breath[:T]
t_breath = t_breath[:T]

# ---- pre-render breath panel background ----------------------------------
# Match dFF panel style: black bg, white thick trace, small right-margin scale bar.
panel = np.zeros((rowH, trace_W, 3), dtype=np.uint8)

y_margin = int(round(rowH * 0.08))
y_lo = rowH - y_margin
y_hi = y_margin + 1
plot_H = y_lo - y_hi

b_min, b_max = float(np.min(breath)), float(np.max(breath))
if b_max == b_min:
    b_max = b_min + 1.0
b_range = b_max - b_min
x_px = np.linspace(0, trace_DrawW - 1, T)
y_px = np.round(y_lo - (breath - b_min) / b_range * plot_H).astype(int)
y_px = np.clip(y_px, 0, rowH - 1)

pts = np.stack([np.round(x_px).astype(int), y_px], axis=1).reshape(-1, 1, 2)
cv2.polylines(panel, [pts], isClosed=False, color=(255, 255, 255),
              thickness=2, lineType=cv2.LINE_AA)

cv2.putText(panel, "breath", (8, 22),
            cv2.FONT_HERSHEY_SIMPLEX, 0.6, (180, 180, 180), 1, cv2.LINE_AA)

# ---- writer ---------------------------------------------------------------
nRows = len(roi_list) + 1
canvasW = rowW
canvasH = rowH * nRows

fourcc = cv2.VideoWriter_fourcc(*"mp4v")
vw = cv2.VideoWriter(str(out_path), fourcc, fps0, (canvasW, canvasH))
assert vw.isOpened(), f"Could not open writer {out_path}"

print(f"  Writing: {out_path} ({canvasW}x{canvasH}, {T} frames @ {fps0:.2f} fps)")

for f in range(T):
    canvas = np.zeros((canvasH, canvasW, 3), dtype=np.uint8)

    # row 0: behavior + breath panel
    ok, bfr = cap_b.read()
    if ok:
        bfr_c = bfr[cy:cy + ch, cx:cx + cw]
        bfr2 = cv2.resize(bfr_c, (bW2, bH2), interpolation=cv2.INTER_AREA)
        canvas[by0:by0 + bH2, bx0:bx0 + bW2] = bfr2

    bp = panel.copy()
    x_line = int(np.clip(round(x_px[f]), 1, trace_W - 2))
    bp[:, x_line - 1:x_line + 2] = (0, 255, 255)   # BGR yellow
    cv2.putText(bp, f"{t_breath[f]:.1f}s", (6, 48),
                cv2.FONT_HERSHEY_SIMPLEX, 0.8, (255, 255, 255), 2, cv2.LINE_AA)
    canvas[0:rowH, patch_sz + 4:patch_sz + 4 + trace_W] = bp

    # rows 1..n: existing ROI frames
    for i, cap in enumerate(caps):
        ok, fr = cap.read()
        if not ok:
            break
        y0 = (i + 1) * rowH
        canvas[y0:y0 + rowH, :rowW] = fr

    vw.write(canvas)

vw.release()
for c in caps:
    c.release()
cap_b.release()
print(f"  Done -> {out_path}")
