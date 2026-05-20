"""Generate ONE test overlay video: vL1 (blue) + vR1 (red) cross markers
drawn on the raw Basler frames. Uses OpenCV which handles FFV1 natively.

Saved into subfolder `representative_whisker_label_video/` next to source.
"""
import os
import glob
import re
import csv
import time
import cv2
import numpy as np

# ===== USER SETTINGS =====
folder_path   = r"C:\260407_KA_electro_dorsal_whisking"
out_subfolder = "representative_whisker_label_video"
run_to_test   = 1              # 1..30

L_regex       = re.compile(r"DLC_Resnet50_whisker_dorsalApr")
R_regex       = re.compile(r"DLC_Resnet50_whisker_dorsal_RApr")
L_point       = "vL1"
R_point       = "vR1"
pmin          = 0.5

half_len      = 8              # cross arm length (px)
thickness     = 2              # cv2.line thickness (px)

# Colors match the analysis plots (blue/red), expressed as BGR for OpenCV
color_L_bgr   = (217, 115,   0)   # vL1 blue
color_R_bgr   = ( 25,  25, 217)   # vR1 red

rotate_deg    = 150            # CCW rotation applied to frame + markers

# Crop applied AFTER rotation, to the rotated canvas. Set None for no crop.
# Use whisker_select_crop.py to pick these interactively.
crop_xywh     = [114, 402, 646, 238]  # x, y, w, h (from whisker_select_crop.py)
# =========================


def load_point(csv_path, point_name):
    """Return np arrays (x, y, p) for the named DLC bodypart."""
    with open(csv_path, "r") as f:
        rows = list(csv.reader(f))
    # rows[0]=scorer, rows[1]=bodyparts, rows[2]=coords(x/y/likelihood)
    bps = []
    for c in rows[1][1:]:
        c = c.strip()
        if not bps or bps[-1] != c:
            bps.append(c)
    assert point_name in bps, f"{point_name} not in {bps} ({csv_path})"
    pid = bps.index(point_name)                # 0-indexed within unique
    cx = 1 + 3 * pid
    cy = 2 + 3 * pid
    cp = 3 + 3 * pid
    data = rows[3:]
    x = np.array([float(r[cx]) for r in data])
    y = np.array([float(r[cy]) for r in data])
    p = np.array([float(r[cp]) for r in data])
    return x, y, p


def draw_cross(frame, x, y, hl, th, color_bgr):
    if np.isnan(x) or np.isnan(y):
        return
    xi = int(round(x))
    yi = int(round(y))
    h, w = frame.shape[:2]
    if xi < 0 or xi >= w or yi < 0 or yi >= h:
        return
    cv2.line(frame, (xi - hl, yi), (xi + hl, yi), color_bgr, thickness=th)
    cv2.line(frame, (xi, yi - hl), (xi, yi + hl), color_bgr, thickness=th)


# ---------- 1. Locate source AVI ----------
raw_avis = [
    a for a in glob.glob(os.path.join(folder_path, f"*run{run_to_test:03d}.avi"))
    if "DLC_Resnet50" not in os.path.basename(a)
]
assert raw_avis, f"No raw AVI for run{run_to_test:03d}"
src_avi = raw_avis[0]
print(f"Source: {src_avi}")

# ---------- 2. Locate CSVs ----------
all_csvs = glob.glob(os.path.join(folder_path, f"*run{run_to_test:03d}*.csv"))
L_csv = next(
    (c for c in all_csvs
     if L_regex.search(os.path.basename(c)) and "dorsal_R" not in os.path.basename(c)),
    None
)
R_csv = next(
    (c for c in all_csvs if R_regex.search(os.path.basename(c))),
    None
)
assert L_csv, f"No L CSV for run{run_to_test:03d}"
assert R_csv, f"No R CSV for run{run_to_test:03d}"
print(f"L CSV:  {L_csv}")
print(f"R CSV:  {R_csv}")

# ---------- 3. Load coordinates ----------
xL, yL, pL = load_point(L_csv, L_point)
xR, yR, pR = load_point(R_csv, R_point)
n_csv = min(len(xL), len(xR))
print(f"CSV frames: L={len(xL)}, R={len(xR)}, common={n_csv}")

# ---------- 4. Output path ----------
out_dir = os.path.join(folder_path, out_subfolder)
os.makedirs(out_dir, exist_ok=True)
stem = os.path.splitext(os.path.basename(src_avi))[0]
out_mp4 = os.path.join(out_dir, f"{stem}_vL1_vR1_cross.mp4")

# ---------- 5. Open input, build rotation, open output ----------
cap = cv2.VideoCapture(src_avi)
assert cap.isOpened(), f"Cannot open {src_avi}"
fps      = cap.get(cv2.CAP_PROP_FPS)
w        = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
h        = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
n_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
print(f"Video (raw): {w}x{h}, {fps:.2f} fps, {n_frames} frames")

# Rotation matrix around the original frame center, CCW by rotate_deg.
# Then expand the canvas to fit the full rotated image and shift the
# matrix so nothing gets clipped.
M = cv2.getRotationMatrix2D((w / 2.0, h / 2.0), rotate_deg, 1.0)
cos_ = abs(M[0, 0])
sin_ = abs(M[0, 1])
new_w = int(h * sin_ + w * cos_)
new_h = int(h * cos_ + w * sin_)
M[0, 2] += (new_w / 2.0) - (w / 2.0)
M[1, 2] += (new_h / 2.0) - (h / 2.0)
print(f"Rotated canvas: {new_w}x{new_h} (+black padding)")

def rotate_xy(x, y):
    """Apply the same affine M to a (x, y) point, return rotated (xr, yr)."""
    xr = M[0, 0] * x + M[0, 1] * y + M[0, 2]
    yr = M[1, 0] * x + M[1, 1] * y + M[1, 2]
    return xr, yr

if crop_xywh is not None:
    cx, cy, cw, ch = [int(v) for v in crop_xywh]
    # Clamp to rotated canvas
    cx = max(0, min(cx, new_w - 1))
    cy = max(0, min(cy, new_h - 1))
    cw = max(1, min(cw, new_w - cx))
    ch = max(1, min(ch, new_h - cy))
    out_w, out_h = cw, ch
    print(f"Crop: x={cx}, y={cy}, w={cw}, h={ch} (output {out_w}x{out_h})")
else:
    cx = cy = 0
    out_w, out_h = new_w, new_h

fourcc = cv2.VideoWriter_fourcc(*"mp4v")
vw = cv2.VideoWriter(out_mp4, fourcc, fps, (out_w, out_h))
assert vw.isOpened(), f"Cannot open writer for {out_mp4}"

# ---------- 6. Frame loop ----------
t0 = time.time()
k = 0
while True:
    ok, frame = cap.read()
    if not ok or k >= n_csv:
        break
    # Rotate frame into expanded canvas, black padding
    frame_rot = cv2.warpAffine(
        frame, M, (new_w, new_h),
        flags=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=(0, 0, 0),
    )
    # Crop (shifts pixel origin to (cx, cy); markers must shift too)
    if crop_xywh is not None:
        frame_rot = frame_rot[cy:cy + ch, cx:cx + cw]

    # Rotate marker coordinates into the same frame, then offset for crop
    if pL[k] >= pmin:
        xLr, yLr = rotate_xy(xL[k], yL[k])
        draw_cross(frame_rot, xLr - cx, yLr - cy, half_len, thickness, color_L_bgr)
    if pR[k] >= pmin:
        xRr, yRr = rotate_xy(xR[k], yR[k])
        draw_cross(frame_rot, xRr - cx, yRr - cy, half_len, thickness, color_R_bgr)
    vw.write(frame_rot)
    k += 1
    if k % 500 == 0:
        print(f"  frame {k}  ({time.time() - t0:.1f}s)")

cap.release()
vw.release()
print(f"\nSaved: {out_mp4}")
print(f"  {k} frames in {time.time() - t0:.1f}s")
