"""Batch overlay all 30 videos with vL1 (blue) + vR1 (red) cross markers,
150° CCW rotation, and the same crop used in the test script.

Outputs saved into subfolder `representative_whisker_label_video/`.
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
run_list      = list(range(1, 31))    # 1..30

L_regex       = re.compile(r"DLC_Resnet50_whisker_dorsalApr")
R_regex       = re.compile(r"DLC_Resnet50_whisker_dorsal_RApr")
L_point       = "vL1"
R_point       = "vR1"
pmin          = 0.5

half_len      = 8
thickness     = 2
color_L_bgr   = (217, 115,   0)
color_R_bgr   = ( 25,  25, 217)

rotate_deg    = 150
crop_xywh     = [114, 402, 646, 238]  # must match test script
# =========================


def load_point(csv_path, point_name):
    with open(csv_path, "r") as f:
        rows = list(csv.reader(f))
    bps = []
    for c in rows[1][1:]:
        c = c.strip()
        if not bps or bps[-1] != c:
            bps.append(c)
    pid = bps.index(point_name)
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
    xi, yi = int(round(x)), int(round(y))
    h, w = frame.shape[:2]
    if xi < 0 or xi >= w or yi < 0 or yi >= h:
        return
    cv2.line(frame, (xi - hl, yi), (xi + hl, yi), color_bgr, thickness=th)
    cv2.line(frame, (xi, yi - hl), (xi, yi + hl), color_bgr, thickness=th)


def process_one(run_id, out_dir):
    raw_avis = [
        a for a in glob.glob(os.path.join(folder_path, f"*run{run_id:03d}.avi"))
        if "DLC_Resnet50" not in os.path.basename(a)
    ]
    if not raw_avis:
        print(f"[run{run_id:03d}] SKIP: no raw AVI")
        return
    src_avi = raw_avis[0]

    all_csvs = glob.glob(os.path.join(folder_path, f"*run{run_id:03d}*.csv"))
    L_csv = next(
        (c for c in all_csvs
         if L_regex.search(os.path.basename(c)) and "dorsal_R" not in os.path.basename(c)),
        None
    )
    R_csv = next(
        (c for c in all_csvs if R_regex.search(os.path.basename(c))),
        None
    )
    if L_csv is None or R_csv is None:
        print(f"[run{run_id:03d}] SKIP: L={L_csv is not None}, R={R_csv is not None}")
        return

    xL, yL, pL = load_point(L_csv, L_point)
    xR, yR, pR = load_point(R_csv, R_point)
    n_csv = min(len(xL), len(xR))

    stem = os.path.splitext(os.path.basename(src_avi))[0]
    out_mp4 = os.path.join(out_dir, f"{stem}_vL1_vR1_cross.mp4")

    cap = cv2.VideoCapture(src_avi)
    if not cap.isOpened():
        print(f"[run{run_id:03d}] SKIP: cannot open {src_avi}")
        return
    fps = cap.get(cv2.CAP_PROP_FPS)
    w   = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    h   = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    # Rotation matrix + expanded canvas
    M = cv2.getRotationMatrix2D((w / 2.0, h / 2.0), rotate_deg, 1.0)
    cos_ = abs(M[0, 0]); sin_ = abs(M[0, 1])
    new_w = int(h * sin_ + w * cos_)
    new_h = int(h * cos_ + w * sin_)
    M[0, 2] += (new_w / 2.0) - (w / 2.0)
    M[1, 2] += (new_h / 2.0) - (h / 2.0)

    if crop_xywh is not None:
        cx, cy, cw, ch = [int(v) for v in crop_xywh]
        cx = max(0, min(cx, new_w - 1))
        cy = max(0, min(cy, new_h - 1))
        cw = max(1, min(cw, new_w - cx))
        ch = max(1, min(ch, new_h - cy))
        out_w, out_h = cw, ch
    else:
        cx = cy = 0
        out_w, out_h = new_w, new_h

    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    vw = cv2.VideoWriter(out_mp4, fourcc, fps, (out_w, out_h))
    if not vw.isOpened():
        print(f"[run{run_id:03d}] SKIP: writer failed for {out_mp4}")
        cap.release()
        return

    t0 = time.time()
    k = 0
    while True:
        ok, frame = cap.read()
        if not ok or k >= n_csv:
            break
        fr = cv2.warpAffine(
            frame, M, (new_w, new_h),
            flags=cv2.INTER_LINEAR,
            borderMode=cv2.BORDER_CONSTANT, borderValue=(0, 0, 0),
        )
        if crop_xywh is not None:
            fr = fr[cy:cy + ch, cx:cx + cw]

        if pL[k] >= pmin:
            xLr = M[0, 0] * xL[k] + M[0, 1] * yL[k] + M[0, 2]
            yLr = M[1, 0] * xL[k] + M[1, 1] * yL[k] + M[1, 2]
            draw_cross(fr, xLr - cx, yLr - cy, half_len, thickness, color_L_bgr)
        if pR[k] >= pmin:
            xRr = M[0, 0] * xR[k] + M[0, 1] * yR[k] + M[0, 2]
            yRr = M[1, 0] * xR[k] + M[1, 1] * yR[k] + M[1, 2]
            draw_cross(fr, xRr - cx, yRr - cy, half_len, thickness, color_R_bgr)

        vw.write(fr)
        k += 1

    cap.release()
    vw.release()
    print(f"[run{run_id:03d}] saved {k} frames in {time.time() - t0:.1f}s  ->  {os.path.basename(out_mp4)}")


# ---- main ----
out_dir = os.path.join(folder_path, out_subfolder)
os.makedirs(out_dir, exist_ok=True)

t_total = time.time()
for rid in run_list:
    process_one(rid, out_dir)

print(f"\nDone. Total: {time.time() - t_total:.1f}s, {len(run_list)} runs.")
