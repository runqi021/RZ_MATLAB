#!/usr/bin/env python
"""
thermal_nostril_breath_single.py  --  SINGLE-video thermal nostril breathing extractor.

Reads a FLIR .ats as TRUE temperature (deg C), follows the DLC-tracked nostrils
(L1 = LEFT, L2 = RIGHT) frame-by-frame, and pulls a breathing signal in actual
degrees C from a small ROI that moves with each nostril. Also builds, per nostril,
a tracking-aligned AVERAGE projection and a BREATH-BAND POWER MAP so you can later
draw a precise ellipse ROI in MATLAB (this script uses a simple default disk ROI so
you can judge the signal right away).

Alignment chain (verified from the sidecar _dlc.json):
  - thermal_dlc_video.py wrote ONE video frame per .ats frame (decim) and upscaled
    by `upscale`. So:  ats_frame = video_frame * decim ;  ats_px = dlc_px / upscale.

Run in the flir env (FileSDK; no matplotlib/pandas there -> we save a .mat and a
MATLAB script does the viewing/ROI drawing):
    C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe thermal_nostril_breath_single.py

Output: <ats_stem>_nostrilC.mat  (next to the .ats), consumed by the MATLAB viewer.
"""
import glob
import json
import os
import re
import sys
import numpy as np
import fnv
import fnv.file
from scipy.io import savemat

# ============================ USER-EDITABLE ============================
# Give ONLY the DLC csv (here, or as argv[1]). The .ats + _dlc.json are auto-
# resolved: DATA_ROOT/<animal>/cam1_*/ , picking the run by the _n# in the csv
# name (animal = csv prefix; n1 = 1st cam1 run chronologically, n2 = 2nd, ...).
DATA_ROOT = r"D:\260615_thermalNbasler"
DLC_CSV   = r"D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos\5916296_nose_n1DLC_Resnet50_260615_thermal_noseJun17shuffle2_snapshot_best-100.csv"

WIN_NATIVE   = 24      # crop window (native .ats px) around each nostril
STACK_FPS    = 400.0   # saved aligned-stack rate; 400 = native (no decimation, full fps)
LIK_THRESH   = 0.6     # below this likelihood -> drop nostril pos & linterp over time
ANCHOR       = 'own'   # crop centering: 'own' (each dot; DEFAULT/best for a nostril's own
                       #   breath) | 'contra' | 'midpoint'. NOTE: contra/midpoint were tested
                       #   and they DEGRADE the breath SNR (inject inter-dot jitter) -- keep 'own'.
OUT          = None    # default: <ats_stem>_nostrilC.mat
# ======================================================================

BPARTS = [("L", "L1", "LEFT"), ("R", "L2", "RIGHT")]   # key, dlc name, side label


def resolve_ats(csv, data_root):
    """From just the DLC csv, find the source .ats + _dlc.json.
    animal = csv prefix; _n# = run index into chronologically-sorted cam1_* folders."""
    base = os.path.basename(csv)
    animal = base.split("_")[0]
    nums = re.findall(r"_n(\d+)", base)
    if not nums:
        sys.exit(f"can't parse run index (_n#) from csv name: {base}")
    k = int(nums[-1])
    adir = os.path.join(data_root, animal)
    runs = sorted(glob.glob(os.path.join(adir, "cam1_*")))
    if not runs:
        sys.exit(f"no cam1_* run folders under {adir}")
    if k > len(runs):
        sys.exit(f"csv says _n{k} but only {len(runs)} cam1 runs in {adir}")
    folder = runs[k - 1]
    ats = sorted(glob.glob(os.path.join(folder, "Rec-*.ats")))
    if not ats:
        sys.exit(f"no Rec-*.ats in {folder}")
    ats = ats[0]
    jp = os.path.splitext(ats)[0] + "_dlc.json"
    if not os.path.isfile(jp):
        sys.exit(f"no sidecar json: {jp}")
    print(f"resolved: animal={animal} n{k} -> {os.path.basename(folder)}")
    return ats, jp


def load_json(jp):
    with open(jp) as f:
        j = json.load(f)
    return j


def load_dlc(csv, n_bodyparts_expected=2):
    """Read a DLC analyzed CSV with numpy: [frame, (x,y,lik) per bodypart]."""
    M = np.loadtxt(csv, delimiter=",", skiprows=3)
    ncol = 1 + 3 * n_bodyparts_expected
    if M.shape[1] != ncol:
        sys.exit(f"DLC csv has {M.shape[1]} cols, expected {ncol} for "
                 f"{n_bodyparts_expected} bodyparts (frame + x,y,lik each)")
    out = {}
    for i, (key, name, _side) in enumerate(BPARTS):
        x = M[:, 1 + 3 * i]
        y = M[:, 2 + 3 * i]
        lik = M[:, 3 + 3 * i]
        out[key] = dict(name=name, x=x, y=y, lik=lik)
    return out, M.shape[0]


def gate_interp(xy, lik, thresh):
    """NaN low-likelihood positions, linearly interpolate over frame index."""
    x = xy.copy()
    bad = lik < thresh
    n = len(x)
    idx = np.arange(n)
    if bad.all():
        return x
    if bad.any():
        x[bad] = np.interp(idx[bad], idx[~bad], x[~bad])
    return x


def extract_one(csv):
    if not os.path.isfile(csv):
        sys.exit(f"no DLC csv: {csv}")
    ats, jp = resolve_ats(csv, DATA_ROOT)
    j = load_json(jp)
    upscale = int(j["upscale"])
    decim = int(j["decim"])
    print(f"json {os.path.basename(jp)}: upscale={upscale} decim={decim} "
          f"ats_frames={j['n_ats_frames']} {j['H']}x{j['W']}")

    dlc, n_rows = load_dlc(csv)
    print(f"DLC: {n_rows} rows, bodyparts L1(LEFT)/L2(RIGHT)")
    if n_rows != int(j["n_video_frames"]):
        print(f"  WARN: csv rows {n_rows} != json video frames {j['n_video_frames']} "
              f"-- run-index mapping may be wrong!")

    im = fnv.file.ImagerFile(ats)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    H, W, N = im.height, im.width, im.num_frames
    if n_rows != N // decim and n_rows != N:
        print(f"  WARN: DLC rows {n_rows} vs ats frames {N} (decim {decim})")
    nframes = min(n_rows, N // decim if decim > 1 else N)

    # nostril centers in NATIVE .ats px (gate + interp low-confidence)
    centers = {}
    for key, _name, _side in [(k, n, s) for k, n, s in BPARTS]:
        cx = gate_interp(dlc[key]["x"] / upscale, dlc[key]["lik"], LIK_THRESH)
        cy = gate_interp(dlc[key]["y"] / upscale, dlc[key]["lik"], LIK_THRESH)
        centers[key] = np.stack([cx, cy], axis=1)  # [N,2] native px

    # --- re-anchor crop centers to kill per-dot jitter (rigid nostril pair) ---
    if ANCHOR != 'own':
        mid = 0.5 * (centers['L'] + centers['R'])              # per-frame midpoint of both dots
        recentered = {}
        for key in ('L', 'R'):
            other = 'R' if key == 'L' else 'L'
            ref = mid if ANCHOR == 'midpoint' else centers[other]
            if ANCHOR not in ('midpoint', 'contra'):
                sys.exit(f"unknown ANCHOR {ANCHOR!r} (use own|contra|midpoint)")
            off = np.nanmedian(centers[key] - ref, axis=0)     # CONSTANT offset -> follows head motion, drops jitter
            recentered[key] = ref + off
        centers = recentered
        print(f"  re-anchored crops: ANCHOR={ANCHOR}")

    win = int(WIN_NATIVE)
    hw = win // 2
    crops = {k: np.full((nframes, win, win), np.nan, np.float32) for k, _, _ in BPARTS}
    times = np.empty(nframes, np.float64)
    t0 = None

    print(f"reading {N} frames @ {H}x{W} (extracting {win}x{win} crops x2)...")
    for k in range(nframes):
        ai = k * decim
        im.get_frame(ai)
        fr = np.array(im.final, dtype=np.float32, copy=True).reshape(H, W)  # COPY (buffer reused)
        ti = im.frame_info.time
        if t0 is None:
            t0 = ti
        times[k] = (ti - t0).total_seconds()
        for key, _name, _side in BPARTS:
            cx, cy = centers[key][k]
            icx, icy = int(round(cx)), int(round(cy))
            x0, x1 = icx - hw, icx - hw + win
            y0, y1 = icy - hw, icy - hw + win
            sx0, sy0 = max(0, x0), max(0, y0)
            sx1, sy1 = min(W, x1), min(H, y1)
            if sx1 > sx0 and sy1 > sy0:
                crops[key][k, sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = fr[sy0:sy1, sx0:sx1]
        if k % 4000 == 0:
            print(f"  {k}/{nframes}", flush=True)

    # --- timing ---
    dt = np.diff(times)
    fps_nom = float(round(1.0 / np.median(dt)))
    fps_mean = (nframes - 1) / (times[-1] - times[0])
    n_exp = int(round(times[-1] * fps_nom)) + 1
    n_drop = n_exp - nframes
    print(f"  fps nominal {fps_nom:.0f}, mean {fps_mean:.2f}, dropped {n_drop} "
          f"({100*n_drop/max(n_exp,1):.1f}%)")

    out = dict(
        src_ats=os.path.abspath(ats), src_dlc=os.path.abspath(csv),
        upscale=upscale, decim=decim, win=win, anchor=ANCHOR,
        times=times, fps_nominal=fps_nom, fps_mean=float(fps_mean),
        n_frames=int(nframes), n_dropped=int(n_drop), H=int(H), W=int(W),
    )

    stack_q = max(1, int(round(fps_mean / STACK_FPS)))
    for key, name, side in BPARTS:
        C = crops[key]                                   # [N,win,win] deg C
        avg = np.nanmean(C, axis=0)                      # average projection
        out[f"{key}_side"] = side
        out[f"{key}_name"] = name
        out[f"{key}_center"] = centers[key][:nframes].astype(np.float32)
        out[f"{key}_avg"] = avg.astype(np.float32)
        out[f"{key}_lik"] = dlc[key]["lik"][:nframes].astype(np.float32)
        out[f"{key}_stack"] = C[::stack_q].astype(np.float32)        # for ROI drawing
        out[f"{key}_stack_fps"] = float(fps_mean / stack_q)
        print(f"  {side} ({name}): {C.shape[0]} frames, avg {np.nanmean(avg):.2f} C")

    outpath = OUT or (os.path.splitext(ats)[0] + "_nostrilC.mat")
    savemat(outpath, out, do_compression=True)
    print(f"saved {outpath}  ({os.path.getsize(outpath)/1e6:.1f} MB)")
    return outpath


def main():
    csv = sys.argv[1] if len(sys.argv) > 1 else DLC_CSV
    extract_one(csv)


if __name__ == "__main__":
    main()
