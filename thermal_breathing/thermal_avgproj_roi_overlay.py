#!/usr/bin/env python
"""
thermal_avgproj_roi_overlay.py  --  full-FOV thermal AVG-PROJECTION (grayscale)
with the manually-drawn nostril ROIs overlaid, for QC of where the ROIs sit on
the actual thermal structure. One (or all) session(s).

For each session it reads the FLIR .ats as temperature, stabilizes the whole frame
by the DLC-tracked nostril MIDPOINT (head-motion follow-the-point, integer-pixel
roll) so fine spatial structure is not motion-blurred, then makes TWO grayscale
projections of the full nose FOV (upscaled, like the DLC-tracking video look):

  RAW    : time-mean of the stabilized temperature  (anatomy / thermal structure)
  LPSUB  : per-pixel RMS of the sub-LP-1Hz signal (raw - lowpass@LP_CUT) = where
           the breath-band modulation lives spatially  (NOT motion-blurred baseline)

The L/R manual ellipse/circle ROIs (drawn in crop-local coords on the tracking-
aligned crop) are mapped back to full-FOV px via each nostril's median center and
drawn on both panels (L = cyan, R = yellow).

Outputs, next to the .ats:
  <stem>_avgproj_raw.png        <stem>_avgproj_lpsub.png
Run in the flir env (fnv + cv2 + scipy + numpy; no matplotlib needed).
"""
import glob
import os
import numpy as np
import cv2
import fnv
import fnv.file
from scipy.io import loadmat
from scipy.signal import butter, filtfilt

# ============================ USER-EDITABLE ============================
DATA_ROOT  = r"D:\260615_thermalNbasler"
BATCH      = True      # False = just RUN_FOLDER; True = every session w/ a _nostrilROI.mat
RUN_FOLDER = r"D:\260615_thermalNbasler\5916296\cam1_20260615_201637_run001"
LP_CUT     = 1.0       # Hz, baseline low-pass; LPSUB map = RMS of (raw - lowpass)
UPSCALE    = 8         # display zoom of the tiny thermal frame
INTERP     = cv2.INTER_CUBIC
CLAHE_RAW  = True      # local-contrast equalize the RAW avg proj (crisper structure)
# ======================================================================


def find_session_files(run_folder):
    ats = glob.glob(os.path.join(run_folder, "Rec-*.ats"))[0]
    stem = os.path.splitext(ats)[0]
    nc = stem + "_nostrilC.mat"
    roi = stem + "_nostrilROI.mat"
    return ats, nc, roi, stem


def roi_params(side_struct, win, center_med):
    """Map a crop-local ROI (center/semiaxes/angle or radius) to full-FOV native px.
    Crop is win x win, centered on the nostril; crop center cc = (win+1)/2 (1-based)
    corresponds to the nostril's native center -> native = center_med + (roi - cc)."""
    cc = (win + 1) / 2.0
    rc = np.asarray(side_struct.center, float)         # [x,y] crop-local, 1-based
    nat = center_med + (rc - cc)                        # full-FOV native px [x,y]
    if getattr(side_struct, "shape", "ellipse") == "circle" or hasattr(side_struct, "radius"):
        r = float(getattr(side_struct, "radius", 2.0))
        axes = (r, r); ang = 0.0
    else:
        sa = np.asarray(side_struct.semiaxes, float)    # semi-axes (radii), native px
        axes = (float(sa[0]), float(sa[1]))
        ang = float(getattr(side_struct, "angle", 0.0))
    return nat, axes, ang


def draw_roi(img, nat, axes, ang, up, color, label):
    cx = nat[0] * up + up / 2.0
    cy = nat[1] * up + up / 2.0
    ax = (max(1, int(round(axes[0] * up))), max(1, int(round(axes[1] * up))))
    cv2.ellipse(img, (int(round(cx)), int(round(cy))), ax, ang, 0, 360, color, 2, cv2.LINE_AA)
    cv2.putText(img, label, (int(round(cx)) + ax[0] + 4, int(round(cy))),
                cv2.FONT_HERSHEY_SIMPLEX, 0.6, color, 2, cv2.LINE_AA)


def gray_to_bgr_up(g2d, lo, hi, up, clahe=False):
    g8 = np.clip((g2d - lo) / (hi - lo + 1e-9) * 255, 0, 255).astype(np.uint8)
    if clahe:
        g8 = cv2.createCLAHE(clipLimit=3.0, tileGridSize=(4, 4)).apply(g8)
    big = cv2.resize(g8, (g8.shape[1] * up, g8.shape[0] * up), interpolation=INTERP)
    return cv2.cvtColor(big, cv2.COLOR_GRAY2BGR)


def process_session(run_folder):
    ats, ncp, roip, stem = find_session_files(run_folder)
    if not (os.path.isfile(ncp) and os.path.isfile(roip)):
        print(f"  SKIP {os.path.basename(run_folder)}: missing nostrilC/ROI"); return None
    animal = os.path.basename(os.path.dirname(run_folder))
    runs = sorted(glob.glob(os.path.join(DATA_ROOT, animal, "cam1_*")))
    k = runs.index(run_folder) + 1
    date = os.path.basename(DATA_ROOT).split("_")[0]
    tag = f"{date} {animal} n{k}"

    S = loadmat(ncp, squeeze_me=True, struct_as_record=False)
    win = int(S["win"])
    Lc = np.asarray(S["L_center"], float); Rc = np.asarray(S["R_center"], float)
    Lmed = np.median(Lc, axis=0); Rmed = np.median(Rc, axis=0)
    Rm = loadmat(roip, squeeze_me=True, struct_as_record=False)
    Lnat, Lax, Lang = roi_params(Rm["L"], win, Lmed)
    Rnat, Rax, Rang = roi_params(Rm["R"], win, Rmed)

    im = fnv.file.ImagerFile(ats)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    H, W, nT = im.height, im.width, im.num_frames
    nT = min(nT, Lc.shape[0])
    print(f"{tag}: thermal {W}x{H}  {nT} frames -> stabilizing + projecting")

    # load full stack + per-frame time
    vol = np.empty((nT, H, W), np.float32); tt = np.empty(nT); t0 = None
    for i in range(nT):
        im.get_frame(i); ti = im.frame_info.time
        if t0 is None: t0 = ti
        tt[i] = (ti - t0).total_seconds()
        vol[i] = np.array(im.final, np.float32, copy=True).reshape(H, W)
        if i % 8000 == 0: print(f"    load {i}/{nT}", flush=True)
    fs = (nT - 1) / (tt[-1] - tt[0])

    # stabilize full frame by the tracked nostril midpoint (integer roll)
    mc = 0.5 * (Lc[:nT] + Rc[:nT]); mc_med = np.median(mc, axis=0)
    for i in range(nT):
        dx, dy = np.round(mc_med - mc[i]).astype(int)
        if dx or dy:
            vol[i] = np.roll(np.roll(vol[i], dy, axis=0), dx, axis=1)

    raw = vol.mean(axis=0)                                # RAW avg proj
    b, a = butter(2, LP_CUT / (fs / 2.0), btype="low")
    vol -= filtfilt(b, a, vol, axis=0).astype(np.float32)  # in-place -> sub-LP-1Hz
    lpsub = vol.std(axis=0)                               # breath-band RMS map
    print(f"    fs={fs:.1f}Hz  raw[{np.percentile(raw,2):.2f},{np.percentile(raw,98):.2f}]C "
          f"lpsub RMS max {lpsub.max():.3f}C")

    out = {}
    for name, g2d, clahe in [("raw", raw, CLAHE_RAW), ("lpsub", lpsub, False)]:
        lo, hi = np.percentile(g2d, [1, 99])
        bgr = gray_to_bgr_up(g2d, lo, hi, UPSCALE, clahe=clahe)
        draw_roi(bgr, Lnat, Lax, Lang, UPSCALE, (255, 255, 0), "L")  # cyan
        draw_roi(bgr, Rnat, Rax, Rang, UPSCALE, (0, 255, 255), "R")  # yellow
        ttl = f"{tag}  {'RAW avg proj (C)' if name=='raw' else f'sub-LP{LP_CUT:g}Hz RMS (breath band)'}"
        cv2.putText(bgr, ttl, (8, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (60, 220, 60), 2, cv2.LINE_AA)
        p = f"{stem}_avgproj_{name}.png"
        cv2.imwrite(p, bgr); out[name] = p
        print(f"    saved {os.path.basename(p)}")
    return out


def main():
    if BATCH:
        rois = glob.glob(os.path.join(DATA_ROOT, "*", "cam1_*", "Rec-*_nostrilROI.mat"))
        folders = sorted({os.path.dirname(r) for r in rois})
        print(f"BATCH: {len(folders)} sessions with a manual ROI")
        for f in folders:
            try:
                process_session(f)
            except Exception as e:
                print(f"  ERROR {f}: {e}")
    else:
        process_session(RUN_FOLDER)


if __name__ == "__main__":
    main()
