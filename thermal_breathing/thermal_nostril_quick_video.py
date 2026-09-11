#!/usr/bin/env python
"""
thermal_nostril_quick_video.py  --  quick thermal-only nostril video for ONE session.

Renders TWO side-by-side THERMAL crops, an 8x8 px window centered on EACH DLC-tracked
nostril (left panel = L nostril, right panel = R nostril, each on its OWN center),
following them frame-by-frame, drawn with a blue-white-red colormap. This is the
breathing readout (warms/cools per breath).

The crop is follow-the-point tracking only (integer-pixel translation); the raw
thermal frame is never warped. One output frame per thermal frame, played back at
OUT_FPS (slow-mo relative to the ~400 Hz acquisition).

Output: <date>_<animal>_n<k>_nostril_thermal.mp4 in the run folder.
Run in the flir env (fnv + cv2 + scipy + numpy).
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
RUN_FOLDER = r"D:\260615_thermalNbasler\5916296\cam1_20260615_201637_run001"
DATA_ROOT  = r"D:\260615_thermalNbasler"
WINDOW     = None          # (t0,t1) seconds to render; None = whole clip
OUT_FPS    = 60.0          # playback fps of the output mp4
FRAME_STEP = 1             # render every Nth thermal frame
CROP_WIN   = 5             # px crop around each nostril (own center)
CROP_PANEL = 320           # display size of each crop panel (px, square)
LP_HZ      = 60.0          # temporal low-pass cutoff (Hz) for denoise; None = off
BASELINE_HZ = 0.5          # baseline = per-pixel low-pass at this cutoff (Hz)
                           #   plotted signal = dT/T = (T - baseline)/baseline
FLIR_EXP   = 0.97648e-3    # s (only used for the on-frame time stamp)
# ======================================================================


def bwr_lut(n=256):
    stops = np.array([0, .25, .5, .75, 1.0])
    cols  = np.array([[0,0,.5],[0,.5,1.],[1,1,1.],[1.,0,0],[.5,0,0]])
    xq = np.linspace(0, 1, n)
    rgb = np.stack([np.interp(xq, stops, cols[:, c]) for c in range(3)], axis=1)
    return (rgb[:, ::-1] * 255).astype(np.uint8)


def crop_patch(fr, cx, cy, win, H, W):
    h = win // 2
    icx, icy = int(round(cx)), int(round(cy))
    out = np.full((win, win), np.nan, np.float32)
    x0, x1, y0, y1 = icx - h, icx - h + win, icy - h, icy - h + win
    sx0, sy0, sx1, sy1 = max(0, x0), max(0, y0), min(W, x1), min(H, y1)
    if sx1 > sx0 and sy1 > sy0:
        out[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = fr[sy0:sy1, sx0:sx1]
    return out


def make_video(run_folder, data_root=DATA_ROOT, window=WINDOW,
               out_fps=OUT_FPS, frame_step=FRAME_STEP):
    ats  = glob.glob(os.path.join(run_folder, "Rec-*.ats"))[0]
    nmat = glob.glob(os.path.join(run_folder, "Rec-*_nostrilC.mat"))[0]
    animal = os.path.basename(os.path.dirname(run_folder))
    runs = sorted(glob.glob(os.path.join(data_root, animal, "cam1_*")))
    k = runs.index(run_folder) + 1
    date = os.path.basename(data_root).split("_")[0]
    out = os.path.join(run_folder, f"{date}_{animal}_n{k}_nostril_thermal.mp4")
    win = CROP_WIN

    S = loadmat(nmat)
    Lc = np.asarray(S["L_center"], float); Rc = np.asarray(S["R_center"], float)

    im = fnv.file.ImagerFile(ats)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    H, W, nT = im.height, im.width, im.num_frames
    print(f"animal={animal} n{k}  thermal {W}x{H}  {nT} frames")

    # load the whole stack (T,H,W) + per-frame time
    vol = np.empty((nT, H, W), np.float32); tt = np.empty(nT); t0 = None
    for i in range(nT):
        im.get_frame(i); ti = im.frame_info.time
        if t0 is None: t0 = ti
        tt[i] = (ti - t0).total_seconds()
        vol[i] = np.array(im.final, np.float32, copy=True).reshape(H, W)
        if i % 5000 == 0: print(f"    load {i}/{nT}", flush=True)
    tt += FLIR_EXP / 2.0

    fs = (nT - 1) / (tt[-1] - tt[0])

    # temporal low-pass (per-pixel, zero-phase) at LP_HZ to denoise
    if LP_HZ is not None:
        b, a = butter(4, LP_HZ / (fs / 2.0), btype="low")
        vol = filtfilt(b, a, vol, axis=0).astype(np.float32)
        print(f"  low-passed @ {LP_HZ:.0f} Hz (fs={fs:.1f} Hz)")

    # dT/T : remove slow baseline (per-pixel low-pass), normalize like dF/F
    b0, a0 = butter(2, BASELINE_HZ / (fs / 2.0), btype="low")
    base = filtfilt(b0, a0, vol, axis=0).astype(np.float32)
    vol = (vol - base) / np.maximum(base, 1e-3)
    print(f"  dT/T computed (baseline low-pass @ {BASELINE_HZ:.2f} Hz)")

    # symmetric colormap range from sampled crops (centered on 0)
    samp = []
    for i in range(0, nT, max(1, nT // 80)):
        j = min(i, Lc.shape[0] - 1)
        samp.append(crop_patch(vol[i], Lc[j, 0], Lc[j, 1], win, H, W))
        samp.append(crop_patch(vol[i], Rc[j, 0], Rc[j, 1], win, H, W))
    sv = np.concatenate([s[~np.isnan(s)].ravel() for s in samp])
    amp = np.percentile(np.abs(sv), 98)
    lo, hi = -amp, amp

    w0, w1 = (tt[0], tt[-1]) if window is None else window
    idx = np.where((tt >= w0) & (tt <= w1))[0][::frame_step]
    print(f"  rendering {len(idx)} frames {w0:.1f}-{w1:.1f}s @ {out_fps:.0f} fps "
          f"-> nostril window [{lo:.2f},{hi:.2f}]C, crop {win}x{win}")

    lut = bwr_lut(); PW = CROP_PANEL

    def panel(fr, cx, cy):
        cr = crop_patch(fr, cx, cy, win, H, W)
        cr = np.where(np.isnan(cr), lo, cr)
        g8 = np.clip((cr - lo) / (hi - lo) * 255, 0, 255).astype(np.uint8)
        return cv2.resize(lut[g8], (PW, PW), interpolation=cv2.INTER_NEAREST)

    vw = cv2.VideoWriter(out, cv2.VideoWriter_fourcc(*"mp4v"), out_fps, (2 * PW, PW), True)
    assert vw.isOpened(), "VideoWriter failed"

    for n, ti in enumerate(idx):
        fr = vol[int(ti)]
        j = min(int(ti), Lc.shape[0] - 1)
        frame = np.hstack([panel(fr, Lc[j, 0], Lc[j, 1]),
                           panel(fr, Rc[j, 0], Rc[j, 1])])
        cv2.putText(frame, "L", (8, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 255, 255), 2, cv2.LINE_AA)
        cv2.putText(frame, "R", (PW + 8, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 255, 255), 2, cv2.LINE_AA)
        cv2.putText(frame, f"t={tt[ti]:.3f}s  T#{ti}", (8, PW - 12),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 1, cv2.LINE_AA)
        vw.write(frame)
        if n % 1000 == 0: print(f"    {n}/{len(idx)}", flush=True)

    vw.release()
    print(f"  saved {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
    return out


if __name__ == "__main__":
    make_video(RUN_FOLDER)
