#!/usr/bin/env python
"""
thermal_basler_epoch_video.py  --  20x-slow sync reels of WHISK and QUIET epochs.

For every session (cam1_* run folder with a whisk DLC csv + _nostrilC.mat) this
auto-detects whisking and quiet epochs from the whisk envelope (same recipe as
whisk_epoch_analyze_auto_RZ.m / _quiet_RZ.m), then writes ONE concatenated reel
per epoch-set per session: only the in-epoch frames, all 400 fps frames kept,
played at OUT_FPS=20 -> 20x slow. A black separator card with the epoch label is
inserted between epochs.

Each frame is the same panel as thermal_basler_sync_video.py:
  LEFT  basler whisker video + DLC whisker points
  RIGHT 30x30 bwr thermal crop at the nostril midpoint (breathing readout)
plus an epoch banner (set, epoch i/N, real time, slow factor).

Output per session:
  <date>_<animal>_n<k>_whiskEpochs_20x.mp4
  <date>_<animal>_n<k>_quietEpochs_20x.mp4
Run in the flir env (fnv + cv2 + scipy + numpy).
Concatenation is PER SESSION (each animal's basler avi is a different resolution,
so cross-session concat would need rescaling).
"""
import glob
import os
import numpy as np
import cv2
import fnv
import fnv.file
from scipy.io import loadmat
from scipy.signal import butter, filtfilt, hilbert

from thermal_basler_sync_video import bwr_lut, crop_patch, find_inputs

# ============================ USER-EDITABLE ============================
DATA_ROOT  = r"D:\260615_thermalNbasler"
WHISK_DIR  = r"D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos"
SETS       = ("whisk", "quiet")  # which epoch sets to render
SLOW       = 20.0                 # slow factor
FPS_W      = 400.0               # whisk/basler acquisition rate (Hz)
OUT_FPS    = FPS_W / SLOW         # 20 fps -> 20x slow (every frame kept)
FRAME_STEP = 1                   # 1 = keep all frames (full 400 fps)
# epoch detection (matches the MATLAB analysis scripts)
BP         = (5.0, 30.0)         # whisk bandpass (Hz)
THR_FRAC   = 0.10               # threshold = THR_FRAC * 95th-pct(envelope)
MIN_DUR    = 1.0               # s, min epoch duration
MERGE_GAP  = 0.20             # s, merge epochs closer than this
# panel / overlay
CROP_WIN   = 30
CROP_PANEL = 486
BASLER_EXP = 0.25e-3
FLIR_EXP   = 0.97648e-3
LIK_MIN    = 0.6
SEP_FRAMES = 10              # black separator frames between epochs
# ======================================================================


def _fill_lin(x):
    x = np.asarray(x, float)
    bad = ~np.isfinite(x)
    if bad.any():
        idx = np.arange(x.size)
        x[bad] = np.interp(idx[bad], idx[~bad], x[~bad])
    return x


def whisk_angles(Wm):
    """Left/Right whisker angle (deg, unwrapped), mirror-left so protraction=+ both."""
    vL1 = Wm[:, [1, 2]]; vL2 = Wm[:, [4, 5]]
    vR1 = Wm[:, [7, 8]]; vR2 = Wm[:, [10, 11]]
    La = np.degrees(np.unwrap(np.arctan2(-(vL2[:, 1] - vL1[:, 1]), -(vL2[:, 0] - vL1[:, 0]))))
    Ra = np.degrees(np.unwrap(np.arctan2(-(vR2[:, 1] - vR1[:, 1]),  (vR2[:, 0] - vR1[:, 0]))))
    return _fill_lin(La), _fill_lin(Ra)


def detect_epochs(Wm):
    """Return {'whisk':[(s0,s1)..], 'quiet':[..]} as inclusive SAMPLE-index ranges."""
    La, Ra = whisk_angles(Wm)
    b, a = butter(3, np.array(BP) / (FPS_W / 2.0), btype="bandpass")
    xL = filtfilt(b, a, La - np.nanmean(La))
    xR = filtfilt(b, a, Ra - np.nanmean(Ra))
    env = (np.abs(hilbert(xL)) + np.abs(hilbert(xR))) / 2.0
    thr = THR_FRAC * np.percentile(env, 95)
    out = {}
    for name, mask in (("whisk", env > thr), ("quiet", env < thr)):
        out[name] = _runs(mask, MIN_DUR, MERGE_GAP)
    return out


def _runs(mask, min_dur, merge_gap):
    m = mask.astype(int)
    d = np.diff(np.r_[0, m, 0])
    s = np.where(d == 1)[0]; e = np.where(d == -1)[0] - 1
    ep = [[int(a), int(b)] for a, b in zip(s, e)]
    if not ep:
        return []
    merged = [ep[0]]
    for a, b in ep[1:]:
        if (a - merged[-1][1]) <= merge_gap * FPS_W:
            merged[-1][1] = b
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged if (b - a + 1) >= min_dur * FPS_W]


def render_set(run_folder, which, ats_im=None, cache=None):
    ats, avi, tsf, nmat, whisk, animal, k = find_inputs(run_folder, DATA_ROOT, WHISK_DIR)
    if whisk is None:
        print(f"  {animal} n{k}: no whisk csv -> skip"); return None
    date = os.path.basename(DATA_ROOT).split("_")[0]
    out = os.path.join(run_folder, f"{date}_{animal}_n{k}_{which}Epochs_20x.mp4")
    win = CROP_WIN

    Wm = np.loadtxt(whisk, delimiter=",", skiprows=3)
    eps = detect_epochs(Wm)[which]
    if not eps:
        print(f"  {animal} n{k} [{which}]: 0 epochs -> skip"); return None

    # basler per-frame time (exposure compensated) + nostril centers
    bts = np.loadtxt(tsf, delimiter=",", skiprows=1)
    bt = (bts[:, 1] - bts[0, 1]) / 1e9 + BASLER_EXP / 2.0
    S = loadmat(nmat); Lc = np.asarray(S["L_center"], float); Rc = np.asarray(S["R_center"], float)

    # thermal timestamps + color window (sampled over the whole clip)
    im = fnv.file.ImagerFile(ats)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    H, W, nT = im.height, im.width, im.num_frames
    tt = np.empty(nT); t0 = None; samp = []
    for i in range(nT):
        im.get_frame(i); ti = im.frame_info.time
        if t0 is None: t0 = ti
        tt[i] = (ti - t0).total_seconds()
        if i % max(1, nT // 80) == 0:
            fr = np.array(im.final, np.float32, copy=True).reshape(H, W)
            j = min(i, Lc.shape[0] - 1); mc = (Lc[j] + Rc[j]) / 2.0
            samp.append(crop_patch(fr, mc[0], mc[1], win, H, W))
    tt += FLIR_EXP / 2.0
    sv = np.concatenate([s[~np.isnan(s)].ravel() for s in samp])
    lo, hi = np.percentile(sv, [2, 98])

    lut = bwr_lut(); PW = CROP_PANEL
    cap = cv2.VideoCapture(avi)
    basW = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)); basH = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    vw = cv2.VideoWriter(out, cv2.VideoWriter_fourcc(*"mp4v"), OUT_FPS, (basW + PW, basH), True)
    assert vw.isOpened(), "VideoWriter failed"
    nf = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    tot = sum(min(b, nf - 1) - a + 1 for a, b in eps)
    print(f"  {animal} n{k} [{which}]: {len(eps)} epochs, {tot} frames "
          f"({tot/FPS_W:.0f}s real -> {tot/OUT_FPS/60:.1f} min @ {OUT_FPS:.0f}fps) basler {basW}x{basH}")

    def patch_tile(fr, cx, cy):
        cr = crop_patch(fr, cx, cy, win, H, W)
        cr = np.where(np.isnan(cr), lo, cr)
        g8 = np.clip((cr - lo) / (hi - lo) * 255, 0, 255).astype(np.uint8)
        return cv2.resize(lut[g8], (PW, basH), interpolation=cv2.INTER_NEAREST)

    col_set = (0, 255, 0) if which == "whisk" else (180, 180, 180)
    cur = -1
    for ei, (a, b) in enumerate(eps):
        b = min(b, nf - 1)
        # separator card
        card = np.zeros((basH, basW + PW, 3), np.uint8)
        cv2.putText(card, f"{which.upper()}  epoch {ei+1}/{len(eps)}", (30, basH // 2 - 10),
                    cv2.FONT_HERSHEY_SIMPLEX, 1.1, col_set, 2, cv2.LINE_AA)
        cv2.putText(card, f"{a/FPS_W:.2f}-{b/FPS_W:.2f}s  ({(b-a+1)/FPS_W:.2f}s, {SLOW:.0f}x slow)",
                    (30, basH // 2 + 30), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 255, 255), 1, cv2.LINE_AA)
        for _ in range(SEP_FRAMES):
            vw.write(card)
        # epoch frames
        idxs = range(a, b + 1, FRAME_STEP)
        if cur < 0 or a < cur:
            cap.set(cv2.CAP_PROP_POS_FRAMES, int(a)); cur = int(a)
        for bi in idxs:
            while cur < bi:
                cap.grab(); cur += 1
            ok, fb = cap.read(); cur += 1
            if not ok:
                break
            if fb.ndim == 2 or fb.shape[2] == 1:
                fb = cv2.cvtColor(fb, cv2.COLOR_GRAY2BGR)
            if fb.shape[1] != basW or fb.shape[0] != basH:
                fb = cv2.resize(fb, (basW, basH))
            if bi < Wm.shape[0]:
                r = Wm[bi]
                for (bx, by, bl, tx, ty, tl, c0) in [
                        (r[1], r[2], r[3], r[4], r[5], r[6], (0, 255, 0)),
                        (r[7], r[8], r[9], r[10], r[11], r[12], (255, 255, 0))]:
                    c = c0 if min(bl, tl) >= LIK_MIN else tuple(int(x * 0.4) for x in c0)
                    p0 = (int(bx), int(by)); p1 = (int(tx), int(ty))
                    cv2.line(fb, p0, p1, c, 1, cv2.LINE_AA)
                    cv2.circle(fb, p0, 3, c, -1, cv2.LINE_AA)
                    cv2.circle(fb, p1, 4, c, 1, cv2.LINE_AA)
            tj = int(np.argmin(np.abs(tt - bt[min(bi, len(bt) - 1)])))
            im.get_frame(tj)
            fr = np.array(im.final, np.float32, copy=True).reshape(H, W)
            j = min(tj, Lc.shape[0] - 1); mc = (Lc[j] + Rc[j]) / 2.0
            frame = np.hstack([fb, patch_tile(fr, mc[0], mc[1])])
            cv2.putText(frame, f"{which} ep {ei+1}/{len(eps)}  t={bi/FPS_W:.3f}s  {SLOW:.0f}x",
                        (8, 22), cv2.FONT_HERSHEY_SIMPLEX, 0.6, col_set, 1, cv2.LINE_AA)
            vw.write(frame)

    cap.release(); vw.release()
    print(f"    saved {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
    return out


def list_sessions():
    sess = []
    for animal in sorted(glob.glob(os.path.join(DATA_ROOT, "*"))):
        if not os.path.isdir(animal) or not os.path.basename(animal).isdigit():
            continue
        for rf in sorted(glob.glob(os.path.join(animal, "cam1_*"))):
            if glob.glob(os.path.join(rf, "Rec-*_nostrilC.mat")):
                sess.append(rf)
    return sess


def main():
    sess = list_sessions()
    print(f"{len(sess)} sessions; sets={SETS}; {SLOW:.0f}x slow @ {OUT_FPS:.0f} fps")
    total = 0.0
    for rf in sess:
        for which in SETS:
            o = render_set(rf, which)
            if o:
                total += os.path.getsize(o) / 1e9
    print(f"DONE. total written {total:.1f} GB")


if __name__ == "__main__":
    main()
