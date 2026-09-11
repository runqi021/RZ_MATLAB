#!/usr/bin/env python
"""
thermal_basler_sync_video.py  --  synchronized dual-cam video for ONE session.

LEFT  : BASLER whisker video with the DLC whisker points overlaid.
RIGHT : a 30x30 thermal crop centered on the MIDPOINT between the two tracked
        nostrils (follows them frame-by-frame), blue-white-red colormap -- the
        breathing readout (warms/cools with each breath).

Frames are paired by EXPOSURE-COMPENSATED relative time (cameras share a WFG
trigger but clocks aren't comparable: basler wall=2026, FLIR internal=1976 ->
align by time-from-first-frame using each cam's own per-frame timestamps, so
thermal dropped frames are handled). Effective sample = trigger + exposure/2:
    BASLER 0.25 ms (+0.125),  FLIR 0.97648 ms (+0.488).

Output: <date>_<animal>_n<k>_sync_video.mp4 in the run folder.
Run in the flir env (fnv + cv2 + scipy + numpy).
"""
import glob
import os
import re
import numpy as np
import cv2
import fnv
import fnv.file
from scipy.io import loadmat

# ============================ USER-EDITABLE ============================
RUN_FOLDER  = r"D:\260615_thermalNbasler\5916297\cam1_20260615_195531_run001"
DATA_ROOT   = r"D:\260615_thermalNbasler"
WHISK_DIR   = r"D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos"
WINDOW      = (20.0, 30.0)   # seconds to render (test window); None = whole clip
OUT_FPS     = 60.0           # playback fps of the output mp4 (400 Hz -> slow-mo)
FRAME_STEP  = 1              # render every Nth basler frame
CROP_WIN    = 30            # px crop around the nostril midpoint
CROP_PANEL  = 486           # display size of the crop (px, square)
BASLER_EXP  = 0.25e-3       # s
FLIR_EXP    = 0.97648e-3    # s
LIK_MIN     = 0.6           # whisker points below this drawn dim
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


def find_inputs(run_folder, data_root, whisk_dir):
    ats = glob.glob(os.path.join(run_folder, "Rec-*.ats"))[0]
    avi = [a for a in glob.glob(os.path.join(run_folder, "*.avi")) if "_dlc" not in os.path.basename(a)][0]
    tsf = os.path.join(run_folder, "timestamps.csv")
    nmat = glob.glob(os.path.join(run_folder, "Rec-*_nostrilC.mat"))
    assert nmat, f"no _nostrilC.mat in {run_folder}"
    nmat = nmat[0]
    animal = os.path.basename(os.path.dirname(run_folder))
    runs = sorted(glob.glob(os.path.join(data_root, animal, "cam1_*")))
    k = runs.index(run_folder) + 1
    cand = glob.glob(os.path.join(whisk_dir, f"{animal}_whisk_n{k}*DLC*.csv"))
    cand.sort(key=lambda c: int(re.search(r"best-(\d+)", c).group(1)) if re.search(r"best-(\d+)", c) else 0)
    whisk = cand[-1] if cand else None
    return ats, avi, tsf, nmat, whisk, animal, k


def make_sync(run_folder, data_root=DATA_ROOT, whisk_dir=WHISK_DIR,
              window=WINDOW, out_fps=OUT_FPS, frame_step=FRAME_STEP):
    ats, avi, tsf, nmat, whisk, animal, k = find_inputs(run_folder, data_root, whisk_dir)
    date = os.path.basename(data_root).split("_")[0]
    out = os.path.join(run_folder, f"{date}_{animal}_n{k}_sync_video.mp4")
    print(f"animal={animal} n{k}  whisk={'yes' if whisk else 'NO'}")
    win = CROP_WIN

    bts = np.loadtxt(tsf, delimiter=",", skiprows=1)
    bt = (bts[:, 1] - bts[0, 1]) / 1e9 + BASLER_EXP / 2.0

    S = loadmat(nmat); Lc = np.asarray(S["L_center"], float); Rc = np.asarray(S["R_center"], float)

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

    Wm = np.loadtxt(whisk, delimiter=",", skiprows=3) if whisk else None

    w0, w1 = (0.0, bt[-1]) if window is None else window
    bidx = np.where((bt >= w0) & (bt <= w1))[0][::frame_step]
    print(f"  rendering {len(bidx)} frames {w0:.1f}-{w1:.1f}s @ {out_fps:.0f} fps -> nostril window [{lo:.2f},{hi:.2f}]C")

    lut = bwr_lut(); PW = CROP_PANEL
    cap = cv2.VideoCapture(avi)
    basW = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)); basH = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    print(f"  basler {basW}x{basH} (whisker coords overlaid at native res)")
    vw = cv2.VideoWriter(out, cv2.VideoWriter_fourcc(*"mp4v"), out_fps, (basW + PW, basH), True)
    assert vw.isOpened(), "VideoWriter failed"

    def patch_tile(fr, cx, cy):
        cr = crop_patch(fr, cx, cy, win, H, W)
        cr = np.where(np.isnan(cr), lo, cr)
        g8 = np.clip((cr - lo) / (hi - lo) * 255, 0, 255).astype(np.uint8)
        return cv2.resize(lut[g8], (PW, basH), interpolation=cv2.INTER_NEAREST)

    cap.set(cv2.CAP_PROP_POS_FRAMES, int(bidx[0])); cur = int(bidx[0])
    for n, bi in enumerate(bidx):
        while cur < bi:
            cap.grab(); cur += 1
        ok, fb = cap.read(); cur += 1
        if not ok: break
        if fb.ndim == 2 or fb.shape[2] == 1:
            fb = cv2.cvtColor(fb, cv2.COLOR_GRAY2BGR)
        if fb.shape[1] != basW or fb.shape[0] != basH:
            fb = cv2.resize(fb, (basW, basH))    # keep native so DLC coords align
        if Wm is not None and bi < Wm.shape[0]:
            r = Wm[bi]
            for (bx, by, bl, tx, ty, tl, col) in [
                (r[1], r[2], r[3], r[4], r[5], r[6], (0, 255, 0)),
                (r[7], r[8], r[9], r[10], r[11], r[12], (255, 255, 0))]:
                c = col if min(bl, tl) >= LIK_MIN else tuple(int(x*0.4) for x in col)
                p0 = (int(bx), int(by)); p1 = (int(tx), int(ty))
                cv2.line(fb, p0, p1, c, 1, cv2.LINE_AA)
                cv2.circle(fb, p0, 3, c, -1, cv2.LINE_AA)
                cv2.circle(fb, p1, 4, c, 1, cv2.LINE_AA)
        tj = int(np.argmin(np.abs(tt - bt[bi])))
        im.get_frame(tj)
        fr = np.array(im.final, np.float32, copy=True).reshape(H, W)
        j = min(tj, Lc.shape[0] - 1); mc = (Lc[j] + Rc[j]) / 2.0
        frame = np.hstack([fb, patch_tile(fr, mc[0], mc[1])])
        cv2.putText(frame, f"t={bt[bi]:.3f}s  B#{bi} T#{tj}", (8, 22),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 1, cv2.LINE_AA)
        vw.write(frame)
        if n % 1000 == 0: print(f"    {n}/{len(bidx)}", flush=True)

    cap.release(); vw.release()
    print(f"  saved {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
    return out


def main():
    make_sync(RUN_FOLDER)


if __name__ == "__main__":
    main()
