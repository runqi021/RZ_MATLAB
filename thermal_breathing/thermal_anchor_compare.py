"""thermal_anchor_compare.py  (flir env)
Quick side-by-side of crop anchoring (own / contra / midpoint) for ONE video.
Reads the first NFRAMES of the .ats, builds the 6 crop stacks, computes the
near-center breath SNR = power(2-10 Hz)/power(10-40 Hz) for each, and saves
<ats_stem>_anchorcompare.mat (traces + SNR) for the MATLAB viewer.
    C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe thermal_anchor_compare.py [DLC_CSV]
"""
import os, sys
import numpy as np
from scipy.signal import butter, filtfilt
from scipy.io import savemat
import fnv, fnv.file
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import thermal_nostril_breath_single as T

NFRAMES = 6000          # quick chunk (~15 s @ 400 Hz)
WIN     = T.WIN_NATIVE
HW      = WIN // 2
ANCHORS = ['own', 'contra', 'midpoint']

def anchor_centers(centers, mode):
    if mode == 'own':
        return centers
    mid = 0.5 * (centers['L'] + centers['R'])
    out = {}
    for key in ('L', 'R'):
        other = 'R' if key == 'L' else 'L'
        ref = mid if mode == 'midpoint' else centers[other]
        out[key] = ref + np.nanmedian(centers[key] - ref, axis=0)
    return out

def snr_map(F, fps):                      # F: [T, win*win]
    nyq = fps/2
    bb,ab = butter(2,[2,10]/np.array(nyq),'bandpass')
    hi = min(40,0.95*nyq); bn,an = butter(2,[10,hi]/np.array(nyq),'bandpass')
    Fc = np.where(np.isfinite(F), F, 0.0)
    Pin = filtfilt(bb,ab,Fc,axis=0).var(0)
    Pno = filtfilt(bn,an,Fc,axis=0).var(0)
    snr = Pin/(Pno+1e-12); snr[Pin < 0.10*Pin.max()] = 0
    return snr, Pin

def best_trace(C, fps):                   # C: [T, win, win]
    T_, w, _ = C.shape
    F = C.reshape(T_, -1)
    snr, _ = snr_map(F, fps)
    S = snr.reshape(w, w)
    yy, xx = np.mgrid[0:w, 0:w]; cy = (w-1)/2; cx = (w-1)/2
    S = S.copy(); S[np.hypot(xx-cx, yy-cy) > 4] = 0      # near-center only
    pk = np.unravel_index(np.argmax(S), S.shape)         # (row,col)
    r0,r1 = max(0,pk[0]-1), min(w,pk[0]+2); c0,c1 = max(0,pk[1]-1), min(w,pk[1]+2)
    tr = np.nanmean(C[:, r0:r1, c0:c1].reshape(T_, -1), axis=1)
    tr = np.where(np.isfinite(tr), tr, np.nanmean(tr))
    bb,ab = butter(2, 1.0/(fps/2), 'low'); detr = -(tr - filtfilt(bb,ab,tr))  # invert (inhale up)
    bb2,ab2 = butter(2,[2,10]/np.array(fps/2),'bandpass')
    bn,an  = butter(2,[10, min(40,0.95*fps/2)]/np.array(fps/2),'bandpass')
    snrv = filtfilt(bb2,ab2,detr).var()/(filtfilt(bn,an,detr).var()+1e-12)
    return detr.astype(np.float32), float(snrv), float(S.max())

def main():
    csv = sys.argv[1] if len(sys.argv) > 1 else T.DLC_CSV
    ats, jp = T.resolve_ats(csv, T.DATA_ROOT)
    j = T.load_json(jp); upscale = int(j["upscale"]); decim = int(j["decim"])
    dlc, n_rows = T.load_dlc(csv)
    im = fnv.file.ImagerFile(ats); im.unit = fnv.Unit.TEMPERATURE_FACTORY; im.temp_type = fnv.TempType.CELSIUS
    H, W, N = im.height, im.width, im.num_frames
    nf = min(NFRAMES, n_rows, N // decim if decim > 1 else N)
    centers = {}
    for key, _n, _s in T.BPARTS:
        cx = T.gate_interp(dlc[key]["x"]/upscale, dlc[key]["lik"], T.LIK_THRESH)[:nf]
        cy = T.gate_interp(dlc[key]["y"]/upscale, dlc[key]["lik"], T.LIK_THRESH)[:nf]
        centers[key] = np.stack([cx, cy], 1)
    cen = {a: anchor_centers(centers, a) for a in ANCHORS}
    crops = {(a,k): np.full((nf,WIN,WIN), np.nan, np.float32) for a in ANCHORS for k in ('L','R')}
    times = np.empty(nf); t0 = None
    print(f"reading {nf} frames @ {H}x{W}...")
    for f in range(nf):
        im.get_frame(f*decim)
        fr = np.array(im.final, np.float32, copy=True).reshape(H, W)
        ti = im.frame_info.time; t0 = ti if t0 is None else t0; times[f] = (ti-t0).total_seconds()
        for a in ANCHORS:
            for k in ('L','R'):
                cx, cy = cen[a][k][f]; icx, icy = int(round(cx)), int(round(cy))
                x0,y0 = icx-HW, icy-HW; x1,y1 = x0+WIN, y0+WIN
                sx0,sy0 = max(0,x0),max(0,y0); sx1,sy1 = min(W,x1),min(H,y1)
                if sx1>sx0 and sy1>sy0:
                    crops[(a,k)][f, sy0-y0:sy1-y0, sx0-x0:sx1-x0] = fr[sy0:sy1, sx0:sx1]
        if f % 2000 == 0: print(f"  {f}/{nf}", flush=True)
    fps = (nf-1)/(times[-1]-times[0])
    print(f"fps {fps:.2f}\n\n=== near-center breath SNR (trace) | (pixel) ===")
    out = dict(fps=float(fps), nframes=int(nf), anchors=ANCHORS, src=os.path.basename(ats))
    print(f"{'anchor':>9} | {'L trace':>8} {'L pix':>6} | {'R trace':>8} {'R pix':>6}")
    for a in ANCHORS:
        row = {}
        for k in ('L','R'):
            tr, snrv, pks = best_trace(crops[(a,k)], fps)
            out[f"{a}_{k}_trace"] = tr; out[f"{a}_{k}_snr"] = snrv; out[f"{a}_{k}_pixsnr"] = pks
            row[k] = (snrv, pks)
        print(f"{a:>9} | {row['L'][0]:8.2f} {row['L'][1]:6.1f} | {row['R'][0]:8.2f} {row['R'][1]:6.1f}")
    outp = os.path.splitext(ats)[0] + "_anchorcompare.mat"
    savemat(outp, out, do_compression=True)
    print(f"\nsaved {outp}")

if __name__ == "__main__":
    main()
