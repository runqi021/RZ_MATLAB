"""pmd_denoise_traces.py -- patch-wise low-rank denoising, then ROI traces.

RUN WITH dlc310 (needs working BLAS + tifffile):
    %USERPROFILE%\\.conda\\envs\\dlc310\\python.exe pmd_denoise_traces.py

THE IDEA (Buchanan et al. 2019, penalised matrix decomposition). Inside a small
spatial patch, calcium signal is LOW-RANK (few cells, smooth footprints) and
temporally CORRELATED (tau ~ 0.5-0.9 s here, so lag-1 autocorrelation ~0.95).
Photon shot noise is full-rank and WHITE (lag-1 autocorrelation ~0). So take the
SVD of each patch and keep only the components whose TEMPORAL singular vector is
autocorrelated; discard the rest. That removes noise using redundancy ACROSS
PIXELS AND TIME TOGETHER -- the axis that neither within-mask reweighting
(measured gain 0.99-1.05x, i.e. none) nor temporal deconvolution exploits.

Patches overlap by half and are recombined with a raised-cosine window so no
patch seams appear in the result.

Reports, per ROI, the transient SNR of the trace extracted from the denoised
movie against the same trace from the raw movie. That ratio is the whole point:
if it is ~1 this avenue is closed too, and the honest answer is a detectability
gate rather than a better detector.
"""
import os
import sys
import time

import numpy as np
import tifffile

# ---------------------------- USER PARAMETERS ----------------------------
SESSIONS = [
    r"D:\260721_Sert_soma_G8s\phys\roi7_1.8x_x-1000-1000_z0_16lp_3000f_00001",
    r"D:\Ventral_surface_summary\Vgat\0323\cell3\-500-990-z-30_8x_3000f_00001",
]
PATCH = 32          # px, patch side
STRIDE = 16         # px, 50% overlap
MAX_T = 3000        # frames processed (memory cap)
ACORR_MIN = 0.25    # keep components with lag-1 temporal autocorrelation above this
MAX_RANK = 20       # cap components tested per patch
# -------------------------------------------------------------------------


def robust_sn(x):
    return np.median(np.abs(np.diff(x))) * 1.4826 / np.sqrt(2)


def tsnr(x):
    return (np.percentile(x, 99) - np.median(x)) / max(robust_sn(x), 1e-9)


def find_movie(folder):
    tifs = [f for f in os.listdir(folder) if f.endswith(".tif")]
    mc = [f for f in tifs if f.endswith("_MC_MC.tif")] or \
         [f for f in tifs if f.endswith("_MC.tif")]
    return os.path.join(folder, sorted(mc)[-1]) if mc else None


def lag1(v):
    v = v - v.mean()
    d = (v * v).sum()
    return float((v[:-1] * v[1:]).sum() / d) if d > 0 else 0.0


def denoise(mov, rng):
    """mov [T,H,W] float32 -> denoised, by patch-wise rank selection."""
    T, H, W = mov.shape
    out = np.zeros((T, H, W), dtype=np.float32)
    wsum = np.zeros((H, W), dtype=np.float32)

    win1 = np.hanning(PATCH + 2)[1:-1].astype(np.float32)
    win2d = np.outer(win1, win1)

    ys = list(range(0, max(H - PATCH, 0) + 1, STRIDE))
    xs = list(range(0, max(W - PATCH, 0) + 1, STRIDE))
    if ys[-1] + PATCH < H:
        ys.append(H - PATCH)
    if xs[-1] + PATCH < W:
        xs.append(W - PATCH)

    kept = []
    for y in ys:
        for x in xs:
            P = mov[:, y:y + PATCH, x:x + PATCH].reshape(T, -1)   # [T, px]
            mu = P.mean(axis=0, keepdims=True)
            Pc = P - mu
            k = min(MAX_RANK, min(Pc.shape) - 1)
            if k < 1:
                continue
            # Randomized SVD: only the top k components are needed, and a full
            # economy SVD here costs O(T*px*min(T,px)) -- ~30 min per session
            # over ~700 patches. The range finder is ~100x cheaper.
            om = rng.standard_normal((Pc.shape[1], k + 5)).astype(np.float32)
            Y = Pc @ om
            # POWER ITERATIONS are mandatory here. Shot noise makes the singular
            # spectrum decay slowly, and a single-pass range finder then returns
            # directions dominated by noise -- every component looks white, the
            # autocorrelation test rejects all of them, and the "denoised" movie
            # collapses to the per-pixel mean (measured: mean kept rank 0.1/20,
            # tSNR -> 0). Two iterations sharpen the spectrum enough to separate.
            for _ in range(2):
                Y = Pc @ (Pc.T @ Y)
            Q, _ = np.linalg.qr(Y)                    # [T, k+5]
            B = Q.T @ Pc                              # [k+5, px]
            Ub, S, Vt = np.linalg.svd(B, full_matrices=False)
            U = Q @ Ub                                # [T, k+5]
            keep = [i for i in range(min(k, len(S))) if lag1(U[:, i]) >= ACORR_MIN]
            kept.append(len(keep))
            if keep:
                rec = (U[:, keep] * S[keep]) @ Vt[keep, :] + mu
            else:
                rec = np.repeat(mu, T, axis=0)
            rec = rec.reshape(T, PATCH, PATCH) * win2d
            out[:, y:y + PATCH, x:x + PATCH] += rec
            wsum[y:y + PATCH, x:x + PATCH] += win2d
    wsum[wsum == 0] = 1.0
    out /= wsum
    return out, (np.mean(kept) if kept else 0.0)


def main():
    for folder in SESSIONS:
        rel = os.path.basename(folder)
        mpath = find_movie(folder)
        mask_p = os.path.join(folder, "roi_maskL.npy")
        if mpath is None or not os.path.isfile(mask_p):
            print(f"{rel}: missing movie or roi_maskL.npy "
                  f"(run export_masks.py) -- skip")
            continue
        maskL = np.load(mask_p)

        t0 = time.time()
        with tifffile.TiffFile(mpath) as tf:
            T = min(MAX_T, len(tf.pages))
            H, W = tf.pages[0].shape
            mov = np.empty((T, H, W), dtype=np.float32)
            for i in range(T):
                mov[i] = tf.pages[i].asarray()
        if maskL.shape != (H, W):
            print(f"{rel}: mask {maskL.shape} != frame {(H, W)} -- skip")
            continue
        print(f"{rel}: {T} x {H} x {W} loaded in {time.time()-t0:.0f}s", flush=True)

        t0 = time.time()
        den, mean_rank = denoise(mov, np.random.default_rng(0))
        print(f"  denoised in {time.time()-t0:.0f}s, mean kept rank "
              f"{mean_rank:.1f} of {MAX_RANK}", flush=True)

        labels = np.unique(maskL)
        labels = labels[labels != 0]
        gains, raws, dens = [], [], []
        for lb in labels:
            m = (maskL == lb).ravel()
            raw = mov.reshape(T, -1)[:, m].mean(axis=1)
            dnz = den.reshape(T, -1)[:, m].mean(axis=1)
            r, d = tsnr(raw), tsnr(dnz)
            raws.append(r)
            dens.append(d)
            gains.append(d / max(r, 1e-9))
        g = np.array(gains)
        print(f"  {len(labels)} ROI   tSNR raw {np.median(raws):.1f} -> "
              f"denoised {np.median(dens):.1f}   GAIN median {np.median(g):.2f}x "
              f"(p25 {np.percentile(g,25):.2f}, p75 {np.percentile(g,75):.2f})")
        np.savez_compressed(os.path.join(folder, "pmd_traces.npz"),
                            raw=np.array(raws), den=np.array(dens),
                            gain=g, labels=labels)
        sys.stdout.flush()


if __name__ == "__main__":
    main()
