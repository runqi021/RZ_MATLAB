"""split_half_validate.py -- label-free validation, detectability gate, and the
objective function for automated parameter search.

RUN WITH dlc310 (working BLAS + tifffile):
    %USERPROFILE%\\.conda\\envs\\dlc310\\python.exe split_half_validate.py [ROOT ...]

THE CONSTRUCTION. Split each ROI's pixels into two INTERLEAVED halves
(checkerboard on x+y, so both halves cover the same territory and see the same
neuropil). Both halves observe the SAME calcium signal; their shot noise is
INDEPENDENT because they are different photons. Therefore:

  * an event present in both halves is real; an event in one only is noise
  * corr(A,B) estimates SNR directly: for s+n_A and s+n_B with independent
    noise, corr = var(s)/(var(s)+var(n))
  * precision and recall are estimable WITH NO LABELS, on every ROI, including
    the ~1500 that were never curated

This yields three things at once: the validation metric, the per-ROI
detectability gate, and the objective to optimise detector parameters against.

WHAT IT DOES NOT COVER. Both halves move together, so a MOTION artifact appears
in both and passes split-half validation. Split-half certifies against shot
noise only. Motion needs the spatial tests (ring / image-gradient projection) in
roi_weighted_traces.py. The two are complementary, not redundant.

Per ROI it reports:
  rel      corr(A,B) on the high-passed traces, Spearman-Brown corrected to the
           full-pixel trace -- the detectability index
  th*      threshold (in units of each trace's own noise) maximising split-half
           agreement
  F1_sh    split-half F1 at th*, a label-free estimate of detection quality
  n        events on the FULL trace at th*
"""
import os
import sys
import time

import numpy as np
import tifffile

# ---------------------------- USER PARAMETERS ----------------------------
ROOTS = sys.argv[1:] or [
    r"D:\260721_Sert_soma_G8s\phys",
    r"D:\Ventral_surface_summary",
]
FPS = 30.0
BASE_WIN_S = 20.0        # rolling baseline window for dF/F
TH_GRID = np.arange(1.0, 6.01, 0.25)
TOL_FR = 6               # +-0.2 s matching tolerance
CONSENSUS_TH = 2.0       # threshold at which both halves must agree to seed the reference
MIN_DIST_S = 0.15
MAX_T = 6000             # frame cap per session
REL_MIN = 0.10           # below this an ROI is gated out as undetectable
OUT_NAME = "split_half.npz"
# -------------------------------------------------------------------------


def robust_sn(x):
    return np.median(np.abs(np.diff(x))) * 1.4826 / np.sqrt(2)


def rolling_baseline(x, w):
    """Rolling 10th percentile: tracks the floor without the post-event
    undershoot a rolling MEDIAN digs after a real transient."""
    n = x.size
    step = max(1, w // 4)
    idx = np.arange(0, n, step)
    vals = np.empty(idx.size)
    for j, i in enumerate(idx):
        lo = max(0, i - w // 2)
        hi = min(n, i + w // 2)
        vals[j] = np.percentile(x[lo:hi], 10)
    return np.interp(np.arange(n), idx, vals)


def dff(f, w):
    b = rolling_baseline(f, w)
    b = np.maximum(b, 1e-6)
    return (f - b) / b


def ac_tau(x, fps):
    """AR(1) tau from autocovariance at lags >= 1 (lag 0 carries white noise)."""
    x = x - np.mean(x)
    L = int(min(max(3, round(0.6 * fps)), len(x) // 4))
    c = np.correlate(x, x, mode="full")[len(x) - 1:len(x) - 1 + L + 1] / len(x)
    c1 = c[1:]
    ok = c1 > 0
    if ok.sum() < 3:
        return 0.5
    lags = np.arange(1, len(c1) + 1)[ok]
    p = np.polyfit(lags, np.log(c1[ok]), 1)
    g = np.exp(p[0])
    if not np.isfinite(g) or g <= 0 or g >= 1:
        return 0.5
    tau = -1.0 / (fps * np.log(g))
    return float(np.clip(tau, 0.1, 3.0))


def kernel(tau, rise_s, fps):
    """Zero-mean matched kernel: the DC component is removed so the filter is
    blind to slow drift by construction. Without this the statistic inherits the
    low-frequency content of dF/F, positive and negative excursions become
    near-symmetric, and the sign-flip null reports FDR ~0.75 on clean traces.
    Subtracting a rolling percentile instead is NOT a fix -- it is envelope
    subtraction and carves negative lobes out of the transients themselves."""
    t = np.arange(0, int(round(2.5 * tau * fps)) + 1) / fps
    tr = max(rise_s, 1.0 / fps) / 2.2
    k = (1 - np.exp(-t / tr)) * np.exp(-t / tau)
    k = k - k.mean()
    return k / max(np.linalg.norm(k), 1e-9)


def matched(x, k):
    """Correlate with the kernel; peaks land on the event ONSET.
    Explicit 'full' indexing -- 'same' would recentre by floor(K/2)."""
    K = k.size
    yf = np.convolve(x - np.median(x), k[::-1])
    return yf[K - 1:K - 1 + x.size]


def find_peaks(z, thr, min_dist):
    """Local maxima above thr, greedily enforcing a refractory."""
    cand = np.flatnonzero((z[1:-1] >= z[:-2]) & (z[1:-1] > z[2:]) &
                          (z[1:-1] >= thr)) + 1
    if cand.size == 0:
        return cand
    order = cand[np.argsort(-z[cand])]
    keep = []
    taken = np.zeros(z.size, bool)
    for c in order:
        lo, hi = max(0, c - min_dist), min(z.size, c + min_dist + 1)
        if not taken[lo:hi].any():
            keep.append(c)
            taken[c] = True
    return np.sort(np.array(keep, dtype=int))


def match_count(a, b, tol):
    if a.size == 0 or b.size == 0:
        return 0
    d = np.abs(a[:, None] - b[None, :])
    return int((d.min(axis=1) <= tol).sum())


def find_movie(folder):
    tifs = [f for f in os.listdir(folder) if f.endswith(".tif")]
    mc = [f for f in tifs if f.endswith("_MC_MC.tif")] or \
         [f for f in tifs if f.endswith("_MC.tif")]
    return os.path.join(folder, sorted(mc)[-1]) if mc else None


def process(folder, rel_name):
    mpath = find_movie(folder)
    mp = os.path.join(folder, "roi_maskL.npy")
    if mpath is None or not os.path.isfile(mp):
        return None
    maskL = np.load(mp)
    labels = np.unique(maskL)
    labels = labels[labels != 0]
    if labels.size == 0:
        return None

    yy, xx = np.indices(maskL.shape)
    parity = ((yy + xx) % 2).ravel()
    flatmask = maskL.ravel()

    idxA, idxB, idxF = [], [], []
    for lb in labels:
        px = np.flatnonzero(flatmask == lb)
        a = px[parity[px] == 0]
        b = px[parity[px] == 1]
        if a.size < 4 or b.size < 4:
            a = px[: px.size // 2]
            b = px[px.size // 2:]
        idxA.append(a)
        idxB.append(b)
        idxF.append(px)

    allpx = np.unique(np.concatenate(idxF))
    posmap = -np.ones(flatmask.size, dtype=np.int64)
    posmap[allpx] = np.arange(allpx.size)

    with tifffile.TiffFile(mpath) as tf:
        T = min(MAX_T, len(tf.pages))
        H, W = tf.pages[0].shape
        if maskL.shape != (H, W):
            return None
        buf = np.empty((allpx.size, T), dtype=np.float32)
        for i in range(T):
            buf[:, i] = tf.pages[i].asarray().ravel()[allpx]

    w = int(BASE_WIN_S * FPS)
    md = max(1, int(round(MIN_DIST_S * FPS)))
    out = []
    for k, lb in enumerate(labels):
        fA = buf[posmap[idxA[k]], :].mean(axis=0)
        fB = buf[posmap[idxB[k]], :].mean(axis=0)
        fF = buf[posmap[idxF[k]], :].mean(axis=0)
        xA, xB, xF = dff(fA, w), dff(fB, w), dff(fF, w)

        # reliability on high-passed traces, so slow drift shared by both
        # halves cannot masquerade as signal agreement
        hp = lambda v: v - rolling_baseline(v, int(2.0 * FPS))  # noqa: E731
        hA, hB = hp(xA), hp(xB)
        sa, sb = hA.std(), hB.std()
        r = float((hA - hA.mean()).dot(hB - hB.mean()) / (len(hA) * sa * sb)) \
            if sa > 0 and sb > 0 else 0.0
        # Spearman-Brown: each half has half the pixels, so correct the
        # reliability up to the full-pixel trace
        rel = 2 * r / (1 + r) if r > -1 else 0.0

        tau = ac_tau(xF, FPS)
        ker = kernel(tau, 0.1, FPS)
        # Normalise by the noise of the INPUT trace, not of the filtered one.
        # robust_sn uses frame-to-frame differences, which is valid only for
        # WHITE noise; matched filtering makes the noise strongly correlated, so
        # applied to the output it underestimates sigma severely and a nominal
        # "3 sigma" lands deep inside the noise (measured: ~85 detections per
        # ROI-minute and sign-flip FDR ~0.85 even on traces where detection is
        # known to work). For a unit-norm kernel, filtered white-noise sigma
        # equals the input sigma, so this is the correct scale.
        snA, snB, snF = robust_sn(xA), robust_sn(xB), robust_sn(xF)
        zA = matched(xA, ker) / max(snA, 1e-9)
        zB = matched(xB, ker) / max(snB, 1e-9)
        zF = matched(xF, ker) / max(snF, 1e-9)

        # CONSENSUS reference: events the two independent halves BOTH find at a
        # permissive threshold. Independent shot noise almost never coincides
        # within +-tol, so this set is high-precision (though incomplete -- it
        # misses weak real events, so recall against it reads high in absolute
        # terms but is a valid basis for COMPARING thresholds).
        pa0 = find_peaks(zA, CONSENSUS_TH, md)
        pb0 = find_peaks(zB, CONSENSUS_TH, md)
        ref = np.array([], dtype=int)
        if pa0.size and pb0.size:
            d = np.abs(pa0[:, None] - pb0[None, :])
            ref = pa0[d.min(axis=1) <= TOL_FR]

        # Two nulls, two roles: consensus gives RECALL, the sign-flipped trace
        # gives FDR (calcium cannot go down, so anything found on -zF is false).
        # Agreement alone is not an objective -- uncorrected it pins at the grid
        # floor, chance-corrected it pins at the ceiling. Yield must be in it.
        best = (-1.0, np.nan, 0.0, 0.0)
        for th in TH_GRID:
            pf = find_peaks(zF, th, md)
            nf = find_peaks(-zF, th, md)
            if pf.size == 0:
                continue
            fdr = min(nf.size / pf.size, 1.0)
            rec = (match_count(ref, pf, TOL_FR) / ref.size) if ref.size else 0.0
            score = rec * (1.0 - fdr)
            if score > best[0]:
                best = (score, th, rec, fdr)
        score, th, rec, fdr = best
        n_full = find_peaks(zF, th, md).size if np.isfinite(th) else 0
        out.append((lb, rel, tau, th, score, n_full, T, rec, fdr, ref.size))

    return np.array(out, dtype=float)


def main():
    t0 = time.time()
    rows, names = [], []
    for root in ROOTS:
        sess = []
        for dp, _dn, fn in os.walk(root):
            if os.sep + "_" in dp:
                continue
            if "roi_maskL.npy" in fn and any(f.endswith(".tif") for f in fn):
                sess.append(dp)
        sess.sort()
        print(f"\n=== {root}  ({len(sess)} sessions) ===")
        for s in sess:
            rel_name = os.path.relpath(s, root)
            try:
                r = process(s, rel_name)
            except Exception as ex:
                print(f"  {rel_name[:52]:52s} FAILED {type(ex).__name__}: {ex}")
                continue
            if r is None or r.size == 0:
                continue
            np.savez_compressed(os.path.join(s, OUT_NAME), rows=r)
            good = r[:, 1] >= REL_MIN
            gm = good & np.isfinite(r[:, 3])
            print(f"  {rel_name[:52]:52s} {r.shape[0]:4d} ROI  "
                  f"rel {np.median(r[:,1]):.2f}  usable {good.sum():4d}"
                  f"({100*good.mean():3.0f}%)  "
                  f"th* {np.median(r[gm,3]) if gm.any() else float('nan'):5.2f}  "
                  f"recall {np.median(r[gm,7]) if gm.any() else float('nan'):.2f}  "
                  f"FDR {np.median(r[gm,8]) if gm.any() else float('nan'):.2f}  "
                  f"ref/min {np.median(r[gm,9]/(r[gm,6]/FPS/60)) if gm.any() else float('nan'):.1f}")
            rows.append(r)
            names.append(rel_name)
            sys.stdout.flush()

    if not rows:
        print("nothing processed")
        return
    R = np.vstack(rows)
    good = R[:, 1] >= REL_MIN
    print(f"\n================ POOLED  ({R.shape[0]} ROIs, {time.time()-t0:.0f}s) ================")
    print(f"reliability corr(A,B): pct 10/25/50/75/90 = "
          f"{np.round(np.percentile(R[:,1], [10,25,50,75,90]), 2)}")
    print(f"usable (rel >= {REL_MIN}): {good.sum()} of {R.shape[0]} "
          f"({100*good.mean():.0f}%)")
    gm = good & np.isfinite(R[:, 3])
    print(f"objective recall*(1-FDR): median {np.median(R[gm,4]):.2f}")
    print(f"  recall vs consensus: median {np.median(R[gm,7]):.2f}   "
          f"FDR (sign-flip): median {np.median(R[gm,8]):.2f}")
    print(f"chosen threshold th*: median {np.median(R[gm,3]):.2f}, "
          f"pct 25/75 {np.percentile(R[gm,3],25):.2f}/{np.percentile(R[gm,3],75):.2f}")
    print(f"tau: median {np.median(R[gm,2]):.2f} s")
    print("\nsplit-half F1 by reliability band:")
    for lo, hi in [(0, .1), (.1, .2), (.2, .4), (.4, .6), (.6, 1.01)]:
        m = (R[:, 1] >= lo) & (R[:, 1] < hi)
        if m.sum() == 0:
            continue
        mm = m & np.isfinite(R[:, 3])
        if mm.sum() == 0: continue
        print(f"  rel {lo:.1f}-{hi:.1f}: n={m.sum():4d}  obj {np.median(R[mm,4]):.2f}"
              f"   th* {np.median(R[mm,3]):5.2f}   recall {np.median(R[mm,7]):.2f}"
              f"   FDR {np.median(R[mm,8]):.2f}   events/min "
              f"{np.median(R[mm,5]/(R[mm,6]/FPS/60)):5.1f}")


if __name__ == "__main__":
    main()
