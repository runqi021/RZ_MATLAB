"""roi_weighted_traces.py -- SNR-optimal ROI traces from the registered movie.

RUN WITH cellpose-gpu (tifffile + h5py):
    %USERPROFILE%\\.conda\\envs\\cellpose-gpu\\python.exe roi_weighted_traces.py

WHY. The existing traces are a FLAT MEAN over every pixel in the cpSAM mask, so a
dim edge pixel counts as much as a bright somatic one. Signal adds coherently
across pixels while noise does not, so weighting each pixel by how much signal it
actually carries raises SNR before detection ever starts. This is what fixes
"misses AND false positives at the same time" -- both are symptoms of a trace
whose SNR is too low for ANY threshold to separate, and no amount of detector
tuning addresses that.

Weights come from the first spatial principal component of the ROI's
pixel x time matrix (mean-subtracted), sign-fixed to be positive and clipped at
zero: the standard CNMF-style spatial footprint. Pixels that co-vary with the
ROI's dominant temporal mode get weight; pixels carrying only noise get ~0.

Also computed, from the same pass, the spatial features for artifact rejection
(the original motion question, which sx/sy could not answer -- motion AUC 0.50):
  core   eroded mask         -- a real event fills the soma
  ring   dilated minus mask  -- an optical artifact moves the surround too
  grad   projection of the per-frame image change on the local image gradient;
         a small translation changes F by grad(F).delta, so a motion artifact's
         spatial signature IS the gradient, by construction.

Outputs <session>/roi_weighted.npz : trace_w [N x T], trace_flat, core, ring,
grad, weights (sparse per ROI), snr_gain.
"""
import json
import os
import sys

import numpy as np
import tifffile

ROOT = r"D:\Ventral_surface_summary"
OUT_NAME = "roi_weighted.npz"
DILATE = 3          # px, ring thickness outside the mask
ERODE = 1           # px, core erosion
MAX_FRAMES = None   # cap for a quick test; None = whole movie
ONLY = "Vgat"        # substring filter; None = all sessions


def robust_sn(x):
    return np.median(np.abs(np.diff(x))) * 1.4826 / np.sqrt(2)


def find_files(folder):
    """Prefer the doubly-registered movie when present, and its matching mask."""
    tifs = [f for f in os.listdir(folder) if f.endswith(".tif")]
    mc = [f for f in tifs if f.endswith("_MC_MC.tif")] or \
         [f for f in tifs if f.endswith("_MC.tif")]
    sam = [f for f in os.listdir(folder) if f.endswith("cpSAM_output.mat")]
    if not mc or not sam:
        return None, None
    # pair the mask to the same registration stage as the movie
    mc_stem = sorted(mc)[-1].replace(".tif", "")
    same = [s for s in sam if s.startswith(mc_stem)]
    return os.path.join(folder, sorted(mc)[-1]), \
        os.path.join(folder, (same or sorted(sam))[-1])


def read_mask(folder):
    p = os.path.join(folder, "roi_maskL.npy")
    if not os.path.isfile(p):
        return None
    return np.load(p)


def binary_dilate(mask, k):
    out = mask.copy()
    for _ in range(k):
        p = np.zeros_like(out)
        p[1:, :] |= out[:-1, :]
        p[:-1, :] |= out[1:, :]
        p[:, 1:] |= out[:, :-1]
        p[:, :-1] |= out[:, 1:]
        out = out | p
    return out


def binary_erode(mask, k):
    out = mask.copy()
    for _ in range(k):
        p = np.ones_like(out)
        p[1:, :] &= out[:-1, :]
        p[:-1, :] &= out[1:, :]
        p[:, 1:] &= out[:, :-1]
        p[:, :-1] &= out[:, 1:]
        out = out & p
    return out


def process(folder, rel):
    tif_path, sam_path = find_files(folder)
    if tif_path is None:
        print(f"  {rel}: no MC tif / cpSAM -- skip")
        return
    maskL = read_mask(folder)
    if maskL is None:
        print(f"  {rel}: no maskL -- skip")
        return

    labels = np.unique(maskL)
    labels = labels[labels != 0]
    N = len(labels)

    with tifffile.TiffFile(tif_path) as tf:
        n_frames = len(tf.pages)
        H, W = tf.pages[0].shape
        T = n_frames if MAX_FRAMES is None else min(MAX_FRAMES, n_frames)
        if maskL.shape != (H, W):
            print(f"  {rel}: mask {maskL.shape} != frame {(H, W)} -- skip")
            return

        # pixel index sets per ROI, plus core / ring
        idx, cidx, ridx = [], [], []
        for lb in labels:
            m = (maskL == lb)
            core = binary_erode(m, ERODE)
            if core.sum() < 4:
                core = m
            ring = binary_dilate(m, DILATE) & ~binary_dilate(m, 1)
            ring &= (maskL == 0)                 # exclude other cells
            idx.append(np.flatnonzero(m.ravel()))
            cidx.append(np.flatnonzero(core.ravel()))
            ridx.append(np.flatnonzero(ring.ravel()))

        allpix = np.unique(np.concatenate(idx))
        pos = {p: i for i, p in enumerate(allpix)}
        buf = np.zeros((len(allpix), T), dtype=np.float32)
        core_t = np.zeros((N, T), dtype=np.float32)
        ring_t = np.zeros((N, T), dtype=np.float32)
        grad_t = np.zeros((N, T), dtype=np.float32)

        mean_img = np.zeros((H, W), dtype=np.float64)
        for i in range(T):
            fr = tf.pages[i].asarray().astype(np.float32)
            mean_img += fr
            flat = fr.ravel()
            buf[:, i] = flat[allpix]
            for k in range(N):
                core_t[k, i] = flat[cidx[k]].mean() if cidx[k].size else 0.0
                ring_t[k, i] = flat[ridx[k]].mean() if ridx[k].size else 0.0
        mean_img /= T

    gy, gx = np.gradient(mean_img)
    gmag = (gx ** 2 + gy ** 2).ravel()

    trace_w = np.zeros((N, T), dtype=np.float32)
    trace_f = np.zeros((N, T), dtype=np.float32)
    weights = []
    gain = np.zeros(N)

    for k in range(N):
        rows = np.array([pos[p] for p in idx[k]])
        Y = buf[rows, :]                       # [px x T]
        flat = Y.mean(axis=0)
        trace_f[k] = flat

        # Matched-filter weights: w_i ~ cov(y_i, m) / sigma_i^2, the SNR-optimal
        # linear combination when per-pixel noise is independent. Anchored to the
        # ROI's own mean trace m, so it cannot run away onto a noise mode.
        #
        # NOT the first spatial PC: on a noise-dominated ROI, PC1 IS the noise
        # mode, and clipping it concentrates weight on a few pixels, which cuts
        # effective averaging and makes the trace noisier than the flat mean
        # (measured: 0.14-0.78x, i.e. worse).
        m = flat - flat.mean()
        vm = float((m * m).mean())
        if vm <= 0:
            u = np.ones(Y.shape[0]) / Y.shape[0]
        else:
            Yc = Y - Y.mean(axis=1, keepdims=True)
            cov = (Yc * m).mean(axis=1)                       # cov(y_i, m)
            sig = np.median(np.abs(np.diff(Y, axis=1)), axis=1) * 1.4826 / np.sqrt(2)
            sig = np.maximum(sig, np.median(sig) * 1e-3 + 1e-9)
            u = np.clip(cov, 0, None) / (sig ** 2)
            if not np.isfinite(u).all() or u.sum() <= 0:
                u = np.ones(Y.shape[0])
            u = u / u.sum()
        wt = u @ Y
        trace_w[k] = wt
        weights.append(u.astype(np.float32))

        # gradient projection: how much of this ROI's frame-to-frame change
        # looks like a translation of the mean image
        gp = gmag[idx[k]]
        gp = gp / max(gp.sum(), 1e-9)
        grad_t[k] = gp @ np.abs(np.diff(Y, axis=1, prepend=Y[:, :1]))

        # transient SNR: how far the active tail rises above baseline, in noise
        # units. Scale-free, so flat and weighted traces are comparable.
        def tsnr(z):
            return (np.percentile(z, 99) - np.median(z)) / max(robust_sn(z), 1e-9)
        gain[k] = tsnr(wt) / max(tsnr(flat), 1e-9)

    np.savez_compressed(os.path.join(folder, OUT_NAME),
                        trace_w=trace_w, trace_flat=trace_f,
                        core=core_t, ring=ring_t, grad=grad_t,
                        snr_gain=gain, labels=labels, T=T)
    print(f"  {rel[:54]:54s} {N:4d} ROI x {T:5d} fr   median SNR gain "
          f"{np.median(gain):.2f}x")


def main():
    sessions = []
    for dp, _dn, fn in os.walk(ROOT):
        if os.sep + "_" in dp:
            continue
        if any(f.endswith("cpSAM_output.mat") for f in fn) and \
           any(f.endswith(".tif") for f in fn):
            sessions.append(dp)
    sessions.sort()
    if ONLY:
        sessions = [s for s in sessions if ONLY in s]
    print(f"{len(sessions)} session(s)")
    for s in sessions:
        rel = os.path.relpath(s, ROOT)
        try:
            process(s, rel)
        except Exception as ex:
            print(f"  {rel}: FAILED {type(ex).__name__}: {ex}")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
