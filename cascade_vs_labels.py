"""cascade_vs_labels.py -- external validation of CASCADE against the hand labels.

RUN WITH cellpose-gpu (needs h5py for the v7.3 .mat; no torch needed):
    %USERPROFILE%\\.conda\\envs\\cellpose-gpu\\python.exe cascade_vs_labels.py

CASCADE was trained on electrophysiology ground truth from other animals and has
never seen this dataset, so agreement with the hand clicks is an EXTERNAL check,
not a fitted score. Nothing here feeds back into the detector.

Read the numbers with the same caveat as before: a "false positive" here only
means CASCADE found something not on the click list, which in an ROI that was
only partly curated may well be a real event. Recall is the more trustworthy
direction.
"""
import json
import os

import h5py
import numpy as np

ROOT = r"D:\Ventral_surface_summary"
TRACE_DIR = os.path.join(ROOT, "_cascade_260722", "traces")
PRED_DIR = os.path.join(ROOT, "_cascade_260722", "predictions")

TOL_FR = 6      # +-0.2 s, as established: rise is 2-4 frames and the GUI snaps
                # clicks to a local max within +-3 frames
FPS = 30.0

# CASCADE infers ACTION POTENTIALS; the hand labels mark CALCIUM EVENTS, and one
# calcium transient contains several spikes. Comparing them one-to-one charges
# CASCADE a false positive for every extra spike inside an event it got right.
# Collapse spikes closer together than this into a single event, timed at the
# first spike -- which is the onset, matching how the events were labelled.
MERGE_S = 0.35


def burst_merge(spk, merge_fr):
    if spk.size == 0:
        return spk
    spk = np.sort(spk)
    keep = np.concatenate(([True], np.diff(spk) > merge_fr))
    return spk[keep]


def read_labels(folder, n_frames=None):
    """roi_spikes(k).spike_idx out of ca_spike_data.mat (v7.3 / HDF5).

    MATLAB v7.3 stores an EMPTY array as a dataset containing its dimensions,
    tagged with the attribute MATLAB_empty=1. Reading it naively yields a
    2-element array of zeros, i.e. a fabricated pair of labels on every ROI that
    was never curated -- which silently inflates the label count and destroys
    any precision estimate. Empties must be detected by the attribute.
    """
    p = os.path.join(folder, "ca_spike_data.mat")
    if not os.path.isfile(p):
        return None
    out = {}
    with h5py.File(p, "r") as f:
        if "roi_spikes" not in f or "spike_idx" not in f["roi_spikes"]:
            return None
        flat = np.array(f["roi_spikes"]["spike_idx"]).ravel()
        for k in range(len(flat)):
            try:
                dset = f[flat[k]]
                if dset.attrs.get("MATLAB_empty", 0) == 1:
                    continue                      # genuinely uncurated / silent
                v = np.array(dset).ravel()
                v = v[np.isfinite(v)].astype(np.int64) - 1     # MATLAB 1-based
                if n_frames is not None:
                    v = v[(v >= 0) & (v < n_frames)]
                if v.size:
                    out[k] = v
            except Exception:
                continue
    return out


def main():
    with open(os.path.join(TRACE_DIR, "index.json")) as f:
        index = json.load(f)

    print(f"{'session':52s} {'ROI':>4s} {'lab':>6s} {'cas':>6s} "
          f"{'P':>5s} {'R':>5s} {'F1':>5s}")
    TP = FP = FN = 0
    per_roi_f1 = []

    for e in index:
        pp = os.path.join(PRED_DIR, e["tag"] + ".npz")
        if not os.path.isfile(pp):
            continue
        labels = read_labels(os.path.join(ROOT, e["session"]), e["n_frames"])
        if not labels:
            continue
        P = np.load(pp, allow_pickle=True)
        disc = P["discrete"]
        if disc.size == 0:
            continue

        tp = fp = fn = 0
        nlab = ncas = 0
        for k, lab in labels.items():
            if k >= len(disc):
                continue          # only ROIs that were actually curated
            det = burst_merge(np.asarray(disc[k], dtype=np.int64),
                              MERGE_S * FPS)
            nlab += lab.size
            ncas += det.size
            if det.size == 0:
                fn += lab.size
                per_roi_f1.append(0.0)
                continue
            D = np.abs(lab[:, None] - det[None, :])
            hit_l = (D.min(axis=1) <= TOL_FR)
            hit_d = (D.min(axis=0) <= TOL_FR)
            a, b, c = int(hit_l.sum()), int((~hit_d).sum()), int((~hit_l).sum())
            tp += a
            fp += b
            fn += c
            p_ = a / max(a + b, 1)
            r_ = a / max(a + c, 1)
            per_roi_f1.append(2 * p_ * r_ / max(p_ + r_, 1e-9))

        if nlab == 0:
            continue
        p_ = tp / max(tp + fp, 1)
        r_ = tp / max(tp + fn, 1)
        f_ = 2 * p_ * r_ / max(p_ + r_, 1e-9)
        print(f"{e['session'][:52]:52s} {len(labels):4d} {nlab:6d} {ncas:6d} "
              f"{p_:5.2f} {r_:5.2f} {f_:5.2f}")
        TP += tp
        FP += fp
        FN += fn

    p_ = TP / max(TP + FP, 1)
    r_ = TP / max(TP + FN, 1)
    print(f"\nPOOLED  TP {TP}  FP {FP}  FN {FN}   precision {p_:.3f}  "
          f"recall {r_:.3f}  F1 {2*p_*r_/max(p_+r_,1e-9):.3f}")
    if per_roi_f1:
        a = np.array(per_roi_f1)
        print(f"per-ROI F1: median {np.median(a):.2f}, "
              f"pct 25/75 {np.percentile(a,25):.2f}/{np.percentile(a,75):.2f}, "
              f">=0.8 in {100*np.mean(a>=0.8):.0f}% of {a.size} ROIs")


if __name__ == "__main__":
    main()
