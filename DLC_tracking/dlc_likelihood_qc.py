"""DLC likelihood QC — pooled likelihood distribution histogram.

A standard quality pass for any DLC analyzed output: pool EVERY keypoint
likelihood (across all bodyparts and all videos) and plot the distribution as a
PERCENTAGE histogram, annotated with total frames and recording time
(time = frames / fps). Reference lines mark common thresholds and report the
fraction of detections below each. A right-shifted mass = good tracking; a left
bump near 0 = frames the model is unsure about (where to refine/re-label).

Run in the dlc310 env:
    conda activate dlc310
    python dlc_likelihood_qc.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ========================= USER-EDITABLE =========================
# A single DLC CSV, or a FOLDER to scan for all DLC analyzed CSVs.
CSV_OR_DIR  = r"D:\260615_thermalNbasler\260615_thermalNbasler-RZ-2026-06-15\videos"
FPS_DEFAULT = 400.0          # used only if no timestamps.csv sits next to a CSV
BINS        = 50
THRESHOLDS  = [0.6, 0.9]     # dashed reference lines + reported % below each
PER_BODYPART = True          # print per-bodypart median + %<0.6 to console
SAVE_PNG    = True           # save figure next to the data
# ================================================================


def find_csvs(path):
    """A single CSV, or all DLC analyzed CSVs in a folder (skip timestamps etc.)."""
    if os.path.isfile(path):
        return [path]
    cands = glob.glob(os.path.join(path, "*.csv"))
    return sorted(c for c in cands if "DLC" in os.path.basename(c))


def fps_for(csv):
    """fps from a sibling timestamps.csv (ns col), else FPS_DEFAULT. Never the header."""
    ts = os.path.join(os.path.dirname(csv), "timestamps.csv")
    if os.path.isfile(ts):
        try:
            a = pd.read_csv(ts, header=None).to_numpy()
            return (a.shape[0] - 1) / ((a[-1, 1] - a[0, 1]) / 1e9)
        except Exception:
            pass
    return FPS_DEFAULT


def load_likelihood(csv):
    """Return a [frames x bodyparts] DataFrame of likelihood columns."""
    df = pd.read_csv(csv, header=[0, 1, 2], index_col=0)
    mask = df.columns.get_level_values(-1) == "likelihood"
    lik = df.loc[:, mask]
    lik.columns = lik.columns.get_level_values(-2)   # bodypart names
    return lik


def main():
    csvs = find_csvs(CSV_OR_DIR)
    assert csvs, f"no DLC CSVs found in {CSV_OR_DIR}"

    pooled = []
    per_bp = {}
    total_frames = 0
    total_time = 0.0
    for csv in csvs:
        lik = load_likelihood(csv)
        fps = fps_for(csv)
        nfr = lik.shape[0]
        total_frames += nfr
        total_time += nfr / fps
        v = lik.to_numpy().ravel()
        pooled.append(v[~np.isnan(v)])
        for bp in lik.columns:
            per_bp.setdefault(bp, []).append(lik[bp].to_numpy())
        print(f"{os.path.basename(csv)}: {nfr} fr @ {fps:.1f} fps -> {nfr/fps:.1f} s")

    pooled = np.concatenate(pooled)
    n = pooled.size
    med = float(np.median(pooled))

    # ---- percentage histogram ----
    fig, ax = plt.subplots(figsize=(8.5, 5))
    w = np.ones(n) * 100.0 / n                   # weights -> bars in % of detections
    ax.hist(pooled, bins=BINS, range=(0, 1), weights=w,
            color=[0.30, 0.50, 0.90], edgecolor="white", linewidth=0.5)
    ax.axvline(med, color="k", ls="-", lw=1.5)
    ax.text(med, ax.get_ylim()[1] * 0.98, f" median {med:.2f}",
            va="top", fontsize=9)
    for th in THRESHOLDS:
        frac = 100.0 * np.mean(pooled < th)
        ax.axvline(th, color="r", ls="--", lw=1)
        ax.text(th, ax.get_ylim()[1] * 0.88, f" <{th:g}: {frac:.1f}%",
                color="r", va="top", fontsize=9)

    ax.set_xlim(0, 1)
    ax.set_xlabel("likelihood")
    ax.set_ylabel("% of keypoint detections")
    ax.grid(alpha=0.3)
    ax.set_title(f"DLC likelihood — {len(csvs)} video(s), {total_frames} frames, "
                 f"{total_time:.1f} s\n{n} keypoint detections "
                 f"({len(per_bp)} bodyparts)")
    fig.tight_layout()

    # ---- console summary ----
    print(f"\nPOOLED: {n} detections | median={med:.3f}")
    for th in THRESHOLDS:
        print(f"  % < {th:g}: {100.0*np.mean(pooled < th):.2f}%")
    if PER_BODYPART:
        print("\nper-bodypart  median  %<0.6:")
        for bp, arrs in per_bp.items():
            vb = np.concatenate(arrs)
            vb = vb[~np.isnan(vb)]
            print(f"  {bp:10s} {np.median(vb):.2f}   {100.0*np.mean(vb < 0.6):.1f}%")

    if SAVE_PNG:
        base = CSV_OR_DIR if os.path.isdir(CSV_OR_DIR) else os.path.dirname(CSV_OR_DIR)
        out = os.path.join(base, "dlc_likelihood_qc.png")
        fig.savefig(out, dpi=150)
        print(f"\nsaved {out}")

    plt.show()


if __name__ == "__main__":
    main()
