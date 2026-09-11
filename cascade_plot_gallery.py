"""cascade_plot_gallery.py -- visualise CASCADE's GCaMP8s event detection.

RUN WITH dlc310:
    %USERPROFILE%\\.conda\\envs\\dlc310\\python.exe cascade_plot_gallery.py

Draws a random sample of ROIs across all sessions. Per panel:
    grey line      dF/F
    filled orange  CASCADE inferred spike rate (spikes/s), on its own scale
    black ticks    discrete inferred spikes
    title          session, ROI, noise level nu, inferred rate

The sample is drawn from ROIs CASCADE assigns any activity to, with a fixed
seed so reruns show the same ROIs.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

ROOT = r"D:\Ventral_surface_summary"
TRACE_DIR = os.path.join(ROOT, "_cascade_260722", "traces")
PRED_DIR = os.path.join(ROOT, "_cascade_260722", "predictions")
OUT_DIR = os.path.join(ROOT, "_cascade_260722", "gallery")

N_TRACES = 100
PER_PAGE = 20
N_COLS = 2
WIN_SEC = 30
SEED = 260722
MIN_SPIKES = 1          # only sample ROIs with at least this many inferred spikes

# One mark per CALCIUM TRANSIENT, at its onset. CASCADE emits many APs inside a
# single transient, so they have to be grouped -- but a fixed time window is the
# wrong rule: it splits long transients and merges genuinely separate fast ones.
# Segment on the trace instead. Consecutive APs belong to the same transient
# while dF/F has NOT fallen back toward baseline between them; a drop below
# RETURN_FRAC of the intervening peak (measured from the event onset) starts a
# new one. No time constant, no tuned window.
RETURN_FRAC = 0.5
MIN_GAP_S = 0.1      # hard refractory, below which nothing counts as separate


def to_events(spk, x, min_gap_fr):
    """AP times -> onsets of distinct calcium transients."""
    if spk.size == 0:
        return spk
    spk = np.sort(spk)
    ons = [spk[0]]
    for a, b in zip(spk[:-1], spk[1:]):
        if b - a < min_gap_fr:
            continue
        seg = x[a:b + 1]
        if seg.size < 3:
            continue
        base = x[ons[-1]]
        peak = np.max(seg)
        amp = peak - base
        # trace returned toward baseline between the two APs -> new transient
        if amp <= 0 or np.min(seg[np.argmax(seg):]) < base + RETURN_FRAC * amp:
            ons.append(b)
    return np.array(ons, dtype=np.int64)


def load_all():
    with open(os.path.join(TRACE_DIR, "index.json")) as f:
        index = json.load(f)
    items = []
    for e in index:
        pp = os.path.join(PRED_DIR, e["tag"] + ".npz")
        if not os.path.isfile(pp):
            continue
        items.append(e)
    return items


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    index = load_all()
    print(f"{len(index)} sessions with predictions")

    # build the ROI pool
    pool = []
    cache = {}
    for e in index:
        P = np.load(os.path.join(PRED_DIR, e["tag"] + ".npz"), allow_pickle=True)
        sp = P["spike_prob"]
        disc = P["discrete"] if P["discrete"].size else None
        cache[e["tag"]] = (sp, P["noise"], disc)
        for k in range(sp.shape[0]):
            n = len(disc[k]) if disc is not None else int(np.nansum(sp[k]) > 0)
            if n >= MIN_SPIKES:
                pool.append((e["tag"], e["session"], k))
    print(f"{len(pool)} ROIs with >= {MIN_SPIKES} inferred spike(s)")

    rng = np.random.default_rng(SEED)
    pick = [pool[i] for i in rng.permutation(len(pool))[:min(N_TRACES, len(pool))]]

    traces = {}
    n_page = int(np.ceil(len(pick) / PER_PAGE))
    pdf = PdfPages(os.path.join(OUT_DIR, "cascade_gallery.pdf"))
    tot_sp = 0

    for p in range(n_page):
        sel = pick[p * PER_PAGE:(p + 1) * PER_PAGE]
        n_row = int(np.ceil(len(sel) / N_COLS))
        fig, axes = plt.subplots(n_row, N_COLS, figsize=(19, 12))
        axes = np.atleast_1d(axes).ravel()

        for j, (tag, sess, k) in enumerate(sel):
            if tag not in traces:
                traces[tag] = np.load(os.path.join(TRACE_DIR, tag + ".npz"))["dFF"]
            x = traces[tag][k]
            sp, noise, disc = cache[tag]
            r = sp[k]
            # discrete spike positions come back as floats; indices must be int
            aps = (np.rint(np.asarray(disc[k], dtype=np.float64)).astype(np.int64)
                   if disc is not None else np.array([], dtype=np.int64))
            spk = to_events(aps, x, MIN_GAP_S * 30.0)
            tot_sp += len(spk)

            T = x.size
            t = np.arange(T) / 30.0
            if WIN_SEC > 0 and T / 30.0 > WIN_SEC:
                w = int(WIN_SEC * 30)
                c = int(spk[len(spk) // 2]) if len(spk) else T // 2
                lo = max(0, min(T - w, c - w // 2))
                hi = lo + w
            else:
                lo, hi = 0, T

            ax = axes[j]
            ax.plot(t[lo:hi], x[lo:hi], color="0.32", lw=0.6, zorder=3)
            ax.set_ylabel("dF/F", fontsize=7)
            ax.tick_params(labelsize=7)

            # spike rate on a twin axis, drawn behind the trace
            ax2 = ax.twinx()
            rr = np.nan_to_num(r[lo:hi])
            ax2.fill_between(t[lo:hi], 0, rr, color="#f0921e", alpha=0.45,
                             lw=0, zorder=1)
            ax2.set_ylim(0, max(np.nanmax(rr) * 1.05, 1e-3))
            ax2.set_ylabel("spikes/s", fontsize=7, color="#b06a10")
            ax2.tick_params(labelsize=6, colors="#b06a10")

            s_in = spk[(spk >= lo) & (spk < hi)]
            if len(s_in):
                ax.plot(t[s_in], x[s_in], "v", color="#1b7f3b", ms=5,
                        mec="none", zorder=5)

            rate = len(spk) / (T / 30.0) * 60.0
            short = sess.split(os.sep)[-1][:24]
            grp = sess.split(os.sep)[0]
            ax.set_title(f"{grp} | {short} ROI{k}   nu={noise[k]:.2f}   "
                         f"{len(spk)} events ({rate:.1f}/min, {len(aps)} APs)",
                         fontsize=7.5)
            ax.grid(alpha=0.25)
            if j < len(sel) - N_COLS:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("time (s)", fontsize=7)

        for j in range(len(sel), len(axes)):
            axes[j].axis("off")

        fig.suptitle("CASCADE GC8s_EXC_30Hz — automatic GCaMP8s EVENT detection "
                     f"(no hand labels).  green = event onset, orange = spike rate"
                     f"   |   page {p+1}/{n_page}",
                     fontsize=11, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.975])
        fig.savefig(os.path.join(OUT_DIR, f"cascade_gallery_p{p+1:02d}.png"),
                    dpi=140)
        pdf.savefig(fig)
        plt.close(fig)
        print(f"page {p+1}/{n_page}")

    pdf.close()
    print(f"\n{tot_sp} spikes across the {len(pick)} sampled ROIs -> {OUT_DIR}")


if __name__ == "__main__":
    main()
