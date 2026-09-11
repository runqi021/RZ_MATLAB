"""ca_detect_gallery100.py -- 100 random detected ROIs pooled across every session.

RUN WITH dlc310:
    %USERPROFILE%\\.conda\\envs\\dlc310\\python.exe ca_detect_gallery100.py

Draws a uniform random sample of USABLE ROIs from both datasets combined, so the
picture is representative of the whole corpus rather than of one session. Each
panel carries the numbers needed to judge it without going back to the data:

    grey        dF/F
    green v     detected event onset
    title       group | session | ROI, reliability, tau, event count, rate,
                and the median spike-triggered rise in units of that ROI's own
                noise -- the single number that says whether the detections sit
                on real transients or on noise

A summary line at the end reports that rise distribution over the whole sample,
which is the honest headline: it is a label-free measure of detection quality.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

ROOTS = [
    (r"D:\260721_Sert_soma_G8s\phys", "Sert"),
    (r"D:\Ventral_surface_summary", "Ventral"),
]
OUT = r"D:\Ventral_surface_summary\_detect_260723\gallery100"
N_PICK, PER_PAGE, N_COLS, WIN_S, FPS, SEED = 100, 10, 2, 40, 30.0, 7
MIN_EV = 1


def robust_sn(x):
    return max(np.median(np.abs(np.diff(x))) * 1.4826 / np.sqrt(2), 1e-9)


def rise_sigma(x, ev, sn):
    """Median event amplitude in noise units.

    The baseline MUST be taken well before the mark. Events here are timed at
    the transient PEAK (the detector shifts +4 frames from the OASIS onset), so
    a window like x[t-3:t+1] sits ON the rising phase and reports roughly a
    third of the true amplitude -- measured 1.1 sigma where the correct value is
    5.1 sigma. Take the baseline at t-14..t-8 and search the peak around the
    mark, which is valid whether the mark lands on the onset or the peak."""
    out = []
    for t in ev:
        if t < 15 or t > len(x) - 16:
            continue
        base = np.median(x[t - 14:t - 8])
        out.append((np.max(x[t - 6:t + 12]) - base) / sn)
    return float(np.median(out)) if out else np.nan


def main():
    os.makedirs(OUT, exist_ok=True)
    pool = []
    for root, tag in ROOTS:
        tdir = os.path.join(root, "_cascade_260722", "traces")
        for dp, _dn, fn in os.walk(root):
            if "ca_events_auto.npz" not in fn:
                continue
            rel_s = os.path.relpath(dp, root)
            d = np.load(os.path.join(dp, "ca_events_auto.npz"), allow_pickle=True)
            usable, ev = d["usable"], d["events"]
            ttag = rel_s.replace(os.sep, "__")
            tp = os.path.join(tdir, ttag + ".npz")
            if not os.path.isfile(tp):
                continue
            for k in range(len(usable)):
                if usable[k] and np.asarray(ev[k]).size >= MIN_EV:
                    pool.append((tag, root, rel_s, ttag, k))
    print(f"{len(pool)} usable ROIs with >= {MIN_EV} event(s) across both datasets")

    rng = np.random.default_rng(SEED)
    pick = [pool[i] for i in rng.permutation(len(pool))[:min(N_PICK, len(pool))]]

    tcache, dcache, rises = {}, {}, []
    npage = int(np.ceil(len(pick) / PER_PAGE))
    pdf = PdfPages(os.path.join(OUT, "gallery100.pdf"))

    for p in range(npage):
        sel = pick[p * PER_PAGE:(p + 1) * PER_PAGE]
        nrow = int(np.ceil(len(sel) / N_COLS))
        fig, axes = plt.subplots(nrow, N_COLS, figsize=(19, 11.5))
        axes = np.atleast_1d(axes).ravel()
        for j, (tag, root, rel_s, ttag, k) in enumerate(sel):
            if ttag not in tcache:
                tcache[ttag] = np.load(
                    os.path.join(root, "_cascade_260722", "traces",
                                 ttag + ".npz"))["dFF"]
            if (root, rel_s) not in dcache:
                dcache[(root, rel_s)] = np.load(
                    os.path.join(root, rel_s, "ca_events_auto.npz"),
                    allow_pickle=True)
            x = tcache[ttag][k]
            d = dcache[(root, rel_s)]
            ev = np.asarray(d["events"][k], int)
            sn = robust_sn(x)
            rs = rise_sigma(x, ev, sn)
            rises.append(rs)

            T = x.size
            t = np.arange(T) / FPS
            if T / FPS > WIN_S:
                w = int(WIN_S * FPS)
                c = int(ev[len(ev) // 2]) if ev.size else T // 2
                lo = max(0, min(T - w, c - w // 2))
                hi = lo + w
            else:
                lo, hi = 0, T

            ax = axes[j]
            ax.plot(t[lo:hi], x[lo:hi], color="0.30", lw=0.7, zorder=2)
            e = ev[(ev >= lo) & (ev < hi)]
            if e.size:
                ax.plot(t[e], x[e], "v", ms=5.5, color="#1b7f3b",
                        mec="none", zorder=4)
            rate = ev.size / (T / FPS / 60.0)
            short = rel_s.split(os.sep)[-1][:26]
            grp = rel_s.split(os.sep)[0] if os.sep in rel_s else tag
            ax.set_title(f"{tag}/{grp} | {short} ROI{k}   rel={float(d['rel'][k]):.2f}  "
                         f"tau={float(d['tau'][k]):.2f}s   {ev.size} ev "
                         f"({rate:.1f}/min)   rise={rs:.1f}\u03c3",
                         fontsize=7.5)
            ax.grid(alpha=.22)
            ax.tick_params(labelsize=7)
            ax.set_ylabel("dF/F", fontsize=7)
            if j >= len(sel) - N_COLS:
                ax.set_xlabel("time (s)", fontsize=7)
        for j in range(len(sel), len(axes)):
            axes[j].axis("off")
        fig.suptitle("Automatic calcium event detection, label-free — 100 random "
                     f"usable ROIs across 49 sessions   |   page {p+1}/{npage}",
                     fontsize=11, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.972])
        fig.savefig(os.path.join(OUT, f"gallery100_p{p+1:02d}.png"), dpi=140)
        pdf.savefig(fig)
        plt.close(fig)
        print(f"page {p+1}/{npage}")
    pdf.close()

    r = np.array([v for v in rises if np.isfinite(v)])
    print(f"\nspike-triggered rise over the {r.size} sampled ROIs (noise units):")
    print(f"  pct 10/25/50/75/90 = {np.round(np.percentile(r, [10,25,50,75,90]), 1)}")
    print(f"  ROIs whose median rise is < 2 sigma: {100*np.mean(r < 2):.0f}%")
    print(f"  ROIs whose median rise is >= 4 sigma: {100*np.mean(r >= 4):.0f}%")
    print(f"\n-> {OUT}")


if __name__ == "__main__":
    main()
