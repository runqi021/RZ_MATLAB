"""ca_detect_ship.py -- the deliverable: automatic calcium event detection with
an honest per-ROI coverage gate, across every session.

RUN WITH dlc310 (needs working LAPACK for the tau fit; shells out to the `oasis`
env for the deconvolution itself, whose Cython core does not need BLAS):

    %USERPROFILE%\\.conda\\envs\\dlc310\\python.exe ca_detect_ship.py [ROOT ...]

DETECTOR. The configuration that reached per-ROI median F1 0.96 against the hand
labels, with every parameter now ESTIMATED rather than tuned:
  tau   per ROI from the autocovariance at lags >= 1 (lag 0 is the only lag that
        white noise touches; including it returns tau ~ 1 frame regardless of
        truth). Reproduces the supervised kernel fit without using any label.
  g     exp(-1/(fps*tau))
  sn    robust noise from frame-to-frame differences of the RAW dF/F -- valid
        there because that noise is white. Never apply it to a filtered trace.
  th    1.5 * sn on the OASIS spike train
  lag   +4 frames: OASIS marks the transient ONSET, events are timed at the peak
  event one per calcium transient, segmented by the trace returning toward
        baseline -- not by a fixed merge window, which splits long transients
        and merges genuinely separate fast ones

GATE. split_half_validate.py gives corr(A,B) between two interleaved pixel
halves: same signal, independent shot noise. ROIs below REL_MIN have no
recoverable signal and are marked unusable rather than emitting events nobody
should trust. Reported as coverage, not hidden.

Writes per session:  ca_events_auto.mat (v7, for MATLAB) and .npz
Writes per root:     _detect_260723/summary.csv and a QC gallery.
"""
import json
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402
from scipy.io import savemat                        # noqa: E402

ROOTS = sys.argv[1:] or [
    r"D:\260721_Sert_soma_G8s\phys",
    r"D:\Ventral_surface_summary",
]
REPO = os.path.dirname(os.path.abspath(__file__))
OASIS_PY = os.path.join(os.environ["USERPROFILE"], ".conda", "envs", "oasis",
                        "python.exe")
FPS = 30.0
TH_SN = 1.5
LAG_FR = 4
MIN_DIST_S = 0.15
RETURN_FRAC = 0.5
REL_MIN = 0.10
# GCaMP8s somatic decay is physically bounded. An ROI whose autocorrelation
# returns tau outside this range has no calcium kinetics in it -- it is noise or
# a shared non-neural fluctuation. This gate is COMPLEMENTARY to reliability,
# not redundant: split-half correlation certifies against shot noise only, and
# noise-only ROIs routinely score rel > 0.9 because both halves see the same
# artifact (measured: ro5 ROI29, rel 0.91, tau 0.13 s, 96 spurious events).
TAU_RANGE = (0.30, 1.50)
# tau also feeds g, so a bad estimate corrupts the deconvolution itself:
# tau 0.13 -> g 0.77 makes OASIS fire on every noise bump. Clamp what is used.
TAU_CLAMP = (0.35, 1.20)
GAL_N, GAL_PER_PAGE, GAL_COLS, GAL_WIN = 100, 20, 2, 30
SEED = 260723


def robust_sn(x):
    return np.median(np.abs(np.diff(x))) * 1.4826 / np.sqrt(2)


def ac_tau(x, fps):
    x = np.asarray(x, float)
    x = x - x.mean()
    L = int(min(max(3, round(0.6 * fps)), len(x) // 4))
    c = np.correlate(x, x, "full")[len(x) - 1:len(x) - 1 + L + 1] / len(x)
    c1 = c[1:]
    ok = c1 > 0
    if ok.sum() < 3:
        return 0.5
    lags = np.arange(1, len(c1) + 1)[ok]
    g = np.exp(np.polyfit(lags, np.log(c1[ok]), 1)[0])
    if not np.isfinite(g) or g <= 0 or g >= 1:
        return 0.5
    return float(np.clip(-1.0 / (fps * np.log(g)), 0.1, 3.0))


def peaks(z, thr, md):
    c = np.flatnonzero((z[1:-1] >= z[:-2]) & (z[1:-1] > z[2:]) &
                       (z[1:-1] >= thr)) + 1
    if c.size == 0:
        return c
    keep, taken = [], np.zeros(z.size, bool)
    for i in c[np.argsort(-z[c])]:
        if not taken[max(0, i - md):min(z.size, i + md + 1)].any():
            keep.append(i)
            taken[i] = True
    return np.sort(np.array(keep, int))


def to_events(spk, x, min_gap):
    """AP times -> onsets of distinct transients (trace-return segmentation)."""
    if spk.size == 0:
        return spk
    spk = np.sort(spk)
    ons = [spk[0]]
    for a, b in zip(spk[:-1], spk[1:]):
        if b - a < min_gap:
            continue
        seg = x[a:b + 1]
        if seg.size < 3:
            continue
        base, peak = x[ons[-1]], np.max(seg)
        amp = peak - base
        if amp <= 0 or np.min(seg[np.argmax(seg):]) < base + RETURN_FRAC * amp:
            ons.append(b)
    return np.array(ons, int)


def run_oasis(Y, g, tag):
    fin = os.path.join(os.environ.get("TEMP", "."), f"ship_in_{tag}.mat")
    fout = os.path.join(os.environ.get("TEMP", "."), f"ship_out_{tag}.mat")
    savemat(fin, {"Y": Y, "g": float(g)}, do_compression=True)
    r = subprocess.run([OASIS_PY, os.path.join(REPO, "ca_oasis_run.py"),
                        "--in", fin, "--out", fout, "--g", f"{g:.6f}",
                        "--lam", "0"], capture_output=True, text=True)
    if r.returncode != 0 or not os.path.isfile(fout):
        for f in (fin, fout):
            if os.path.isfile(f):
                os.remove(f)
        raise RuntimeError(f"oasis failed ({r.returncode}): {r.stderr[-200:]}")
    from scipy.io import loadmat
    S = loadmat(fout)["S"]
    for f in (fin, fout):
        if os.path.isfile(f):
            os.remove(f)
    return S


def load_reliability(folder, n):
    p = os.path.join(folder, "split_half.npz")
    rel = np.full(n, np.nan)
    if os.path.isfile(p):
        try:
            rows = np.load(p)["rows"]
            m = min(n, rows.shape[0])
            rel[:m] = rows[:m, 1]
        except Exception:
            pass
    return rel


def main():
    for root in ROOTS:
        tdir = os.path.join(root, "_cascade_260722", "traces")
        idxf = os.path.join(tdir, "index.json")
        if not os.path.isfile(idxf):
            print(f"{root}: no exported traces -- run cascade_export_dff.py")
            continue
        index = json.load(open(idxf))
        outdir = os.path.join(root, "_detect_260723")
        os.makedirs(outdir, exist_ok=True)
        print(f"\n=== {root}  ({len(index)} sessions) ===")

        summary, gal_pool = [], []
        for si, e in enumerate(index, 1):
            dff = np.load(os.path.join(tdir, e["tag"] + ".npz"))["dFF"]
            N, T = dff.shape
            folder = os.path.join(root, e["session"])

            sn = np.array([max(robust_sn(dff[k]), 1e-9) for k in range(N)])
            tau = np.array([ac_tau(dff[k], FPS) for k in range(N)])
            tau_used = np.clip(tau, *TAU_CLAMP)      # what the deconvolution sees
            g = np.exp(-(1.0 / FPS) / tau_used)
            rel = load_reliability(folder, N)

            S = np.zeros((T, N))
            gq = np.round(g * 100) / 100
            for gv in np.unique(gq):
                sel = np.flatnonzero(gq == gv)
                try:
                    S[:, sel] = run_oasis(dff[sel].T, gv, f"{si:03d}_{int(gv*100)}")
                except Exception as ex:
                    print(f"    g={gv:.2f}: {ex}")

            md = max(1, int(round(MIN_DIST_S * FPS)))
            events, nev = [], np.zeros(N, int)
            for k in range(N):
                z = S[:, k] / sn[k]
                ap = peaks(z, TH_SN, md) + LAG_FR
                ap = ap[(ap >= 0) & (ap < T)]
                ev = to_events(ap, dff[k], md)
                events.append(ev)
                nev[k] = ev.size

            # Two INDEPENDENT gates. Missing reliability means UNKNOWN, not
            # unusable -- defaulting it to fail silently excluded 807 of 898
            # Ventral ROIs on the first run purely because split_half had not
            # been computed there yet.
            tau_ok = (tau >= TAU_RANGE[0]) & (tau <= TAU_RANGE[1])
            rel_ok = np.isnan(rel) | (rel >= REL_MIN)
            usable = tau_ok & rel_ok
            rate = nev / (T / FPS / 60.0)
            savemat(os.path.join(folder, "ca_events_auto.mat"),
                    {"events": np.array([ev + 1 for ev in events], dtype=object),
                     "n_events": nev, "rate_per_min": rate, "tau": tau,
                     "sn": sn, "reliability": rel, "usable": usable.astype(np.uint8),
                     "tau_used": tau_used, "tau_ok": tau_ok.astype(np.uint8),
                     "th_sn": TH_SN, "lag_fr": LAG_FR, "fps": FPS}, do_compression=True)
            np.savez_compressed(os.path.join(folder, "ca_events_auto.npz"),
                                events=np.array(events, dtype=object), n=nev,
                                rate=rate, tau=tau, sn=sn, rel=rel, usable=usable)

            for k in range(N):
                if usable[k] and nev[k] > 0:
                    gal_pool.append((e["tag"], e["session"], k))
            summary.append((e["session"], N, int(usable.sum()),
                            float(np.nanmedian(rel)), float(np.median(tau)),
                            float(np.median(rate[usable])) if usable.any() else np.nan,
                            int(nev.sum())))
            print(f"  [{si:2d}/{len(index)}] {e['session'][:46]:46s} {N:4d} ROI  "
                  f"usable {usable.sum():4d} ({100*usable.mean():3.0f}%)  "
                  f"tau {np.median(tau):.2f}s (tau_ok {tau_ok.sum():3d})  {nev[usable].sum():5d} ev  "
                  f"{np.median(rate[usable]) if usable.any() else float('nan'):5.1f}/min",
                  flush=True)

        with open(os.path.join(outdir, "summary.csv"), "w") as f:
            f.write("session,n_roi,n_usable,median_reliability,median_tau_s,"
                    "median_rate_per_min,total_events\n")
            for r in summary:
                f.write(",".join(str(x) for x in r) + "\n")

        tot = sum(r[1] for r in summary)
        use = sum(r[2] for r in summary)
        ev = sum(r[6] for r in summary)
        print(f"  -> {tot} ROI, {use} usable ({100*use/max(tot,1):.0f}%), "
              f"{ev} events   summary.csv written")
        gallery(root, tdir, gal_pool, outdir)


def gallery(root, tdir, pool, outdir):
    if not pool:
        return
    rng = np.random.default_rng(SEED)
    pick = [pool[i] for i in rng.permutation(len(pool))[:min(GAL_N, len(pool))]]
    cache = {}
    npage = int(np.ceil(len(pick) / GAL_PER_PAGE))
    for p in range(npage):
        sel = pick[p * GAL_PER_PAGE:(p + 1) * GAL_PER_PAGE]
        nrow = int(np.ceil(len(sel) / GAL_COLS))
        fig, axes = plt.subplots(nrow, GAL_COLS, figsize=(19, 12))
        axes = np.atleast_1d(axes).ravel()
        for j, (tag, sess, k) in enumerate(sel):
            if tag not in cache:
                cache[tag] = np.load(os.path.join(tdir, tag + ".npz"))["dFF"]
            x = cache[tag][k]
            d = np.load(os.path.join(root, sess, "ca_events_auto.npz"),
                        allow_pickle=True)
            ev = np.asarray(d["events"][k], int)
            T = x.size
            t = np.arange(T) / FPS
            if T / FPS > GAL_WIN:
                w = int(GAL_WIN * FPS)
                c = int(ev[len(ev) // 2]) if ev.size else T // 2
                lo = max(0, min(T - w, c - w // 2))
                hi = lo + w
            else:
                lo, hi = 0, T
            ax = axes[j]
            ax.plot(t[lo:hi], x[lo:hi], color="0.32", lw=0.6)
            e_in = ev[(ev >= lo) & (ev < hi)]
            if e_in.size:
                ax.plot(t[e_in], x[e_in], "v", ms=5, color="#1b7f3b", mec="none")
            ax.set_title(f"{sess.split(os.sep)[0]} | {sess.split(os.sep)[-1][:22]} "
                         f"ROI{k}  rel={float(d['rel'][k]):.2f}  tau={float(d['tau'][k]):.2f}s  "
                         f"{ev.size} ev", fontsize=7.5)
            ax.grid(alpha=.25)
            ax.tick_params(labelsize=7)
            ax.set_ylabel("dF/F", fontsize=7)
            if j >= len(sel) - GAL_COLS:
                ax.set_xlabel("time (s)", fontsize=7)
        for j in range(len(sel), len(axes)):
            axes[j].axis("off")
        fig.suptitle(f"Automatic calcium event detection (label-free) — "
                     f"green = event onset   |   page {p+1}/{npage}",
                     fontsize=11, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.975])
        fig.savefig(os.path.join(outdir, f"gallery_p{p+1:02d}.png"), dpi=140)
        plt.close(fig)
    print(f"  -> gallery: {npage} page(s) in {outdir}")


if __name__ == "__main__":
    main()
