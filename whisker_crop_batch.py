#!/usr/bin/env python3
"""
whisker_crop_batch.py  --  crop the two whisker ROIs out of EVERY run (active +
archived) under a root, for the whisker_crop_track.m pipeline.

The camera frame differs between scopes but is constant within a scope, so ROIs
are drawn ONCE per scope and applied to all of that scope's runs.

Two phases (run --draw first, then the batch):
  1) DRAW (interactive, fast):   pick the 2 ROIs for each scope, cached to
     <root>\<scope>\whisker_rois.json  (a scope already cropped is seeded from
     its existing *_whiskcrop.mat, no redraw needed).
       C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe whisker_crop_batch.py <root> --draw
  2) CROP (heavy, ~8-10 min/run, background-able): crop every run using the
     cached ROIs; skips runs that already have a *_whiskcrop.mat.
       C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe whisker_crop_batch.py <root>

Options: --names whisk_L,whisk_R  --preview blend|max|mean  --projframes 3000
         --max_bytes 1.5e9  --redraw (force redraw in --draw)
"""
import os, sys, glob, json, argparse
import numpy as np
import cv2
from scipy.io import savemat, loadmat


def load_timestamps(ts):
    arr = np.genfromtxt(ts, delimiter=",", skip_header=1)
    cam_ns = arr[:, 1]
    return (cam_ns - cam_ns[0]) / 1e9


def frame_size(avi):
    cap = cv2.VideoCapture(avi)
    wh = (int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)), int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT)))
    cap.release()
    return wh


def projection(avi, nframes, preview):
    cap = cv2.VideoCapture(avi)
    acc = None; mx = None; c = 0
    while c < nframes:
        ok, fr = cap.read()
        if not ok:
            break
        g = fr[:, :, 0]
        acc = g.astype(np.float64) if acc is None else acc + g
        mx = g.copy() if mx is None else np.maximum(mx, g)
        c += 1
    cap.release()
    mean = acc / c
    mean8 = np.clip((mean - np.percentile(mean, 1)) / max(np.percentile(mean, 99) - np.percentile(mean, 1), 1e-6), 0, 1) * 255
    mean8 = mean8.astype(np.uint8)
    if preview == "mean":
        return mean8
    if preview == "max":
        return mx
    return np.clip(0.6 * mean8.astype(np.float32) + 0.7 * mx.astype(np.float32), 0, 255).astype(np.uint8)


def avg_one(avi, nframes):
    """average projection (first nframes) of one run, contrast-stretched uint8 + raw mean."""
    cap = cv2.VideoCapture(avi); acc = None; c = 0
    while c < nframes:
        ok, fr = cap.read()
        if not ok:
            break
        g = fr[:, :, 0].astype(np.float64)
        acc = g if acc is None else acc + g
        c += 1
    cap.release()
    return None if acc is None else acc / c


def composite_avg(runs, nframes):
    """mean of each run's average projection -> one composite across ALL runs."""
    acc = None; n = 0
    for folder in runs:
        avi = glob.glob(os.path.join(folder, "*.avi"))[0]
        m = avg_one(avi, nframes)
        if m is None:
            continue
        if acc is not None and m.shape != acc.shape:
            print(f"    SKIP {os.path.basename(folder)} (size {m.shape} != {acc.shape})", flush=True)
            continue
        acc = m if acc is None else acc + m
        n += 1
        print(f"    + {os.path.basename(folder)} ({n})", flush=True)
    if acc is None:
        return None
    m = acc / n
    lo, hi = np.percentile(m, 1), np.percentile(m, 99)
    return (np.clip((m - lo) / max(hi - lo, 1e-6), 0, 1) * 255).astype(np.uint8)


def discover(root):
    """run folders (have a *.avi + timestamps.csv) grouped by scope (1st dir under root)."""
    runs = {}
    for avi in glob.glob(os.path.join(root, "**", "*.avi"), recursive=True):
        folder = os.path.dirname(avi)
        if not os.path.isfile(os.path.join(folder, "timestamps.csv")):
            continue
        rel = os.path.relpath(folder, root)
        scope = rel.split(os.sep)[0]
        runs.setdefault(scope, []).append(folder)
    return {s: sorted(v) for s, v in sorted(runs.items())}


def seed_from_existing(runs, names):
    """if some run in this scope already has a *_whiskcrop.mat, reuse its ROIs."""
    for folder in runs:
        m = glob.glob(os.path.join(folder, "*_whiskcrop.mat"))
        if not m:
            continue
        try:
            d = loadmat(m[0])
            rois = {nm: [int(v) for v in d[f"roi_{nm}"].ravel()] for nm in names}
            return rois, os.path.basename(m[0])
        except Exception:
            continue
    return None, None


def draw_rois_img(base, names):
    """run selectROI for each name on a ready grayscale base image."""
    disp = cv2.cvtColor(base, cv2.COLOR_GRAY2BGR)
    rois = {}
    for nm in names:
        x, y, w, h = cv2.selectROI(f"Draw {nm}  (drag, ENTER=ok)", disp, showCrosshair=True, fromCenter=False)
        cv2.destroyAllWindows()
        if w == 0 or h == 0:
            sys.exit(f"no ROI for {nm}")
        rois[nm] = [int(x), int(y), int(w), int(h)]
        cv2.rectangle(disp, (x, y), (x + w, y + h), (0, 255, 0), 2)
        cv2.putText(disp, nm, (x, max(12, y - 4)), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 0), 2)
    return rois, disp


def draw_rois(avi, names, preview, projframes):
    return draw_rois_img(projection(avi, projframes, preview), names)


def crop_run(folder, rois, names, max_bytes):
    avi = glob.glob(os.path.join(folder, "*.avi"))[0]
    ts = os.path.join(folder, "timestamps.csv")
    run = os.path.basename(folder)
    t_all = load_timestamps(ts); N = len(t_all)
    fps = (N - 1) / (t_all[-1] - t_all[0])

    binf = {}
    for nm in names:
        x, y, w, h = rois[nm]; b = 1
        while N * (h // b) * (w // b) > max_bytes:
            b += 1
        binf[nm] = b
    cubes = {nm: np.empty((rois[nm][3] // binf[nm], rois[nm][2] // binf[nm], N), np.uint8) for nm in names}

    cap = cv2.VideoCapture(avi); col = 0
    while col < N:
        ok, fr = cap.read()
        if not ok:
            break
        g = fr[:, :, 0]
        for nm in names:
            x, y, w, h = rois[nm]; b = binf[nm]
            c = g[y:y + h, x:x + w]
            if b > 1:
                c = cv2.resize(c, (w // b, h // b), interpolation=cv2.INTER_AREA)
            cubes[nm][:, :, col] = c
        col += 1
        if col % 10000 == 0:
            print(f"    {run} {col}/{N}", flush=True)
    cap.release()
    for nm in names:
        cubes[nm] = cubes[nm][:, :, :col]
    t_keep = (t_all[:col] - t_all[0]).astype(np.float64)

    md = {"names": np.array(names, dtype=object), "t_s": t_keep, "fps": float(fps),
          "src_avi": os.path.abspath(avi)}
    for nm in names:
        md[f"mov_{nm}"] = cubes[nm]; md[f"roi_{nm}"] = np.array(rois[nm], np.float64)
        md[f"bin_{nm}"] = float(binf[nm])
    savemat(os.path.join(folder, run + "_whiskcrop.mat"), md, do_compression=True)
    return col


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--composite-build", dest="composite_build", action="store_true",
                    help="phase 0: build per-scope average composite (across ALL runs) to draw on")
    ap.add_argument("--draw", action="store_true", help="phase 1: pick ROIs per scope")
    ap.add_argument("--redraw", action="store_true", help="force redraw even if cached")
    ap.add_argument("--names", default="whisk_L,whisk_R")
    ap.add_argument("--preview", choices=["max", "mean", "blend"], default="blend")
    ap.add_argument("--projframes", type=int, default=3000)
    ap.add_argument("--max_bytes", type=float, default=1.5e9)
    args = ap.parse_args()
    names = [s.strip() for s in args.names.split(",")]

    scopes = discover(args.root)
    print(f"{len(scopes)} scopes, {sum(len(v) for v in scopes.values())} runs total", flush=True)

    if args.composite_build:
        for scope, runs in scopes.items():
            out = os.path.join(args.root, scope, "whisker_composite.png")
            if os.path.isfile(out) and not args.redraw:
                print(f"[{scope}] composite exists (skip; --redraw to rebuild)", flush=True); continue
            print(f"[{scope}] building avg composite across {len(runs)} runs ...", flush=True)
            comp = composite_avg(runs, args.projframes)
            if comp is None:
                print(f"[{scope}] no frames"); continue
            cv2.imwrite(out, comp)
            print(f"[{scope}] saved {out}  ({comp.shape[1]}x{comp.shape[0]})", flush=True)
        print("\nComposites built. Now run with --draw to pick one ROI per scope on each composite.")
        return

    if args.draw:
        for scope, runs in scopes.items():
            cache = os.path.join(args.root, scope, "whisker_rois.json")
            comp_png = os.path.join(args.root, scope, "whisker_composite.png")
            if os.path.isfile(cache) and not args.redraw:
                print(f"[{scope}] cached ROIs -> {cache} (skip; --redraw to redo)")
                continue
            ref = list(frame_size(glob.glob(os.path.join(runs[0], '*.avi'))[0]))
            if os.path.isfile(comp_png):
                print(f"[{scope}] draw 2 ROIs on the all-runs avg composite ...")
                base = cv2.imread(comp_png, cv2.IMREAD_GRAYSCALE)
                rois, disp = draw_rois_img(base, names)
                cv2.imwrite(os.path.join(args.root, scope, "whisker_rois.png"), disp)
            else:
                seeded, src = seed_from_existing(runs, names)
                if seeded and not args.redraw:
                    rois = seeded
                    print(f"[{scope}] seeded ROIs from {src} (no composite; --composite-build for one)")
                else:
                    avi = glob.glob(os.path.join(runs[0], "*.avi"))[0]
                    print(f"[{scope}] draw 2 ROIs on {os.path.basename(runs[0])} (no composite) ...")
                    rois, disp = draw_rois(avi, names, args.preview, args.projframes)
                    cv2.imwrite(os.path.join(args.root, scope, "whisker_rois.png"), disp)
            with open(cache, "w") as f:
                json.dump({"ref_wh": ref, "rois": rois, "names": names}, f, indent=2)
            print(f"[{scope}] saved {cache}: {rois}")
        print("\nDraw phase done. Now run WITHOUT --draw to crop all runs.")
        return

    # crop phase
    done = skipped = failed = 0
    for scope, runs in scopes.items():
        cache = os.path.join(args.root, scope, "whisker_rois.json")
        if not os.path.isfile(cache):
            print(f"[{scope}] NO ROI cache -- run with --draw first. Skipping scope.")
            continue
        with open(cache) as f:
            cc = json.load(f)
        rois = {nm: cc["rois"][nm] for nm in names}; ref = tuple(cc["ref_wh"])
        for folder in runs:
            run = os.path.basename(folder)
            if glob.glob(os.path.join(folder, "*_whiskcrop.mat")):
                print(f"  [{scope}/{run}] SKIP (done)"); skipped += 1; continue
            wh = frame_size(glob.glob(os.path.join(folder, "*.avi"))[0])
            if wh != ref:
                print(f"  [{scope}/{run}] SKIP frame {wh} != scope ref {ref} (draw separately)"); skipped += 1; continue
            print(f"  [{scope}/{run}] cropping ...", flush=True)
            try:
                n = crop_run(folder, rois, names, args.max_bytes)
                print(f"  [{scope}/{run}] OK ({n} frames)"); done += 1
            except Exception as e:
                print(f"  [{scope}/{run}] FAILED: {e}"); failed += 1
    print(f"\n=== BATCH: cropped {done}, skipped {skipped}, failed {failed} ===")


if __name__ == "__main__":
    main()
