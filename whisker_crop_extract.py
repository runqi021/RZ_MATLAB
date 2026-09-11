#!/usr/bin/env python3
"""
whisker_crop_extract.py
=======================
Interactively crop TWO whisker ROIs (e.g. left + right cheek) out of a Basler
FFV1 .avi and save them as compact cubes for MATLAB whisker tracking.

You draw each box on a MAX-projection preview (the whisker marker sweeps a fan, so
the max image shows the full range to enclose). Each ROI is then cropped across
every frame, size-budget binned to fit the MAT 2 GB/var limit, and saved; a small
CLAHE-enhanced cropped mp4 is also written per ROI so you can watch it.

Run it yourself (needs the interactive window). In the Claude prompt:
    ! C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe C:\\Users\\Admin\\Desktop\\RZ_MATLAB\\whisker_crop_extract.py "<run_folder>"

  optional: --names whisk_L,whisk_R   --projframes 3000   --max_bytes 1.5e9
            --rois x,y,w,h;x,y,w,h     (skip interactive, pass both boxes)

Output (in the run folder):
  <run>_whiskcrop.mat   mov_<name> uint8 [H,W,T], roi_<name>, bin_<name>, t_s, fps
  <run>_<name>_crop.mp4  cropped CLAHE preview per ROI
  <run>_whisk_rois.png   max projection with the chosen boxes
"""
import os, sys, glob, argparse
import numpy as np
import cv2
from scipy.io import savemat


def find_one(folder, pat):
    h = glob.glob(os.path.join(folder, pat))
    if not h:
        sys.exit(f"no {pat} in {folder}")
    return h[0]


def load_timestamps(ts_path):
    arr = np.genfromtxt(ts_path, delimiter=",", skip_header=1)
    cam_ns = arr[:, 1]
    return (cam_ns - cam_ns[0]) / 1e9


def projection(avi, nframes):
    """mean + max over the first nframes mono frames (max shows the whisker fan)."""
    cap = cv2.VideoCapture(avi)
    if not cap.isOpened():
        sys.exit(f"cannot open {avi}")
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
    if acc is None:
        sys.exit("no frames read for preview")
    mean8 = np.clip((acc / c - np.percentile(acc / c, 1)) /
                    max(np.ptp(np.percentile(acc / c, [1, 99])), 1e-6), 0, 1) * 255
    return mean8.astype(np.uint8), mx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--names", default="whisk_L,whisk_R")
    ap.add_argument("--projframes", type=int, default=3000)
    ap.add_argument("--preview", choices=["max", "mean", "blend"], default="blend",
                    help="ROI-drawing image: max=whisker sweep envelope, mean=static face, "
                         "blend=both (default)")
    ap.add_argument("--max_bytes", type=float, default=1.5e9)
    ap.add_argument("--rois", default=None, help="'x,y,w,h;x,y,w,h' to skip interactive")
    args = ap.parse_args()

    folder = args.folder
    avi = find_one(folder, "*.avi")
    ts = find_one(folder, "timestamps.csv")
    run = os.path.basename(os.path.abspath(os.path.normpath(folder)))
    names = [s.strip() for s in args.names.split(",")]

    t_all = load_timestamps(ts)
    fps = (len(t_all) - 1) / (t_all[-1] - t_all[0])
    print(f"{run}: {len(t_all)} frames, {t_all[-1]:.2f}s, {fps:.2f} fps", flush=True)

    print(f"building {args.preview} projection for ROI drawing ...", flush=True)
    mean8, mx = projection(avi, args.projframes)
    if args.preview == "mean":
        base = mean8
    elif args.preview == "max":
        base = mx                                  # whisker sweep envelope
    else:  # blend: static face (mean) brightened by the sweep envelope (max)
        base = np.clip(0.6 * mean8.astype(np.float32) + 0.7 * mx.astype(np.float32), 0, 255).astype(np.uint8)
    disp = cv2.cvtColor(base, cv2.COLOR_GRAY2BGR)

    # --- pick ROIs ---
    rois = {}
    if args.rois:
        for nm, box in zip(names, args.rois.split(";")):
            rois[nm] = tuple(int(v) for v in box.split(","))
    else:
        for nm in names:
            x, y, w, h = cv2.selectROI(f"Draw {nm}  (drag, ENTER=ok)", disp,
                                       showCrosshair=True, fromCenter=False)
            cv2.destroyAllWindows()
            if w == 0 or h == 0:
                sys.exit(f"no ROI for {nm}")
            rois[nm] = (int(x), int(y), int(w), int(h))
            cv2.rectangle(disp, (x, y), (x + w, y + h), (0, 255, 0), 2)
            cv2.putText(disp, nm, (x, max(12, y - 4)), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 255, 0), 2)
    for nm in names:
        print(f"  {nm} ROI = {rois[nm]}", flush=True)
    cv2.imwrite(os.path.join(folder, run + "_whisk_rois.png"), disp)

    # --- size-budget bin per ROI ---
    N = len(t_all)
    binf = {}
    for nm in names:
        x, y, w, h = rois[nm]
        b = 1
        while N * (h // b) * (w // b) > args.max_bytes:
            b += 1
        binf[nm] = b
        print(f"  {nm}: bin {b} -> {w//b}x{h//b}, ~{N*(h//b)*(w//b)/1e6:.0f} MB", flush=True)

    cubes = {nm: np.empty((rois[nm][3] // binf[nm], rois[nm][2] // binf[nm], N), np.uint8) for nm in names}

    # --- single pass: crop every frame into both cubes ---
    cap = cv2.VideoCapture(avi)
    col = 0
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
        if col % 5000 == 0:
            print(f"  {col}/{N}", flush=True)
    cap.release()
    for nm in names:
        cubes[nm] = cubes[nm][:, :, :col]
    t_keep = (t_all[:col] - t_all[0]).astype(np.float64)
    print(f"cropped {col} frames", flush=True)

    # --- cropped CLAHE preview mp4 per ROI (decim 4, 60 fps) ---
    clahe = cv2.createCLAHE(2.5, (8, 8))
    for nm in names:
        H, W, T = cubes[nm].shape
        up = max(1, int(round(240 / max(H, W))))     # upscale small crops for visibility
        vw = cv2.VideoWriter(os.path.join(folder, f"{run}_{nm}_crop.mp4"),
                             cv2.VideoWriter_fourcc(*"mp4v"), 60.0, (W * up, H * up), isColor=False)
        for k in range(0, T, 4):
            im = clahe.apply(cubes[nm][:, :, k])
            vw.write(cv2.resize(im, (W * up, H * up), interpolation=cv2.INTER_NEAREST))
        vw.release()

    out = os.path.join(folder, run + "_whiskcrop.mat")
    md = {"names": np.array(names, dtype=object), "t_s": t_keep, "fps": float(fps),
          "src_avi": os.path.abspath(avi)}
    for nm in names:
        md[f"mov_{nm}"] = cubes[nm]
        md[f"roi_{nm}"] = np.array(rois[nm], np.float64)
        md[f"bin_{nm}"] = float(binf[nm])
    savemat(out, md, do_compression=True)
    print(f"saved {out}")
    print(f"saved per-ROI crop mp4s + {run}_whisk_rois.png")


if __name__ == "__main__":
    main()
