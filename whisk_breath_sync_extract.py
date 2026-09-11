#!/usr/bin/env python
"""
whisk_breath_sync_extract.py  --  pull the timing + a breathing proxy from BOTH
cameras of one run folder, so the frame-level thermal<->Basler alignment can be
measured (they share the 400 Hz WFG trigger but their absolute clocks are NOT
comparable, so alignment must come from a shared physical signal = breathing).

For one run folder it reads:
  - thermal .ats (fnv)  : per-frame IRIG time (every frame) + a breathing proxy
                          (whole-frame mean temp, and mean of the highest-variance
                          pixels = nostril airflow region)
  - Basler  .avi (cv2)  : nostril-ROI mean intensity (breathing proxy) per frame
  - timestamps.csv      : Basler per-frame camera-clock + wall time

Saves <run>_sync.mat for whisk_breath_sync_report.m to align + plot.

  C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe whisk_breath_sync_extract.py \
      <run_folder> --nostril x,y,w,h  [--bdecim 1] [--tvar-frac 0.05]

--nostril is the Basler nostril ROI (pixels in the FULL frame); read it off
<stem>_proj.png from basler_preview.py. If omitted, whole-frame mean is used.
"""
import argparse, os, sys, glob, csv
import numpy as np
import cv2
import fnv, fnv.file
from scipy.io import savemat


def find_one(folder, pattern):
    hits = glob.glob(os.path.join(folder, pattern))
    if not hits:
        sys.exit(f"no {pattern} in {folder}")
    return hits[0]


def read_thermal(ats_path, tvar_frac):
    im = fnv.file.ImagerFile(ats_path)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    H, W, N = im.height, im.width, im.num_frames
    P = H * W
    X = np.empty((N, P), dtype=np.float32)
    t = np.empty(N, dtype=np.float64)
    t0 = None
    for i in range(N):
        im.get_frame(i)
        X[i] = np.asarray(im.final, dtype=np.float32)
        ti = im.frame_info.time
        if t0 is None:
            t0 = ti
        t[i] = (ti - t0).total_seconds()
        if (i % 5000) == 0:
            print(f"  thermal {i}/{N}", flush=True)
    mean_trace = X.mean(axis=1)                       # whole-frame mean temp
    v = X.var(axis=0)                                 # per-pixel temporal variance
    k = max(1, int(round(tvar_frac * P)))
    sel = np.argsort(v)[-k:]                          # nostril/airflow pixels
    nostril_trace = X[:, sel].mean(axis=1)
    var_img = v.reshape(H, W)
    return t, mean_trace, nostril_trace, var_img, H, W, N


def read_basler(avi_path, roi, bdecim):
    cap = cv2.VideoCapture(avi_path)
    if not cap.isOpened():
        sys.exit("OpenCV could not open the Basler AVI")
    x, y, w, h = roi if roi is not None else (None, None, None, None)
    trace = []
    n = 0
    while True:
        ok = cap.grab()
        if not ok:
            break
        if (n % bdecim) == 0:
            ok, fr = cap.retrieve()
            if not ok:
                break
            g = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY)
            if roi is None:
                trace.append(float(g.mean()))
            else:
                trace.append(float(g[y:y + h, x:x + w].mean()))
        n += 1
        if (n % 5000) == 0:
            print(f"  basler {n}", flush=True)
    cap.release()
    return np.asarray(trace, dtype=np.float64), n


def read_csv_times(csv_path):
    fidx, cam_ns, wall_s = [], [], []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            fidx.append(int(row["frame_idx"]))
            cam_ns.append(int(row["camera_timestamp_ns"]))
            wall_s.append(float(row["wall_time_s"]))
    return (np.asarray(cam_ns, dtype=np.int64), np.asarray(wall_s, dtype=np.float64))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", help="run folder containing the .ats, .avi and timestamps.csv")
    ap.add_argument("--nostril", default=None, help="Basler nostril ROI 'x,y,w,h' (full-frame px)")
    ap.add_argument("--bdecim", type=int, default=1, help="Basler frame decimation (default 1 = all)")
    ap.add_argument("--tvar-frac", type=float, default=0.05,
                    help="fraction of thermal pixels (highest variance) used as nostril proxy")
    args = ap.parse_args()

    folder = args.folder
    ats = find_one(folder, "*.ats")
    avi = find_one(folder, "*.avi")
    csvp = find_one(folder, "timestamps.csv")
    roi = tuple(int(v) for v in args.nostril.split(",")) if args.nostril else None

    print("reading thermal ...", flush=True)
    t_t, t_mean, t_nostril, var_img, Ht, Wt, Nt = read_thermal(ats, args.tvar_frac)
    print("reading basler nostril ROI ...", flush=True)
    b_breath, Nb_read = read_basler(avi, roi, args.bdecim)
    cam_ns, wall_s = read_csv_times(csvp)

    run_name = os.path.basename(os.path.abspath(os.path.normpath(folder)))
    out = os.path.join(folder, run_name + "_sync.mat")
    savemat(out, {
        "thermal_t": t_t,                         # s, IRIG relative
        "thermal_mean": t_mean,                   # whole-frame mean temp
        "thermal_nostril": t_nostril,             # high-variance pixel mean (breathing)
        "thermal_var_img": var_img,               # [Ht x Wt] per-pixel temporal variance
        "basler_breath": b_breath,                # nostril-ROI mean intensity
        "basler_bdecim": int(args.bdecim),
        "basler_cam_ns": cam_ns,                  # camera hardware clock (ns)
        "basler_wall_s": wall_s,                  # PC epoch (s)
        "basler_roi": np.asarray(roi if roi else [0, 0, 0, 0]),
        "Nt": int(Nt), "Ht": int(Ht), "Wt": int(Wt),
        "ats": os.path.abspath(ats), "avi": os.path.abspath(avi),
    }, do_compression=True)
    print(f"saved {out}")
    print(f"  thermal frames {Nt}, basler frames read {Nb_read}, csv frames {len(cam_ns)}")


if __name__ == "__main__":
    main()
