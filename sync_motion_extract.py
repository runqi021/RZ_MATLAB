#!/usr/bin/env python
"""
sync_motion_extract.py  --  waveform-agnostic thermal<->Basler frame alignment.

Breathing waveforms differ between thermal (airflow temperature) and Basler
(marker motion), so they don't cross-correlate cleanly. Abrupt body/head
MOVEMENTS, however, appear on the SAME trigger pulse in BOTH cameras. This reads
both videos once and stores per-frame MOTION ENERGY (mean |frame_t - frame_t-1|)
for each, plus a proper bottom-center thermal nostril-ROI temperature and the
Basler nostril-ROI mean, so sync_motion_report.m can lock the integer lag and
also keep breathing proxies.

  C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe sync_motion_extract.py <run_folder> \
      --bnostril x,y,w,h  --tnostril r0,r1,c0,c1  [--bscale 0.5]
"""
import argparse, os, sys, glob, csv
import numpy as np
import cv2
import fnv, fnv.file
from scipy.io import savemat


def find_one(folder, pat):
    h = glob.glob(os.path.join(folder, pat))
    if not h:
        sys.exit(f"no {pat} in {folder}")
    return h[0]


def read_thermal(ats, troi):
    im = fnv.file.ImagerFile(ats)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    H, W, N = im.height, im.width, im.num_frames
    r0, r1, c0, c1 = troi
    mot = np.zeros(N); nos = np.zeros(N); t = np.zeros(N)
    prev = None; t0 = None
    for i in range(N):
        im.get_frame(i)
        # COPY: im.final is a view onto a buffer fnv overwrites on the next
        # get_frame(); without copy, prev would alias it and the diff is 0.
        fr = np.array(im.final, dtype=np.float32, copy=True).reshape(H, W)
        if prev is not None:
            mot[i] = np.abs(fr - prev).mean()
        prev = fr
        nos[i] = fr[r0:r1, c0:c1].mean()
        ti = im.frame_info.time
        if t0 is None:
            t0 = ti
        t[i] = (ti - t0).total_seconds()
        if i % 5000 == 0:
            print(f"  thermal {i}/{N}", flush=True)
    return mot, nos, t, H, W, N


def read_basler(avi, broi, scale):
    cap = cv2.VideoCapture(avi)
    if not cap.isOpened():
        sys.exit("cannot open basler avi")
    x, y, w, h = broi
    mot = []; nos = []
    prev = None; n = 0
    while True:
        ok = cap.grab()
        if not ok:
            break
        ok, fr = cap.retrieve()
        if not ok:
            break
        g = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY)
        nos.append(float(g[y:y + h, x:x + w].mean()))
        gs = cv2.resize(g, None, fx=scale, fy=scale, interpolation=cv2.INTER_AREA).astype(np.float32)
        mot.append(0.0 if prev is None else float(np.abs(gs - prev).mean()))
        prev = gs
        n += 1
        if n % 5000 == 0:
            print(f"  basler {n}", flush=True)
    cap.release()
    return np.asarray(mot), np.asarray(nos), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--bnostril", required=True, help="basler nostril ROI 'x,y,w,h'")
    ap.add_argument("--tnostril", required=True, help="thermal nostril ROI 'r0,r1,c0,c1' (native px)")
    ap.add_argument("--bscale", type=float, default=0.5)
    args = ap.parse_args()

    folder = args.folder
    ats = find_one(folder, "*.ats"); avi = find_one(folder, "*.avi")
    csvp = find_one(folder, "timestamps.csv")
    broi = tuple(int(v) for v in args.bnostril.split(","))
    troi = tuple(int(v) for v in args.tnostril.split(","))

    print("thermal ...", flush=True)
    tmot, tnos, tt, H, W, Nt = read_thermal(ats, troi)
    print("basler ...", flush=True)
    bmot, bnos, Nb = read_basler(avi, broi, args.bscale)

    cam_ns, wall_s = [], []
    with open(csvp) as f:
        for row in csv.DictReader(f):
            cam_ns.append(int(row["camera_timestamp_ns"]))
            wall_s.append(float(row["wall_time_s"]))

    run = os.path.basename(os.path.abspath(os.path.normpath(folder)))
    out = os.path.join(folder, run + "_syncmot.mat")
    savemat(out, {
        "thermal_motion": tmot, "thermal_nostril": tnos, "thermal_t": tt,
        "basler_motion": bmot, "basler_nostril": bnos,
        "basler_cam_ns": np.asarray(cam_ns, np.int64),
        "basler_wall_s": np.asarray(wall_s, np.float64),
        "bnostril": np.asarray(broi), "tnostril": np.asarray(troi),
        "Nt": int(Nt), "Nb": int(Nb), "Ht": int(H), "Wt": int(W),
    }, do_compression=True)
    print(f"saved {out}  (thermal {Nt}, basler {Nb})")


if __name__ == "__main__":
    main()
