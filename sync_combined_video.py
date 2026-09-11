#!/usr/bin/env python
"""
sync_combined_video.py  --  one synchronized mp4 with the Basler face (top) and the
thermal snout (bottom), frame-paired off the shared 400 Hz WFG trigger so you can
watch both cameras together and eyeball the alignment.

Pairing:  basler_frame = thermal_frame + LAG   (default LAG=0: both cameras caught
the same first pulse; thermal just stopped earlier). A wrong LAG shows up as shared
movements happening at different moments in the two panels -> this video is also a
sync check.

Reads thermal with an explicit COPY (im.final is a reused buffer view; without the
copy every frame aliases the last one -> a static "single image" video).

  C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe sync_combined_video.py <run_folder> \
      [--decim 4] [--play-fps 60] [--lag 0] [--width 640]
"""
import argparse, os, sys, glob
import numpy as np
import cv2
import fnv, fnv.file


def find_one(folder, pat):
    h = glob.glob(os.path.join(folder, pat))
    if not h:
        sys.exit(f"no {pat} in {folder}")
    return h[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder")
    ap.add_argument("--decim", type=int, default=4)
    ap.add_argument("--play-fps", type=float, default=60.0)
    ap.add_argument("--lag", type=int, default=0, help="basler_frame = thermal_frame + LAG")
    ap.add_argument("--width", type=int, default=640, help="output width (px)")
    ap.add_argument("--clip-pct", type=float, nargs=2, default=[0.5, 99.7])
    args = ap.parse_args()

    folder = args.folder
    ats = find_one(folder, "*.ats")
    avi = find_one(folder, "*.avi")
    run = os.path.basename(os.path.abspath(os.path.normpath(folder)))

    # ---- thermal: keep decimated frames (COPY!) in a dict keyed by frame index ----
    im = fnv.file.ImagerFile(ats)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    Ht, Wt, Nt = im.height, im.width, im.num_frames
    keep = set(range(0, Nt, args.decim))
    ther = {}
    for i in range(Nt):
        im.get_frame(i)
        if i in keep:
            ther[i] = np.array(im.final, dtype=np.float32, copy=True).reshape(Ht, Wt)
        if i % 5000 == 0:
            print(f"  thermal {i}/{Nt}", flush=True)
    allt = np.stack(list(ther.values()))
    lo = np.percentile(allt, args.clip_pct[0]); hi = np.percentile(allt, args.clip_pct[1])
    print(f"thermal window [{lo:.2f} {hi:.2f}] C, {len(ther)} frames", flush=True)

    # ---- panel geometry ----
    Wc = args.width
    bh = int(round(Wc * 504 / 992))           # basler panel height (aspect of full frame)
    th = int(round(Wc * Ht / Wt))             # thermal panel height
    bar = 26                                   # label bar
    Hc = bar + bh + th
    clahe = cv2.createCLAHE(clipLimit=2.5, tileGridSize=(8, 8))

    out = os.path.join(folder, run + f"_combined_lag{args.lag}.mp4")
    vw = cv2.VideoWriter(out, cv2.VideoWriter_fourcc(*"mp4v"), args.play_fps, (Wc, Hc), isColor=True)

    cap = cv2.VideoCapture(avi)
    if not cap.isOpened():
        sys.exit("cannot open basler avi")
    Nb = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    n = 0; written = 0
    while True:
        ok = cap.grab()
        if not ok:
            break
        j = n - args.lag                       # thermal index paired with this basler frame
        if j in ther:
            ok, fr = cap.retrieve()
            if not ok:
                break
            # basler panel
            g = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY)
            g = clahe.apply(g)
            gb = cv2.resize(g, (Wc, bh), interpolation=cv2.INTER_AREA)
            gb = cv2.cvtColor(gb, cv2.COLOR_GRAY2BGR)
            # thermal panel
            tt = np.clip((ther[j] - lo) / (hi - lo), 0, 1)
            t8 = (tt * 255).astype(np.uint8)
            t8 = cv2.resize(t8, (Wc, th), interpolation=cv2.INTER_NEAREST)
            tc = cv2.applyColorMap(t8, cv2.COLORMAP_INFERNO)
            # compose
            canvas = np.zeros((Hc, Wc, 3), np.uint8)
            canvas[bar:bar + bh] = gb
            canvas[bar + bh:bar + bh + th] = tc
            tsec = j / 400.0
            cv2.putText(canvas, f"t={tsec:6.2f}s  basler#{n}  thermal#{j}  lag={args.lag}",
                        (6, 18), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 1, cv2.LINE_AA)
            cv2.putText(canvas, "BASLER", (6, bar + 16), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 0), 1, cv2.LINE_AA)
            cv2.putText(canvas, "THERMAL", (6, bar + bh + 16), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 0), 1, cv2.LINE_AA)
            vw.write(canvas)
            written += 1
        n += 1
        if n % 5000 == 0:
            print(f"  basler {n}/{Nb}", flush=True)
    cap.release(); vw.release()
    print(f"saved {out}  ({written} frames, {written/args.play_fps:.1f}s at {args.play_fps:.0f} fps)")


if __name__ == "__main__":
    main()
