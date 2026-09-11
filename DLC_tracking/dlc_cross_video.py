#!/usr/bin/env python3
"""
dlc_cross_video.py  --  overlay DLC tracked dots as CROSS markers on each video.

Reads each video's DLC .csv (save_as_csv output; parsed without pandas so it runs
in the flir env too), draws a colored cross at each bodypart whose likelihood >
pmin, and writes <stem>_cross.mp4 next to the video. Decimated for a viewable size.

  C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe dlc_cross_video.py "<...config.yaml>" \
      [--pmin 0.5] [--decim 4] [--play-fps 60] [--halflen 9] [--thick 2] [--loop]

--loop keeps re-checking and makes a cross video as soon as each video's CSV
appears (useful while analyze_videos is still running).
"""
import os, sys, csv, glob, time, argparse
import numpy as np
import cv2

# distinct BGR per dot (L = greens/cyan, R = reds/orange)
COLORS = {"vL0": (0, 255, 0), "vL1": (255, 255, 0), "vR0": (0, 0, 255), "vR1": (0, 165, 255)}
DEFCOL = (0, 255, 255)


def read_dlc_csv(path):
    """returns (bodyparts list, data [T, nbp, 3] = x,y,likelihood)."""
    with open(path) as f:
        rows = list(csv.reader(f))
    # row0=scorer, row1=bodyparts, row2=coords(x,y,likelihood), row3+=data
    bps = rows[1][1:]                      # repeated x3
    parts = bps[0::3]
    data = np.array([[float(v) for v in r[1:]] for r in rows[3:]], dtype=np.float64)
    T = data.shape[0]
    data = data.reshape(T, len(parts), 3)
    return parts, data


def find_csv(video):
    d = os.path.dirname(video)
    stem = os.path.splitext(os.path.basename(video))[0]
    hits = [c for c in glob.glob(os.path.join(d, stem + "*DLC*.csv"))]
    return hits[0] if hits else None


def make_cross_video(video, csvpath, pmin, decim, play_fps, halflen, thick):
    parts, data = read_dlc_csv(csvpath)
    cap = cv2.VideoCapture(video)
    W = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)); H = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    out = os.path.splitext(video)[0] + "_cross.mp4"
    vw = cv2.VideoWriter(out, cv2.VideoWriter_fourcc(*"mp4v"), play_fps, (W, H), isColor=True)
    n = 0; T = data.shape[0]
    while True:
        ok = cap.grab()
        if not ok:
            break
        if (n % decim) == 0:
            ok, fr = cap.retrieve()
            if not ok:
                break
            if n < T:
                for j, bp in enumerate(parts):
                    x, y, p = data[n, j]
                    if p < pmin or not np.isfinite(x):
                        continue
                    x, y = int(round(x)), int(round(y))
                    c = COLORS.get(bp, DEFCOL)
                    cv2.line(fr, (x - halflen, y), (x + halflen, y), c, thick, cv2.LINE_AA)
                    cv2.line(fr, (x, y - halflen), (x, y + halflen), c, thick, cv2.LINE_AA)
                vw.write(fr)
        n += 1
    cap.release(); vw.release()
    return out


def videos_of(config):
    scope = os.path.dirname(os.path.dirname(config))
    return sorted(glob.glob(os.path.join(scope, "**", "*.avi"), recursive=True))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--pmin", type=float, default=0.5)
    ap.add_argument("--decim", type=int, default=4)
    ap.add_argument("--play-fps", type=float, default=60.0)
    ap.add_argument("--halflen", type=int, default=9)
    ap.add_argument("--thick", type=int, default=2)
    ap.add_argument("--loop", action="store_true")
    ap.add_argument("--sleep", type=int, default=120)
    args = ap.parse_args()

    vids = videos_of(args.config)
    print(f"{len(vids)} videos for {args.config}", flush=True)
    while True:
        pending = 0
        for v in vids:
            out = os.path.splitext(v)[0] + "_cross.mp4"
            if os.path.isfile(out):
                continue
            c = find_csv(v)
            if not c:
                pending += 1; continue
            print(f"  cross video: {os.path.basename(v)}", flush=True)
            try:
                o = make_cross_video(v, c, args.pmin, args.decim, args.play_fps, args.halflen, args.thick)
                print(f"    saved {o}", flush=True)
            except Exception as e:
                print(f"    FAILED {v}: {e}", flush=True)
        if not args.loop or pending == 0:
            break
        print(f"  {pending} videos awaiting CSV; sleeping {args.sleep}s", flush=True)
        time.sleep(args.sleep)
    print("=== cross videos done ===", flush=True)


if __name__ == "__main__":
    main()
