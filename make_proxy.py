#!/usr/bin/env python3
"""
make_proxy.py
=============
Make a small, smooth, correctly-timed H.264 .mp4 PROXY of Basler FFV1 recordings
for VISUAL INSPECTION (scrub / drag / frame-step in any player).

Why: the lossless FFV1 .avi is huge (~10 GB) and VLC's FFV1 decoder can glitch
("snow"). An H.264 proxy plays cleanly everywhere, is ~50-100x smaller, scrubs
instantly, and -- because we read the TRUE fps from timestamps.csv -- plays at
the correct speed/duration (the AVI header fps is wrong; see fix_avi_timing.py).

This is for EYEBALLING ONLY. Quantitative analysis uses the original FFV1 +
timestamps.csv, never this proxy.

Usage:
  python make_proxy.py <run_folder | file.avi | parent_dir>   [--overwrite]
        [--crf 20]      # quality: lower=better/bigger (18 visually lossless, 23 small)
        [--fps F]       # override (else measured from timestamps.csv)
        [--scale 1.0]   # downscale factor for even smaller files (e.g. 0.5)
        [--suffix _proxy]
"""

import os
import sys
import glob
import argparse
import subprocess
import numpy as np


def measure_fps(ts_path):
    arr = np.genfromtxt(ts_path, delimiter=",", skip_header=1)  # idx, cam_ns, wall_s
    cam_ns = arr[:, 1]
    return (len(cam_ns) - 1) / ((cam_ns[-1] - cam_ns[0]) / 1e9)


def find_jobs(path):
    if os.path.isfile(path) and path.lower().endswith(".avi"):
        return [(path, os.path.join(os.path.dirname(path), "timestamps.csv"))]
    jobs = []
    for avi in glob.glob(os.path.join(path, "**", "*.avi"), recursive=True):
        ts = os.path.join(os.path.dirname(avi), "timestamps.csv")
        jobs.append((avi, ts if os.path.isfile(ts) else None))
    return jobs


def make_one(avi, ts, args):
    stem, _ = os.path.splitext(avi)
    out = stem + args.suffix + ".mp4"
    if os.path.isfile(out) and not args.overwrite:
        print(f"SKIP (exists): {out}   [--overwrite to redo]")
        return True

    if args.fps is not None:
        fps = args.fps
    elif ts and os.path.isfile(ts):
        fps = measure_fps(ts)
    else:
        print(f"SKIP: no timestamps.csv beside {os.path.basename(avi)} and no --fps; "
              f"pass --fps to proxy this one (timing unknown).")
        return False
    print(f"{os.path.basename(avi)}: fps={fps if fps else 'header'} crf={args.crf}")

    # even dimensions required by yuv420p; optional downscale
    vf = "scale=trunc(iw*{s}/2)*2:trunc(ih*{s}/2)*2".format(s=args.scale) \
         if args.scale != 1.0 else "scale=trunc(iw/2)*2:trunc(ih/2)*2"

    cmd = ["ffmpeg", "-y"]
    if fps is not None:
        cmd += ["-r", f"{fps:.6f}"]
    cmd += ["-i", avi, "-an", "-vf", vf,
            "-c:v", "libx264", "-preset", "veryfast", "-crf", str(args.crf),
            "-pix_fmt", "yuv420p", "-movflags", "+faststart", out]
    print("  encoding H.264 proxy...")
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        print("  FFMPEG FAILED:\n" + (r.stderr or "")[-1500:])
        return False
    mb = os.path.getsize(out) / 1e6
    print(f"  -> {out}  ({mb:.0f} MB)")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="run folder, a .avi file, or a parent dir (batch)")
    ap.add_argument("--crf", type=int, default=20)
    ap.add_argument("--fps", type=float, default=None)
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--suffix", default="_proxy")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    jobs = find_jobs(args.path)
    if not jobs:
        sys.exit(f"No .avi found under: {args.path}")
    print(f"Found {len(jobs)} video(s).\n")
    ok = 0
    for avi, ts in jobs:
        ok += make_one(avi, ts, args)
        print()
    print(f"Done: {ok}/{len(jobs)} succeeded.")


if __name__ == "__main__":
    main()
