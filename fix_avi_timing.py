#!/usr/bin/env python3
"""
fix_avi_timing.py
=================
Fix the WRONG playback timing of Basler FFV1 recordings.

The problem: basler_dual_acq.py stamps the AVI with a NOMINAL frame rate
(e.g. 200 fps) but the camera actually runs slower (~111 fps here). So players
(VLC etc.) play the video ~1.8x too fast and show the wrong duration. The real
per-frame times live in the sibling `timestamps.csv` -- not in the AVI header.

The fix: rewrite the video PTS at the TRUE constant frame rate (measured from
timestamps.csv) using ffmpeg's `setts` bitstream filter. This is a LOSSLESS
STREAM COPY -- the FFV1 frames are not re-encoded, only re-timed -- so it is
fast (~15-20x realtime) and bit-faithful. Output is an .mkv next to the source;
the original .avi is never touched.

Note: this is a VIEWING fix. The analysis pipeline already uses timestamps.csv
directly, so it does not need this.

Usage:
  python fix_avi_timing.py <run_folder | file.avi | parent_dir>   [--overwrite]
        [--fps F]   # override measured rate (else computed from timestamps.csv)
        [--suffix _fixedfps]

  # one run:
  python fix_avi_timing.py "C:\\...\\cam1_20260602_180811_run001"
  # batch every run under a parent folder:
  python fix_avi_timing.py "C:\\...\\260602_orofacial_ChAT_ChroME"
"""

import os
import sys
import glob
import argparse
import subprocess
import numpy as np


def measure_fps(ts_path):
    """True constant fps = (N-1) / (t_last - t_first), from camera timestamps."""
    arr = np.genfromtxt(ts_path, delimiter=",", skip_header=1)  # idx, cam_ns, wall_s
    cam_ns = arr[:, 1]
    span_s = (cam_ns[-1] - cam_ns[0]) / 1e9
    n = len(cam_ns)
    return (n - 1) / span_s, n, span_s


def find_jobs(path):
    """Return list of (avi, timestamps_csv) pairs to process."""
    if os.path.isfile(path) and path.lower().endswith(".avi"):
        ts = os.path.join(os.path.dirname(path), "timestamps.csv")
        return [(path, ts)]
    # directory: find every .avi that has a sibling timestamps.csv
    jobs = []
    for avi in glob.glob(os.path.join(path, "**", "*.avi"), recursive=True):
        ts = os.path.join(os.path.dirname(avi), "timestamps.csv")
        if os.path.isfile(ts):
            jobs.append((avi, ts))
    return jobs


def fix_one(avi, ts, fps_override, suffix, overwrite):
    stem, _ = os.path.splitext(avi)
    out = stem + suffix + ".mkv"
    if os.path.isfile(out) and not overwrite:
        print(f"SKIP (exists): {out}   [--overwrite to redo]")
        return True

    if fps_override is not None:
        fps, n, span = fps_override, None, None
        print(f"{os.path.basename(avi)}: using override fps={fps:.5f}")
    elif os.path.isfile(ts):
        fps, n, span = measure_fps(ts)
        print(f"{os.path.basename(avi)}: {n} frames / {span:.3f} s -> {fps:.5f} fps")
    else:
        print(f"ERROR: no timestamps.csv beside {avi} and no --fps given. Skipping.")
        return False

    # setts rewrites PTS/DTS to frame_index / fps (in the stream timebase TB),
    # copying the FFV1 stream untouched (lossless, no re-encode).
    expr = f"setts=pts=N/({fps:.6f}*TB):dts=N/({fps:.6f}*TB)"
    cmd = ["ffmpeg", "-y", "-i", avi, "-map", "0", "-c", "copy",
           "-bsf:v", expr, out]
    print("  ffmpeg re-timing (lossless copy)...")
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        print("  FFMPEG FAILED:\n" + (r.stderr or "")[-1500:])
        return False
    print(f"  -> {out}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="run folder, a .avi file, or a parent dir (batch)")
    ap.add_argument("--fps", type=float, default=None,
                    help="override fps (default: measured from timestamps.csv)")
    ap.add_argument("--suffix", default="_fixedfps")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    jobs = find_jobs(a.path)
    if not jobs:
        sys.exit(f"No .avi (with sibling timestamps.csv) found under: {a.path}")
    print(f"Found {len(jobs)} video(s) to re-time.\n")

    ok = 0
    for avi, ts in jobs:
        ok += fix_one(avi, ts, a.fps, a.suffix, a.overwrite)
        print()
    print(f"Done: {ok}/{len(jobs)} succeeded.")


if __name__ == "__main__":
    main()
