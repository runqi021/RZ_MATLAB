#!/usr/bin/env python
"""
thermal_inspect_videos.py  --  make two viewable mp4s from a FLIR .ats so you can
EXAMINE the thermal recording (and decide how to track from it) without the FLIR GUI.

Reads the .ats as TRUE temperature (deg C, no AGC/gamma), reports timing / dropped
frames, then writes:
  <stem>_track.mp4    grayscale, robust LINEAR window  -> best for SEEING structure
                      (nostril edges / airflow) to judge a tracking approach
  <stem>_heat.mp4     inferno heatmap, fixed temperature window -> the radiometric view

Timestamps are read for EVERY frame (cheap) so the drop report is exact; pixels are
only decoded every --decim-th frame (fast) for the videos.

Run in the FileSDK env:
  C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe thermal_inspect_videos.py <in.ats> \
      [--decim 4] [--play-fps 100] [--upscale 6]

Companion to thermal_ats_to_mat.py / thermal_breath_svd.m (those are NOT touched).
"""
import argparse, os, sys
import numpy as np
import cv2
import fnv, fnv.file


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ats", help="input .ats file")
    ap.add_argument("--decim", type=int, default=4,
                    help="decode every Nth frame for the video (default 4 -> 400Hz becomes 100Hz)")
    ap.add_argument("--play-fps", type=float, default=100.0,
                    help="playback frame rate of the mp4 (default 100)")
    ap.add_argument("--upscale", type=int, default=6,
                    help="nearest-neighbour zoom (56x176 is tiny); default 6")
    ap.add_argument("--clip-pct", type=float, nargs=2, default=[0.5, 99.7],
                    help="robust percentile window over the sampled frames")
    ap.add_argument("--outdir", default=None, help="output dir (default: alongside input)")
    args = ap.parse_args()

    if not os.path.isfile(args.ats):
        sys.exit(f"no such file: {args.ats}")
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.ats))
    stem = os.path.splitext(os.path.basename(args.ats))[0]

    im = fnv.file.ImagerFile(args.ats)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    H, W, N = im.height, im.width, im.num_frames
    print(f"reading {args.ats}\n  {N} frames, {H}x{W} px", flush=True)

    times = np.empty(N, dtype=np.float64)
    t0 = None
    frames = []          # sampled frames (every decim-th) for the video
    fidx_kept = []
    for i in range(N):
        im.get_frame(i)
        ti = im.frame_info.time
        if t0 is None:
            t0 = ti
        times[i] = (ti - t0).total_seconds()
        if (i % args.decim) == 0:
            # COPY: im.final is a view onto a buffer fnv overwrites on the next
            # get_frame(); without copy every stored frame aliases the LAST one
            # (-> a static "single image" video).
            frames.append(np.array(im.final, dtype=np.float32, copy=True).reshape(H, W))
            fidx_kept.append(i)
        if (i % 5000) == 0:
            print(f"  {i}/{N}", flush=True)

    # ---- timing / dropped-frame report (every frame) ----
    dt = np.diff(times)
    med = np.median(dt)
    fps_nom = round(1.0 / med)
    fps_mean = (N - 1) / (times[-1] - times[0])
    n_expected = int(round(times[-1] * fps_nom)) + 1
    n_drop = n_expected - N
    gaps = int((dt > 1.5 * med).sum())
    print("=== THERMAL TIMING ===")
    print(f"  frames present : {N}")
    print(f"  fps            : {fps_mean:.3f} (nominal {fps_nom}, median dt {med*1e3:.3f} ms)")
    print(f"  duration       : {times[-1]:.4f} s")
    print(f"  expected@nom    : {n_expected}   DROPPED {n_drop}  (gaps>1.5x median: {gaps})")

    stack = np.stack(frames, axis=0)                       # [t x H x W]
    lo = np.percentile(stack, args.clip_pct[0])
    hi = np.percentile(stack, args.clip_pct[1])
    print(f"  temp window    : [{lo:.2f} {hi:.2f}] C  ({len(frames)} frames in videos)")

    up = args.upscale
    sz = (W * up, H * up)
    track_path = os.path.join(outdir, stem + "_track.mp4")
    heat_path = os.path.join(outdir, stem + "_heat.mp4")
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    vt = cv2.VideoWriter(track_path, fourcc, args.play_fps, sz, isColor=True)
    vh = cv2.VideoWriter(heat_path, fourcc, args.play_fps, sz, isColor=True)

    for fr in stack:
        g = np.clip((fr - lo) / (hi - lo), 0, 1)
        g8 = (g * 255).astype(np.uint8)
        g8u = cv2.resize(g8, sz, interpolation=cv2.INTER_NEAREST)
        vt.write(cv2.cvtColor(g8u, cv2.COLOR_GRAY2BGR))           # grayscale tracking view
        vh.write(cv2.applyColorMap(g8u, cv2.COLORMAP_INFERNO))    # heatmap view
    vt.release()
    vh.release()
    dur = len(frames) / args.play_fps
    print(f"saved {track_path}")
    print(f"saved {heat_path}   ({dur:.1f}s at {args.play_fps:.0f} fps, decim {args.decim})")


if __name__ == "__main__":
    main()
