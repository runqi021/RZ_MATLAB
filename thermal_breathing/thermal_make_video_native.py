#!/usr/bin/env python
"""
thermal_make_video_native.py  --  watch a FLIR .ats with NOTHING thrown away.

Reads every native frame as true temperature, places each frame on its REAL
time using the per-frame timestamps (dropped-frame gaps filled by interpolation
onto a uniform grid at the nominal rate -- NO decimation, full temporal detail),
and writes an .mp4 with a FIXED linear temperature window (no gamma, no per-frame
auto-gain). What you see is true relative temperature on a true clock.

Run in the flir env:
  C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe thermal_make_video_native.py <in.ats> [--play-fps 120]

Output: <stem>_view_native.mp4 next to the input.
"""
import argparse, os, sys
import numpy as np
import fnv, fnv.file
from scipy.interpolate import interp1d
import cv2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ats")
    ap.add_argument("--play-fps", type=float, default=120.0,
                    help="playback fps in the mp4 (native is ~400 Hz; 120 = ~3.3x slow-mo to inspect sniffs)")
    ap.add_argument("--upscale", type=int, default=6)
    ap.add_argument("--clip", type=float, nargs=2, default=[0.5, 99.7],
                    help="percentiles for the fixed global temperature window")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if not os.path.isfile(args.ats):
        sys.exit(f"no such file: {args.ats}")
    out = args.out or os.path.splitext(args.ats)[0] + "_view_native.mp4"

    im = fnv.file.ImagerFile(args.ats)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    H, W, N = im.height, im.width, im.num_frames
    print(f"reading {N} frames {H}x{W}", flush=True)

    data = np.empty((N, H * W), dtype=np.float32)
    t = np.empty(N, dtype=np.float64)
    t0 = None
    for i in range(N):
        im.get_frame(i)
        data[i, :] = np.asarray(im.final, dtype=np.float32)
        ti = im.frame_info.time
        if t0 is None:
            t0 = ti
        t[i] = (ti - t0).total_seconds()
        if i % 4000 == 0:
            print(f"  read {i}/{N}", flush=True)

    # --- true time axis: uniform grid at nominal rate, gaps interpolated ---
    dt = np.diff(t)
    fps_nom = float(round(1.0 / np.median(dt)))
    n_exp = int(round(t[-1] * fps_nom)) + 1
    n_drop = n_exp - N
    print(f"  nominal {fps_nom:.0f} Hz, {n_drop} dropped of {n_exp} "
          f"({100*n_drop/n_exp:.1f}%) -- filling on real timestamps", flush=True)
    t_uni = np.arange(n_exp, dtype=np.float64) / fps_nom
    data = interp1d(t, data, axis=0, kind='linear', bounds_error=False,
                    fill_value='extrapolate', assume_sorted=True)(t_uni).astype(np.float32)

    # --- fixed global linear window (no gamma, no per-frame scaling) ---
    lo, hi = np.percentile(data, args.clip)
    print(f"  window [{lo:.2f}, {hi:.2f}] C ; writing {data.shape[0]} frames @ {args.play_fps:.0f} fps", flush=True)

    Hs, Ws = H * args.upscale, W * args.upscale
    vw = cv2.VideoWriter(out, cv2.VideoWriter_fourcc(*'mp4v'), args.play_fps, (Ws, Hs))
    for k in range(data.shape[0]):
        g = (data[k].reshape(H, W) - lo) / (hi - lo)
        g8 = np.clip(g * 255.0, 0, 255).astype(np.uint8)
        rgb = cv2.applyColorMap(g8, cv2.COLORMAP_INFERNO)
        rgb = cv2.resize(rgb, (Ws, Hs), interpolation=cv2.INTER_NEAREST)
        vw.write(rgb)
    vw.release()
    dur = data.shape[0] / args.play_fps
    print(f"saved {out}  ({dur:.0f} s video, {dur*fps_nom/data.shape[0]:.2f}x slow-mo vs real time)", flush=True)


if __name__ == "__main__":
    main()
