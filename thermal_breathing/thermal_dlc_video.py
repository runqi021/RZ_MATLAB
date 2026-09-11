#!/usr/bin/env python
"""
thermal_dlc_video.py  --  turn a FLIR .ats into a video you can DLC-track the
NOSTRIL in. The point is a crisp, stable landmark, NOT radiometry: absolute
temperature (for breathing) is read straight from the .ats later in the tracked
ROI, so here we are free to use whatever "unit" makes the nostril easiest to see.

DEFAULT render = per-frame CLAHE local-contrast grayscale (mode 'clahe'):
  the raw fixed-temperature window wastes most of the dynamic range on the warm
  face, so the nostril edge is low contrast; CLAHE equalizes local contrast every
  frame so the nostril/face structure stays crisp and stable while the head moves.
Other modes: 'temp' (fixed global temperature window, == the old _track look),
             'stretch' (per-frame robust percentile stretch, no CLAHE).

CRITICAL — frame mapping: this writes ONE video frame per .ats frame (decim=1),
i.e. video frame k == .ats frame k. NO time-resampling / interpolation (that would
invent frames with no .ats index and break the temperature read-back). With
--decim D, video frame k == .ats frame k*D (written to the sidecar .json).

Run in the flir env:
  C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe thermal_dlc_video.py <in.ats> \
      [--mode clahe|temp|stretch] [--upscale 6] [--decim 1] [--play-fps 100]

Output: <stem>_dlc.mp4  +  <stem>_dlc.json (mapping / params), next to the input.
"""
import argparse, json, os, sys
import numpy as np
import cv2
import fnv, fnv.file


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ats")
    ap.add_argument("--mode", choices=["clahe", "temp", "stretch"], default="clahe")
    ap.add_argument("--upscale", type=int, default=6, help="zoom; thermal frames are tiny")
    ap.add_argument("--decim", type=int, default=1,
                    help="write every Dth .ats frame (default 1 = all; keeps trivial mapping)")
    ap.add_argument("--play-fps", type=float, default=100.0, help="mp4 playback fps")
    ap.add_argument("--interp", choices=["nearest", "cubic", "lanczos"], default="cubic",
                    help="upscale interpolation (cubic = smoother gradients for the CNN)")
    ap.add_argument("--clahe-clip", type=float, default=3.0)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    if not os.path.isfile(args.ats):
        sys.exit(f"no such file: {args.ats}")
    out = args.out or os.path.splitext(args.ats)[0] + "_dlc.mp4"
    jout = os.path.splitext(out)[0] + ".json"

    im = fnv.file.ImagerFile(args.ats)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS
    H, W, N = im.height, im.width, im.num_frames
    print(f"{os.path.basename(args.ats)}: {N} frames, {H}x{W} px  mode={args.mode}", flush=True)

    interp = {"nearest": cv2.INTER_NEAREST, "cubic": cv2.INTER_CUBIC,
              "lanczos": cv2.INTER_LANCZOS4}[args.interp]
    up = args.upscale
    sz = (W * up, H * up)

    # for mode 'temp' we need a fixed global window -> one decimated pre-pass
    lo = hi = None
    if args.mode == "temp":
        samp = []
        for i in range(0, N, max(1, N // 600)):
            im.get_frame(i)
            samp.append(np.array(im.final, dtype=np.float32, copy=True).reshape(H, W))
        lo, hi = (float(x) for x in np.percentile(np.stack(samp, 0), [0.5, 99.7]))
        print(f"  fixed temperature window [{lo:.2f}, {hi:.2f}] C", flush=True)

    clahe = cv2.createCLAHE(clipLimit=args.clahe_clip, tileGridSize=(4, 4))

    def render(fr):
        if args.mode == "temp":
            g = np.clip((fr - lo) / (hi - lo), 0, 1)
            g8 = (g * 255).astype(np.uint8)
        elif args.mode == "stretch":
            a, b = np.percentile(fr, [1, 99])
            g8 = np.clip((fr - a) / (b - a + 1e-6) * 255, 0, 255).astype(np.uint8)
        else:  # clahe
            a, b = np.percentile(fr, [1, 99])
            g8 = np.clip((fr - a) / (b - a + 1e-6) * 255, 0, 255).astype(np.uint8)
            g8 = clahe.apply(g8)
        return cv2.resize(g8, sz, interpolation=interp)

    vw = cv2.VideoWriter(out, cv2.VideoWriter_fourcc(*"mp4v"), args.play_fps, sz, isColor=True)
    if not vw.isOpened():
        sys.exit("VideoWriter failed to open")
    kept = list(range(0, N, args.decim))
    for k, i in enumerate(kept):
        im.get_frame(i)
        fr = np.array(im.final, dtype=np.float32, copy=True).reshape(H, W)
        vw.write(cv2.cvtColor(render(fr), cv2.COLOR_GRAY2BGR))
        if k % 4000 == 0:
            print(f"  {k}/{len(kept)}", flush=True)
    vw.release()

    meta = {"ats": os.path.abspath(args.ats), "mode": args.mode, "upscale": up,
            "interp": args.interp, "decim": args.decim, "play_fps": args.play_fps,
            "n_ats_frames": N, "n_video_frames": len(kept), "H": H, "W": W,
            "ats_frame_of_video_frame": "k*decim",
            "note": "absolute temperature for breathing comes from the .ats, NOT this video"}
    with open(jout, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"saved {out}  ({len(kept)} frames @ {args.play_fps:.0f} fps)\nsaved {jout}", flush=True)


if __name__ == "__main__":
    main()
