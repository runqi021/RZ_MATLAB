#!/usr/bin/env python
"""
basler_preview.py  --  make a small, viewable mp4 from a (large, lossless FFV1)
Basler run .avi so you can EXAMINE the whisker / nostril bright-dot markers and
decide a tracking approach. Also writes a mean + max-projection PNG (the max image
makes the bright dots pop, useful for choosing crop ROIs / thresholds).

MATLAB's VideoReader cannot decode FFV1, so this uses OpenCV (ffmpeg backend).

  python basler_preview.py <run.avi> [--decim 4] [--scale 0.5] [--play-fps 100]

Output (next to the .avi):
  <stem>_preview.mp4      decimated, downscaled, real-time-ish playback
  <stem>_proj.png         left: mean image, right: max image (dots are brightest)

Reuses the FFV1->OpenCV convention from orofacial_crop_extract.py.
"""
import argparse, os, sys
import numpy as np
import cv2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("avi", help="input Basler run .avi (FFV1)")
    ap.add_argument("--decim", type=int, default=4, help="keep every Nth frame (default 4)")
    ap.add_argument("--scale", type=float, default=0.5, help="spatial downscale (default 0.5)")
    ap.add_argument("--play-fps", type=float, default=100.0, help="playback fps (default 100)")
    ap.add_argument("--max-frames", type=int, default=0, help="stop after N source frames (0=all)")
    args = ap.parse_args()

    if not os.path.isfile(args.avi):
        sys.exit(f"no such file: {args.avi}")
    outdir = os.path.dirname(os.path.abspath(args.avi))
    stem = os.path.splitext(os.path.basename(args.avi))[0]

    cap = cv2.VideoCapture(args.avi)
    if not cap.isOpened():
        sys.exit("OpenCV could not open the AVI (FFV1/ffmpeg backend missing?)")
    W = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    H = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    Nhdr = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    fps_hdr = cap.get(cv2.CAP_PROP_FPS)
    print(f"basler {W}x{H}, header {Nhdr} frames @ {fps_hdr:.3f} fps", flush=True)

    ow = int(round(W * args.scale)); oh = int(round(H * args.scale))
    out_mp4 = os.path.join(outdir, stem + "_preview.mp4")
    vw = cv2.VideoWriter(out_mp4, cv2.VideoWriter_fourcc(*"mp4v"),
                         args.play_fps, (ow, oh), isColor=True)

    sum_img = np.zeros((H, W), dtype=np.float64)
    max_img = np.zeros((H, W), dtype=np.uint8)
    n_read = 0; n_kept = 0
    while True:
        ok = cap.grab()
        if not ok:
            break
        if (n_read % args.decim) == 0:
            ok, fr = cap.retrieve()
            if not ok:
                break
            g = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY)
            sum_img += g
            np.maximum(max_img, g, out=max_img)
            small = cv2.resize(fr, (ow, oh), interpolation=cv2.INTER_AREA)
            vw.write(small)
            n_kept += 1
        n_read += 1
        if (n_read % 5000) == 0:
            print(f"  {n_read} frames", flush=True)
        if args.max_frames and n_read >= args.max_frames:
            break
    cap.release(); vw.release()

    mean_img = (sum_img / max(1, n_kept))
    mean8 = np.clip(255 * (mean_img - mean_img.min()) / (np.ptp(mean_img) + 1e-9), 0, 255).astype(np.uint8)
    proj = np.hstack([mean8, max_img])
    proj_path = os.path.join(outdir, stem + "_proj.png")
    cv2.imwrite(proj_path, proj)

    print(f"read {n_read} frames, kept {n_kept}")
    print(f"saved {out_mp4}  ({n_kept/args.play_fps:.1f}s at {args.play_fps:.0f} fps)")
    print(f"saved {proj_path}  (left=mean, right=max; dots brightest in max)")


if __name__ == "__main__":
    main()
