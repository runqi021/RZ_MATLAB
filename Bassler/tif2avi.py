"""Convert multi-frame TIFF stacks to uncompressed grayscale AVI.

Usage:
    python tif2avi.py "C:\path\to\folder" --fps 90
    python tif2avi.py "C:\path\to\folder" --fps 90 --delete
"""

import argparse
import glob
import os
import sys

import cv2
import tifffile
import numpy as np


def convert_one(tif_path, fps):
    """Convert a single TIFF stack to AVI. Returns output path."""
    avi_path = os.path.splitext(tif_path)[0] + '.avi'
    if os.path.exists(avi_path):
        print(f"  SKIP (exists): {avi_path}")
        return avi_path

    print(f"  Reading {os.path.basename(tif_path)}...", end='', flush=True)
    stack = tifffile.imread(tif_path)
    nframes = stack.shape[0]
    h, w = stack.shape[1], stack.shape[2]
    print(f" {nframes} frames, {w}x{h}")

    fourcc = cv2.VideoWriter_fourcc(*'FFV1')  # lossless
    writer = cv2.VideoWriter(avi_path, fourcc, fps, (w, h), False)
    if not writer.isOpened():
        # Fallback to raw
        fourcc = cv2.VideoWriter_fourcc(*'GREY')
        writer = cv2.VideoWriter(avi_path, fourcc, fps, (w, h), False)
    if not writer.isOpened():
        # Last resort — MJPG (lossy but always works)
        fourcc = cv2.VideoWriter_fourcc(*'MJPG')
        writer = cv2.VideoWriter(avi_path, fourcc, fps, (w, h), False)

    for i in range(nframes):
        frame = stack[i]
        if frame.dtype != np.uint8:
            frame = (frame / frame.max() * 255).astype(np.uint8)
        writer.write(frame)

    writer.release()
    avi_mb = os.path.getsize(avi_path) / 1e6
    print(f"  -> {os.path.basename(avi_path)}  ({avi_mb:.0f} MB)")
    return avi_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', help='Folder containing .tif files')
    parser.add_argument('--fps', type=float, default=90)
    parser.add_argument('--delete', action='store_true',
                        help='Delete original .tif after successful conversion')
    args = parser.parse_args()

    tifs = sorted(glob.glob(os.path.join(args.folder, '*.tif')))
    tifs += sorted(glob.glob(os.path.join(args.folder, '*.tiff')))
    if not tifs:
        print(f"No .tif files in {args.folder}")
        sys.exit(1)

    print(f"Converting {len(tifs)} TIFF stacks -> AVI (fps={args.fps})\n")

    for tif_path in tifs:
        avi_path = convert_one(tif_path, args.fps)
        if args.delete and os.path.exists(avi_path) and os.path.getsize(avi_path) > 0:
            os.remove(tif_path)
            print(f"  Deleted {os.path.basename(tif_path)}")

    print(f"\nDone. {len(tifs)} files converted.")


if __name__ == '__main__':
    main()
