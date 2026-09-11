#!/usr/bin/env python3
"""
orofacial_crop_extract.py
=========================
Reads a cropped, time-clipped ROI out of a Basler FFV1 .avi and saves it as a
compact .mat. Normally CALLED BY orofacial_crop_analyze.m (one-click MATLAB
entry point), but can also be run standalone.

Why Python: Basler recordings are FFV1-encoded (lossless); MATLAB's VideoReader
cannot decode them, OpenCV (ffmpeg backend) can. Reading happens here; all
analysis (per-pixel dF/F, spatial FFT, time maps, SVD) happens in MATLAB.

All settings can be overridden from the command line (MATLAB passes them); the
constants below are the standalone defaults.

  python orofacial_crop_extract.py [--avi ...] [--out ...] [--t_start S]
        [--t_end S] [--temporal_stride N] [--spatial_bin N] [--roi x,y,w,h]
"""

import os
import sys
import argparse
import numpy as np
import cv2
from scipy.io import savemat

# ------------------------- STANDALONE DEFAULTS -------------------------
RUN_DIR  = r"C:\Users\Admin\Desktop\260602_orofacial_ChAT_ChroME\cam1\cam1_20260602_180811_run001"
AVI_PATH = os.path.join(RUN_DIR, "cam1_20260602_180811_run001.avi")
TS_PATH  = os.path.join(RUN_DIR, "timestamps.csv")
OUT_MAT  = os.path.join(RUN_DIR, "crop_for_analysis.mat")
T_START_S, T_END_S = 0.0, None     # manual time-window clip (s); None=full
TEMPORAL_STRIDE, SPATIAL_BIN = 1, 1
ROI_XYWH = None                    # None -> interactive; else (x,y,w,h)
PREVIEW_NFRAMES = 60
MAX_BYTES = 1.6e9                  # cube size budget (MAT v5 2 GB/var limit)
# -----------------------------------------------------------------------


def load_timestamps(ts_path):
    """Per-frame time in seconds, relative to the first frame (camera clock)."""
    arr = np.genfromtxt(ts_path, delimiter=",", skip_header=1)  # idx, cam_ns, wall_s
    cam_ns = arr[:, 1]
    return (cam_ns - cam_ns[0]) / 1e9


def video_props(avi_path):
    """(n_frames, header_fps, width, height) from OpenCV."""
    cap = cv2.VideoCapture(avi_path)
    if not cap.isOpened():
        sys.exit(f"ERROR: OpenCV could not open {avi_path}")
    n = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    fps = cap.get(cv2.CAP_PROP_FPS)
    w = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    h = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    cap.release()
    return n, fps, w, h


def build_timeaxis(avi_path, ts_path, fps_override):
    """Return (t_all_seconds, fps_eff, source_str).

    Priority: explicit --fps  >  timestamps.csv  >  AVI header fps.
    For --fps / header (no per-frame clock) a uniform time axis is synthesized.
    Needed because Pylon-Viewer recordings have no timestamps.csv and a header
    fps that can be bogus (e.g. a 2P-frame-triggered breath cam stamped 98.6 Hz
    but actually delivered at the 30 Hz imaging rate -- pass --fps 30).
    """
    n_hdr, fps_hdr, _, _ = video_props(avi_path)
    if fps_override and fps_override > 0:
        return np.arange(n_hdr) / fps_override, float(fps_override), \
            f"--fps override ({fps_override:g} Hz)"
    if ts_path and os.path.isfile(ts_path):
        t = load_timestamps(ts_path)
        fps = (len(t) - 1) / (t[-1] - t[0]) if len(t) > 1 else fps_hdr
        return t, float(fps), "timestamps.csv (camera clock)"
    if fps_hdr and fps_hdr > 0 and n_hdr > 0:
        return np.arange(n_hdr) / fps_hdr, float(fps_hdr), \
            f"AVI header ({fps_hdr:.2f} Hz) -- no timestamps.csv"
    sys.exit("No timestamps.csv, no --fps, and unreadable header fps. Pass --fps.")


def build_preview(avi_path, nframes):
    """Average the first `nframes` mono frames and contrast-stretch to uint8."""
    cap = cv2.VideoCapture(avi_path)
    if not cap.isOpened():
        sys.exit(f"ERROR: OpenCV could not open {avi_path}")
    acc, count = None, 0
    while count < nframes:
        ok, frame = cap.read()
        if not ok:
            break
        g = frame[:, :, 0].astype(np.float64)   # Basler mono -> 3 equal channels
        acc = g if acc is None else acc + g
        count += 1
    cap.release()
    if acc is None:
        sys.exit("ERROR: could not read any frames for preview.")
    mean_img = acc / count
    lo, hi = np.percentile(mean_img, [1, 99])
    return (np.clip((mean_img - lo) / max(hi - lo, 1e-6), 0, 1) * 255).astype(np.uint8)


def select_roi(preview):
    """Drag a box; ENTER/SPACE confirm. Returns (x, y, w, h)."""
    disp = cv2.cvtColor(preview, cv2.COLOR_GRAY2BGR)
    x, y, w, h = cv2.selectROI("Select ROI  (drag box, ENTER=ok, C=cancel)",
                               disp, showCrosshair=True, fromCenter=False)
    cv2.destroyAllWindows()
    if w == 0 or h == 0:
        sys.exit("No ROI selected -- aborting.")
    return int(x), int(y), int(w), int(h)


def parse_args():
    none = lambda s: None if str(s).lower() in ("", "none") else float(s)
    roi  = lambda s: None if s is None else tuple(int(v) for v in s.split(","))
    ap = argparse.ArgumentParser()
    ap.add_argument("--avi", default=AVI_PATH)
    ap.add_argument("--ts", default=TS_PATH)
    ap.add_argument("--out", default=OUT_MAT)
    ap.add_argument("--t_start", type=float, default=T_START_S)
    ap.add_argument("--t_end", type=none, default=T_END_S)
    ap.add_argument("--temporal_stride", type=int, default=TEMPORAL_STRIDE)
    ap.add_argument("--spatial_bin", type=int, default=SPATIAL_BIN)
    ap.add_argument("--roi", type=roi, default=ROI_XYWH)
    ap.add_argument("--fps", type=float, default=0.0,
                    help="force fps (overrides timestamps/header); 0 = auto")
    ap.add_argument("--fullframe", action="store_true",
                    help="use the entire frame (no crop, no interactive ROI)")
    ap.add_argument("--max_bytes", type=float, default=MAX_BYTES)
    return ap.parse_args()


def main():
    a = parse_args()
    if not os.path.isfile(a.avi):
        sys.exit(f"ERROR: video not found: {a.avi}")

    t_all, fps_eff, src = build_timeaxis(a.avi, a.ts, a.fps)
    print(f"Recording: {len(t_all)} frames, {t_all[-1]:.2f} s, "
          f"{fps_eff:.2f} fps  [time source: {src}]")

    # --- manual time window -> frame indices ---
    t_end = t_all[-1] if a.t_end is None else a.t_end
    f0 = int(np.searchsorted(t_all, a.t_start, side="left"))
    f1 = int(np.searchsorted(t_all, t_end, side="right"))
    sel_idx = np.arange(f0, f1, a.temporal_stride)
    n_keep = len(sel_idx)
    if n_keep == 0:
        sys.exit("Time window selects 0 frames -- check --t_start / --t_end.")
    print(f"Clip [{a.t_start}, {t_end:.2f}] s -> frames [{f0}, {f1}) stride "
          f"{a.temporal_stride} -> {n_keep} frames")

    # --- ROI ---
    if a.fullframe:
        _, _, fw, fh = video_props(a.avi)
        x, y, w, h = 0, 0, fw, fh
        print(f"Full frame (no crop): {fw}x{fh}")
    elif a.roi is None:
        print("Building preview for ROI selection...")
        x, y, w, h = select_roi(build_preview(a.avi, PREVIEW_NFRAMES))
    else:
        x, y, w, h = a.roi
    print(f"ROI (x,y,w,h) = ({x},{y},{w},{h})")

    # --- pick spatial bin so the cube fits the size budget ---
    sbin = max(1, int(a.spatial_bin))
    while n_keep * (h // sbin) * (w // sbin) > a.max_bytes:
        sbin += 1
    h_b, w_b = h // sbin, w // sbin
    print(f"Spatial bin = {sbin} -> ROI {w_b}x{h_b}, cube ~{n_keep*h_b*w_b/1e6:.0f} MB")

    # --- read sequentially, crop, bin, store ---
    cap = cv2.VideoCapture(a.avi)
    if not cap.isOpened():
        sys.exit(f"ERROR: OpenCV could not open {a.avi}")
    mov = np.empty((h_b, w_b, n_keep), dtype=np.uint8)   # [H, W, T] for MATLAB
    sel_set = set(int(i) for i in sel_idx)
    write_col, idx = 0, 0
    while idx < f1:
        ok, frame = cap.read()
        if not ok:
            break
        if idx in sel_set:
            g = frame[y:y + h, x:x + w, 0]
            if sbin > 1:
                g = cv2.resize(g, (w_b, h_b), interpolation=cv2.INTER_AREA)
            mov[:, :, write_col] = g
            write_col += 1
            if write_col % 2000 == 0:
                print(f"  ...{write_col}/{n_keep} frames")
        idx += 1
    cap.release()
    mov = mov[:, :, :write_col]
    t_keep = t_all[sel_idx[:write_col]] - t_all[sel_idx[0]]
    print(f"Stored {write_col} frames.")

    savemat(a.out, {
        "mov": mov,                          # uint8 [H, W, T]
        "t_s": t_keep.astype(np.float64),    # [T] seconds, from camera timestamps
        "fps_eff": float(fps_eff),
        "roi_xywh": np.array([x, y, w, h], dtype=np.float64),
        "spatial_bin": float(sbin),
        "temporal_stride": float(a.temporal_stride),
        "t_start_s": float(a.t_start),
        "src_avi": a.avi,
    }, do_compression=True)
    print(f"Saved -> {a.out}")


if __name__ == "__main__":
    main()
