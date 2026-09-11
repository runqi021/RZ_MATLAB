#!/usr/bin/env python3
"""
make_viewable_video_copies.py
=============================
Make EASILY VIEWABLE copies of Basler behaviour videos, played back at the
CORRECT frame rate.

Two separate problems with the raw Basler recordings, both fixed here:

  1. CODEC. The raw runs are FFV1-in-AVI. FFV1 is lossless and great for
     storage, but VLC / ImageJ / MATLAB VideoReader all choke on it (see
     `feedback_ffv1_video_codec`). The copy is re-encoded to H.264/mp4, which
     every player, ImageJ (with the standard mp4 reader), MATLAB VideoReader,
     PowerPoint and browsers can open.

  2. FRAME RATE. The AVI header carries a NOMINAL rate (30 fps here) that is
     not the rate the camera actually ran at, so playback speed and duration
     are wrong. The true rate is recovered per-run (see below) and baked into
     the copy, so it plays in real time.

The ORIGINAL .avi files are never touched. These copies are for VIEWING only --
the analysis pipeline reads the raw frames and its own timing source.

True-fps priority chain
-----------------------
  1. --fps                    explicit override
  2. sibling timestamps.csv   measured (N-1)/(t_last-t_first)  [Basler free-run
                              recordings; see fix_avi_timing.py]
  3. sibling *_meta.mat       ScanImage scanFrameRate_raw -- correct when the
                              camera is 2P-triggered frame-for-frame, i.e. one
                              video frame per imaging frame. The script CHECKS
                              this (video frame count vs framesPerSlice) and
                              warns loudly when they disagree.
  4. otherwise                skipped, with a message

Usage
-----
  python make_viewable_video_copies.py <parent_dir | run_folder | file.avi>
        [--fps F]        force this frame rate for every video
        [--crf N]        H.264 quality, lower = better (default 17, visually
                         near-lossless; 0 = truly lossless, much bigger)
        [--preset P]     x264 speed/size preset (default veryfast)
        [--imagej]       ALSO write an MJPEG .avi, for ImageJ installs whose
                         mp4 reader is missing
        [--outdir DIR]   collect all copies here instead of beside each source
        [--suffix S]     default "_view"
        [--overwrite]    redo copies that already exist
        [--dry-run]      report what would be done, encode nothing

  # one experiment, every run under it:
  python make_viewable_video_copies.py "D:\\260721_Sert_soma_G8s\\phys"
"""

import os
import sys
import glob
import shutil
import argparse
import subprocess

import numpy as np

try:
    import scipy.io as sio
except ImportError:
    sio = None


# =========================================================================
#  ffmpeg discovery
# =========================================================================

def find_ffmpeg():
    """Return a WORKING ffmpeg path.

    The conda ffmpeg on this machine is present but crashes on launch
    (0xC0000139, bad DLL entry point), so a binary being on PATH is not
    enough -- each candidate is actually executed before being accepted.
    """
    candidates = []
    on_path = shutil.which("ffmpeg")
    if on_path:
        candidates.append(on_path)
    try:
        import imageio_ffmpeg
        candidates.append(imageio_ffmpeg.get_ffmpeg_exe())
    except Exception:
        pass

    for exe in candidates:
        try:
            r = subprocess.run([exe, "-version"], capture_output=True, timeout=30)
            if r.returncode == 0:
                return exe
        except Exception:
            continue
    return None


# =========================================================================
#  frame rate recovery
# =========================================================================

def fps_from_timestamps(ts_path):
    """True constant fps from the camera timestamp log."""
    arr = np.genfromtxt(ts_path, delimiter=",", skip_header=1)  # idx, cam_ns, wall_s
    cam_ns = arr[:, 1]
    span_s = (cam_ns[-1] - cam_ns[0]) / 1e9
    n = len(cam_ns)
    return (n - 1) / span_s, f"timestamps.csv ({n} fr / {span_s:.3f} s)"


def _pick_meta(folder):
    """Prefer the _meta.mat of the RAW tif over any _preproc/_MC variant."""
    hits = glob.glob(os.path.join(folder, "*_meta.mat"))
    if not hits:
        return None
    raw = [h for h in hits
           if not any(t in os.path.basename(h) for t in ("_preproc", "_MC", "_AVG"))]
    return (raw or hits)[0]


def fps_from_meta(folder, n_video):
    """ScanImage frame rate, valid because the camera is 2P frame-triggered.

    Returns (fps, note, n_imaging) -- n_imaging is None when it could not be
    read, so the caller can skip the consistency check.
    """
    if sio is None:
        return None, "scipy not available", None
    meta_path = _pick_meta(folder)
    if meta_path is None:
        return None, "no _meta.mat", None
    try:
        m = sio.loadmat(meta_path, squeeze_me=True)
    except Exception as e:
        return None, f"_meta.mat unreadable ({e})", None

    fps = m.get("scanFrameRate_raw")
    src = "scanFrameRate_raw"
    if fps is None or not np.isfinite(float(fps)) or float(fps) <= 0:
        fps = m.get("fps")
        src = "fps (rounded)"
    if fps is None or not np.isfinite(float(fps)) or float(fps) <= 0:
        return None, "_meta.mat has no usable frame rate", None

    n_img = None
    try:
        fps_slice = m.get("framesPerSlice")
        n_slices = m.get("numSlices")
        if fps_slice is not None and n_slices is not None:
            n_img = int(float(fps_slice) * float(n_slices))
    except Exception:
        pass

    return float(fps), f"{os.path.basename(meta_path)} [{src}]", n_img


def count_frames(avi):
    """Frame count of the source video (OpenCV's bundled ffmpeg reads FFV1)."""
    try:
        import cv2
    except ImportError:
        return None
    c = cv2.VideoCapture(avi)
    if not c.isOpened():
        return None
    n = int(c.get(cv2.CAP_PROP_FRAME_COUNT))
    c.release()
    return n or None


def resolve_fps(avi, fps_override):
    """Work down the priority chain. Returns (fps, note, warning_or_None)."""
    folder = os.path.dirname(avi)
    n_video = count_frames(avi)

    if fps_override is not None:
        return fps_override, "--fps override", None, n_video

    ts = os.path.join(folder, "timestamps.csv")
    if os.path.isfile(ts):
        try:
            fps, note = fps_from_timestamps(ts)
            return fps, note, None, n_video
        except Exception as e:
            print(f"    timestamps.csv unusable ({e}), falling back to _meta.mat")

    fps, note, n_img = fps_from_meta(folder, n_video)
    if fps is None:
        return None, note, None, n_video

    warn = None
    if n_video and n_img:
        if n_video != n_img:
            ratio = n_video / n_img
            warn = (f"video has {n_video} frames but imaging has {n_img} "
                    f"(ratio {ratio:.3f}) -- camera and 2P did not run "
                    f"frame-for-frame; fps may be off")
    return fps, note, warn, n_video


# =========================================================================
#  encoding
# =========================================================================

def encode(ffmpeg, src, out, fps, crf, preset, mjpg=False, dry=False):
    if mjpg:
        # MJPEG in AVI: what ImageJ's built-in AVI reader understands.
        vcodec = ["-c:v", "mjpeg", "-q:v", "3", "-pix_fmt", "yuvj420p"]
    else:
        vcodec = ["-c:v", "libx264", "-preset", preset, "-crf", str(crf),
                  "-pix_fmt", "yuv420p", "-movflags", "+faststart"]

    # -r BEFORE -i sets the INPUT frame rate. Combined with a re-encode this
    # is what actually re-times the video. (On a stream COPY it only relabels
    # the header and players ignore it -- that is why fix_avi_timing.py has to
    # use the `setts` bitstream filter instead.)
    cmd = ([ffmpeg, "-hide_banner", "-loglevel", "error", "-y",
            "-r", f"{fps:.6f}", "-i", src] + vcodec + ["-an", out])

    if dry:
        print(f"    [dry-run] would write {out}")
        return True

    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("    FFMPEG FAILED:\n" + (r.stderr or "")[-1500:])
        return False
    return True


# =========================================================================
#  driver
# =========================================================================

def find_videos(path):
    if os.path.isfile(path):
        return [path] if path.lower().endswith((".avi", ".mkv")) else []
    vids = sorted(glob.glob(os.path.join(path, "**", "*.avi"), recursive=True))
    # never pick up copies this script (or fix_avi_timing.py) already made
    return [v for v in vids
            if not any(t in os.path.basename(v).lower()
                       for t in ("_view", "_fixedfps"))]


def main():
    ap = argparse.ArgumentParser(
        description="Make viewable, correctly-timed copies of Basler videos.")
    ap.add_argument("path", help="parent dir (batch), a run folder, or a .avi")
    ap.add_argument("--fps", type=float, default=None)
    ap.add_argument("--crf", type=int, default=17)
    ap.add_argument("--preset", default="veryfast")
    ap.add_argument("--imagej", action="store_true",
                    help="also write an MJPEG .avi for ImageJ")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--suffix", default="_view")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    ffmpeg = find_ffmpeg()
    if ffmpeg is None and not a.dry_run:
        sys.exit("No working ffmpeg found. Install one, or run:\n"
                 "    pip install imageio-ffmpeg")
    print(f"ffmpeg: {ffmpeg}\n")

    vids = find_videos(a.path)
    if not vids:
        sys.exit(f"No source .avi found under: {a.path}")
    print(f"Found {len(vids)} video(s).\n")

    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)

    done = skipped = failed = 0
    warnings = []

    for i, src in enumerate(vids, 1):
        run = os.path.basename(os.path.dirname(src))
        stem = os.path.splitext(os.path.basename(src))[0]
        print(f"[{i}/{len(vids)}] {run}/{os.path.basename(src)}")

        fps, note, warn, n_video = resolve_fps(src, a.fps)
        if fps is None:
            print(f"    SKIP -- could not determine fps: {note}\n")
            skipped += 1
            continue

        dur = (n_video / fps) if n_video else float("nan")
        print(f"    fps = {fps:.4f}   [{note}]")
        print(f"    {n_video} frames -> {dur:.1f} s real time")
        if warn:
            print(f"    WARNING: {warn}")
            warnings.append(f"{run}: {warn}")

        outdir = a.outdir or os.path.dirname(src)
        # in a shared --outdir, prefix with the run name to keep names unique
        base = (f"{run}__{stem}" if a.outdir else stem) + a.suffix

        targets = [(os.path.join(outdir, base + ".mp4"), False)]
        if a.imagej:
            targets.append((os.path.join(outdir, base + "_mjpg.avi"), True))

        for out, is_mjpg in targets:
            if os.path.isfile(out) and not a.overwrite:
                print(f"    SKIP (exists): {os.path.basename(out)}")
                skipped += 1
                continue
            ok = encode(ffmpeg, src, out, fps, a.crf, a.preset,
                        mjpg=is_mjpg, dry=a.dry_run)
            if ok:
                done += 1
                if not a.dry_run:
                    mb_in = os.path.getsize(src) / 1e6
                    mb_out = os.path.getsize(out) / 1e6
                    print(f"    -> {os.path.basename(out)}"
                          f"   {mb_in:.0f} MB -> {mb_out:.1f} MB")
            else:
                failed += 1
        print()

    print(f"Done: {done} written, {skipped} skipped, {failed} failed.")
    if warnings:
        print("\nFrame-count mismatches (check these runs):")
        for w in warnings:
            print(f"  - {w}")


if __name__ == "__main__":
    main()
