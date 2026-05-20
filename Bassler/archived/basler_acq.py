"""Basler camera acquisition helper — called from MATLAB or command line.

Usage from MATLAB:
    system('python basler_acq.py preview --binH 2 --binV 2 --fps 90 --pixel Mono8')
    system('python basler_acq.py acquire --saveDir ... --fps 90 --dur 60 --period 300 --nRuns 60 ...')

Usage standalone:
    python basler_acq.py preview --binH 2 --binV 2 --fps 90
    python basler_acq.py acquire --saveDir C:/data --baseName basler --fps 90 --dur 60 --period 300 --nRuns 60
"""

import argparse
import sys
import time
import os
from datetime import datetime

import numpy as np
from pypylon import pylon


def open_camera(fps, pixel_format, bin_h, bin_v, roi=None):
    """Open and configure the first available Basler camera.
    roi: optional (x, y, w, h) tuple for sensor crop.
    """
    tlf = pylon.TlFactory.GetInstance()
    camera = pylon.InstantCamera(tlf.CreateFirstDevice())
    camera.Open()

    print(f"Connected: {camera.GetDeviceInfo().GetFriendlyName()}")

    camera.PixelFormat.Value = pixel_format
    camera.BinningHorizontal.Value = bin_h
    camera.BinningVertical.Value = bin_v

    # Apply ROI if specified
    if roi is not None:
        rx, ry, rw, rh = roi
        camera.OffsetX.Value = 0
        camera.OffsetY.Value = 0
        camera.Width.Value = rw
        camera.Height.Value = rh
        camera.OffsetX.Value = rx
        camera.OffsetY.Value = ry

    camera.AcquisitionFrameRateEnable.Value = True
    camera.AcquisitionFrameRate.Value = float(fps)

    # Cap exposure to frame period
    max_exp = int(1e6 / fps)
    if camera.ExposureTime.Value > max_exp:
        camera.ExposureTime.Value = float(max_exp)

    w, h = camera.Width.Value, camera.Height.Value
    print(f"Resolution: {w}x{h}  |  {fps} fps  |  ~{w*h*fps/1e6:.1f} MB/s")

    return camera


def cmd_preview(args):
    """Live preview — shows camera feed, press 'q' to quit."""
    camera = open_camera(args.fps, args.pixel, args.binH, args.binV)

    try:
        import cv2
        use_cv2 = True
    except ImportError:
        use_cv2 = False

    if not use_cv2:
        # Grab one frame and report — no GUI without cv2
        camera.StartGrabbing(pylon.GrabStrategy_LatestImageOnly)
        grab = camera.RetrieveResult(5000, pylon.TimeoutHandling_ThrowException)
        if grab.GrabSucceeded():
            img = grab.Array
            print(f"Preview frame OK: {img.shape}, dtype={img.dtype}, "
                  f"min={img.min()} max={img.max()}")
            print("Install opencv-python for live preview window: pip install opencv-python")
        grab.Release()
        camera.StopGrabbing()
        camera.Close()
        return

    # Live preview with OpenCV
    import cv2

    camera.StartGrabbing(pylon.GrabStrategy_LatestImageOnly)
    win = f"Basler Preview | {args.pixel} bin{args.binH}x{args.binV} {args.fps}fps"
    cv2.namedWindow(win, cv2.WINDOW_NORMAL)

    print("Preview running:")
    print("  'c' = draw crop ROI (drag rectangle, then Enter/Space to confirm, Esc to cancel)")
    print("  'r' = reset to full frame")
    print("  'q' = accept and proceed to acquisition")

    roi = None  # (x, y, w, h) on sensor

    while camera.IsGrabbing():
        grab = camera.RetrieveResult(1000, pylon.TimeoutHandling_Return)
        if grab is None or not grab.GrabSucceeded():
            continue
        frame = grab.Array.copy()
        grab.Release()

        # Draw ROI overlay if set
        display = cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR) if len(frame.shape) == 2 else frame.copy()
        if roi is not None:
            rx, ry, rw, rh = roi
            cv2.rectangle(display, (rx, ry), (rx+rw, ry+rh), (0, 255, 255), 2)
            label = f"ROI: {rw}x{rh} @ ({rx},{ry})"
            cv2.putText(display, label, (rx, max(ry-10, 20)),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 255), 2)

        h, w = frame.shape[:2]
        info = f"{w}x{h} | {args.fps}fps"
        if roi:
            info += f" | crop: {roi[2]}x{roi[3]}"
        cv2.putText(display, info, (10, 30),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 255, 255), 2)
        cv2.imshow(win, display)

        key = cv2.waitKey(30) & 0xFF
        if key == ord('q'):
            break
        elif key == ord('c'):
            # Pause grabbing, let user draw ROI
            print("  Draw crop rectangle, then press Enter/Space to confirm (Esc to cancel)...")
            rect = cv2.selectROI(win, frame, fromCenter=False, showCrosshair=True)
            if rect[2] > 0 and rect[3] > 0:
                # Snap to multiples of 2 (Basler sensor requirement)
                rx = int(rect[0]) & ~1
                ry = int(rect[1]) & ~1
                rw = int(rect[2]) & ~1
                rh = int(rect[3]) & ~1
                roi = (rx, ry, rw, rh)
                print(f"  ROI set: {rw}x{rh} at offset ({rx}, {ry})")
            else:
                print("  ROI selection cancelled.")
        elif key == ord('r'):
            roi = None
            print("  ROI reset to full frame.")
        if cv2.getWindowProperty(win, cv2.WND_PROP_VISIBLE) < 1:
            break

    camera.StopGrabbing()
    cv2.destroyAllWindows()

    # Apply ROI to camera if set
    if roi is not None:
        camera.StopGrabbing()
        rx, ry, rw, rh = roi
        camera.OffsetX.Value = 0  # reset first to avoid range errors
        camera.OffsetY.Value = 0
        camera.Width.Value = rw
        camera.Height.Value = rh
        camera.OffsetX.Value = rx
        camera.OffsetY.Value = ry
        w, h = camera.Width.Value, camera.Height.Value
        max_fps = camera.AcquisitionFrameRate.Max
        print(f"\nCrop applied: {w}x{h} at offset ({rx},{ry})")
        print(f"Max fps at this ROI: {max_fps:.0f}")
        print(f"Current fps: {args.fps}")
        # Output ROI as comma-separated for MATLAB to parse
        print(f"ROI={rx},{ry},{rw},{rh}")
    else:
        print(f"\nFull frame: {camera.Width.Value}x{camera.Height.Value}")
        print("ROI=none")

    camera.Close()
    print("Preview closed.")


def cmd_acquire(args):
    """Run N acquisitions, saving lossless AVI (one per run).
    Grabs frames into a queue, writes to disk in a background thread.
    """
    import cv2
    import threading
    import queue

    roi = getattr(args, '_roi', None)
    camera = open_camera(args.fps, args.pixel, args.binH, args.binV, roi=roi)

    os.makedirs(args.saveDir, exist_ok=True)
    frames_per_run = round(args.fps * args.dur)

    for run_idx in range(1, args.nRuns + 1):
        run_start = time.time()
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        fname = os.path.join(args.saveDir,
                             f"{args.baseName}_{ts}_run{run_idx:03d}.avi")

        # Background writer thread
        frame_q = queue.Queue(maxsize=500)
        write_done = threading.Event()
        written_count = [0]

        def disk_writer():
            writer = None
            while True:
                item = frame_q.get()
                if item is None:  # sentinel
                    break
                if writer is None:
                    h, w = item.shape[:2]
                    fourcc = cv2.VideoWriter_fourcc(*'FFV1')  # lossless
                    writer = cv2.VideoWriter(fname, fourcc, args.fps, (w, h), False)
                writer.write(item)
                written_count[0] += 1
            if writer is not None:
                writer.release()
            write_done.set()

        writer_thread = threading.Thread(target=disk_writer, daemon=True)
        writer_thread.start()

        # Live preview window (non-blocking)
        win_name = f"Acquiring run {run_idx}/{args.nRuns}"
        cv2.namedWindow(win_name, cv2.WINDOW_NORMAL)
        latest_frame = [None]

        camera.StartGrabbing(pylon.GrabStrategy_OneByOne)
        acq_start = time.time()
        print(f"[Run {run_idx}/{args.nRuns}] START  {fname}")

        grabbed = 0
        dropped = 0
        for fi in range(frames_per_run):
            grab = camera.RetrieveResult(5000, pylon.TimeoutHandling_ThrowException)
            if grab.GrabSucceeded():
                frame = grab.Array.copy()  # copy before release
                grabbed += 1
                latest_frame[0] = frame
                try:
                    frame_q.put_nowait(frame)
                except queue.Full:
                    dropped += 1  # queue full, drop frame
            grab.Release()

            # Update preview + progress every 0.5s worth of frames
            if fi % max(1, int(args.fps / 2)) == 0:
                if latest_frame[0] is not None:
                    cv2.imshow(win_name, latest_frame[0])
                    cv2.waitKey(1)
                elapsed_acq = time.time() - acq_start
                pct = (fi + 1) / frames_per_run * 100
                print(f"\r  {fi+1}/{frames_per_run} frames ({pct:.0f}%)  "
                      f"elapsed {elapsed_acq:.0f}s  queue={frame_q.qsize()}", end='', flush=True)

        camera.StopGrabbing()
        cv2.destroyWindow(win_name)
        print()  # newline after progress

        # Signal writer to finish and wait
        frame_q.put(None)
        print(f"  Flushing {frame_q.qsize()} queued frames to disk...", end='', flush=True)
        writer_thread.join(timeout=120)
        print(" done.")

        fsize_mb = os.path.getsize(fname) / 1e6
        elapsed_total = time.time() - run_start
        msg = f"[Run {run_idx}/{args.nRuns}] DONE   grabbed={grabbed} written={written_count[0]}"
        if dropped > 0:
            msg += f" DROPPED={dropped}"
        msg += f"  {fsize_mb:.0f} MB  {elapsed_total:.0f}s"
        print(msg)

        # Wait for next period
        if run_idx < args.nRuns:
            wait = args.period - elapsed_total
            if wait > 0:
                next_time = datetime.now().strftime("%H:%M:%S")
                print(f"  Waiting {wait:.0f}s until next run (next ~{next_time})...")
                # Show countdown in preview window during wait
                while wait > 0:
                    cv2.waitKey(1)
                    time.sleep(min(wait, 10))
                    wait = args.period - (time.time() - run_start)
            else:
                print(f"  WARNING: run took {elapsed_total:.1f}s, exceeds period {args.period}s")

    camera.Close()
    print(f"\n=== Completed {args.nRuns} runs ===")


def main():
    parser = argparse.ArgumentParser(description="Basler camera acquisition")
    parser.add_argument('--fps', type=float, default=90)
    parser.add_argument('--pixel', type=str, default='Mono8')
    parser.add_argument('--binH', type=int, default=1)
    parser.add_argument('--binV', type=int, default=1)

    sub = parser.add_subparsers(dest='cmd')

    sub.add_parser('preview')

    acq = sub.add_parser('acquire')
    acq.add_argument('--saveDir', required=True)
    acq.add_argument('--baseName', default='basler')
    acq.add_argument('--dur', type=float, default=60)
    acq.add_argument('--period', type=float, default=300)
    acq.add_argument('--nRuns', type=int, default=50)
    acq.add_argument('--roi', type=str, default=None,
                     help='x,y,w,h sensor crop (from preview)')

    args = parser.parse_args()

    # Parse ROI string into tuple
    roi = None
    if hasattr(args, 'roi') and args.roi and args.roi != 'none':
        roi = tuple(int(v) for v in args.roi.split(','))

    if args.cmd == 'preview':
        cmd_preview(args)
    elif args.cmd == 'acquire':
        args._roi = roi
        cmd_acquire(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
