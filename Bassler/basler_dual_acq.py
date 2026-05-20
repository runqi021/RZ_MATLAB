"""Dual Basler camera acquisition with hardware trigger support.

Designed for 2-photon synchronization:
  - Breathing camera: hardware-triggered, 1 frame per 2P frame (Line3 falling edge)
  - Whisking camera:  free-running at higher fps (e.g. 150 Hz)
Both save lossless FFV1 AVI + per-frame timestamp CSV for post-hoc alignment.

Commands:
  list     -- enumerate connected cameras (find serial numbers)
  preview  -- live preview with ROI selection (one camera at a time)
  acquire  -- multi-run dual-camera acquisition

Usage:
  python basler_dual_acq.py list
  python basler_dual_acq.py preview --serial <SN> --binH 3 --binV 3
  python basler_dual_acq.py acquire --saveDir <path> --dur 60 --period 300 --nRuns 60 \
      --cam1-serial <SN> --cam1-name breathing --cam1-trigger hardware --cam1-fps 30 \
      --cam1-binH 3 --cam1-binV 3 \
      --cam2-serial <SN> --cam2-name whisking --cam2-trigger freerun --cam2-fps 150 \
      --cam2-binH 1 --cam2-binV 1
"""

import argparse
import os
import sys
import time
import threading
import queue
from datetime import datetime

import numpy as np
from pypylon import pylon

try:
    import cv2
    HAS_CV2 = True
except ImportError:
    HAS_CV2 = False


# ── Helpers ──────────────────────────────────────────────────────────────

def parse_roi(s):
    """Parse 'x,y,w,h' string to (x, y, w, h) int tuple."""
    if not s or s == "none":
        return None
    return tuple(int(v) for v in s.split(","))


def list_cameras():
    """Enumerate and print all connected Basler cameras."""
    tlf = pylon.TlFactory.GetInstance()
    devices = tlf.EnumerateDevices()
    if not devices:
        print("No Basler cameras found.")
        return []
    print(f"Found {len(devices)} camera(s):")
    for i, d in enumerate(devices):
        print(f"  [{i}] {d.GetFriendlyName()}")
        print(f"       SN={d.GetSerialNumber()}  Model={d.GetModelName()}")
    return devices


def find_device(devices, serial):
    """Find DeviceInfo by serial number."""
    for d in devices:
        if d.GetSerialNumber() == serial:
            return d
    return None


def resolve_devices(devices, serials):
    """Resolve a list of serial numbers to DeviceInfo objects.

    *serials* is a list where each entry is either:
      - None  -> camera slot not requested (returns None)
      - ""    -> auto-assign from remaining devices (returns None if none left)
      - "AUTO" (case-insensitive) -> same as ""
      - an explicit serial number string

    Returns a list of DeviceInfo|None of the same length as *serials*.
    Exits with error if an explicit serial is not found or duplicates an
    already-claimed device.
    """
    resolved = [None] * len(serials)
    used = set()

    for i, sn in enumerate(serials):
        if sn is None:
            continue
        tag = f"cam{i+1}"
        if sn and sn.upper() != "AUTO":
            dev = find_device(devices, sn)
            if dev is None:
                print(f"ERROR: {tag} serial {sn} not found")
                sys.exit(1)
            if dev.GetSerialNumber() in used:
                print(f"ERROR: {tag} duplicates a serial already assigned to another slot")
                sys.exit(1)
            resolved[i] = dev
            used.add(dev.GetSerialNumber())
        else:
            # Auto-assign from remaining devices
            if not devices:
                if i == 0:
                    print("ERROR: No cameras found")
                    sys.exit(1)
                print(f"WARNING: no device available for {tag} -- slot unused")
                continue
            picked = None
            for d in devices:
                if d.GetSerialNumber() not in used:
                    picked = d
                    break
            if picked is None:
                print(f"WARNING: no free device for {tag} -- slot unused")
                continue
            resolved[i] = picked
            used.add(picked.GetSerialNumber())
            print(f"[{tag}] Auto-assigned: SN={picked.GetSerialNumber()}")

    return resolved


# ── Camera setup ─────────────────────────────────────────────────────────

def open_camera(device_info, name, trigger, fps, pixel, bin_h, bin_v,
                roi=None, burst_count=1, exposure_us=0, edge="falling"):
    """Open, configure, and return a pylon.InstantCamera.

    trigger:
      "hardware"  — FrameStart on Line3 falling edge, 1 frame per trigger
      "burst"     — FrameBurstStart on Line3 falling edge, *burst_count*
                    frames at *fps* Hz per trigger, then idle until next
      "freerun"   — internal clock at *fps* Hz (no trigger)
    """
    tlf = pylon.TlFactory.GetInstance()
    cam = pylon.InstantCamera(tlf.CreateDevice(device_info))
    cam.Open()

    sn = device_info.GetSerialNumber()
    print(f"[{name}] Opened: {device_info.GetFriendlyName()} (SN={sn})")

    # Size transport-layer buffer: up to 2 GB worth of frames, capped at 500.
    # Too many buffers crashes the USB transport on large-sensor cameras.
    frame_bytes = cam.Width.Value * cam.Height.Value  # after binning + ROI
    max_buf = max(200, min(500, int(2e9 / max(frame_bytes, 1))))
    try:
        cam.MaxNumBuffer = max_buf
    except Exception:
        pass

    cam.PixelFormat.Value = pixel
    cam.BinningHorizontal.Value = bin_h
    cam.BinningVertical.Value = bin_v

    # Configure Line3 as a digital input.  On Basler ace USB cameras Line3
    # is a bidirectional GPIO; it only receives triggers when LineMode=Input.
    # Use try/except to stay compatible with cameras that expose Line3 as
    # a fixed-direction opto-input (no LineMode setter).
    try:
        cam.LineSelector.Value = "Line3"
        try:
            cam.LineMode.Value = "Input"
        except Exception:
            pass  # fixed-direction line -- nothing to do
        cam.LineDebouncerTime.Value = 5.0  # µs; ignore contact bounce / ringing
    except Exception as e:
        print(f"[{name}] Warning: could not configure Line3: {e}")

    # Clear any stale TriggerMode on every selector before configuring the
    # one we want.  Without this, a prior session that left FrameBurstStart
    # on could starve a fresh FrameStart trigger.
    for sel in ("AcquisitionStart", "FrameBurstStart", "FrameStart"):
        try:
            cam.TriggerSelector.Value = sel
            cam.TriggerMode.Value = "Off"
        except Exception:
            pass

    # Sensor crop (snap to camera's actual increment requirements)
    # Always reset to full sensor first — Basler cameras persist ROI
    # across open/close cycles, so a stale crop from a previous session
    # would stick if we don't clear it.
    cam.OffsetX.Value = 0
    cam.OffsetY.Value = 0
    cam.Width.Value = cam.WidthMax.Value
    cam.Height.Value = cam.HeightMax.Value
    if roi:
        xi = cam.OffsetX.Inc
        yi = cam.OffsetY.Inc
        wi = cam.Width.Inc
        hi = cam.Height.Inc
        rx = (roi[0] // xi) * xi
        ry = (roi[1] // yi) * yi
        rw = (roi[2] // wi) * wi
        rh = (roi[3] // hi) * hi
        cam.Width.Value = rw
        cam.Height.Value = rh
        cam.OffsetX.Value = rx
        cam.OffsetY.Value = ry
        print(f"[{name}] ROI: {rw}x{rh} at ({rx},{ry})  "
              f"(inc: w={wi} h={hi} x={xi} y={yi})")

    activation = "RisingEdge" if edge == "rising" else "FallingEdge"

    if trigger == "hardware":
        # Each edge on Line3 -> one frame
        cam.TriggerSelector.Value = "FrameStart"
        cam.TriggerMode.Value = "On"
        cam.TriggerSource.Value = "Line3"
        cam.TriggerActivation.Value = activation
        cam.AcquisitionFrameRateEnable.Value = False
        print(f"[{name}] Trigger: Line3 {activation}  (1 frame / pulse)")

    elif trigger == "burst":
        # Matches verified Pylon Viewer settings for acA1300-200um:
        try:
            cam.SensorReadoutMode.Value = "Normal"
        except Exception:
            pass
        cam.TriggerSelector.Value = "FrameBurstStart"
        cam.TriggerMode.Value = "On"
        cam.TriggerSource.Value = "Line3"
        cam.TriggerActivation.Value = activation
        cam.TriggerSelector.Value = "FrameStart"
        cam.TriggerMode.Value = "Off"
        cam.AcquisitionBurstFrameCount.Value = burst_count
        cam.AcquisitionFrameRateEnable.Value = True
        cam.AcquisitionFrameRate.Value = float(fps)
        print(f"[{name}] Burst trigger: Line3 {activation}  "
              f"({burst_count} frames @ {fps} Hz per pulse)")

    else:  # freerun
        cam.TriggerMode.Value = "Off"
        cam.AcquisitionFrameRateEnable.Value = True
        cam.AcquisitionFrameRate.Value = float(fps)
        print(f"[{name}] Free-run: {fps} Hz")

    # Exposure: use explicit value if provided, otherwise auto-cap
    if exposure_us > 0:
        cam.ExposureTime.Value = float(exposure_us)
    elif fps > 0:
        max_exp = int(0.9e6 / fps)
        if cam.ExposureTime.Value > max_exp:
            cam.ExposureTime.Value = float(max_exp)

    w, h = cam.Width.Value, cam.Height.Value
    exp_ms = cam.ExposureTime.Value / 1000
    print(f"[{name}] {w}x{h}  bin={bin_h}x{bin_v}  exp={exp_ms:.1f} ms")

    # Report actual achievable frame rate vs requested
    try:
        actual_fps = cam.ResultingFrameRate.Value
        print(f"[{name}] Max achievable FPS: {actual_fps:.1f} Hz"
              f"  (requested: {fps} Hz)")
        if trigger == "burst" and actual_fps < fps * 0.95:
            burst_dur_ms = (burst_count - 1) / actual_fps * 1000
            print(f"[{name}] WARNING: camera caps at {actual_fps:.0f} Hz, "
                  f"burst of {burst_count} will take {burst_dur_ms:.1f} ms — "
                  f"need ROI crop or more binning for {fps} Hz")
    except Exception:
        pass

    # Diagnostic: report instantaneous Line3 state + active trigger selectors.
    # If Line3 reads the same value for 2+ seconds while 2P is triggering,
    # the problem is upstream (wiring / voltage / TTL polarity), not software.
    try:
        cam.LineSelector.Value = "Line3"
        lvl = cam.LineStatus.Value
        print(f"[{name}] Line3 input state: {'HIGH' if lvl else 'LOW'}")
    except Exception:
        pass
    try:
        active = []
        for sel in ("FrameStart", "FrameBurstStart", "AcquisitionStart"):
            cam.TriggerSelector.Value = sel
            if cam.TriggerMode.Value == "On":
                active.append(f"{sel}<-{cam.TriggerSource.Value}/"
                              f"{cam.TriggerActivation.Value}")
        # Leave TriggerSelector on the user-selected trigger type for runtime.
        if trigger == "hardware":
            cam.TriggerSelector.Value = "FrameStart"
        elif trigger == "burst":
            cam.TriggerSelector.Value = "FrameBurstStart"
        print(f"[{name}] Active triggers: {active if active else 'none'}")
    except Exception:
        pass
    return cam


# ── Per-camera session acquisition (thread target) ───────────────────────

def acquire_camera_session(cam, name, save_dir, session_s, fps_hint, trigger,
                           split_gap_s, single_run, results, latest_frame):
    """Continuous triggered acquisition for an entire session.

    Cameras stay armed for *session_s* wall-clock seconds.  Frames are
    recorded whenever triggers arrive.  When triggers pause for longer
    than *split_gap_s*, the current AVI is closed.

    If *single_run* is True, exits after the first completed recording
    (i.e. after the first trigger gap).  Otherwise starts a new file
    on the next trigger.

    Output per run:
        {name}_{YYYYMMDD_HHMMSS}_run{NNN}.avi
        {name}_{YYYYMMDD_HHMMSS}_run{NNN}_timestamps.csv

    Populates *results* dict with summary statistics.
    *latest_frame* is a single-element list written by this thread and
    read by the main-thread display loop (no lock needed under the GIL).
    """
    is_triggered = trigger in ("hardware", "burst")
    timeout_ms = 5000 if is_triggered else 2000
    last_print_t = 0.0

    cam.StartGrabbing(pylon.GrabStrategy_OneByOne)
    t_session = time.time()

    run_idx = 0
    total_grabbed = 0

    # Current run state
    run_dir = None
    run_avi_path = None
    run_timestamps = []     # [(frame_idx, cam_ns, wall_s), ...]
    run_grabbed = 0
    run_active = False
    last_frame_wall = None

    # Background writer thread — unlimited RAM queue, writes AVI on disk
    # Grab loop never blocks. Queue absorbs bursts, writer catches up.
    write_q = queue.Queue()
    vw = [None]              # [VideoWriter handle]

    def writer_thread():
        """Pulls frames from queue, writes to AVI."""
        while True:
            item = write_q.get()
            if item is None:      # sentinel → stop
                break
            if item == "CLOSE":   # close current AVI
                if vw[0] is not None:
                    vw[0].release()
                    vw[0] = None
                write_q.task_done()
                continue
            avi_path, frame_data = item
            if vw[0] is None:
                h, w = frame_data.shape[:2]
                vw[0] = cv2.VideoWriter(avi_path, cv2.VideoWriter_fourcc(*"FFV1"),
                                        float(fps_hint), (w, h), False)
                if not vw[0].isOpened():
                    raise RuntimeError(
                        f"cv2.VideoWriter failed to open {avi_path} with FFV1. "
                        "Verify your OpenCV build has ffmpeg + FFV1 support "
                        "(check cv2.getBuildInformation() for 'FFMPEG: YES')."
                    )
            vw[0].write(frame_data)
            write_q.task_done()

    writer = threading.Thread(target=writer_thread, daemon=True)
    writer.start()

    if is_triggered:
        print(f"[{name}] Armed — waiting for triggers on Line3 ...")

    def close_run():
        """Flush writer, save timestamps, print summary."""
        nonlocal run_timestamps, run_dir, run_avi_path, run_grabbed, run_active
        # Signal writer to close AVI, wait for all queued frames to flush
        write_q.put("CLOSE")
        write_q.join()
        if run_timestamps and run_dir:
            csv_path = os.path.join(run_dir, "timestamps.csv")
            with open(csv_path, "w") as f:
                f.write("frame_idx,camera_timestamp_ns,wall_time_s\n")
                for idx, cts, wts in run_timestamps:
                    f.write(f"{idx},{cts},{wts:.6f}\n")
            fsize = os.path.getsize(run_avi_path) / 1e6 if os.path.exists(run_avi_path) else 0
            actual_fps = 0.0
            if len(run_timestamps) >= 2:
                dt = run_timestamps[-1][2] - run_timestamps[0][2]
                if dt > 0:
                    actual_fps = (len(run_timestamps) - 1) / dt
            print(f"[{name}] Run {run_idx} saved: {run_grabbed} frames  "
                  f"{fsize:.0f} MB  {actual_fps:.1f} Hz")
        run_timestamps = []
        run_dir = None
        run_avi_path = None
        run_grabbed = 0
        run_active = False

    def open_run():
        """Create folder + AVI path for this run."""
        nonlocal run_idx, run_timestamps, run_dir, run_avi_path, run_grabbed, run_active
        run_idx += 1
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        cam_dir = os.path.join(save_dir, name)
        run_dir = os.path.join(cam_dir, f"{name}_{ts}_run{run_idx:03d}")
        os.makedirs(run_dir, exist_ok=True)
        run_avi_path = os.path.join(run_dir, f"{name}_{ts}_run{run_idx:03d}.avi")
        run_timestamps = []
        run_grabbed = 0
        run_active = True
        print(f"[{name}] Run {run_idx} started  ({name}/{os.path.basename(run_dir)}/)")

    # ── main grab loop (zero display overhead — main thread handles that) ──
    while True:
        now = time.time()
        if now - t_session >= session_s:
            break

        # Detect trigger gap → close run (and exit if single_run)
        if last_frame_wall is not None and run_active:
            if now - last_frame_wall > split_gap_s:
                close_run()
                if single_run:
                    break

        try:
            grab = cam.RetrieveResult(timeout_ms, pylon.TimeoutHandling_Return)
            if grab is None or not grab.GrabSucceeded():
                if grab is not None:
                    grab.Release()
                continue
            frame = grab.Array.copy()
            cam_ts = grab.TimeStamp
            grab.Release()
        except Exception:
            continue

        now = time.time()          # refresh after grab (may have blocked)
        last_frame_wall = now

        # Lazy open: start new folder on first frame / after gap
        if not run_active:
            open_run()

        # Queue frame for background AVI writing (never blocks grab loop)
        write_q.put((run_avi_path, frame))
        run_timestamps.append((run_grabbed, cam_ts, now))
        run_grabbed += 1
        total_grabbed += 1

        # Share latest frame with main-thread display (GIL-safe)
        latest_frame[0] = frame

        # Console progress (~1 Hz)
        if now - last_print_t >= 1.0:
            last_print_t = now
            print(f"[{name}]  run {run_idx}: {run_grabbed} frames  "
                  f"q={write_q.qsize()}  elapsed {now - t_session:.0f}/{session_s:.0f}s")

    # Close final run
    if run_active:
        close_run()

    # Stop writer thread
    write_q.put(None)
    writer.join(timeout=120)

    # Ask pylon how many frames it had to drop at the transport layer
    # (i.e. could not be delivered to our grab loop).  Non-zero here means
    # MaxNumBuffer was too small OR the USB link couldn't keep up.
    n_underrun = 0
    try:
        n_underrun = cam.NumReadyBuffers.Value  # still-queued in pylon, harmless
    except Exception:
        pass
    try:
        n_lost = cam.GetStreamGrabberParams().Statistic_Failed_Buffer_Count.Value
        if n_lost > 0:
            print(f"[{name}] WARNING: pylon dropped {n_lost} frame(s) "
                  f"(transport-layer; increase MaxNumBuffer if recurring)")
    except Exception:
        pass

    cam.StopGrabbing()

    results.update({
        "name": name,
        "total_grabbed": total_grabbed,
        "total_runs": run_idx,
    })
    print(f"[{name}] Session done: {run_idx} run(s), {total_grabbed} total frames")


# ── Per-camera N-runs free-run acquisition (thread target) ───────────────

def acquire_camera_nruns(cam, name, save_dir, dur_s, period_s, n_runs,
                         fps_hint, t_global_start, results, latest_frame):
    """Free-run periodic acquisition: dur_s record + idle, repeated n_runs times.

    All cameras share *t_global_start* so their run boundaries align across
    cameras (within thread-scheduling jitter, typically <50 ms).  Within
    each run the camera is StartGrabbing'd for *dur_s* wall-clock seconds,
    then StopGrabbing'd until the next scheduled start.

    Output per run mirrors session-mode layout:
        <save_dir>/<name>/<name>_<YYYYMMDD_HHMMSS>_run<NNN>/
            <name>_<YYYYMMDD_HHMMSS>_run<NNN>.avi
            timestamps.csv
    """
    # Background writer thread — same RAM-first pattern as session mode.
    # Unbounded queue: grab loop never blocks; writer drains in parallel.
    write_q = queue.Queue()
    vw = [None]

    def writer_thread():
        while True:
            item = write_q.get()
            if item is None:           # sentinel → terminate writer
                break
            if item == "CLOSE":        # close current AVI (between runs)
                if vw[0] is not None:
                    vw[0].release()
                    vw[0] = None
                write_q.task_done()
                continue
            avi_path, frame_data = item
            if vw[0] is None:
                h, w = frame_data.shape[:2]
                vw[0] = cv2.VideoWriter(avi_path, cv2.VideoWriter_fourcc(*"FFV1"),
                                        float(fps_hint), (w, h), False)
                if not vw[0].isOpened():
                    raise RuntimeError(
                        f"cv2.VideoWriter failed to open {avi_path} with FFV1. "
                        "Verify your OpenCV build has ffmpeg + FFV1 support "
                        "(check cv2.getBuildInformation() for 'FFMPEG: YES')."
                    )
            vw[0].write(frame_data)
            write_q.task_done()

    writer = threading.Thread(target=writer_thread, daemon=True)
    writer.start()

    total_grabbed = 0
    last_print_t = 0.0
    last_idle_print_t = 0.0

    # Keep camera grabbing continuously across the whole session — during
    # idle gaps we still drain frames (otherwise the pylon buffer overflows)
    # but only push them to the latest_frame preview, not to the AVI writer.
    # This lets the user monitor the animal between runs.
    cam.StartGrabbing(pylon.GrabStrategy_OneByOne)

    def drain_until(deadline_s, status_prefix):
        """Pull frames + update preview until wall-clock reaches deadline_s.
        Frames are NOT recorded — pure preview-feed during idle gaps."""
        nonlocal last_idle_print_t
        while True:
            remain = deadline_s - time.time()
            if remain <= 0:
                return
            try:
                grab = cam.RetrieveResult(min(500, max(10, int(remain * 1000))),
                                          pylon.TimeoutHandling_Return)
                if grab is None or not grab.GrabSucceeded():
                    if grab is not None:
                        grab.Release()
                    continue
                frame = grab.Array.copy()
                grab.Release()
                latest_frame[0] = frame
            except Exception:
                continue
            now = time.time()
            if now - last_idle_print_t >= 5.0:
                last_idle_print_t = now
                print(f"[{name}] {status_prefix}  preview live  "
                      f"({deadline_s - now:.0f}s remaining)")

    for run_idx in range(1, n_runs + 1):
        # Drain frames into the preview until this run's scheduled start
        run_start_target = t_global_start + (run_idx - 1) * period_s
        wait_s = run_start_target - time.time()
        if wait_s > 0:
            print(f"[{name}] Idle {wait_s:.0f}s until run {run_idx}/{n_runs} "
                  f"(preview stays on) ...")
            drain_until(run_start_target, f"idle before run {run_idx}/{n_runs}")
        elif run_idx > 1:
            print(f"[{name}] WARNING: run {run_idx} starting {-wait_s:.1f}s late "
                  f"(dur > period, or previous run flush slow)")

        # Open AVI path + folder for this run
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        cam_dir = os.path.join(save_dir, name)
        run_dir = os.path.join(cam_dir, f"{name}_{ts}_run{run_idx:03d}")
        os.makedirs(run_dir, exist_ok=True)
        avi_path = os.path.join(run_dir, f"{name}_{ts}_run{run_idx:03d}.avi")

        timestamps = []
        run_grabbed = 0

        run_t0 = time.time()
        print(f"[{name}] Run {run_idx}/{n_runs} START  "
              f"({name}/{os.path.basename(run_dir)}/)")

        while time.time() - run_t0 < dur_s:
            try:
                grab = cam.RetrieveResult(2000, pylon.TimeoutHandling_Return)
                if grab is None or not grab.GrabSucceeded():
                    if grab is not None:
                        grab.Release()
                    continue
                frame = grab.Array.copy()
                cam_ts = grab.TimeStamp
                grab.Release()
            except Exception:
                continue

            now = time.time()
            write_q.put((avi_path, frame))
            timestamps.append((run_grabbed, cam_ts, now))
            run_grabbed += 1
            total_grabbed += 1
            latest_frame[0] = frame

            if now - last_print_t >= 1.0:
                last_print_t = now
                print(f"[{name}]  run {run_idx}: {run_grabbed} frames  "
                      f"q={write_q.qsize()}  elapsed {now - run_t0:.0f}/{dur_s:.0f}s")

        # Flush any queued frames + close AVI handle BEFORE writing CSV,
        # so timestamps.csv reflects only frames actually committed to disk.
        # NOTE: we do NOT StopGrabbing — the camera keeps running so the
        # preview stays live during the idle gap before the next run.
        write_q.put("CLOSE")
        write_q.join()

        if timestamps:
            csv_path = os.path.join(run_dir, "timestamps.csv")
            with open(csv_path, "w") as f:
                f.write("frame_idx,camera_timestamp_ns,wall_time_s\n")
                for idx, cts, wts in timestamps:
                    f.write(f"{idx},{cts},{wts:.6f}\n")
            fsize = os.path.getsize(avi_path) / 1e6 if os.path.exists(avi_path) else 0
            actual_fps = 0.0
            if len(timestamps) >= 2:
                dt = timestamps[-1][2] - timestamps[0][2]
                if dt > 0:
                    actual_fps = (len(timestamps) - 1) / dt
            print(f"[{name}] Run {run_idx} saved: {run_grabbed} frames  "
                  f"{fsize:.0f} MB  {actual_fps:.1f} Hz")

    # Stop the continuous-grab session that was started before run 1
    try:
        cam.StopGrabbing()
    except Exception:
        pass

    # Report transport-layer drops at session end
    try:
        n_lost = cam.GetStreamGrabberParams().Statistic_Failed_Buffer_Count.Value
        if n_lost > 0:
            print(f"[{name}] WARNING: pylon dropped {n_lost} frame(s) "
                  f"(transport-layer; increase MaxNumBuffer if recurring)")
    except Exception:
        pass

    write_q.put(None)
    writer.join(timeout=120)

    results.update({
        "name": name,
        "total_grabbed": total_grabbed,
        "total_runs": n_runs,
    })
    print(f"[{name}] All runs done: {n_runs} runs, {total_grabbed} total frames")


# ── Commands ─────────────────────────────────────────────────────────────

def cmd_list(_args):
    list_cameras()


def cmd_preview(args):
    """Live preview for one camera -- identify it and optionally set ROI."""
    if not HAS_CV2:
        print("ERROR: opencv-python required for preview:  pip install opencv-python")
        sys.exit(1)

    devices = list_cameras()
    if not devices:
        return

    dev = None
    if args.serial and args.serial.upper() != "AUTO":
        dev = find_device(devices, args.serial)
        if dev is None:
            print(f"Camera SN={args.serial} not found.")
            return
    else:
        dev = devices[0]
        print(f"No serial specified -- using first camera: SN={dev.GetSerialNumber()}")

    cam = open_camera(dev, "preview", "freerun", args.fps, args.pixel,
                      args.binH, args.binV)

    cam.StartGrabbing(pylon.GrabStrategy_LatestImageOnly)
    sn = dev.GetSerialNumber()
    w_sensor, h_sensor = cam.Width.Value, cam.Height.Value
    win = f"Preview SN={sn} | {w_sensor}x{h_sensor} {args.pixel} bin{args.binH}x{args.binV}"
    cv2.namedWindow(win, cv2.WINDOW_NORMAL)

    # Scale window to fit screen while keeping correct aspect ratio
    max_dim = 800
    scale = min(max_dim / w_sensor, max_dim / h_sensor, 1.0)
    cv2.resizeWindow(win, int(w_sensor * scale), int(h_sensor * scale))

    print("  'c' = draw crop ROI    'r' = reset to full frame    'q' = done")

    roi = None

    while cam.IsGrabbing():
        try:
            grab = cam.RetrieveResult(1000, pylon.TimeoutHandling_Return)
            if grab is None or not grab.GrabSucceeded():
                continue
            frame = grab.Array.copy()
            grab.Release()
        except Exception:
            continue

        disp = cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR) if frame.ndim == 2 else frame.copy()
        if roi:
            rx, ry, rw, rh = roi
            cv2.rectangle(disp, (rx, ry), (rx + rw, ry + rh), (0, 255, 255), 2)
            cv2.putText(disp, f"ROI {rw}x{rh}@({rx},{ry})", (rx, max(ry - 10, 20)),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 255), 2)

        h, w = frame.shape[:2]
        cv2.putText(disp, f"SN={sn}  {w}x{h}  {args.fps}fps",
                    (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 255, 255), 2)
        cv2.imshow(win, disp)

        key = cv2.waitKey(30) & 0xFF
        if key == ord("q"):
            break
        elif key == ord("c"):
            print("  Draw rectangle, then Enter/Space to confirm (Esc to cancel) ...")
            rect = cv2.selectROI(win, frame, fromCenter=False, showCrosshair=True)
            # Reset window zoom (selectROI can leave it zoomed to selection)
            h, w = frame.shape[:2]
            cv2.resizeWindow(win, min(w, 1280), min(h, 960))
            if rect[2] > 0 and rect[3] > 0:
                xi = cam.OffsetX.Inc
                yi = cam.OffsetY.Inc
                wi = cam.Width.Inc
                hi = cam.Height.Inc
                roi = ((int(rect[0]) // xi) * xi, (int(rect[1]) // yi) * yi,
                       (int(rect[2]) // wi) * wi, (int(rect[3]) // hi) * hi)
                print(f"  ROI set: {roi[2]}x{roi[3]} at ({roi[0]},{roi[1]})")
                print(f"  Press 'c' to redraw, 'r' to clear, 'q' to accept")
            else:
                print("  ROI cancelled.")
        elif key == ord("r"):
            roi = None
            print("  ROI reset.")
        if cv2.getWindowProperty(win, cv2.WND_PROP_VISIBLE) < 1:
            break

    cam.StopGrabbing()
    cv2.destroyAllWindows()

    if roi:
        print(f"ROI={roi[0]},{roi[1]},{roi[2]},{roi[3]}")
    else:
        print("ROI=none")
    cam.Close()
    print("Preview closed.")


def cmd_acquire(args):
    """Session-based dual-camera acquisition.

    Cameras stay armed for --session seconds.  Files are auto-split
    when triggers pause for > --split-gap seconds (i.e. between 2P runs).
    """
    if not HAS_CV2:
        print("ERROR: opencv-python required:  pip install opencv-python")
        sys.exit(1)

    # Write PID file so MATLAB GUI can monitor / kill this process
    pidfile = os.path.join(os.environ.get('TEMP', '.'), 'basler_acq.pid')
    with open(pidfile, 'w') as f:
        f.write(str(os.getpid()))

    devices = list_cameras()
    if not devices:
        sys.exit(1)

    # Build per-slot config list from args (cam1 is required; cam2/cam3
    # opt-in when their --camN-serial was supplied).
    cam_cfgs = [{
        "serial":   args.cam1_serial,
        "name":     args.cam1_name,
        "trigger":  args.cam1_trigger,
        "fps":      args.cam1_fps,
        "pixel":    args.cam1_pixel,
        "binH":     args.cam1_binH,
        "binV":     args.cam1_binV,
        "roi":      args.cam1_roi,
        "edge":     args.cam1_edge,
        "burst":    args.cam1_burst,
        "exposure": args.cam1_exposure,
    }]
    if args.cam2_serial is not None:
        cam_cfgs.append({
            "serial":   args.cam2_serial,
            "name":     args.cam2_name,
            "trigger":  args.cam2_trigger,
            "fps":      args.cam2_fps,
            "pixel":    args.cam2_pixel,
            "binH":     args.cam2_binH,
            "binV":     args.cam2_binV,
            "roi":      args.cam2_roi,
            "edge":     args.cam2_edge,
            "burst":    args.cam2_burst,
            "exposure": args.cam2_exposure,
        })
    if args.cam3_serial is not None:
        cam_cfgs.append({
            "serial":   args.cam3_serial,
            "name":     args.cam3_name,
            "trigger":  args.cam3_trigger,
            "fps":      args.cam3_fps,
            "pixel":    args.cam3_pixel,
            "binH":     args.cam3_binH,
            "binV":     args.cam3_binV,
            "roi":      args.cam3_roi,
            "edge":     args.cam3_edge,
            "burst":    args.cam3_burst,
            "exposure": args.cam3_exposure,
        })

    devs = resolve_devices(devices, [c["serial"] for c in cam_cfgs])

    # Filter to slots that resolved to a real device
    active = [(dev, cfg) for dev, cfg in zip(devs, cam_cfgs) if dev is not None]
    if not active:
        print("ERROR: no cameras resolved")
        sys.exit(1)

    cams = []
    for dev, cfg in active:
        roi = parse_roi(cfg["roi"])
        cams.append(open_camera(dev, cfg["name"], cfg["trigger"],
                                cfg["fps"], cfg["pixel"],
                                cfg["binH"], cfg["binV"], roi,
                                burst_count=cfg["burst"],
                                exposure_us=cfg["exposure"],
                                edge=cfg["edge"]))

    os.makedirs(args.saveDir, exist_ok=True)

    mode = "single-run" if args.single_run else f"session {args.session/3600:.1f}h"
    print(f"\n{'=' * 60}")
    print(f"  {len(cams)} camera(s)  |  {mode}  |  "
          f"trigger gap = {args.split_gap}s")
    print(f"{'=' * 60}")

    # Launch one thread per camera — grab loops have zero display overhead
    PREVIEW_MAX_DIM = 400

    threads = []
    displays = []
    results = []
    for cam, (_, cfg) in zip(cams, active):
        holder = [None]
        res = {}
        th = threading.Thread(
            target=acquire_camera_session,
            args=(cam, cfg["name"], args.saveDir, args.session,
                  cfg["fps"], cfg["trigger"], args.split_gap,
                  args.single_run, res, holder))
        threads.append(th)
        displays.append((cfg["name"], holder))
        results.append(res)

    # Create preview windows on the main thread (required by OpenCV)
    for cam_name, _ in displays:
        cv2.namedWindow(cam_name, cv2.WINDOW_NORMAL)

    # Track which windows have been sized to their first frame — WINDOW_NORMAL
    # keeps whatever size the window was created at and stretches frames to
    # fit, so we need to resize once the true frame dimensions are known.
    window_sized = {name: False for name, _ in displays}

    for t in threads:
        t.start()

    # Main-thread display loop — ~60 Hz, no impact on grab threads
    while any(t.is_alive() for t in threads):
        for cam_name, holder in displays:
            f = holder[0]
            if f is not None:
                h, w = f.shape[:2]
                if max(h, w) > PREVIEW_MAX_DIM:
                    s = PREVIEW_MAX_DIM / max(h, w)
                    f = cv2.resize(f, (int(w * s), int(h * s)),
                                   interpolation=cv2.INTER_NEAREST)
                if not window_sized[cam_name]:
                    dh, dw = f.shape[:2]
                    cv2.resizeWindow(cam_name, dw, dh)
                    window_sized[cam_name] = True
                cv2.imshow(cam_name, f)
        cv2.waitKey(16)  # ~60 Hz

    for t in threads:
        t.join(timeout=10)

    cv2.destroyAllWindows()
    for cam in cams:
        try:
            cam.Close()
        except Exception:
            pass
    print(f"\n=== Session complete ===")

    # Clean up PID file
    try:
        os.remove(pidfile)
    except Exception:
        pass


def cmd_acquire_nruns(args):
    """N-runs free-run periodic acquisition (untriggered).

    All cameras run free-run at their requested fps.  Each camera grabs
    for *--dur* seconds, then idles until the next *--period* boundary,
    repeated *--nRuns* times.  Run boundaries are aligned across cameras
    via a shared wall-clock schedule.
    """
    if not HAS_CV2:
        print("ERROR: opencv-python required:  pip install opencv-python")
        sys.exit(1)

    pidfile = os.path.join(os.environ.get('TEMP', '.'), 'basler_acq.pid')
    with open(pidfile, 'w') as f:
        f.write(str(os.getpid()))

    devices = list_cameras()
    if not devices:
        sys.exit(1)

    # Build per-slot config.  Trigger is hardcoded to "freerun" — this
    # subcommand exists specifically for the untriggered N-runs use case;
    # use `acquire` if you need hardware/burst triggers.
    cam_cfgs = [{
        "serial":   args.cam1_serial,
        "name":     args.cam1_name,
        "fps":      args.cam1_fps,
        "pixel":    args.cam1_pixel,
        "binH":     args.cam1_binH,
        "binV":     args.cam1_binV,
        "roi":      args.cam1_roi,
        "exposure": args.cam1_exposure,
    }]
    if args.cam2_serial is not None:
        cam_cfgs.append({
            "serial":   args.cam2_serial,
            "name":     args.cam2_name,
            "fps":      args.cam2_fps,
            "pixel":    args.cam2_pixel,
            "binH":     args.cam2_binH,
            "binV":     args.cam2_binV,
            "roi":      args.cam2_roi,
            "exposure": args.cam2_exposure,
        })
    if args.cam3_serial is not None:
        cam_cfgs.append({
            "serial":   args.cam3_serial,
            "name":     args.cam3_name,
            "fps":      args.cam3_fps,
            "pixel":    args.cam3_pixel,
            "binH":     args.cam3_binH,
            "binV":     args.cam3_binV,
            "roi":      args.cam3_roi,
            "exposure": args.cam3_exposure,
        })

    devs = resolve_devices(devices, [c["serial"] for c in cam_cfgs])

    active = [(dev, cfg) for dev, cfg in zip(devs, cam_cfgs) if dev is not None]
    if not active:
        print("ERROR: no cameras resolved")
        sys.exit(1)

    cams = []
    for dev, cfg in active:
        roi = parse_roi(cfg["roi"])
        cams.append(open_camera(dev, cfg["name"], "freerun",
                                cfg["fps"], cfg["pixel"],
                                cfg["binH"], cfg["binV"], roi,
                                exposure_us=cfg["exposure"]))

    os.makedirs(args.saveDir, exist_ok=True)

    total_session_s = args.period * (args.nRuns - 1) + args.dur
    print(f"\n{'=' * 60}")
    print(f"  {len(cams)} camera(s)  |  N-runs free-run mode")
    print(f"  dur={args.dur}s  period={args.period}s  nRuns={args.nRuns}")
    print(f"  total session ~{total_session_s:.0f}s "
          f"({total_session_s/60:.1f} min)")
    print(f"{'=' * 60}")

    PREVIEW_MAX_DIM = 400

    # Shared run-start schedule: 2 s offset gives all camera threads time
    # to spin up before run 1 begins, so all cameras see run 1's StartGrab
    # at roughly the same wall time.
    t_global_start = time.time() + 2.0

    threads = []
    displays = []
    results = []
    for cam, (_, cfg) in zip(cams, active):
        holder = [None]
        res = {}
        th = threading.Thread(
            target=acquire_camera_nruns,
            args=(cam, cfg["name"], args.saveDir, args.dur, args.period,
                  args.nRuns, cfg["fps"], t_global_start, res, holder))
        threads.append(th)
        displays.append((cfg["name"], holder))
        results.append(res)

    for cam_name, _ in displays:
        cv2.namedWindow(cam_name, cv2.WINDOW_NORMAL)
    window_sized = {name: False for name, _ in displays}

    for t in threads:
        t.start()

    # Main-thread display loop (~60 Hz) — required by OpenCV for cv2.imshow
    while any(t.is_alive() for t in threads):
        for cam_name, holder in displays:
            f = holder[0]
            if f is not None:
                h, w = f.shape[:2]
                if max(h, w) > PREVIEW_MAX_DIM:
                    s = PREVIEW_MAX_DIM / max(h, w)
                    f = cv2.resize(f, (int(w * s), int(h * s)),
                                   interpolation=cv2.INTER_NEAREST)
                if not window_sized[cam_name]:
                    dh, dw = f.shape[:2]
                    cv2.resizeWindow(cam_name, dw, dh)
                    window_sized[cam_name] = True
                cv2.imshow(cam_name, f)
        cv2.waitKey(16)

    for t in threads:
        t.join(timeout=10)

    cv2.destroyAllWindows()
    for cam in cams:
        try:
            cam.Close()
        except Exception:
            pass
    print(f"\n=== N-runs acquisition complete ===")

    try:
        os.remove(pidfile)
    except Exception:
        pass


# ── CLI ──────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Dual Basler camera acquisition with hardware trigger")
    sub = p.add_subparsers(dest="cmd")

    # list
    sub.add_parser("list", help="Enumerate connected cameras")

    # preview
    pv = sub.add_parser("preview", help="Live preview for one camera")
    pv.add_argument("--serial",  default="")
    pv.add_argument("--fps",     type=float, default=30)
    pv.add_argument("--pixel",   default="Mono8")
    pv.add_argument("--binH",    type=int, default=1)
    pv.add_argument("--binV",    type=int, default=1)

    # acquire
    a = sub.add_parser("acquire", help="Dual-camera acquisition")
    a.add_argument("--saveDir",    required=True)
    a.add_argument("--session",    type=float, default=4*3600,
                   help="total session time in seconds (default: 4 hours)")
    a.add_argument("--split-gap",  type=float, default=10,
                   help="seconds of no triggers before closing file (default: 10)")
    a.add_argument("--single-run", action="store_true",
                   help="exit after first completed recording (one trigger gap)")

    # cam 1 -- breathing (hardware-triggered by default)
    a.add_argument("--cam1-serial",  default="AUTO")
    a.add_argument("--cam1-name",    default="breathing")
    a.add_argument("--cam1-trigger", default="hardware",
                   choices=["hardware", "burst", "freerun"])
    a.add_argument("--cam1-fps",     type=float, default=30,
                   help="expected trigger rate (exposure cap + AVI metadata)")
    a.add_argument("--cam1-binH",    type=int,   default=3)
    a.add_argument("--cam1-binV",    type=int,   default=3)
    a.add_argument("--cam1-pixel",   default="Mono8")
    a.add_argument("--cam1-roi",     default=None, help="x,y,w,h from preview")
    a.add_argument("--cam1-edge",    default="falling",
                   choices=["falling", "rising"],
                   help="trigger edge polarity (default: falling)")
    a.add_argument("--cam1-burst",   type=int,   default=1,
                   help="frames per trigger (only used when trigger=burst)")
    a.add_argument("--cam1-exposure", type=float, default=0,
                   help="exposure time in microseconds (0 = auto)")

    # cam 2 -- whisking (burst-triggered by default)
    a.add_argument("--cam2-serial",  default=None,
                   help="omit to run single-camera; 'AUTO' for auto-assign")
    a.add_argument("--cam2-name",    default="whisking")
    a.add_argument("--cam2-trigger", default="burst",
                   choices=["hardware", "burst", "freerun"])
    a.add_argument("--cam2-fps",     type=float, default=150,
                   help="internal burst frame rate (Hz)")
    a.add_argument("--cam2-binH",    type=int,   default=1)
    a.add_argument("--cam2-binV",    type=int,   default=1)
    a.add_argument("--cam2-pixel",   default="Mono8")
    a.add_argument("--cam2-roi",     default=None, help="x,y,w,h from preview")
    a.add_argument("--cam2-edge",    default="falling",
                   choices=["falling", "rising"],
                   help="trigger edge polarity (default: falling)")
    a.add_argument("--cam2-burst",   type=int,   default=5,
                   help="frames per trigger (only used when trigger=burst)")
    a.add_argument("--cam2-exposure", type=float, default=0,
                   help="exposure time in microseconds (0 = auto)")

    # cam 3 -- optional third camera (same pattern as cam2)
    a.add_argument("--cam3-serial",  default=None,
                   help="omit to skip cam3; 'AUTO' for auto-assign")
    a.add_argument("--cam3-name",    default="cam3")
    a.add_argument("--cam3-trigger", default="hardware",
                   choices=["hardware", "burst", "freerun"])
    a.add_argument("--cam3-fps",     type=float, default=30)
    a.add_argument("--cam3-binH",    type=int,   default=1)
    a.add_argument("--cam3-binV",    type=int,   default=1)
    a.add_argument("--cam3-pixel",   default="Mono8")
    a.add_argument("--cam3-roi",     default=None, help="x,y,w,h from preview")
    a.add_argument("--cam3-edge",    default="falling",
                   choices=["falling", "rising"],
                   help="trigger edge polarity (default: falling)")
    a.add_argument("--cam3-burst",   type=int,   default=1,
                   help="frames per trigger (only used when trigger=burst)")
    a.add_argument("--cam3-exposure", type=float, default=0,
                   help="exposure time in microseconds (0 = auto)")

    # acquire-nruns -- free-run periodic acquisition (untriggered)
    ar = sub.add_parser("acquire-nruns",
                        help="N-runs free-run periodic acquisition (untriggered)")
    ar.add_argument("--saveDir", required=True)
    ar.add_argument("--dur",     type=float, default=60,
                    help="seconds of recording per run (default: 60)")
    ar.add_argument("--period",  type=float, default=300,
                    help="seconds between RUN STARTS (default: 300)")
    ar.add_argument("--nRuns",   type=int,   default=10,
                    help="total number of runs (default: 10)")

    # cam 1
    ar.add_argument("--cam1-serial",   default="AUTO",
                    help="serial number, or 'AUTO' for first available")
    ar.add_argument("--cam1-name",     default="cam1")
    ar.add_argument("--cam1-fps",      type=float, default=90)
    ar.add_argument("--cam1-pixel",    default="Mono8")
    ar.add_argument("--cam1-binH",     type=int,   default=1)
    ar.add_argument("--cam1-binV",     type=int,   default=1)
    ar.add_argument("--cam1-roi",      default=None, help="x,y,w,h from preview")
    ar.add_argument("--cam1-exposure", type=float, default=0,
                    help="exposure time in microseconds (0 = auto)")

    # cam 2 (omit --cam2-serial to skip)
    ar.add_argument("--cam2-serial",   default=None,
                    help="omit to skip cam2; 'AUTO' for auto-assign")
    ar.add_argument("--cam2-name",     default="cam2")
    ar.add_argument("--cam2-fps",      type=float, default=90)
    ar.add_argument("--cam2-pixel",    default="Mono8")
    ar.add_argument("--cam2-binH",     type=int,   default=1)
    ar.add_argument("--cam2-binV",     type=int,   default=1)
    ar.add_argument("--cam2-roi",      default=None)
    ar.add_argument("--cam2-exposure", type=float, default=0)

    # cam 3 (omit --cam3-serial to skip)
    ar.add_argument("--cam3-serial",   default=None,
                    help="omit to skip cam3; 'AUTO' for auto-assign")
    ar.add_argument("--cam3-name",     default="cam3")
    ar.add_argument("--cam3-fps",      type=float, default=90)
    ar.add_argument("--cam3-pixel",    default="Mono8")
    ar.add_argument("--cam3-binH",     type=int,   default=1)
    ar.add_argument("--cam3-binV",     type=int,   default=1)
    ar.add_argument("--cam3-roi",      default=None)
    ar.add_argument("--cam3-exposure", type=float, default=0)

    args = p.parse_args()

    if args.cmd == "list":
        cmd_list(args)
    elif args.cmd == "preview":
        cmd_preview(args)
    elif args.cmd == "acquire":
        cmd_acquire(args)
    elif args.cmd == "acquire-nruns":
        cmd_acquire_nruns(args)
    else:
        p.print_help()


if __name__ == "__main__":
    main()
