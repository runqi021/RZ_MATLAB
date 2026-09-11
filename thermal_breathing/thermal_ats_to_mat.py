#!/usr/bin/env python
"""
thermal_ats_to_mat.py  --  FLIR .ats radiometric reader -> compact .mat for MATLAB

Reads a FLIR Science File SDK (.ats) thermal sequence as TRUE temperature (deg C),
derives the real frame rate from per-frame timestamps (never hardcoded), anti-alias
decimates to a breathing-appropriate rate, and saves a compact [T x H x W] stack as
a .mat that thermal_breath_svd.m consumes.

Run inside the conda env that has the FileSDK wheel installed:
    C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe thermal_ats_to_mat.py <in.ats> [--target-fps 50] [--out out.mat]

Companion: thermal_breath_svd.m  (does the SVD breath readout + QC figure in MATLAB).
breath_svd_pc1.m is a separate script and is NOT touched.
"""
import argparse, os, sys
import numpy as np
import fnv, fnv.file
from scipy.signal import decimate
from scipy.interpolate import interp1d
from scipy.io import savemat


def read_ats(path):
    """Read an .ats as true temperature. Returns (data[N,P] f32, H, W, times[N] s)."""
    im = fnv.file.ImagerFile(path)
    # radiometric temperature in Celsius (no AGC / no gamma -- the whole point)
    im.unit = fnv.Unit.TEMPERATURE_FACTORY
    im.temp_type = fnv.TempType.CELSIUS

    H, W = im.height, im.width
    N = im.num_frames
    P = H * W

    data = np.empty((N, P), dtype=np.float32)
    times = np.empty(N, dtype=np.float64)  # seconds relative to frame 0
    t0 = None
    for i in range(N):
        im.get_frame(i)
        data[i, :] = np.asarray(im.final, dtype=np.float32)
        ti = im.frame_info.time
        if t0 is None:
            t0 = ti
        times[i] = (ti - t0).total_seconds()
        if (i % 2000) == 0:
            print(f"  read {i}/{N}", flush=True)
    return data, H, W, times


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ats", help="input .ats file")
    ap.add_argument("--target-fps", type=float, default=50.0,
                    help="decimate down to ~this rate (default 50 Hz; Nyquist 25 Hz >> breathing/sniffing)")
    ap.add_argument("--out", default=None, help="output .mat (default: alongside input)")
    args = ap.parse_args()

    if not os.path.isfile(args.ats):
        sys.exit(f"no such file: {args.ats}")
    out = args.out or os.path.splitext(args.ats)[0] + "_thermal.mat"

    print(f"reading {args.ats}", flush=True)
    data, H, W, t = read_ats(args.ats)
    N = data.shape[0]

    # --- TRUE time axis from the timestamps. The camera nominally runs at a
    # fixed rate but DROPS frames (hundreds to thousands per 60 s file), so the
    # frames present are NOT uniformly spaced. We must put them back on real
    # wall-clock time before any spectral analysis. ---
    dt = np.diff(t)
    fps_nominal = float(round(1.0 / np.median(dt)))      # intended camera rate (e.g. 400)
    fps_mean = (N - 1) / (t[-1] - t[0])                  # actual achieved mean rate
    n_expected = int(round(t[-1] * fps_nominal)) + 1     # frames if none were dropped
    n_dropped = n_expected - N
    print(f"  {N} frames read, {H}x{W} px", flush=True)
    print(f"  nominal {fps_nominal:.0f} Hz, mean {fps_mean:.2f} Hz, "
          f"expected {n_expected}, DROPPED {n_dropped} "
          f"({100*n_dropped/n_expected:.1f}%)", flush=True)

    # resample onto a uniform grid at the nominal rate using the real times
    # (dropped-frame gaps filled by linear interpolation in time)
    t_uni = np.arange(n_expected, dtype=np.float64) / fps_nominal
    f_interp = interp1d(t, data, axis=0, kind='linear',
                        bounds_error=False, fill_value='extrapolate',
                        assume_sorted=True)
    data = f_interp(t_uni).astype(np.float32)
    del f_interp
    print(f"  resampled to uniform {fps_nominal:.0f} Hz: {data.shape[0]} frames", flush=True)

    # anti-alias temporal decimation to ~target fps (input is now truly uniform)
    q = max(1, int(round(fps_nominal / args.target_fps)))
    if q > 1:
        # decimate works along an axis with an FIR zero-phase filter; q<=13 per stage
        stages = []
        qq = q
        while qq > 13:
            stages.append(2)
            qq = int(round(qq / 2))
        stages.append(qq)
        dec = data
        for s in stages:
            if s > 1:
                dec = decimate(dec, s, axis=0, ftype='fir', zero_phase=True).astype(np.float32)
        fps_dec = fps_nominal / int(np.prod(stages))
        t_dec = np.arange(dec.shape[0], dtype=np.float64) / fps_dec
    else:
        dec = data
        t_dec = t_uni
        fps_dec = fps_nominal
    print(f"  decimated factor {q}: {dec.shape[0]} frames @ {fps_dec:.3f} fps", flush=True)

    # reshape to [T x H x W] for MATLAB
    Tn = dec.shape[0]
    stack = dec.reshape(Tn, H, W)
    savemat(out, {
        "stack": stack,                 # [T x H x W] single, deg C, uniform time
        "fps": float(fps_dec),          # decimated frame rate (Hz), exact/uniform
        "fps_nominal": float(fps_nominal),
        "fps_mean": float(fps_mean),    # achieved mean rate before resampling
        "decim_factor": int(q),
        "n_frames_read": int(N),
        "n_expected": int(n_expected),
        "n_dropped": int(n_dropped),
        "H": int(H), "W": int(W), "T": int(Tn),
        "t": t_dec.astype(np.float64),  # uniform frame times (s)
        "src": os.path.abspath(args.ats),
    }, do_compression=True)
    sz = os.path.getsize(out) / 1e6
    print(f"saved {out}  ({sz:.1f} MB)", flush=True)


if __name__ == "__main__":
    main()
