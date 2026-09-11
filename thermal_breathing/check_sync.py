#!/usr/bin/env python
"""
check_sync.py  --  evaluate whether a FLIR .ats and a Basler run are synchronized.

Compares the two cameras' OWN clocks:
  - Basler : timestamps.csv  (camera_timestamp_ns = hardware clock; wall_time_s = PC epoch)
  - FLIR   : per-frame frame_info.time (IRIG time-of-day, UTC; fnv stamps a 1976 placeholder year)

Reports frame rate, duration, dropped frames for each, and the absolute start offset
(via Basler PC-epoch UTC vs FLIR IRIG UTC time-of-day). Run in the flir env.
"""
import sys, csv, datetime as dt
import numpy as np
import fnv, fnv.file

ats_path = sys.argv[1]
csv_path = sys.argv[2]

# ---------------- Basler ----------------
fidx, cam_ns, wall_s = [], [], []
with open(csv_path) as f:
    r = csv.DictReader(f)
    for row in r:
        fidx.append(int(row['frame_idx']))
        cam_ns.append(int(row['camera_timestamp_ns']))
        wall_s.append(float(row['wall_time_s']))
cam_ns = np.array(cam_ns, dtype=np.int64)
wall_s = np.array(wall_s, dtype=np.float64)
Nb = len(cam_ns)
cam_s = (cam_ns - cam_ns[0]) / 1e9            # camera clock, seconds from frame 0
dtb = np.diff(cam_s)
fps_b = (Nb - 1) / cam_s[-1]
med_dtb = np.median(dtb)
nom_b = round(1.0 / med_dtb)
ndrop_b = int(round(cam_s[-1] * nom_b)) + 1 - Nb
b_start_utc = dt.datetime.utcfromtimestamp(wall_s[0])
b_end_utc   = dt.datetime.utcfromtimestamp(wall_s[-1])
print("=== BASLER ===")
print(f"  frames        : {Nb}")
print(f"  fps (cam clk) : {fps_b:.3f}  (nominal {nom_b}, median dt {med_dtb*1e3:.3f} ms)")
print(f"  duration      : {cam_s[-1]:.4f} s")
print(f"  dropped       : {ndrop_b}  (gaps>2x median: {(dtb>2*med_dtb).sum()})")
print(f"  start wall UTC: {b_start_utc.isoformat()}  (epoch {wall_s[0]:.6f})")
print(f"  end   wall UTC: {b_end_utc.isoformat()}")

# ---------------- FLIR thermal ----------------
im = fnv.file.ImagerFile(ats_path)
Nt = im.num_frames
tt = np.empty(Nt, dtype=np.float64)
t0 = None
times_dt = []
for i in range(Nt):
    im.get_frame(i)                  # timestamps only (do NOT touch im.final -> fast)
    ti = im.frame_info.time
    if t0 is None:
        t0 = ti
        first_dt = ti
    if i == Nt - 1:
        last_dt = ti
    tt[i] = (ti - t0).total_seconds()
dtt = np.diff(tt)
fps_t = (Nt - 1) / tt[-1]
med_dtt = np.median(dtt)
nom_t = round(1.0 / med_dtt)
ndrop_t = int(round(tt[-1] * nom_t)) + 1 - Nt
print("=== FLIR THERMAL ===")
print(f"  frames        : {Nt}")
print(f"  fps (IRIG)    : {fps_t:.3f}  (nominal {nom_t}, median dt {med_dtt*1e3:.3f} ms)")
print(f"  duration      : {tt[-1]:.4f} s")
print(f"  dropped       : {ndrop_t}  (gaps>2x median: {(dtt>2*med_dtt).sum()})")
print(f"  start IRIG    : {first_dt.isoformat()}  (fnv 1976 placeholder year)")
print(f"  end   IRIG    : {last_dt.isoformat()}")

# ---------------- alignment ----------------
# FLIR IRIG is UTC time-of-day; rebuild absolute UTC on the real date, compare to Basler epoch UTC
real_date = b_start_utc.date()    # both recorded same day
t_first_utc = dt.datetime.combine(real_date, first_dt.time())
# handle midnight wrap if thermal hour < basler hour by a lot
t_first_epoch = t_first_utc.replace(tzinfo=dt.timezone.utc).timestamp()
offset = t_first_epoch - wall_s[0]
print("=== ALIGNMENT ===")
print(f"  thermal start (UTC, rebuilt): {t_first_utc.isoformat()}  epoch {t_first_epoch:.6f}")
print(f"  basler  start (UTC)         : {b_start_utc.isoformat()}  epoch {wall_s[0]:.6f}")
print(f"  start offset (thermal - basler): {offset*1e3:.1f} ms")
print(f"  duration diff (thermal - basler): {(tt[-1]-cam_s[-1])*1e3:.1f} ms")
print(f"  => same rate? {abs(fps_t-fps_b)<1:.0f}  | overlap ~{min(tt[-1],cam_s[-1]):.2f} s")
