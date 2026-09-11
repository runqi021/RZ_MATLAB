#!/usr/bin/env python
"""
thermal_basler_sync_batch.py  --  make_sync() for EVERY session that has the
needed parts. Loops all DATA_ROOT/<animal>/cam1_* run folders that contain a
_nostrilC.mat (and .ats/.avi/timestamps), writing <date>_<animal>_n<k>_sync_video.mp4
into each. Errors on one session are logged and the batch continues.

Defaults to FULL clips. Full 80 s clips are large (~0.5-0.7 GB each) and slow;
set WINDOW / FRAME_STEP below to trim. Run in the flir env.
"""
import glob
import os
import sys
import traceback

import thermal_basler_sync_video as sv

# ============================ USER-EDITABLE ============================
DATA_ROOT  = r"D:\260615_thermalNbasler"
WHISK_DIR  = r"D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos"
WINDOW     = None     # None = full clip; or (t0,t1) seconds
OUT_FPS    = 100.0    # playback fps
FRAME_STEP = 2        # render every Nth basler frame (2 -> half size/time)
# ======================================================================


def main():
    runs = []
    for animal in sorted(glob.glob(os.path.join(DATA_ROOT, "*"))):
        if not os.path.isdir(animal):
            continue
        for rf in sorted(glob.glob(os.path.join(animal, "cam1_*"))):
            if glob.glob(os.path.join(rf, "Rec-*_nostrilC.mat")):
                runs.append(rf)
    if not runs:
        sys.exit(f"no run folders with _nostrilC.mat under {DATA_ROOT}")
    print(f"BATCH: {len(runs)} sessions\n")

    ok, failed = [], []
    for i, rf in enumerate(runs, 1):
        print(f"===== [{i}/{len(runs)}] {os.path.relpath(rf, DATA_ROOT)} =====", flush=True)
        try:
            sv.make_sync(rf, DATA_ROOT, WHISK_DIR, WINDOW, OUT_FPS, FRAME_STEP)
            ok.append(rf)
        except Exception as e:
            print(f"  ERROR: {e}", flush=True); traceback.print_exc()
            failed.append((rf, str(e)))
        print(flush=True)

    print(f"\n==== DONE: {len(ok)} ok, {len(failed)} failed ====")
    for rf, why in failed:
        print(f"  FAIL {os.path.relpath(rf, DATA_ROOT)}: {why}")


if __name__ == "__main__":
    main()
