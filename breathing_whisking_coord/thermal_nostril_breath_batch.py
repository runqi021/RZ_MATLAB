#!/usr/bin/env python
"""
thermal_nostril_breath_batch.py  --  batch the nostril thermal extraction.

Runs thermal_nostril_breath_single.extract_one() on EVERY DLC csv in a videos
folder. Each .ats + _dlc.json is auto-resolved from the csv name (animal prefix
+ _n# run index), and a _nostrilC.mat is written next to each .ats. Analysis is
separate (thermal_nostril_breath_analyze.m), done later.

Run in the flir env:
    C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe thermal_nostril_breath_batch.py
"""
import glob
import os
import sys
import traceback

from thermal_nostril_breath_single import extract_one

# ============================ USER-EDITABLE ============================
VIDEOS_DIR = r"D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos"
# ======================================================================


def main():
    csvs = sorted(glob.glob(os.path.join(VIDEOS_DIR, "*DLC*.csv")))
    if not csvs:
        sys.exit(f"no DLC csvs in {VIDEOS_DIR}")
    print(f"BATCH: {len(csvs)} videos in {VIDEOS_DIR}\n")

    done, failed = [], []
    for i, c in enumerate(csvs, 1):
        name = os.path.basename(c)
        print(f"===== [{i}/{len(csvs)}] {name} =====", flush=True)
        try:
            out = extract_one(c)
            done.append(out)
        except SystemExit as e:                 # resolve/validation bail-out
            print(f"  SKIP: {e}", flush=True)
            failed.append((name, str(e)))
        except Exception as e:                  # anything else: log, keep going
            print(f"  ERROR: {e}", flush=True)
            traceback.print_exc()
            failed.append((name, str(e)))
        print(flush=True)

    print(f"\n==== DONE: {len(done)} ok, {len(failed)} failed ====")
    for name, why in failed:
        print(f"  FAIL {name}: {why}")


if __name__ == "__main__":
    main()
