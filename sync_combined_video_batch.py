#!/usr/bin/env python
"""
sync_combined_video_batch.py  --  run sync_combined_video.py on EVERY run folder
under a root (each folder must have a *.ats, a *.avi and timestamps.csv). Resumable:
skips any folder that already has the matching *_combined_lag{LAG}.mp4.

  C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe sync_combined_video_batch.py \
      C:\\260613_breathing_thermalNbasler  [--decim 4] [--play-fps 60] [--lag 0] [--width 640]
"""
import argparse, os, sys, glob, subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
SINGLE = os.path.join(HERE, "sync_combined_video.py")


def run_folders(root, filt=None):
    # recursive: catches active runs AND archived/<run>/ (one level deeper)
    out = []
    for ats in glob.glob(os.path.join(root, "**", "*.ats"), recursive=True):
        folder = os.path.dirname(ats)
        if filt and filt not in folder:
            continue
        if glob.glob(os.path.join(folder, "*.avi")) and \
           os.path.isfile(os.path.join(folder, "timestamps.csv")):
            out.append(folder)
    return sorted(set(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--decim", type=int, default=4)
    ap.add_argument("--play-fps", type=float, default=60.0)
    ap.add_argument("--lag", type=int, default=0)
    ap.add_argument("--width", type=int, default=640)
    ap.add_argument("--filter", default=None,
                    help="only folders whose path contains this substring (e.g. 'hived' for archived/arhived)")
    args = ap.parse_args()

    folders = run_folders(args.root, args.filter)
    print(f"found {len(folders)} run folders under {args.root}"
          + (f" (filter '{args.filter}')" if args.filter else "") + "\n", flush=True)
    done, skipped, failed = [], [], []
    for k, folder in enumerate(folders, 1):
        run = os.path.basename(folder)
        outmp4 = os.path.join(folder, run + f"_combined_lag{args.lag}.mp4")
        lock = os.path.join(folder, ".combining.lock")
        tag = f"[{k}/{len(folders)}] {run}"
        if os.path.isfile(outmp4):
            print(f"{tag}  SKIP (exists)", flush=True)
            skipped.append(run); continue
        # atomic lock so concurrent batch instances never build the same folder
        try:
            fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, str(os.getpid()).encode()); os.close(fd)
        except FileExistsError:
            print(f"{tag}  SKIP (locked by another run)", flush=True)
            skipped.append(run); continue
        print(f"{tag}  building ...", flush=True)
        cmd = [sys.executable, SINGLE, folder,
               "--decim", str(args.decim), "--play-fps", str(args.play_fps),
               "--lag", str(args.lag), "--width", str(args.width)]
        r = subprocess.run(cmd)
        try:
            os.remove(lock)
        except OSError:
            pass
        if r.returncode == 0 and os.path.isfile(outmp4):
            print(f"{tag}  OK", flush=True); done.append(run)
        else:
            print(f"{tag}  FAILED (rc={r.returncode})", flush=True); failed.append(run)

    print("\n=== BATCH SUMMARY ===")
    print(f"  built  : {len(done)}")
    print(f"  skipped: {len(skipped)}")
    print(f"  failed : {len(failed)}  {failed if failed else ''}")


if __name__ == "__main__":
    main()
