"""Serialize DLC train+analyze jobs on a single GPU.

The A1000 (6 GB) fits one DLC job at a time, so this waits for the GPU to be FREE
(sustained, i.e. the currently-running job has finished) before starting each
config in turn. Pass any number of project config.yaml paths; they run in order.

  C:\\Users\\Admin\\.conda\\envs\\dlc310\\python.exe dlc_queue.py "<cfg1>" "<cfg2>" ...

Safe to launch while another DLC job is already running -- it will wait for that
job to release the GPU before starting the first queued config.
"""
import subprocess, sys, os, time

DLC_PY = r"C:\Users\Admin\.conda\envs\dlc310\python.exe"
TRAIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dlc_train_analyze.py")
FREE_MIB = 800           # GPU considered free below this many MiB used
SUSTAIN = 3              # consecutive free polls (60 s apart) before starting


def gpu_used_mib():
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"]
        ).decode().strip().splitlines()[0]
        return int(out)
    except Exception:
        return 99999


def wait_free():
    free = 0
    while free < SUSTAIN:
        time.sleep(60)
        u = gpu_used_mib()
        free = free + 1 if u < FREE_MIB else 0
        print(f"  GPU {u} MiB (free streak {free}/{SUSTAIN})", flush=True)


def main():
    configs = sys.argv[1:]
    if not configs:
        sys.exit("pass one or more config.yaml paths")
    print(f"DLC queue: {len(configs)} jobs", flush=True)
    for i, cfg in enumerate(configs, 1):
        print(f"\n[{i}/{len(configs)}] waiting for GPU to free before: {cfg}", flush=True)
        wait_free()
        print(f"[{i}/{len(configs)}] GPU free -> training {cfg}", flush=True)
        r = subprocess.run([DLC_PY, TRAIN, cfg])
        print(f"[{i}/{len(configs)}] exit {r.returncode}: {cfg}", flush=True)
    print("\n=== DLC queue done ===", flush=True)


if __name__ == "__main__":
    main()
