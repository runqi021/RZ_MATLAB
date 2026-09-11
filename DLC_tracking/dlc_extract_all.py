"""Extract frames for ALL per-scope DLC projects under a root (automatic, uniform).

Non-interactive (userfeedback=False) so it can run unattended. Frames land in each
project's labeled-data\<video>\ ready for labeling.

Run in the DLC env:
    C:\\Users\\Admin\\.conda\\envs\\dlc310\\python.exe dlc_extract_all.py
"""
import deeplabcut
import glob, os

root = r"C:\260613_breathing_thermalNbasler"
algo = "uniform"     # fast, evenly-spaced; 'kmeans' is slower/diverse

configs = sorted(glob.glob(os.path.join(root, "*", "*_experimental-RZ-*", "config.yaml")))
print(f"{len(configs)} DLC projects found")
for cfg in configs:
    print(f"\n=== extract_frames: {cfg} ===", flush=True)
    try:
        deeplabcut.extract_frames(cfg, mode="automatic", algo=algo, userfeedback=False)
        print(f"OK {cfg}", flush=True)
    except Exception as e:
        print(f"FAILED {cfg}: {e}", flush=True)
print("\n=== frame extraction done for all projects ===")
