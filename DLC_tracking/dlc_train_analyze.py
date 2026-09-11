"""Train + evaluate + analyze a per-scope DLC project (after labeling is done).

Reusable: pass the project's config.yaml. Videos analyzed = all raw Basler .avi
under the scope (the project's videos, incl. archived). PyTorch engine -> let the
training schedule defaults drive (do NOT pass maxiters/saveiters/displayiters).

  C:\\Users\\Admin\\.conda\\envs\\dlc310\\python.exe dlc_train_analyze.py "<...\\config.yaml>"
"""
import deeplabcut
import sys, os, glob

config = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\260613_breathing_thermalNbasler\5916297\5916297_experimental-RZ-2026-06-14\config.yaml"
assert os.path.isfile(config), f"config not found: {config}"
scope_dir = os.path.dirname(os.path.dirname(config))   # <root>\<scope>
print(f"config: {config}\nscope:  {scope_dir}", flush=True)

print("=== create_training_dataset ===", flush=True)
deeplabcut.create_training_dataset(config)

print("=== train_network (defaults) ===", flush=True)
deeplabcut.train_network(config)

print("=== evaluate_network ===", flush=True)
deeplabcut.evaluate_network(config, plotting=True)

vids = sorted(glob.glob(os.path.join(scope_dir, "**", "*.avi"), recursive=True))
print(f"=== analyze_videos: {len(vids)} videos ===", flush=True)
for v in vids:
    print("   ", v, flush=True)
deeplabcut.analyze_videos(config, vids, save_as_csv=True)

print("=== done: CSVs saved alongside each video ===", flush=True)
