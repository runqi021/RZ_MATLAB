"""Analyze-only: evaluate + analyze_videos using existing snapshots.
Use after stopping training early once the model has converged.
"""

import deeplabcut
import glob
import os

# ===== USER SETTINGS =====
video_dir = r"C:\260407_KA_electro_dorsal_whisking"
config    = r"C:\Users\Admin\Desktop\whisker_dorsal-RZ-2026-04-08\config.yaml"

assert os.path.isfile(config), f"Config not found: {config}"
print(f"Using config: {config}")

# ===== EVALUATE =====
print("=== Evaluating ===")
deeplabcut.evaluate_network(config, plotting=True)

# ===== ANALYZE VIDEOS → CSV =====
print("=== Analyzing videos ===")
vids = sorted(glob.glob(os.path.join(video_dir, "*.avi")))
print(f"Found {len(vids)} videos")
deeplabcut.analyze_videos(config, vids, save_as_csv=True)

# ===== LABELED OVERLAY VIDEOS =====
print("=== Creating labeled videos ===")
deeplabcut.create_labeled_video(
    config, vids,
    save_frames=False,
    color_by="bodypart",
    draw_skeleton=False,
    trailpoints=3,
)

print(f"\n=== Done === CSVs + labeled videos saved alongside videos in: {video_dir}")
