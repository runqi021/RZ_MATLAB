"""Step 2: Train network, evaluate, analyze all videos, save CSV.
Run after labeling is complete.
"""

import deeplabcut
import glob
import os

# ===== USER SETTINGS =====
video_dir = r"C:\Users\Admin\Desktop\260505_breathing_wt\cam1"

# Explicit config path — L-side project lives on Desktop, NOT inside video_dir.
# Flip this to the R-project path when retraining R.
config = r"C:\Users\Admin\Desktop\260505_breathing_wt\cam1\260505_breathing_wt_C3H-RZ-2026-05-05\config.yaml"
assert os.path.isfile(config), f"Config not found: {config}"
print(f"Using config: {config}")

# ===== CREATE TRAINING DATASET =====
print("=== Creating training dataset ===")
deeplabcut.create_training_dataset(config)

# ===== TRAIN =====
# PyTorch engine: let defaults from pytorch_config.yaml drive the training schedule.
print("=== Training ===")
deeplabcut.train_network(config)

# ===== EVALUATE =====
print("=== Evaluating ===")
deeplabcut.evaluate_network(config, plotting=True)

# ===== ANALYZE VIDEOS → CSV =====
print("=== Analyzing videos ===")
# Recurse — videos sit in cam1\cam1_<ts>_runNNN\*.avi from run_basler_n_runs_n_cam
vids = sorted(glob.glob(os.path.join(video_dir, "**", "*.avi"), recursive=True))
print(f"Found {len(vids)} videos")
deeplabcut.analyze_videos(config, vids, save_as_csv=True)

print(f"\n=== Done ===")
print(f"CSV files saved alongside videos in: {video_dir}")
