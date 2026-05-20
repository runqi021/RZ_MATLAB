"""Run DLC training pipeline: create dataset → train → evaluate → analyze all videos → CSV.

Usage:
    python dlc_train_and_analyze.py
"""

import deeplabcut
import glob, os

config = r"C:\Users\Admin\Desktop\whisker_dorsal-RZ-2026-04-08\config.yaml"
video_dir = r"C:\260407_KA_electro_dorsal_whisking"

# 1. Create training dataset
print("=== Creating training dataset ===")
deeplabcut.create_training_dataset(config)

# 2. Train
print("=== Training (50k iterations) ===")
deeplabcut.train_network(config, maxiters=50000, saveiters=10000, displayiters=1000)

# 3. Evaluate
print("=== Evaluating ===")
deeplabcut.evaluate_network(config, plotting=True)

# 4. Analyze videos → CSV
print("=== Analyzing videos ===")
vids = sorted(glob.glob(os.path.join(video_dir, "*.avi")))
print(f"Found {len(vids)} videos")
deeplabcut.analyze_videos(config, vids, save_as_csv=True)

# # 5. Create labeled videos for QC (uncomment if needed)
# print("=== Creating labeled videos ===")
# deeplabcut.create_labeled_video(config, vids, draw_skeleton=True)

print("\n=== Done ===")
print(f"CSV files saved alongside videos in: {video_dir}")
