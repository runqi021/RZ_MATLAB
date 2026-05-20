"""Step 3: Create labeled overlay videos for QC.
Run after analyze_videos is complete.
"""

import deeplabcut
import glob
import os

# ===== USER SETTINGS =====
video_dir = r"C:\Users\Admin\Desktop\260505_breathing_wt\cam1"

config = r"C:\Users\Admin\Desktop\260505_breathing_wt\cam1\260505_breathing_wt_C3H-RZ-2026-05-05\config.yaml"
print(f"Using config: {config}")

# ===== CREATE LABELED VIDEOS =====
# Recurse — videos sit in cam1\cam1_<ts>_runNNN\*.avi from run_basler_n_runs_n_cam
vids = sorted(glob.glob(os.path.join(video_dir, "**", "*.avi"), recursive=True))
print(f"Found {len(vids)} videos")

deeplabcut.create_labeled_video(
    config, vids,
    save_frames=False,
    color_by="bodypart",
    draw_skeleton=False,
    trailpoints=3,
)

print(f"\n=== Done === Labeled videos saved in: {video_dir}")
