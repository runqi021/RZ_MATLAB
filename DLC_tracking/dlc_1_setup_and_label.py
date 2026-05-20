"""Step 1: Create DLC project, extract frames, open labeling GUI.
Edit settings below, then run. Close the GUI when done labeling.
"""

import deeplabcut
import os
import glob
import cv2
import ruamel.yaml

# ===== USER SETTINGS =====
video_dir    = r"C:\Users\Admin\Desktop\260505_breathing_wt\cam1"
proj_name    = "260505_breathing_wt_C3H"
experimenter = "RZ"
bodyparts    = ["dot1", "dot2", "dot3", "dot4"]
work_dir     = video_dir                   # DLC project created inside video folder
numframes    = 20                           # frames per video to extract
algo         = "uniform"                    # "uniform" or "kmeans"
net_type     = "resnet_101"                 # "resnet_50", "resnet_101", "hrnet_w32"

# ===== FIND VIDEOS =====
# Recurse — videos sit in cam1\cam1_<ts>_runNNN\*.avi from run_basler_n_runs_n_cam
vids = sorted(glob.glob(os.path.join(video_dir, "**", "*.avi"), recursive=True))
print(f"Found {len(vids)} videos in {video_dir}")
assert vids, f"No .avi files under {video_dir}"

# Read first video's dimensions for the crop field (DLC default is full frame)
cap = cv2.VideoCapture(vids[0])
vw = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
vh = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
cap.release()
crop_str = f"0, {vw}, 0, {vh}"
print(f"Video size: {vw} x {vh}  →  crop = '{crop_str}'")

# ===== CREATE PROJECT =====
print("\n=== Creating project ===")
config_path = deeplabcut.create_new_project(
    proj_name, experimenter, [vids[0]],
    working_directory=work_dir, copy_videos=False, multianimal=False,
)
print(f"Config: {config_path}")

# Patch config: add all videos (no copying), set bodyparts
yaml = ruamel.yaml.YAML()
with open(config_path, 'r') as f:
    cfg = yaml.load(f)

cfg['video_sets'] = {}
for v in vids:
    cfg['video_sets'][v] = {'crop': crop_str}
cfg['bodyparts'] = bodyparts
cfg['numframes2pick'] = numframes
cfg['net_type'] = net_type

with open(config_path, 'w') as f:
    yaml.dump(cfg, f)

# Clean up the single copied video
proj_vid_dir = os.path.join(os.path.dirname(config_path), 'videos')
for f in os.listdir(proj_vid_dir):
    fp = os.path.join(proj_vid_dir, f)
    if os.path.isfile(fp):
        os.remove(fp)

print(f"Added {len(vids)} videos, bodyparts={bodyparts}")

# ===== EXTRACT FRAMES =====
print(f"\n=== Extracting frames (algo={algo}) ===")
deeplabcut.extract_frames(config_path, mode='automatic', algo=algo, userfeedback=False)

# ===== OPEN LABELING GUI =====
print("\n=== Opening labeling GUI ===")
print("Label your frames, then CLOSE the GUI window when done.")
print(f"\nFor next step, run dlc_2_train_and_analyze.py with this config:")
print(f"  {config_path}\n")

deeplabcut.label_frames(config_path)

import napari
napari.run()
