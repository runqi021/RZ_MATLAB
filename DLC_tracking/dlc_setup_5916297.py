"""Create the DLC project for scope 5916297 (full raw Basler videos, 4 whisker dots).

Project creation ONLY — registers all videos + sets bodyparts. Frame extraction
and labeling are deliberately left for later (run them when ready):
    deeplabcut.extract_frames(config_path, mode='automatic', algo='uniform', userfeedback=False)
    deeplabcut.label_frames(config_path)

Run in the DLC env:
    C:\\Users\\Admin\\.conda\\envs\\dlc310\\python.exe dlc_setup_5916297.py
"""
import deeplabcut
import os, glob
import cv2
import ruamel.yaml

# ===== USER SETTINGS =====
video_dir    = r"C:\260613_breathing_thermalNbasler\5916297"
proj_name    = "5916297_experimental"      # DLC task name (path-safe; no spaces)
experimenter = "RZ"
bodyparts    = ["vL0", "vL1", "vR0", "vR1"]
work_dir     = video_dir                    # project created inside the scope folder
numframes    = 20                           # frames/video to pick LATER
net_type     = "resnet_50"                  # matches the existing whisker projects

# ===== FIND VIDEOS (raw Basler .avi only; ignores the _combined/_crop .mp4) =====
vids = sorted(glob.glob(os.path.join(video_dir, "**", "*.avi"), recursive=True))
print(f"Found {len(vids)} raw Basler videos under {video_dir}")
assert vids, "no .avi found"
for v in vids:
    print("  ", v)

cap = cv2.VideoCapture(vids[0])
vw = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)); vh = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
cap.release()
crop_str = f"0, {vw}, 0, {vh}"
print(f"First video {vw}x{vh} -> crop '{crop_str}'")

# ===== CREATE PROJECT =====
config_path = deeplabcut.create_new_project(
    proj_name, experimenter, [vids[0]],
    working_directory=work_dir, copy_videos=False, multianimal=False,
)
print(f"Config: {config_path}")

# ===== PATCH CONFIG: all videos, bodyparts, defaults =====
yaml = ruamel.yaml.YAML()
with open(config_path) as f:
    cfg = yaml.load(f)
cfg['video_sets'] = {v: {'crop': crop_str} for v in vids}
cfg['bodyparts'] = bodyparts
cfg['numframes2pick'] = numframes
cfg['net_type'] = net_type
with open(config_path, 'w') as f:
    yaml.dump(cfg, f)

# remove the single auto-copied video stub (we use copy_videos=False)
pvd = os.path.join(os.path.dirname(config_path), 'videos')
for f in os.listdir(pvd):
    fp = os.path.join(pvd, f)
    if os.path.isfile(fp):
        os.remove(fp)

print(f"\nProject ready: {len(vids)} videos, bodyparts={bodyparts}, net={net_type}")
print("Frames NOT extracted yet. When ready:")
print(f"  deeplabcut.extract_frames(r'{config_path}', mode='automatic', algo='uniform', userfeedback=False)")
print(f"  deeplabcut.label_frames(r'{config_path}')")
