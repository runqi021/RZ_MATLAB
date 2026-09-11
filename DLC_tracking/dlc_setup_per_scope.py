"""Create ONE DLC project per animal/scope under a root (e.g. C:\\260613_breathing_thermalNbasler).

For each scope folder (5840027, 5916296, ...) it makes a project named
"<scope>_experimental" with all of that scope's raw Basler videos (active +
archived) and bodyparts vL0/vL1/vR0/vR1. Project creation ONLY -- frames are
extracted/labeled later. Skips a scope that already has a "<scope>_experimental-*"
project.

Run in the DLC env:
    C:\\Users\\Admin\\.conda\\envs\\dlc310\\python.exe dlc_setup_per_scope.py
"""
import deeplabcut
import os, glob
import cv2
import ruamel.yaml

# ===== USER SETTINGS =====
root         = r"C:\260613_breathing_thermalNbasler"
experimenter = "RZ"
bodyparts    = ["vL0", "vL1", "vR0", "vR1"]
numframes    = 20
net_type     = "resnet_50"


def scopes_under(root):
    out = []
    for name in sorted(os.listdir(root)):
        d = os.path.join(root, name)
        if os.path.isdir(d) and glob.glob(os.path.join(d, "**", "*.avi"), recursive=True):
            out.append(name)
    return out


yaml = ruamel.yaml.YAML()
made, skipped = [], []
for scope in scopes_under(root):
    scope_dir = os.path.join(root, scope)
    proj_name = f"{scope}_experimental"
    if glob.glob(os.path.join(scope_dir, f"{proj_name}-*")):
        print(f"[{scope}] already has a project -- skip")
        skipped.append(scope); continue

    vids = sorted(glob.glob(os.path.join(scope_dir, "**", "*.avi"), recursive=True))
    print(f"[{scope}] {len(vids)} videos")
    # seed with the SMALLEST file to minimise the forced copy (symlink needs admin)
    seed = min(vids, key=os.path.getsize)
    cap = cv2.VideoCapture(vids[0])
    vw = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)); vh = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    cap.release()
    crop_str = f"0, {vw}, 0, {vh}"

    config_path = deeplabcut.create_new_project(
        proj_name, experimenter, [seed],
        working_directory=scope_dir, copy_videos=False, multianimal=False,
    )
    with open(config_path) as f:
        cfg = yaml.load(f)
    cfg['video_sets'] = {v: {'crop': crop_str} for v in vids}
    cfg['bodyparts'] = bodyparts
    cfg['numframes2pick'] = numframes
    cfg['net_type'] = net_type
    with open(config_path, 'w') as f:
        yaml.dump(cfg, f)

    pvd = os.path.join(os.path.dirname(config_path), 'videos')
    for f in os.listdir(pvd):
        fp = os.path.join(pvd, f)
        if os.path.isfile(fp):
            os.remove(fp)

    print(f"[{scope}] created {config_path}  ({len(vids)} videos, {vw}x{vh})")
    made.append((scope, config_path))

print("\n=== DLC PER-SCOPE SUMMARY ===")
for scope, cp in made:
    print(f"  made    {scope}: {cp}")
for scope in skipped:
    print(f"  skipped {scope}")
print("Frames NOT extracted. Per project: deeplabcut.extract_frames(cfg,...) then label_frames(cfg).")
