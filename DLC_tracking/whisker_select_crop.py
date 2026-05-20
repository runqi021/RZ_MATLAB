"""One-shot crop picker.

Opens a middle frame from one video, applies the same rotation used by
whisker_label_video_test.py, and lets you drag a rectangle to select the
crop ROI. Prints the coordinates so you can paste them into
whisker_label_video_test.py as `crop_xywh = [x, y, w, h]`.

Controls:
  - Drag a rectangle with the mouse
  - ENTER or SPACE = confirm
  - C = cancel
"""
import os
import glob
import cv2

# ===== USER SETTINGS =====
folder_path   = r"C:\260407_KA_electro_dorsal_whisking"
run_to_probe  = 1                # which run to preview
frame_index   = 2700             # middle-ish frame (5400 total → 2700)
rotate_deg    = 150              # MUST match whisker_label_video_test.py
# =========================

raw_avis = [
    a for a in glob.glob(os.path.join(folder_path, f"*run{run_to_probe:03d}.avi"))
    if "DLC_Resnet50" not in os.path.basename(a)
]
assert raw_avis, f"No raw AVI for run{run_to_probe:03d}"
src_avi = raw_avis[0]
print(f"Source: {src_avi}")

cap = cv2.VideoCapture(src_avi)
assert cap.isOpened(), f"Cannot open {src_avi}"
w = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
h = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
cap.set(cv2.CAP_PROP_POS_FRAMES, frame_index)
ok, frame = cap.read()
cap.release()
assert ok, f"Could not read frame {frame_index}"

# Same rotation as the labeling script
M = cv2.getRotationMatrix2D((w / 2.0, h / 2.0), rotate_deg, 1.0)
cos_ = abs(M[0, 0]); sin_ = abs(M[0, 1])
new_w = int(h * sin_ + w * cos_)
new_h = int(h * cos_ + w * sin_)
M[0, 2] += (new_w / 2.0) - (w / 2.0)
M[1, 2] += (new_h / 2.0) - (h / 2.0)
frame_rot = cv2.warpAffine(
    frame, M, (new_w, new_h),
    flags=cv2.INTER_LINEAR,
    borderMode=cv2.BORDER_CONSTANT, borderValue=(0, 0, 0),
)

print(f"Rotated frame: {new_w}x{new_h}. Drag ROI, press ENTER to confirm, C to cancel.")
roi = cv2.selectROI("select crop (ENTER=ok, C=cancel)", frame_rot,
                    showCrosshair=True, fromCenter=False)
cv2.destroyAllWindows()

x, y, cw, ch = [int(v) for v in roi]
if cw == 0 or ch == 0:
    print("No ROI selected (cancelled). Nothing to copy.")
else:
    print("\n==== PASTE THIS into whisker_label_video_test.py ====")
    print(f"crop_xywh = [{x}, {y}, {cw}, {ch}]   # x, y, width, height")
    print("=====================================================")
