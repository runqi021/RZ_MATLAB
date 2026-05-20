"""
Pick a crop rectangle on the first frame of the behavior AVI and save it
as behavior_crop.json in the FOV folder. combine_roi_videos.py will pick
it up automatically.

Controls:
  - Drag a rectangle.
  - Press ENTER or SPACE to confirm.
  - Press C to cancel.
"""
import argparse, json, cv2
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--fov", required=True)
parser.add_argument("--frame", type=int, default=100,
                    help="Frame index to preview (post-toss frames work well)")
args = parser.parse_args()

fov = Path(args.fov)
avi = sorted(fov.glob("Basler_*.avi"))[0]
cap = cv2.VideoCapture(str(avi))
cap.set(cv2.CAP_PROP_POS_FRAMES, args.frame)
ok, fr = cap.read()
assert ok, f"Could not read frame {args.frame} from {avi}"
cap.release()

win = "Drag crop box — ENTER to confirm, C to cancel"
r = cv2.selectROI(win, fr, showCrosshair=True, fromCenter=False)
cv2.destroyAllWindows()

x, y, w, h = [int(v) for v in r]
if w == 0 or h == 0:
    print("  Cancelled (empty rect).")
    raise SystemExit(1)

out = {"x": x, "y": y, "w": w, "h": h, "avi": avi.name, "frame_h": fr.shape[0], "frame_w": fr.shape[1]}
out_path = fov / "behavior_crop.json"
out_path.write_text(json.dumps(out, indent=2))
print(f"  Saved crop {out} -> {out_path}")
