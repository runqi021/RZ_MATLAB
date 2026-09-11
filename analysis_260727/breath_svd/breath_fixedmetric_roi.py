"""
breath_fixedmetric_roi.py
=========================
Pick the anatomical ROI for the fixed-metric breathing pipeline.

Opens one preview window, you drag a box, ENTER confirms.  Same interaction as
orofacial_crop_extract.py, with two deliberate differences:

  * ONE ROI FOR THE WHOLE SESSION, not one per video.  A per-video box is a
    per-video spatial weighting, which is exactly the non-comparability the
    fixed-metric pipeline exists to remove.  Runs are registered to a common
    reference, so a single box propagates correctly to all of them and then
    follows the anatomy as the FOV drifts.

  * THE PREVIEW IS OVERLAID WITH BREATHING-BAND POWER.  Chest wall motion is
    NOT spatially uniform, and for breathing that is physiology rather than
    nuisance: ribcage and abdominal compartments move with different amplitudes
    and different phases (the classic two-compartment picture).  On the 260723
    session one region held 0.68x of its amplitude over 75 min while another
    fell to 0.10x, and their motion axes differed by ~70 degrees.  So "breathing
    amplitude" is only defined once you say WHICH COMPARTMENT you measured.
    Hot colour = strong motion in BAND; put the box on one compartment, not
    across both.

USAGE
-----
  C:/Users/Admin/.conda/envs/dlc310/python.exe breath_fixedmetric_roi.py
  (optionally:  ... breath_fixedmetric_roi.py "D:\\some\\other\\session")

Writes <ROOT_DIR>/breath_roi.mat with roi_xywh = [x, y, w, h].
breath_fixedmetric_extract.py picks that file up automatically and wipes its
cache if the box changed.

Press ENTER/SPACE to confirm, C or ESC to cancel.  Selecting nothing keeps any
existing ROI untouched.
"""

import os
import sys
import glob

import cv2
import numpy as np
from scipy.io import savemat, loadmat

# ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR   = r"D:\260724_heteroSHI_nRuns\cam1"
PREVIEW_RUN = 0        # index into the sorted run list; 0 = first run
NFRAMES    = 900       # frames used to build the preview + power map
BAND       = (0.5, 8.0)
FPS_FALLBACK = 60.0
BIN        = 2         # spatial binning for the power map only (speed)
OVERLAY    = 0.45      # 0 = plain grayscale preview, 1 = pure power map
# ---------------------------------------------------------------------------


def log(m):
    print(m, flush=True)


def read_fps(folder, fallback):
    ts = os.path.join(folder, "timestamps.csv")
    if not os.path.isfile(ts):
        return float(fallback)
    try:
        t = np.loadtxt(ts, delimiter=",", skiprows=1, usecols=1) * 1e-9
        return float((t.size - 1) / (t[-1] - t[0])) if t.size >= 10 else float(fallback)
    except Exception:
        return float(fallback)


def build_preview(avi, nframes, fps, band, b):
    """Mean image + per-pixel power fraction inside `band`, both full-res."""
    cap = cv2.VideoCapture(avi)
    frames = []
    while len(frames) < nframes:
        ok, fr = cap.read()
        if not ok:
            break
        frames.append(cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY).astype(np.float32))
    cap.release()
    if not frames:
        sys.exit("ERROR: could not read any frames from %s" % avi)
    Y = np.stack(frames)
    mean_img = Y.mean(0)

    H, W = mean_img.shape
    Yb = Y[:, :H // b * b, :W // b * b].reshape(len(Y), H // b, b, W // b, b).mean((2, 4))
    Yb = Yb - Yb.mean(0)
    T = len(Yb)
    F = np.fft.rfft(Yb * np.hanning(T)[:, None, None], axis=0)
    f = np.fft.rfftfreq(T, 1.0 / fps)
    P = np.abs(F) ** 2
    inb = (f >= band[0]) & (f < band[1])
    pw = P[inb].sum(0)
    pw = cv2.resize(pw, (W, H), interpolation=cv2.INTER_LINEAR)
    return mean_img, pw


def stretch(img, lo_p=1, hi_p=99):
    lo, hi = np.percentile(img, [lo_p, hi_p])
    return (np.clip((img - lo) / max(hi - lo, 1e-6), 0, 1) * 255).astype(np.uint8)


def main():
    root = ROOT_DIR
    if len(sys.argv) > 1:
        root = sys.argv[1]
    if not os.path.isdir(root):
        sys.exit("ROOT_DIR not found: %s" % root)

    runs = []
    for d in sorted(glob.glob(os.path.join(root, "*"))):
        if os.path.isdir(d) and not os.path.basename(d).startswith("."):
            avi = sorted(glob.glob(os.path.join(d, "*.avi")))
            if len(avi) == 1:
                runs.append(avi[0])
    if not runs:
        sys.exit("no run folders with a single .avi under %s" % root)

    k = min(PREVIEW_RUN, len(runs) - 1)
    avi = runs[k]
    fps = read_fps(os.path.dirname(avi), FPS_FALLBACK)
    log("preview from %s  (%.3f fps, %d frames)" % (os.path.basename(avi), fps, NFRAMES))
    log("building breathing-band power map (%.1f-%.1f Hz)..." % BAND)

    mean_img, pw = build_preview(avi, NFRAMES, fps, BAND, BIN)

    gray = cv2.cvtColor(stretch(mean_img), cv2.COLOR_GRAY2BGR)
    heat = cv2.applyColorMap(stretch(np.log10(pw + 1e-6), 5, 99.5), cv2.COLORMAP_INFERNO)
    disp = cv2.addWeighted(gray, 1.0 - OVERLAY, heat, OVERLAY, 0)

    prev = None
    roi_file = os.path.join(root, "breath_roi.mat")
    if os.path.isfile(roi_file):
        try:
            prev = tuple(int(v) for v in loadmat(roi_file)["roi_xywh"].ravel()[:4])
            x, y, w, h = prev
            cv2.rectangle(disp, (x, y), (x + w, y + h), (0, 255, 0), 1)
            log("existing ROI shown in green: x=%d y=%d w=%d h=%d" % prev)
        except Exception:
            prev = None

    log("")
    log("  drag a box over the region whose motion you want to measure")
    log("  hot colour = strong breathing-band motion")
    log("  ENTER/SPACE = confirm,  C/ESC = cancel (keeps any existing ROI)")

    title = "Select breathing ROI  (drag, ENTER=ok, C=cancel)"
    x, y, w, h = cv2.selectROI(title, disp, showCrosshair=True, fromCenter=False)
    cv2.destroyAllWindows()
    cv2.waitKey(1)

    if w == 0 or h == 0:
        log("\nno ROI selected -- nothing written%s"
            % ("; existing ROI kept" if prev else ""))
        return 0

    roi = (int(x), int(y), int(w), int(h))
    frac = float(pw[y:y + h, x:x + w].sum() / max(pw.sum(), 1e-30))
    log("\nROI (x,y,w,h) = (%d, %d, %d, %d)   %d x %d px" % (roi + (w, h)))
    log("  holds %.1f%% of the frame's total breathing-band motion power" % (100 * frac))
    if w < 40 or h < 40:
        log("  WARNING: box under 40 px on a side -- phase correlation gets noisy;")
        log("           this signal is only ~1.2 px peak-to-peak to begin with.")

    savemat(roi_file, dict(roi_xywh=np.array(roi, float),
                           preview_run=os.path.basename(avi),
                           band=np.array(BAND, float),
                           power_frac=frac))
    log("\nwrote %s" % roi_file)
    if prev is not None and prev != roi:
        log("ROI changed from %s -- the extractor will clear its cache and re-decode." % (prev,))
    log("next: python breath_fixedmetric_extract.py   (then breath_fixedmetric_analyze.m)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
