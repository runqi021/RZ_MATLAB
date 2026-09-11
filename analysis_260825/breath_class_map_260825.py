"""breath_class_map_260825.py

Average chest-motion MAPS per breath-cycle class (normal vs GASP vs fail).

WHY THIS IS A SEPARATE PYTHON STEP
The behaviour videos are FFV1, which MATLAB's VideoReader will not open, so any
per-pixel work on them has to go through OpenCV. This writes a .mat and
breath_class_map_260825.m draws it -- the same Python -> .mat -> MATLAB split
the thermal pipeline uses.

WHY IT HAS TO READ THE VIDEO AT ALL
Nothing already on disk holds per-pixel time series. breath_pc1.mat keeps ONE
spatial mode (u1_img) plus the 1-D trace; breath_fixedmetric.mat keeps a 2-D
rigid displacement per frame (D_all) and a single diff image per run. A map that
separates cycle classes needs the frames themselves.

WHAT IT COMPUTES, per recording and per class
    map_foot  = mean over that class's cycles of the binned frame at the FOOT
    map_peak  = mean over that class's cycles of the binned frame at the PEAK
    map_diff  = mean over cycles of (peak frame - foot frame)
                -> the average inspiratory displacement map. Same quantity as
                   the pipeline's diff_imgs, but split by class, so it is
                   directly comparable to u1_img / diffImg.
    map_absdiff = mean over cycles of |peak - foot|
                -> direction-free motion magnitude, for when a class's
                   displacements partly cancel.

GEOMETRY IS COPIED FROM THE EXTRACTOR, NOT REINVENTED
crop_shift() and binned() below are the same operations as
breath_fixedmetric_extract.py:284-306, including the PER-RUN registration shift
(reg_shift). That shift exists because the FOV drifts up to 9 px over a session,
so a fixed box would measure a different piece of chest in the first and last
run. Dropping it would put the classes on subtly different anatomy.

FRAME INDEXING
cyc_*_idx in breath_cycle_class_pc1.mat are RAW pc1 frame indices, 1-based
(MATLAB). stride/sub_every is 1 for this session, so pc1 sample k is video frame
k-1 (0-based). The script asserts stride == 1 rather than assuming it.

USAGE
    python breath_class_map_260825.py --session "C:\\260824_..._vagotomized\\phys"

OUTPUT  <session>\\breath_class_maps.mat

Runqi Zhang / 2026-08-25
"""

import argparse
import glob
import os
import sys

import cv2
import numpy as np
from scipy.io import loadmat, savemat

CLASS_CODES = [-1, 0, 1]                 # fail, normal, gasp
CLASS_NAMES = ["fail", "normal", "gasp"]


def crop_shift(img, roi, shift_yx=(0, 0)):
    """Crop `roi` from `img`, displaced by this run's registration shift.

    Verbatim behaviour of breath_fixedmetric_extract.py:284.
    """
    if roi is None:
        return img
    x, y, w, h = roi
    dy, dx = shift_yx
    H, W = img.shape
    x = int(np.clip(x + dx, 0, max(W - w, 0)))
    y = int(np.clip(y + dy, 0, max(H - h, 0)))
    return img[y:y + h, x:x + w]


def binned(img, b):
    """Block-mean spatial binning; trailing partial blocks dropped."""
    h, w = img.shape
    return img[:h // b * b, :w // b * b].reshape(h // b, b, w // b, b).mean((1, 3))


def find_avi(folder):
    for pat in ("cam*.avi", "Basler_*.avi", "*.avi"):
        hits = sorted(glob.glob(os.path.join(folder, pat)))
        hits = [h for h in hits if "_labeled" not in h and "_mjpeg" not in h]
        if hits:
            return hits[0]
    return None


def read_frames(avi, wanted, roi, shift, b):
    """Binned, cropped frames for the 0-based indices in `wanted`.

    Sequential grab() with retrieve() only on wanted frames: seeking inside FFV1
    is unreliable, and grab-without-retrieve is roughly half the per-frame cost,
    so a full sweep retrieving ~200 of 3000 frames is fast anyway.
    """
    cap = cv2.VideoCapture(avi)
    if not cap.isOpened():
        raise RuntimeError("could not open %s" % avi)
    want = set(int(v) for v in wanted)
    out, i = {}, 0
    while want:
        if not cap.grab():
            break
        if i in want:
            ok, fr = cap.retrieve()
            if not ok:
                break
            g = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY).astype(np.float32)
            out[i] = binned(crop_shift(g, roi, shift), b).astype(np.float32)
            want.discard(i)
        i += 1
    cap.release()
    return out, i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--session", required=True,
                    help=r"session phys folder holding breath_roi.mat and the run folders")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    sess = args.session
    out_file = args.out or os.path.join(sess, "breath_class_maps.mat")

    roi_mat = loadmat(os.path.join(sess, "breath_roi.mat"))
    roi = tuple(int(v) for v in np.ravel(roi_mat["roi_xywh"])[:4])

    fm = loadmat(os.path.join(sess, "breath_fixedmetric.mat"),
                 variable_names=["reg_shift", "run_names", "bin_factor",
                                 "stride", "sub_every", "u1_img"])
    b = int(np.ravel(fm["bin_factor"])[0])
    stride = int(np.ravel(fm["stride"])[0])
    if stride != 1:
        sys.exit("stride is %d, not 1 -- the frame mapping in this script assumes 1" % stride)
    reg = np.array(fm["reg_shift"], float)                  # [nRun x 2] = (dy, dx)
    names = [str(np.ravel(n)[0]) for n in np.ravel(fm["run_names"])]
    u1 = np.array(fm["u1_img"], float)

    print("session : %s" % sess)
    print("roi     : x=%d y=%d w=%d h=%d, bin %d -> map %s" %
          (roi[0], roi[1], roi[2], roi[3], b, u1.shape))

    runs = sorted(d for d in glob.glob(os.path.join(sess, "*"))
                  if os.path.isdir(d)
                  and os.path.isfile(os.path.join(d, "breath_cycle_class_pc1.mat")))
    if not runs:
        sys.exit("no run folder under %s has breath_cycle_class_pc1.mat" % sess)
    print("labelled: %d recordings\n" % len(runs))

    H, W = u1.shape
    nR, nC = len(runs), len(CLASS_CODES)
    m_foot = np.full((H, W, nR, nC), np.nan, np.float32)
    m_peak = np.full((H, W, nR, nC), np.nan, np.float32)
    m_diff = np.full((H, W, nR, nC), np.nan, np.float32)
    m_absd = np.full((H, W, nR, nC), np.nan, np.float32)
    counts = np.zeros((nR, nC), int)
    used, skipped = [], []

    for r, folder in enumerate(runs):
        name = os.path.basename(folder.rstrip("\\/"))
        avi = find_avi(folder)
        if avi is None:
            skipped.append((name, "no avi"))
            print("  %-44s SKIP (no avi)" % name[:44])
            continue
        if name in names:
            shift = tuple(reg[names.index(name)])
        else:
            # A run absent from the fixedmetric list has no measured drift
            # correction. Using (0,0) would silently measure a different piece
            # of chest, so it is skipped instead.
            skipped.append((name, "not in run_names"))
            print("  %-44s SKIP (not in breath_fixedmetric run_names)" % name[:44])
            continue

        cl = loadmat(os.path.join(folder, "breath_cycle_class_pc1.mat"),
                     variable_names=["class_final", "cyc_foot_idx", "cyc_peak_idx"])
        cls = np.ravel(cl["class_final"]).astype(int)
        foot = np.ravel(cl["cyc_foot_idx"]).astype(int) - 1     # 1-based -> 0-based
        peak = np.ravel(cl["cyc_peak_idx"]).astype(int) - 1

        need = np.concatenate([foot, peak])
        frames, nread = read_frames(avi, need, roi, shift, b)
        missing = [i for i in need if i not in frames]
        if missing:
            print("  %-44s WARN %d/%d frames past end of video (%d frames)"
                  % (name[:44], len(missing), len(need), nread))

        for ci, code in enumerate(CLASS_CODES):
            sel = np.where(cls == code)[0]
            sel = [i for i in sel if foot[i] in frames and peak[i] in frames]
            counts[r, ci] = len(sel)
            if not sel:
                continue
            F = np.stack([frames[foot[i]] for i in sel])
            P = np.stack([frames[peak[i]] for i in sel])
            D = P - F
            m_foot[:, :, r, ci] = F.mean(0)
            m_peak[:, :, r, ci] = P.mean(0)
            m_diff[:, :, r, ci] = D.mean(0)
            m_absd[:, :, r, ci] = np.abs(D).mean(0)

        used.append(name)
        print("  %-44s  fail %3d | normal %3d | gasp %3d"
              % (name[:44], counts[r, 0], counts[r, 1], counts[r, 2]))

    savemat(out_file, {
        "map_foot": m_foot, "map_peak": m_peak,
        "map_diff": m_diff, "map_absdiff": m_absd,
        "counts": counts,
        "run_names": np.array(used + [n for n, _ in skipped], dtype=object),
        "run_used": np.array([os.path.basename(r.rstrip("\\/")) for r in runs], dtype=object),
        "class_codes": np.array(CLASS_CODES),
        "class_names": np.array(CLASS_NAMES, dtype=object),
        "roi_xywh": np.array(roi),
        "bin_factor": b,
        "u1_img": u1,
        "session": sess,
    }, do_compression=True)

    tot = counts.sum(0)
    print("\npooled cycles: fail %d, normal %d, gasp %d" % (tot[0], tot[1], tot[2]))
    if skipped:
        print("skipped %d run(s): %s" % (len(skipped), ", ".join(n for n, _ in skipped)))
    print("wrote %s" % out_file)


if __name__ == "__main__":
    main()
