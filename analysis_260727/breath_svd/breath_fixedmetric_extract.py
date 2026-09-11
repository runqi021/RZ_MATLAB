"""
breath_fixedmetric_extract.py
=============================
Acquisition-comparable breathing readout from Basler chest video.

WHY THIS EXISTS
---------------
Per-video SVD (breath_svd_pc1.m) re-derives its spatial basis for every video,
so "PC1" is a different physical quantity in every recording: arbitrary scale,
arbitrary sign, and arbitrary basis.  Z-scoring fixes only the scale.  This
script produces metrics that are comparable across runs and animals BY
CONSTRUCTION rather than by post-hoc normalisation:

  1. disp   -- rigid sub-pixel translation of the chest wall, in PIXELS, by
               phase correlation of each frame against that run's own mean
               frame, then projected onto a motion axis frozen from the pooled
               baseline runs.  No basis ambiguity (translation is translation in
               every video), no sign ambiguity, physical units, and immune to
               slow FOV drift because each run is referenced to its own mean.
               >>> This is the primary metric. <<<

  2. fb     -- frozen-basis SVD.  The spatial mode u1 is computed ONCE from the
               pooled baseline runs and then applied unchanged to every run, so
               all runs share one yardstick and the sign is decided once.

  3. pv     -- per-video SVD PC1, the incumbent method, computed only so the
               other two can be benchmarked against it.

Run means are rigidly registered to a reference before the frozen basis is
built or applied; on real data this FOV drifts several pixels over an hour and
an unregistered fixed basis slowly slides off the target.

OUTPUT
------
<ROOT_DIR>/breath_fixedmetric.mat  -- consumed by breath_fixedmetric_analyze.m
<ROOT_DIR>/.breathcache/*.npy      -- binned cubes; delete the folder to redo

ENVIRONMENT
-----------
Run with the dlc310 interpreter:
  C:/Users/Admin/.conda/envs/dlc310/python.exe breath_fixedmetric_extract.py

NOT cellpose-gpu -- its BLAS/LAPACK segfaults on matmul and svd (verified), so
anything with linear algebra dies silently there.

Basler AVIs are FFV1; MATLAB VideoReader cannot decode them, which is why all
video I/O lives here rather than in the .m file.
"""

import os
import re
import sys
import glob
import time
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

import cv2
import numpy as np
from scipy.io import savemat

# ----------------------------- USER PARAMETERS -----------------------------
ROOT_DIR      = r"D:\260724_heteroSHI_nRuns\cam1"  # folder of cam*_run### folders

# Which runs BUILD THE SHARED SPATIAL BASIS: the frozen mode u1 and the frozen
# motion axis are computed from these runs only, then applied unchanged to every
# run.  Run numbers come from the folder name (run001 -> 1); None = use all runs.
#
# This is a technical choice about which runs best define the spatial mode, and
# it is NOT the amplitude reference.  "1.0 = one baseline breath" is set by
# BASELINE_RUNS in breath_fixedmetric_analyze.m, which you can change in seconds
# without re-extracting.  Pick runs here with clean, representative breathing;
# pick runs there to define the condition you want to measure against.
BASIS_RUNS = list(range(1, 16))             # runs 001-015

FPS_FALLBACK  = 60.0    # used only if timestamps.csv is missing/unusable
BAND          = (0.5, 8.0)   # Hz, breathing band used for QC power fractions
BIN           = 4       # spatial binning factor for the SVD basis

# ONE temporal stride for everything: displacement, cube, and therefore all
# three metrics.  Using a single stride means every metric lands natively on the
# same time base and no interpolation is needed anywhere.
#
# It MUST satisfy Nyquist for BAND -- fps/STRIDE has to stay above 2*BAND[1] --
# and that is asserted at startup rather than left to trust.  At 60 fps with an
# 8 Hz ceiling the largest legal value is 3.  (Breathing here is ~2.2 Hz, so a
# stride of 12 -> 5 Hz sampling would alias it to 2.33 Hz and produce a
# clean-looking, entirely wrong trace.)
#
# It is also the main speed control.  Frames are demuxed either way, but only
# every STRIDE-th is colour-converted and phase-correlated, which is where the
# time goes.  Measured, per 3600-frame run: stride 1 = 66.7 s, stride 2 = 36.4 s,
# stride 3 = 31.9 s.  Stride 2 takes 1.8x off the clock and still leaves Nyquist
# at 15 Hz against an 8 Hz ceiling; stride 3 buys only 12% more and halves that
# margin, so 2 is the recommended setting.
#
# "auto" (the default now) picks the largest stride that still leaves a healthy
# Nyquist margin for BAND at the fps ACTUALLY measured, capped at STRIDE_MAX.
# A hard-wired 2 was written for the 60 fps free-running sessions this pipeline
# was developed on and is WRONG on a 2P-triggered 30 fps session: 30/2 = 15 Hz
# against a 16 Hz requirement for an 8 Hz band, i.e. aliasing. Set an integer to
# force it.
STRIDE        = "auto"
STRIDE_MAX    = 3       # never decimate harder than this, even if fps allows
NYQ_MARGIN    = 2.5     # need fps/STRIDE > NYQ_MARGIN * BAND[1] (2.0 = bare Nyquist)

# Stride used ONLY when assembling the basis matrix, applied on top of the
# cached cube.  The basis is a spatial object, so it needs a good sample of the
# spatial covariance, not a Nyquist-valid time series -- this can be coarse.
BASIS_EVERY   = 6
MAX_BASIS_FR  = 3000    # cap on pooled frames used to build the basis

# Worker processes for the two per-run passes.  Runs are independent, but the
# FFV1 decoder is already internally multithreaded, so one process alone spreads
# across several cores and throughput saturates early.  Measured on this machine
# (20 cores, NVMe at 1.2 GB/s so decoding is CPU-bound, not disk-bound):
# 1 worker 1.0x, 2 -> 2.0x, 4 -> 2.75x, 8 -> 3.0x, 12 -> 3.1x.  Past 4 the extra
# processes buy essentially nothing, so 4 is the default.  Set 1 to disable.
N_WORKERS     = 4

REF_NFRAMES   = 300     # frames averaged to form each run's reference image

# A run is dropped when it delivered less than this fraction of the frames IT
# WAS ASKED FOR.  An aborted recording is usually still decodable, so it has to
# be excluded on length or it drags the common (min) length down and clips every
# other run.  The reference is the run's OWN intended count (2P metadata, else
# the hardware timestamps), never the cohort median -- see the drop block in
# main() for why the median version was wrong.
SHORT_RUN_FRAC = 0.9
SIGN_FLIP     = False   # flip the frozen polarity of all three metrics
FORCE         = False   # True = ignore the cube cache and re-decode every video

# Optional [x, y, w, h] box, full-resolution pixels, restricting BOTH the phase
# correlation and the SVD basis.  None = whole frame.  For chest imaging a whole-
# frame metric averages over compartments that move differently (ribcage vs
# abdomen), so an explicit box is strongly preferred -- see breath_fixedmetric_roi.py.
ROI_XYWH      = None
# ---------------------------------------------------------------------------


def log(msg):
    print(msg, flush=True)


def find_runs(root):
    """Run folders holding exactly one .avi, sorted by name."""
    out = []
    for d in sorted(glob.glob(os.path.join(root, "*"))):
        if not os.path.isdir(d) or os.path.basename(d).startswith("."):
            continue
        avi = sorted(glob.glob(os.path.join(d, "*.avi")))
        if len(avi) != 1:
            continue
        m = re.search(r"run(\d+)", os.path.basename(d))
        # Wall-clock time of the run.  mtime is the file's last write, i.e. the
        # END of the recording; on this rig it lands start+60.0 s for every
        # complete run, so it carries the same information as the folder-name
        # start time plus one run duration.  The folder-name stamp is parsed as
        # a cross-check because a mtime that does NOT sit at start+duration is a
        # good sign the run was cut short.
        ts = re.search(r"_(\d{8})_(\d{6})_run", os.path.basename(d))
        t_name = None
        if ts:
            try:
                t_name = time.mktime(time.strptime(ts.group(1) + ts.group(2),
                                                   "%Y%m%d%H%M%S"))
            except ValueError:
                t_name = None
        out.append(dict(folder=d, avi=avi[0], name=os.path.basename(d),
                        num=int(m.group(1)) if m else len(out) + 1,
                        mtime=os.path.getmtime(avi[0]), t_name=t_name))
    return out


def read_fps_2p(folder):
    """Breath-cam fps from the 2P METADATA -- the primary source.

    The breath camera is 2P-frame-TRIGGERED, so its rate IS the calcium imaging
    frame rate. That is what breath_svd_pc1.m uses (via detect_session_fps) and
    what the rest of this codebase treats as authoritative; this pipeline was the
    one place that did not, and it cost a 2x time-base error on 260728_vglut2.
    Verified on that session: 6000 camera frames for 6000 2P frames, exactly 1:1.

    Reads *_ch1_meta.mat (falling back to *_meta.mat) and takes fps, or
    scanFrameRate_raw. Returns None if there is no usable metadata, so the caller
    can fall through to the timestamps.
    """
    cands = sorted(glob.glob(os.path.join(folder, "*_ch1_meta.mat"))) or \
            sorted(glob.glob(os.path.join(folder, "*_meta.mat")))
    for mf in cands:
        try:
            from scipy.io import loadmat
            M = loadmat(mf, squeeze_me=True, variable_names=("fps", "scanFrameRate_raw"))
            for k in ("fps", "scanFrameRate_raw"):
                if k in M:
                    v = float(np.asarray(M[k]).ravel()[0])
                    if np.isfinite(v) and 1.0 < v < 1000.0:
                        return v
        except Exception:
            continue
    return None


def read_fps(folder, fallback):
    """Breath-cam fps. 2P metadata first, then hardware timestamps, then fallback.

    The AVI header fps on this rig is not trustworthy (see the fix_avi_timing
    work), so it is never consulted.
    """
    v = read_fps_2p(folder)
    if v is not None:
        return v, True
    # GLOB, do not assume the bare name. The legacy Pylon layout wrote
    # "timestamps.csv", but the dual-cam GUI writes
    # "cam1_<YYYYMMDD>_<HHMMSS>_run###_timestamps.csv". Matching only the bare
    # name meant the file was never found on GUI-acquired sessions and fps
    # silently fell back to FPS_FALLBACK -- 60 on a 30 fps camera, i.e. a time
    # axis 2x too fast and every breath frequency doubled, with nothing in the
    # log to say so. Now the fallback is announced by the caller.
    cands = sorted(glob.glob(os.path.join(folder, "*timestamps.csv")))
    exact = os.path.join(folder, "timestamps.csv")
    if os.path.isfile(exact):
        cands.insert(0, exact)
    for ts in cands:
        try:
            t = np.loadtxt(ts, delimiter=",", skiprows=1, usecols=1) * 1e-9
            if t.size < 10 or not np.isfinite(t).all() or t[-1] <= t[0]:
                continue
            return float((t.size - 1) / (t[-1] - t[0])), True
        except Exception:
            continue
    return float(fallback), False


def expected_frames(folder):
    """How many camera frames this run was SUPPOSED to record.

    Needed to tell "aborted" from "deliberately shorter than its neighbours".
    Sessions routinely mix acquisition lengths (a 3000-frame run next to a
    6000-frame one), so the only meaningful completeness test is against the
    run's own intended count.

    Primary source is the 2P metadata: the breath camera is 2P-frame-TRIGGERED,
    one camera frame per 2P frame (verified 1:1 on 260728_vglut2), so the
    intended count is framesPerSlice x numSlices.  Falls back to the number of
    rows in the hardware timestamps, which is one row per frame actually
    written, and is present even for runs whose calcium data has not been
    processed yet (no *_meta.mat).

    Returns (n_frames, source_label), or (None, None) if neither exists.
    """
    cands = sorted(glob.glob(os.path.join(folder, "*_ch1_meta.mat"))) or \
            sorted(glob.glob(os.path.join(folder, "*_meta.mat")))
    for mf in cands:
        try:
            from scipy.io import loadmat
            M = loadmat(mf, squeeze_me=True,
                        variable_names=("framesPerSlice", "numSlices"))
            if "framesPerSlice" not in M:
                continue
            fps_slice = float(np.asarray(M["framesPerSlice"]).ravel()[0])
            nz = float(np.asarray(M["numSlices"]).ravel()[0]) if "numSlices" in M else 1.0
            n = fps_slice * max(nz, 1.0)
            if np.isfinite(n) and n >= 1:
                return int(round(n)), "2P meta"
        except Exception:
            continue

    cands = sorted(glob.glob(os.path.join(folder, "*timestamps.csv")))
    exact = os.path.join(folder, "timestamps.csv")
    if os.path.isfile(exact):
        cands.insert(0, exact)
    for ts in cands:
        try:
            with open(ts, "r") as fh:
                n = sum(1 for ln in fh if ln.strip()) - 1   # minus the header
            if n >= 1:
                return int(n), "timestamps"
        except Exception:
            continue
    return None, None


def crop_shift(img, roi, shift_yx=(0, 0)):
    """Crop `roi` from `img`, displaced by this run's registration shift.

    The box has to follow the anatomy rather than sit at fixed pixel
    coordinates: this FOV drifts up to 9 px over a session, so a fixed box
    slowly slides off the target and the "same" ROI would measure a different
    piece of chest in the first and last run -- silently reintroducing the
    per-video-basis problem the whole pipeline exists to remove.
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
    """Block-mean spatial binning; trailing partial blocks are dropped."""
    h, w = img.shape
    return img[:h // b * b, :w // b * b].reshape(h // b, b, w // b, b).mean((1, 3))


def ref_frame(avi, nref):
    """Full-frame reference image = mean of the first `nref` frames.

    Always full-frame: registration between runs must be measured on the whole
    image, before any ROI is applied, or the ROI's own drift corrupts the shift
    estimate that is supposed to correct it.

    Deliberately not the whole-run mean: that would cost an extra full decode of
    every video, and the choice of reference only sets an additive offset in the
    displacement, which is removed by mean-centring anyway.  Phase correlation
    is accurate over many pixels, so a slow within-run drift away from this
    reference costs nothing.
    """
    cap = cv2.VideoCapture(avi)
    acc, n = None, 0
    while n < nref:
        ok, fr = cap.read()
        if not ok:
            break
        g = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY).astype(np.float32)
        acc = g if acc is None else acc + g
        n += 1
    cap.release()
    if acc is None:
        raise RuntimeError("no frames decoded from %s" % avi)
    return acc / n


def _worker_ref(a):
    """Build one run's reference image (parallel phase 1). Writes, returns status."""
    avi, nref, out = a
    cv2.setNumThreads(1)          # workers must not each grab every core
    try:
        if not os.path.isfile(out):
            np.save(out, ref_frame(avi, nref))
        return out, None
    except Exception as e:
        return out, str(e)


def _worker_disp(a):
    """One run's displacement + cube (parallel phase 3). Writes to cache.

    Results go to disk rather than back through the process pool: a cube is tens
    of MB and shipping it over IPC would cost more than the compute saved.
    """
    avi, ref_path, roi, shift, b, stride, d_path, c_path = a
    cv2.setNumThreads(1)
    try:
        ref = np.load(ref_path)
        D, sub = pass_displacement(avi, ref, roi, tuple(shift), b, stride)
        np.save(d_path, D)
        np.save(c_path, sub)
        return d_path, int(len(D)), float(D[:, 1].std()), None
    except Exception as e:
        return d_path, 0, 0.0, str(e)


def top_pc(X, iters=300, tol=1e-10):
    """Leading PC of a centred n x p matrix.

    Only the FIRST component is ever needed, so this uses power iteration on
    X^T X (never formed) rather than a full eigendecomposition.  Forming the
    n x n Gram matrix and calling eigh computes all n eigenvectors to use one:
    measured on real data that is 11.9 s per run versus 0.19 s here, a 62x
    difference, and the two agree to |v.v_exact| = 1.000000 with an identical
    variance fraction.  Falls back to the exact route if it fails to converge,
    which can happen when the top two eigenvalues are nearly degenerate.

    Returns (unit spatial mode v [p], scores Xv [n], variance fraction).
    """
    # einsum accumulates in float64 without materialising a float64 copy of X
    # (which would be an extra ~120 MB per call at these sizes).
    total = float(np.einsum("ij,ij->", X, X, dtype=np.float64))
    rng = np.random.default_rng(0)          # fixed seed -> reproducible
    v = rng.standard_normal(X.shape[1]).astype(X.dtype)
    v /= max(np.linalg.norm(v), 1e-30)
    lam_prev, converged = 0.0, False
    for _ in range(iters):
        vn = X.T @ (X @ v)
        lam = float(np.linalg.norm(vn))
        if lam <= 0:
            break
        vn /= lam
        if abs(lam - lam_prev) <= tol * max(lam, 1e-30):
            v, lam_prev, converged = vn, lam, True
            break
        v, lam_prev = vn, lam

    if not converged:
        G = X @ X.T
        w, V = np.linalg.eigh(G)
        order = np.argsort(w)[::-1]
        lam_prev = max(float(w[order[0]]), 1e-30)
        v = (X.T @ V[:, order[0]]) / np.sqrt(lam_prev)
        v /= max(np.linalg.norm(v), 1e-30)

    return v, X @ v, lam_prev / max(total, 1e-30)


def reg_shift(mov, ref):
    """Integer (dy, dx) aligning `mov` onto `ref` by FFT cross-correlation."""
    A = np.fft.rfft2(mov - mov.mean())
    B = np.fft.rfft2(ref - ref.mean())
    c = np.fft.irfft2(A * np.conj(B), s=mov.shape)
    iy, ix = np.unravel_index(np.argmax(c), c.shape)
    dy = iy if iy < mov.shape[0] // 2 else iy - mov.shape[0]
    dx = ix if ix < mov.shape[1] // 2 else ix - mov.shape[1]
    return int(dy), int(dx)


def pass_displacement(avi, ref, roi, shift_yx, b, cube_every):
    """Full decode: per-frame sub-pixel displacement + a binned, registered cube.

    Displacement is measured against the run's OWN mean frame, so a slow
    between-run FOV drift cannot leak into it.  The cached cube, by contrast, is
    shifted onto the common reference grid because the frozen basis needs every
    run on the same pixel lattice.
    """
    cap = cv2.VideoCapture(avi)
    # `ref` arrives full-frame; cropping it at this run's shifted box puts the
    # frames and the cube on a common anatomical grid, so no further rolling is
    # needed.  With no ROI there is nothing to track, so the cube is rolled onto
    # the reference grid the old way instead.
    refc = crop_shift(ref, roi, shift_yx)
    win = cv2.createHanningWindow((refc.shape[1], refc.shape[0]), cv2.CV_32F)
    dy, dx = shift_yx
    D, sub, i = [], [], 0
    while True:
        # grab() demuxes and decodes; retrieve() does the colour conversion and
        # copy, which is roughly half the per-frame cost.  Skipping retrieve on
        # frames we do not keep is where the speed-up comes from -- 66.7 s ->
        # 36.4 s per run at stride 2.  Displacement and cube share one stride,
        # so both land on the same time base and nothing needs interpolating.
        if not cap.grab():
            break
        if i % cube_every == 0:
            ok, fr = cap.retrieve()
            if not ok:
                break
            g = cv2.cvtColor(fr, cv2.COLOR_BGR2GRAY).astype(np.float32)
            gc = crop_shift(g, roi, shift_yx)
            (sx, sy), _ = cv2.phaseCorrelate(refc, gc, win)  # OpenCV gives (dx, dy)
            D.append((sx, sy))
            gb = gc if roi is not None else np.roll(g, (-dy, -dx), axis=(0, 1))
            sub.append(binned(gb, b).astype(np.float32))
        i += 1
    cap.release()
    return np.asarray(D, np.float64), np.stack(sub).astype(np.uint8)


def band_frac(x, fps, band):
    """Fraction of 0.1-Nyquist power that falls inside `band`."""
    x = np.asarray(x, float)
    x = x - x.mean()
    P = np.abs(np.fft.rfft(x * np.hanning(x.size))) ** 2
    f = np.fft.rfftfreq(x.size, 1.0 / fps)
    inb = (f >= band[0]) & (f < band[1])
    tot = (f >= 0.1) & (f < 0.98 * fps / 2)
    return float(P[inb].sum() / max(P[tot].sum(), 1e-30))


def dom_freq(x, fps, fmin=0.3):
    x = np.asarray(x, float)
    x = x - x.mean()
    P = np.abs(np.fft.rfft(x * np.hanning(x.size))) ** 2
    f = np.fft.rfftfreq(x.size, 1.0 / fps)
    m = f > fmin
    return float(f[m][np.argmax(P[m])])


def fix_sign(x):
    """Deterministic polarity convention: make the first difference right-skewed.

    Breathing is temporally asymmetric (one phase is sharper than the other), so
    the sign of skew(diff(x)) is a stable, reproducible convention.  It is a
    CONVENTION, not a claim about which phase is inspiration -- set SIGN_FLIP to
    invert it once you know the true polarity for your prep.  What matters for
    comparability is that it is decided once, on pooled baseline data, and then
    applied unchanged to every run.
    """
    d = np.diff(np.asarray(x, float))
    d = d - d.mean()
    s = (d ** 3).mean() / max((d ** 2).mean() ** 1.5, 1e-30)
    return -1.0 if s < 0 else 1.0


def main():
    global STRIDE            # resolved from the measured fps below when "auto"
    root = sys.argv[1] if len(sys.argv) > 1 else ROOT_DIR
    if not os.path.isdir(root):
        log("ROOT_DIR not found: %s" % root)
        return 1

    runs = find_runs(root)
    if not runs:
        log("no run folders with a single .avi under %s" % root)
        return 1
    log("breath_fixedmetric_extract: %d run(s) under %s" % (len(runs), root))

    # ROI: the explicit parameter wins; otherwise pick up whatever
    # breath_fixedmetric_roi.py last wrote for this folder.
    roi = ROI_XYWH
    roi_file = os.path.join(root, "breath_roi.mat")
    if roi is None and os.path.isfile(roi_file):
        try:
            from scipy.io import loadmat
            roi = tuple(int(v) for v in loadmat(roi_file)["roi_xywh"].ravel()[:4])
            log("ROI loaded from breath_roi.mat: x=%d y=%d w=%d h=%d" % roi)
        except Exception as e:
            log("could not read %s (%s) -- using the full frame" % (roi_file, e))
            roi = None
    if roi is None:
        log("ROI: full frame (run breath_fixedmetric_roi.py to pick one)")

    cache = os.path.join(root, ".breathcache")
    # A cache built under a different ROI describes a different piece of chest.
    # Reusing it would silently mix two spatial weightings, so the ROI is
    # stamped into the cache and any change wipes it.
    stamp = os.path.join(cache, "roi_stamp.npy")
    if os.path.isdir(cache) and not FORCE:
        prev = None
        if os.path.isfile(stamp):
            v = np.load(stamp)
            prev = None if v.size == 0 else tuple(int(x) for x in v)
        if prev != (tuple(roi) if roi is not None else None):
            log("ROI changed (%s -> %s) -- clearing cache" % (prev, roi))
            shutil.rmtree(cache)
    if FORCE and os.path.isdir(cache):
        shutil.rmtree(cache)
    os.makedirs(cache, exist_ok=True)
    np.save(stamp, np.array([] if roi is None else list(roi), dtype=np.int64))

    base_nums = set(BASIS_RUNS) if BASIS_RUNS else {r["num"] for r in runs}
    for r in runs:
        r["is_base"] = r["num"] in base_nums
    if not any(r["is_base"] for r in runs):
        log("none of BASIS_RUNS matched the folders present -- aborting")
        return 1

    # ---- phase 1: reference images + fps ---------------------------------
    log("\n[1/4] reference images + fps  (%d worker%s)"
        % (N_WORKERS, "" if N_WORKERS == 1 else "s"))
    for r in runs:
        r["mf_path"] = os.path.join(cache, r["name"] + "_mean.npy")
    jobs = [(r["avi"], REF_NFRAMES, r["mf_path"]) for r in runs]
    errs = {}
    if N_WORKERS > 1:
        with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
            for out, e in ex.map(_worker_ref, jobs):
                if e:
                    errs[out] = e
    else:
        for j in jobs:
            out, e = _worker_ref(j)
            if e:
                errs[out] = e

    good, bad = [], []
    for r in runs:
        e = errs.get(r["mf_path"])
        if e or not os.path.isfile(r["mf_path"]):
            # A truncated / aborted recording must not kill a multi-hour job.
            bad.append((r["name"], e or "no reference image produced"))
            continue
        r["mean"] = np.load(r["mf_path"])
        r["fps"], r["fps_ok"] = read_fps(r["folder"], FPS_FALLBACK)
        good.append(r)
        log("  %-32s %.4f fps%s  %s" % (
            r["name"], r["fps"], "" if r["fps_ok"] else " (FALLBACK)",
            "basis" if r["is_base"] else "-"))

    if bad:
        log("")
        for nm, e in bad:
            log("  SKIPPED (unreadable, likely truncated): %s  [%s]" % (nm, e))
    runs = good
    if not runs:
        log("no readable runs -- aborting")
        return 1

    # mtime should sit one run-duration after the folder-name start time.  A run
    # that breaks that pattern was cut short even if it still decodes.
    off = [r["mtime"] - r["t_name"] for r in runs if r["t_name"]]
    if off:
        med = float(np.median(off))
        odd = [r["name"] for r in runs
               if r["t_name"] and abs((r["mtime"] - r["t_name"]) - med) > 5.0]
        log("\n  mtime sits %.1f s after the folder-name start time (median)" % med)
        for nm in odd:
            log("  WARNING: %s breaks that pattern -- likely a short/aborted run" % nm)

    fps_cam = float(np.median([r["fps"] for r in runs]))
    fps_spread = float(np.max([r["fps"] for r in runs]) - np.min([r["fps"] for r in runs]))
    nFall = sum(1 for r in runs if not r["fps_ok"])
    if nFall:
        log("\n  WARNING: %d of %d runs had NO usable fps source and fell back to %.1f."
            % (nFall, len(runs), FPS_FALLBACK))
        log("           Check that each run folder has 2P metadata (*_ch1_meta.mat) or a")
        log("           *timestamps.csv. A wrong fps silently rescales every trace.")
    fps_min = float(np.min([r["fps"] for r in runs]))
    if fps_spread > 0.5:
        # NOT a fault: the breath cam is 2P-TRIGGERED, so each run inherits its own
        # scan rate, and runs with different scan configs legitimately differ. What
        # would be a fault is pretending they share one rate -- the per-run value is
        # therefore carried through to the output as fps_run, and the exported
        # per-recording trace uses it. fps_cam (the median) is kept only for the
        # cross-run comparability arrays, which share one sample grid by design.
        log("  runs differ in fps by %.2f Hz (%.1f-%.1f). This is expected for a"
            % (fps_spread, fps_min, float(np.max([r["fps"] for r in runs]))))
        log("  2P-triggered camera with mixed scan configs; per-run fps is saved as fps_run.")

    # Resolve STRIDE against the fps we actually measured, not the one this
    # pipeline happened to be written for. Use the SLOWEST run: one global stride
    # has to satisfy Nyquist for every run, not just the typical one.
    stride_req = STRIDE
    if isinstance(STRIDE, str):
        s = 1
        while s + 1 <= STRIDE_MAX and (fps_min / (s + 1)) > NYQ_MARGIN * BAND[1]:
            s += 1
        STRIDE = s
        log("\n  STRIDE auto -> %d  (slowest run %.3f fps / %d = %.2f Hz, need > %.1f Hz for a %.1f Hz band)"
            % (STRIDE, fps_min, STRIDE, fps_min / STRIDE, NYQ_MARGIN * BAND[1], BAND[1]))
    else:
        log("\n  STRIDE fixed at %d by the user" % STRIDE)

    nyq_ok = (fps_min / STRIDE) > 2.0 * BAND[1]
    if not nyq_ok:
        log("\nABORT: STRIDE=%d gives %.2f Hz sampling, which is below the"
            % (STRIDE, fps_cam / STRIDE))
        log("       %.1f Hz needed for a %.1f Hz band ceiling. Breathing would alias."
            % (2 * BAND[1], BAND[1]))
        if not isinstance(stride_req, str):
            log("       Set STRIDE = \"auto\" to have it chosen from the measured fps.")
        return 1

    # Every trace lives on the decimated grid, so `fps` from here on is the
    # EFFECTIVE rate.  It is saved as such, because everything downstream
    # (bandpass, breath rate, time axis) must use the rate the samples have.
    fps = fps_cam / STRIDE
    ref_run = next(r for r in runs if r["is_base"])
    ref_mean = ref_run["mean"]
    log("  camera %.4f fps / stride %d -> %.4f Hz effective; reference = %s"
        % (fps_cam, STRIDE, fps, ref_run["name"]))

    # ---- phase 2: registration ------------------------------------------
    log("\n[2/4] registering run means to reference")
    for r in runs:
        r["shift"] = reg_shift(r["mean"], ref_mean)
    sh = np.array([r["shift"] for r in runs])
    log("  drift range: dy %+d..%+d px, dx %+d..%+d px"
        % (sh[:, 0].min(), sh[:, 0].max(), sh[:, 1].min(), sh[:, 1].max()))

    # ---- phase 3: displacement + cached cubes ---------------------------
    log("\n[3/4] per-frame displacement (the slow pass, %d worker%s, stride %d)"
        % (N_WORKERS, "" if N_WORKERS == 1 else "s", STRIDE))
    todo = []
    for r in runs:
        # STRIDE and BIN are part of the cache KEY, not just the contents. A cube
        # decimated at one stride is a different time series from the same video at
        # another, and the old flat name silently reused it: re-running 260728 after
        # a stride change left some runs at 1500 samples and others at 3000, in one
        # .mat, with nothing to indicate it. Different key -> clean re-decode.
        tag = "_s%d_b%d" % (STRIDE, BIN)
        r["d_path"] = os.path.join(cache, r["name"] + tag + "_D.npy")
        r["cube"] = os.path.join(cache, r["name"] + tag + "_cube.npy")
        if os.path.isfile(r["d_path"]) and os.path.isfile(r["cube"]):
            log("  %-32s cached" % r["name"])
            continue
        todo.append((r["avi"], r["mf_path"], roi, r["shift"], BIN, STRIDE,
                     r["d_path"], r["cube"]))

    if todo:
        t0 = time.time()
        done = 0
        derr = {}
        if N_WORKERS > 1:
            with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
                futs = {ex.submit(_worker_disp, j): j for j in todo}
                for f in as_completed(futs):
                    dp, n, sd, e = f.result()
                    done += 1
                    nm = os.path.basename(dp).replace("_D.npy", "")
                    if e:
                        derr[dp] = e
                        log("  %-32s FAILED: %s" % (nm, e))
                    else:
                        log("  %-32s %5d fr  |dy| sd=%.4f px  [%d/%d]"
                            % (nm, n, sd, done, len(todo)))
        else:
            for j in todo:
                dp, n, sd, e = _worker_disp(j)
                done += 1
                nm = os.path.basename(dp).replace("_D.npy", "")
                if e:
                    derr[dp] = e
                    log("  %-32s FAILED: %s" % (nm, e))
                else:
                    log("  %-32s %5d fr  |dy| sd=%.4f px  [%d/%d]"
                        % (nm, n, sd, done, len(todo)))
        el = time.time() - t0
        log("  %d run(s) in %.1f min (%.1f s/run wall)" % (len(todo), el / 60, el / len(todo)))

    keep2 = []
    for r in runs:
        if not (os.path.isfile(r["d_path"]) and os.path.isfile(r["cube"])):
            log("  SKIPPED (no displacement produced): %s" % r["name"])
            continue
        r["D"] = np.load(r["d_path"])
        keep2.append(r)
    runs = keep2
    if not runs:
        log("no runs produced displacement -- aborting")
        return 1

    # A truncated recording still decodes frame-by-frame even when its container
    # index is broken, so it survives phase 1.  It must be dropped HERE, because
    # the common length below is a min(): one 32 s run among 60 s runs would
    # silently clip every other run to 32 s and nothing downstream would say so.
    #
    # The question is "did this run deliver the frames IT was asked for", NOT
    # "is this run as long as the others".  Comparing each run to the cohort
    # MEDIAN -- what this did until 2026-08-07 -- silently assumed every run in a
    # session shares one intended duration.  On a session that deliberately
    # mixes acquisition lengths the median lands BETWEEN them, and every run of
    # the shorter kind is discarded as "short": on 260807_sst-soma-g8s, six
    # complete 3000-frame (100 s) runs were thrown away against a 154.7 s median
    # set by the six 6000-frame runs.  Nothing was wrong with them.
    #
    # So the test is now per-run completeness against that run's own intended
    # count.  len(D) lives on the strided grid, so the expected count is strided
    # the same way (ceil, matching pass_displacement's `i % stride == 0`).
    lens = np.array([len(r["D"]) for r in runs], float)
    keep, unchecked = [], []
    for r, L in zip(runs, lens):
        n_exp, src = expected_frames(r["folder"])
        r["n_expected"] = float(n_exp) if n_exp else np.nan
        if n_exp is None:
            unchecked.append(r["name"])
            keep.append(True)
            continue
        exp_dec = int(np.ceil(n_exp / float(STRIDE)))
        r["frac_complete"] = L / max(exp_dec, 1)
        ok = L >= SHORT_RUN_FRAC * exp_dec
        keep.append(ok)
        if not ok:
            log("\n  DROPPED (truncated: %d of %d expected fr = %.0f%%, %s): %s"
                % (int(L), exp_dec, 100 * L / max(exp_dec, 1), src, r["name"]))
    if unchecked:
        log("")
        for nm in unchecked:
            log("  NOT CHECKED for truncation (no *_meta.mat and no timestamps.csv): %s"
                % nm)
    runs = [r for r, k in zip(runs, keep) if k]
    if not runs:
        log("every run failed the completeness test -- aborting")
        return 1

    # Runs of different intended lengths are legitimate, but the common window
    # below is still a min(), so say out loud when one run is about to set a
    # much shorter comparability window than the rest.
    lens = np.array([len(r["D"]) for r in runs], float)
    if lens.min() < 0.5 * lens.max():
        short = runs[int(np.argmin(lens))]
        log("\n  NOTE the common window is set by %s (%.0f s) while the longest run is"
            % (short["name"], lens.min() / fps))
        log("       %.0f s. That is fine if the durations were intended; if that run was"
            % (lens.max() / fps))
        log("       cut short, exclude its folder and re-run to widen the window.")

    # T = the COMMON window, set by the shortest surviving run.  It is what makes
    # the amplitude metric comparable between runs, and it is the right length for
    # freezing the axis and the basis, where every run must contribute equally.
    #
    # It is NOT the right length for the exported breath TRACE.  Runs here differ
    # legitimately in duration (100 s vs 200 s recordings), so the shortest one
    # used to clip every other run's trace to its length -- on 260728_vglut2 that
    # meant a 50 s window on 200 s recordings, i.e. 25% of the data, which then
    # silently truncated every downstream phase analysis via T = min(breath, Ca).
    #
    # So each run is ALSO projected over its own FULL length below, into the
    # *_full arrays.  The truncated arrays keep their exact previous numerics.
    T = min(len(r["D"]) for r in runs)
    nR = len(runs)
    run_len = np.array([len(r["D"]) for r in runs], int)
    Tmax = int(run_len.max())
    log("  common length T = %d frames (%.2f s)  [used for axis/basis + comparability]"
        % (T, T / fps))
    log("  run lengths %d-%d frames (%.1f-%.1f s); full-length traces kept to Tmax = %d"
        % (run_len.min(), Tmax, run_len.min() / fps, Tmax / fps, Tmax))
    if Tmax > T:
        log("  NOTE %d of %d runs are longer than the common window and would have been"
            % (int((run_len > T).sum()), nR))
        log("       clipped; their full traces are in fb_full / disp_full (run_len marks the end)")

    # ---- phase 4: frozen axis, frozen basis, projections ----------------
    log("\n[4/4] freezing axis + basis on pooled baseline, projecting all runs")

    # Motion axis, frozen on pooled baseline displacements.  Each run is
    # mean-centred first so a between-run offset cannot tilt the axis.
    pool = np.vstack([r["D"][:T] - r["D"][:T].mean(0) for r in runs if r["is_base"]])
    _, sv, vt = np.linalg.svd(pool, full_matrices=False)
    axis = vt[0]
    axis_var = float(sv[0] ** 2 / (sv ** 2).sum())
    log("  motion axis (dx,dy) = (%+.3f, %+.3f), %.1f%% of displacement variance"
        % (axis[0], axis[1], 100 * axis_var))

    # Frozen spatial basis from pooled baseline frames.  BASIS_EVERY thins the
    # cached cube further: the basis only needs the spatial covariance well
    # sampled, not a Nyquist-valid time series.
    # Each run is centred on ITS OWN mean before pooling.  Centring on the
    # pooled mean instead lets between-run structure -- illumination drift,
    # residual sub-pixel misregistration -- into the basis, and on real data it
    # dominates: the pooled-mean mode came out near-orthogonal to the per-run
    # mode (|dot| = 0.32) and correlated only 0.31 with rigid displacement
    # versus 0.54 for this one.  The basis must describe within-run dynamics,
    # not which run a frame came from.
    base_runs = [r for r in runs if r["is_base"]]
    per_run = max(1, int(MAX_BASIS_FR / len(base_runs)))
    chunks, Hs, Ws = [], None, None
    pooled_sum, pooled_n = None, 0
    for r in base_runs:
        c = np.load(r["cube"])[::BASIS_EVERY]
        if Hs is None:
            Hs, Ws = c.shape[1], c.shape[2]
        c = c[:per_run].reshape(-1, Hs * Ws).astype(np.float32)
        pooled_sum = c.sum(0) if pooled_sum is None else pooled_sum + c.sum(0)
        pooled_n += len(c)
        chunks.append(c - c.mean(0))
    X = np.concatenate(chunks, 0)
    del chunks
    basis_mean = pooled_sum / pooled_n      # kept for display only, not centring
    log("  basis matrix %d frames x %d binned px (bin=%d, %d baseline runs, per-run centred)"
        % (X.shape[0], X.shape[1], BIN, len(base_runs)))
    u1, _, vf1 = top_pc(X)
    var_frac = np.array([vf1])
    log("  frozen mode u1 explains %.1f%% of pooled baseline pixel variance" % (100 * vf1))
    del X

    # Polarity, decided once on the pooled baseline displacement projection.
    sgn_disp = fix_sign(pool @ axis)
    if SIGN_FLIP:
        sgn_disp = -sgn_disp
    axis = axis * sgn_disp
    log("  polarity: axis sign %+d (SIGN_FLIP=%s)" % (int(sgn_disp), SIGN_FLIP))

    disp = np.zeros((T, nR))
    fb = np.zeros((T, nR))
    pv = np.zeros((T, nR))
    diff_imgs = np.zeros((Hs, Ws, nR), np.float32)
    D_all = np.zeros((T, 2, nR))
    pv_corr = np.zeros(nR)
    pv_sign = np.zeros(nR)
    # Full-length traces, NaN-padded to Tmax. NaN rather than zero so a consumer
    # that forgets run_len gets an obvious hole instead of a silent flat tail.
    disp_full = np.full((Tmax, nR), np.nan)
    fb_full = np.full((Tmax, nR), np.nan)

    # Displacement and cube share STRIDE, so all three metrics are already on the
    # same grid at the same length -- no interpolation anywhere, which also
    # removes the ~1-2% amplitude error that resampling used to introduce.
    for j, r in enumerate(runs):
        Lj = int(run_len[j])
        D = r["D"][:T].copy()
        D_all[:, :, j] = D
        disp[:, j] = (D - D.mean(0)) @ axis

        cube_all = np.load(r["cube"]).astype(np.float32)

        # ---- full-length projections, on the SAME frozen axis and basis -------
        # Computed separately rather than by slicing the truncated result, because
        # each is mean-centred over its own window and the two means differ. The
        # truncated arrays therefore keep their exact previous values.
        Df = r["D"][:Lj]
        disp_full[:Lj, j] = (Df - Df.mean(0)) @ axis
        Mf = cube_all[:Lj].reshape(Lj, -1)
        Mf = Mf - Mf.mean(0)
        fbf = Mf @ u1
        fb_full[:Lj, j] = fbf - fbf.mean()
        del Mf, fbf

        cube = cube_all[:T]

        # Inspiration-minus-expiration difference map, for the peak GUI's image
        # panel (breathing_peak_gui_pc1.compute_diff_map).  That function needs
        # either a ready-made diffImg or the full U/V/sv factors; giving it the
        # image directly is both simpler and more faithful, since this is the
        # exact top-minus-bottom frame difference rather than a rank-k SVD
        # reconstruction of it.  Frames are ranked by `disp`, the physical
        # trace, so the map means "where the chest moved between the extremes
        # of the breath" whichever metric is being viewed.
        s = disp[:, j]
        nf = max(1, min(100, len(s) // 2))
        ordr = np.argsort(s)
        diff_imgs[:, :, j] = (cube[ordr[-nf:]].mean(0) - cube[ordr[:nf]].mean(0))

        # Per-run centring here too, to match how the basis was built -- which
        # also makes the frozen-basis trace immune to between-run illumination
        # drift, the same way the displacement metric already is.
        M = cube.reshape(len(cube), -1)
        M = M - M.mean(0)

        fb_sub = M @ u1
        fb[:, j] = fb_sub - fb_sub.mean()

        # Per-video PC1: its own basis, hence its own arbitrary sign and scale.
        Mc = M - M.mean(0)
        v_pv, pv_sub, _ = top_pc(Mc)
        pv[:, j] = pv_sub - pv_sub.mean()

        # How far this run's own PC1 has rotated away from the frozen mode.
        dot = float(np.dot(v_pv, u1))
        pv_corr[j] = abs(dot)
        pv_sign[j] = float(np.sign(dot))
        del cube, cube_all, M, Mc

    # Polarity is anchored to `disp`, the only metric whose direction is
    # physically defined (image coordinates).  Applying the waveform heuristic
    # to each metric independently is NOT enough: on this data it handed fb the
    # opposite polarity to disp in all 21 runs, so "decide the sign once" has to
    # mean once against a shared reference, not once per metric.
    bmask = np.array([r["is_base"] for r in runs])
    bidx = np.where(bmask)[0]

    naive_fb = fix_sign(fb[:, bmask].ravel())
    rsum = float(sum(np.corrcoef(disp[:, j], fb[:, j])[0, 1] for j in bidx))
    s_fb = -1.0 if rsum < 0 else 1.0
    if SIGN_FLIP:
        s_fb = -s_fb
    fb *= s_fb
    fb_full *= s_fb          # same frozen polarity; disp_full already carries it
                             # through the signed axis, which was fixed before the loop
    fb_naive_disagreed = bool(np.sign(naive_fb) != np.sign(s_fb))

    # pv's sign is per-video by nature, so it is anchored per video.  Count how
    # often the unaided heuristic would have got it wrong -- that count is the
    # cost of having no shared reference at all.
    pv_naive_bad = 0
    for j in range(nR):
        naive = fix_sign(pv[:, j])
        s = -1.0 if np.corrcoef(disp[:, j], pv[:, j])[0, 1] < 0 else 1.0
        if np.sign(naive) != np.sign(s):
            pv_naive_bad += 1
        pv[:, j] *= s * (-1.0 if SIGN_FLIP else 1.0)

    log("  polarity anchored to disp: fb global sign %+d (waveform heuristic %s)"
        % (int(s_fb), "disagreed" if fb_naive_disagreed else "agreed"))
    log("  per-video polarity: waveform heuristic would have flipped %d/%d runs wrongly"
        % (pv_naive_bad, nR))

    log("\n  per-run QC (basis runs marked *)")
    log("  %-30s %8s %8s %8s %8s" % ("run", "dom_Hz", "bandfrac", "sd_px", "|u1.pv1|"))
    for j, r in enumerate(runs):
        log("  %-30s %8.3f %8.3f %8.4f %8.3f%s" % (
            r["name"], dom_freq(disp[:, j], fps), band_frac(disp[:, j], fps, BAND),
            disp[:, j].std(), pv_corr[j], " *" if r["is_base"] else ""))

    out = os.path.join(root, "breath_fixedmetric.mat")
    savemat(out, dict(
        root_dir=root,
        run_names=np.array([r["name"] for r in runs], dtype=object),
        run_num=np.array([r["num"] for r in runs], float),
        # Kept under the historical name; it flags the BASIS runs.  The
        # analysis may override which runs define the amplitude reference.
        is_baseline=np.array([1.0 if r["is_base"] else 0.0 for r in runs]),
        # Wall clock, for a real time axis instead of a run-index axis.
        # run_mtime = video last-write (end of recording), epoch seconds.
        # run_tname = folder-name start time, epoch seconds (NaN if unparsable).
        # run_hours = hours elapsed from the first run, derived from mtime.
        run_mtime=np.array([r["mtime"] for r in runs], float),
        run_tname=np.array([r["t_name"] if r["t_name"] else np.nan
                            for r in runs], float),
        run_hours=np.array([(r["mtime"] - runs[0]["mtime"]) / 3600.0
                            for r in runs], float),
        # Intended camera frame count per run (2P framesPerSlice x numSlices, or
        # timestamps rows; NaN if neither existed).  Runs legitimately differ in
        # duration, so this is what "complete" was judged against.
        run_expected=np.array([r.get("n_expected", np.nan) for r in runs], float),
        fps=fps, T=float(T), band=np.array(BAND, float),
        # PER-RUN effective rate. The camera is 2P-triggered, so runs with
        # different scan configs have different rates and no single number is
        # right for all of them. Anything TEMPORAL must use fps_run(j); `fps`
        # above is the median and exists only for the shared-grid arrays.
        fps_run=np.array([r["fps"] / STRIDE for r in runs], float),
        fps_cam_run=np.array([r["fps"] for r in runs], float),
        fps_ok_run=np.array([1.0 if r["fps_ok"] else 0.0 for r in runs], float),
        bin_factor=float(BIN), sub_every=float(STRIDE),
        fps_camera=fps_cam, stride=float(STRIDE),
        disp=disp, fb=fb, pv=pv, D_all=D_all,
        # Full-length traces, NaN-padded to Tmax, with run_len marking each end.
        # Use these for anything TEMPORAL (breath phase, peak/trough marking);
        # use disp/fb/pv above for anything COMPARING runs, since those share the
        # common window T that makes amplitudes comparable.
        disp_full=disp_full, fb_full=fb_full, run_len=run_len.astype(float),
        Tmax=float(Tmax),
        motion_axis=axis, axis_var=axis_var,
        # Rebuilt from the FINAL run list: `sh` was formed in phase 2, before
        # short/unreadable runs were dropped, so using it directly would emit a
        # column of a different length from every other per-run field.
        reg_shift=np.array([r["shift"] for r in runs], float),
        u1_img=u1.reshape(Hs, Ws), basis_mean_img=basis_mean.reshape(Hs, Ws),
        diff_imgs=diff_imgs,
        basis_var_frac=var_frac, pv_align=pv_corr, pv_align_sign=pv_sign,
        ref_mean_img=ref_mean, sign_flip=float(SIGN_FLIP),
        fb_sign=s_fb, pv_naive_bad=float(pv_naive_bad),
        fb_naive_disagreed=float(fb_naive_disagreed),
    ), do_compression=True)
    log("\nwrote %s" % out)
    log("next: run breath_fixedmetric_analyze.m in MATLAB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
