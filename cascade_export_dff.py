"""cascade_export_dff.py -- pull dF/F out of the v7.3 .mat files for CASCADE.

RUN WITH THE cellpose-gpu ENV (it is the only one carrying h5py):

    %USERPROFILE%\\.conda\\envs\\cellpose-gpu\\python.exe cascade_export_dff.py

The pipeline writes `*_dFF.mat` as MATLAB v7.3, which is HDF5 and unreadable by
scipy.io. h5py is deliberately NOT installed into dlc310 -- that env hosts
DeepLabCut and pulling a second HDF5 binary into it is not worth the risk -- so
the export happens here and CascadeTorch consumes plain .npy in dlc310.

Note h5py returns MATLAB arrays transposed: dFF is [T x N] in MATLAB and arrives
as [N x T], which is exactly the [neurons x time] layout CASCADE wants.

Writes <OUT_DIR>/<session>.npz  with  dFF [N x T], and an index.json.
"""
import json
import os
import sys

import h5py
import numpy as np

# root may be given on the command line, e.g. D:\260721_Sert_soma_G8s\phys
ROOT = sys.argv[1] if len(sys.argv) > 1 else \
    r"D:\Ventral_surface_summary"
OUT_DIR = os.path.join(ROOT, "_cascade_260722", "traces")
FPS = 30.0   # every session reports scanFrameRate 30 via detect_session_fps


def find_sessions(root):
    hits = []
    for dirpath, _dirnames, filenames in os.walk(root):
        if os.sep + "_" in dirpath:          # skip our own output folders
            continue
        cand = [f for f in filenames if f.endswith("_dFF.mat")]
        if cand:
            hits.append((dirpath, sorted(cand)[0]))
    return sorted(hits)


def read_dff(path):
    with h5py.File(path, "r") as f:
        if "dFF" not in f:
            return None
        a = np.array(f["dFF"], dtype=np.float64)   # [N x T] after transpose
    return a


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    sessions = find_sessions(ROOT)
    print(f"{len(sessions)} session(s) with _dFF.mat")

    index = []
    for i, (folder, fname) in enumerate(sessions, 1):
        rel = os.path.relpath(folder, ROOT)
        try:
            dff = read_dff(os.path.join(folder, fname))
        except Exception as e:
            print(f"  [{i:2d}] {rel}: READ FAILED ({type(e).__name__}: {e})")
            continue
        if dff is None or dff.ndim != 2:
            print(f"  [{i:2d}] {rel}: no usable dFF")
            continue

        tag = rel.replace(os.sep, "__")
        out = os.path.join(OUT_DIR, tag + ".npz")
        np.savez_compressed(out, dFF=dff, fps=FPS)
        index.append({"session": rel, "tag": tag, "n_roi": int(dff.shape[0]),
                      "n_frames": int(dff.shape[1]), "fps": FPS})
        print(f"  [{i:2d}] {rel[:56]:56s} {dff.shape[0]:4d} ROI x {dff.shape[1]:5d} fr")

    with open(os.path.join(OUT_DIR, "index.json"), "w") as f:
        json.dump(index, f, indent=1)

    tot = sum(e["n_roi"] for e in index)
    print(f"\nexported {len(index)} sessions, {tot} ROIs -> {OUT_DIR}")


if __name__ == "__main__":
    main()
