"""export_masks.py -- dump cpSAM maskL to .npy so the extraction can run elsewhere.

RUN WITH cellpose-gpu (the only env with h5py):
    %USERPROFILE%\\.conda\\envs\\cellpose-gpu\\python.exe export_masks.py

cellpose-gpu is the only env that can READ the v7.3 .mat, but its BLAS/LAPACK is
broken (matmul, lstsq and corrcoef all die with Windows 0xc06d007f), so it cannot
run the extraction itself. dlc310 has working linear algebra and tifffile but no
h5py. This bridges them: masks out as .npy, extraction in dlc310.

Writes <session>/roi_maskL.npy alongside each cpSAM_output.mat.
"""
import os
import sys

import h5py
import numpy as np

# root may be passed on the command line, e.g. D:\260721_Sert_soma_G8s\phys
ROOT = sys.argv[1] if len(sys.argv) > 1 else \
    r"D:\Ventral_surface_summary"


def main():
    n = 0
    for dp, _dn, fn in os.walk(ROOT):
        if os.sep + "_" in dp:
            continue
        sam = sorted(f for f in fn if f.endswith("cpSAM_output.mat"))
        if not sam or not any(f.endswith(".tif") for f in fn):
            continue
        # pair with the most-registered movie, same as the extractor does
        src = os.path.join(dp, sam[-1])
        try:
            with h5py.File(src, "r") as f:
                if "maskL" not in f:
                    print(f"  no maskL: {os.path.relpath(dp, ROOT)}")
                    continue
                mask = np.array(f["maskL"]).T      # h5py transposes MATLAB
        except Exception as e:
            print(f"  FAILED {os.path.relpath(dp, ROOT)}: {e}")
            continue
        np.save(os.path.join(dp, "roi_maskL.npy"), mask.astype(np.int32))
        lbl = np.unique(mask)
        print(f"  {os.path.relpath(dp, ROOT)[:56]:56s} {mask.shape}  "
              f"{len(lbl)-1:4d} ROI  <- {sam[-1][-28:]}")
        n += 1
    print(f"\n{n} mask(s) exported")


if __name__ == "__main__":
    main()
