"""cascade_fix_discrete.py -- recompute discrete spike times from saved rates.

RUN WITH dlc310:
    %USERPROFILE%\\.conda\\envs\\dlc310\\python.exe cascade_fix_discrete.py

The first run stored the wrong half of infer_discrete_spikes' return value
(approximations, i.e. the reconstructed rate, instead of the spike times). The
network inference itself is fine and expensive, so this redoes only the cheap
discrete step in place, reusing the saved spike_prob.
"""
import os
import sys
import time

import numpy as np

CASCADE_DIR = os.path.join(os.environ.get("TEMP", "."), "CascadeTorch")
ROOT = r"D:\Ventral_surface_summary"
PRED_DIR = os.path.join(ROOT, "_cascade_260722", "predictions")
MODEL = "GC8s_EXC_30Hz_smoothing25ms_high_noise"

sys.path.insert(0, CASCADE_DIR)
from cascade2p import utils_discrete_spikes as uds   # noqa: E402


def main():
    model_folder = os.path.join(CASCADE_DIR, "Pretrained_models")
    files = sorted(f for f in os.listdir(PRED_DIR) if f.endswith(".npz"))
    print(f"{len(files)} prediction file(s)")
    t0 = time.time()
    tot_roi = tot_spk = 0

    for i, fn in enumerate(files, 1):
        p = os.path.join(PRED_DIR, fn)
        P = np.load(p, allow_pickle=True)
        sp = P["spike_prob"].astype(np.float64)

        _, spikes = uds.infer_discrete_spikes(sp, MODEL,
                                              model_folder=model_folder,
                                              verbosity=0)
        disc = np.array([np.rint(np.asarray(s, dtype=np.float64)).astype(np.int64)
                         for s in spikes], dtype=object)

        n = int(sum(len(s) for s in disc))
        dur_min = sp.shape[1] / 30.0 / 60.0
        tot_roi += sp.shape[0]
        tot_spk += n
        print(f"[{i:2d}/{len(files)}] {str(P['session'])[:50]:50s} "
              f"{sp.shape[0]:4d} ROI  {n:6d} spikes  "
              f"({n/max(sp.shape[0]*dur_min,1e-9):5.1f}/ROI-min)")

        np.savez_compressed(p, spike_prob=P["spike_prob"], noise=P["noise"],
                            session=P["session"], fps=P["fps"], discrete=disc)

    print(f"\n{tot_spk} spikes over {tot_roi} ROIs in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
