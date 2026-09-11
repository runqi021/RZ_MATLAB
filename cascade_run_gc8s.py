"""cascade_run_gc8s.py -- automatic GCaMP8s soma event detection with CASCADE.

RUN WITH THE dlc310 ENV (torch 2.9 + CUDA, and its LAPACK works):

    %USERPROFILE%\\.conda\\envs\\dlc310\\python.exe cascade_run_gc8s.py

Why CASCADE rather than the threshold/OASIS detectors: it is a deep network
trained on ground truth from SIMULTANEOUS electrophysiology + 2P calcium
imaging, so the supervision comes from real recorded spikes in other animals
rather than from anyone's clicks on this dataset. The chosen model,

    GC8s_EXC_30Hz_smoothing25ms_high_noise

is trained on DS30/31/32 (GCaMP8f/8m/8s, excitatory somata) resampled to 30 Hz,
which matches this data exactly: GCaMP8s somata imaged at 30 Hz.

CASCADE picks the ensemble member matching each trace's own measured noise
level, so there is no threshold to tune per session. Output is a spike RATE per
frame (spikes/s); discrete events come from the companion inference step.

Inputs : <ROOT>/_cascade_260722/traces/*.npz  (from cascade_export_dff.py)
Outputs: <ROOT>/_cascade_260722/predictions/*.npz  with
         spike_prob [N x T] (rate, NaN-padded at the edges), discrete spike
         times per ROI, and the noise level CASCADE measured for each ROI.
"""
import json
import os
import sys
import time

import numpy as np

CASCADE_DIR = os.path.join(os.environ.get("TEMP", "."), "CascadeTorch")
ROOT = r"D:\Ventral_surface_summary"
TRACE_DIR = os.path.join(ROOT, "_cascade_260722", "traces")
PRED_DIR = os.path.join(ROOT, "_cascade_260722", "predictions")

# Model choice follows the cell type. The EXC models are trained on excitatory
# somatic ground truth; applying them to inhibitory interneurons is a genuine
# mismatch -- Sst and Vgat are GABAergic, and CASCADE ships interneuron models
# (mostly GCaMP8 ground truth) precisely for this case.
MODEL_BY_GROUP = {
    "Vglut2": "GC8s_EXC_30Hz_smoothing25ms_high_noise",   # glutamatergic -> EXC
    "ChAT":   "GC8s_EXC_30Hz_smoothing25ms_high_noise",   # cholinergic, no
                                                          # dedicated model; EXC
                                                          # is the closest fit
    # Interneurons_GC8+ was TRIED and REJECTED for Sst/Vgat: it emits a flat
    # ~4 Hz rate for 60-90% of the recording, uncorrelated with the trace
    # (CV 0.2-0.5, corr ~0). These are the cleanest ROIs here (nu 0.35-0.69),
    # far below that model's nu=2 floor, and it falls back to a prior rate.
    # The EXC model produces sparse, onset-locked output on the same data.
    "Sst":    "GC8s_EXC_30Hz_smoothing25ms_high_noise",
    "Vgat":   "GC8s_EXC_30Hz_smoothing25ms_high_noise",
}
DEFAULT_MODEL = "GC8s_EXC_30Hz_smoothing25ms_high_noise"
ONLY_GROUPS = ["Sst", "Vgat"]   # e.g. ['Sst','Vgat'] to redo just those; [] = all
# Discrete spike inference only makes sense at low firing rates. Interneurons
# run at ~4 Hz here (246 spikes/min), and at 30 Hz sampling with GCaMP8s the
# individual APs in such a train are simply not resolvable -- infer_discrete_spikes
# also becomes intractably slow trying to place ~4000 spikes per ROI. For these
# cell types the calibrated RATE is the meaningful output.
DISCRETE_BY_GROUP = {"Vglut2": True, "ChAT": True, "Sst": True, "Vgat": True}
DO_DISCRETE_DEFAULT = True

sys.path.insert(0, CASCADE_DIR)
from cascade2p import cascade                      # noqa: E402
from cascade2p import utils_discrete_spikes as uds  # noqa: E402


def main():
    os.makedirs(PRED_DIR, exist_ok=True)
    model_folder = os.path.join(CASCADE_DIR, "Pretrained_models")

    with open(os.path.join(TRACE_DIR, "index.json")) as f:
        index = json.load(f)
    if ONLY_GROUPS:
        index = [e for e in index
                 if e["session"].split(os.sep)[0].replace("_test", "") in ONLY_GROUPS]
    print(f"{len(index)} sessions, {sum(e['n_roi'] for e in index)} ROIs\n")

    # make sure every needed model is present before starting
    for m in sorted({MODEL_BY_GROUP.get(
            e["session"].split(os.sep)[0].replace("_test", ""), DEFAULT_MODEL)
            for e in index}):
        if not os.path.isdir(os.path.join(model_folder, m)):
            print(f"downloading model {m} ...")
            cascade.download_model(m, model_folder=model_folder, verbose=1)

    t0 = time.time()
    for i, e in enumerate(index, 1):
        grp = e["session"].split(os.sep)[0].replace("_test", "")
        MODEL = MODEL_BY_GROUP.get(grp, DEFAULT_MODEL)
        out_path = os.path.join(PRED_DIR, e["tag"] + ".npz")
        d = np.load(os.path.join(TRACE_DIR, e["tag"] + ".npz"))
        dff = d["dFF"].astype(np.float64)           # [N x T]

        # CASCADE expects dF/F as a RATIO, not percent. Guard against a trace
        # set that is accidentally in percent -- it would land on the wrong
        # noise-level model and quietly degrade every prediction.
        med_amp = np.nanpercentile(np.abs(dff), 99)
        if med_amp > 20:
            print(f"  !! {e['session']}: |dF/F| p99 = {med_amp:.1f}; expected a "
                  f"ratio (~0.1-3). Skipping -- check the units.")
            continue

        print(f"[{i:2d}/{len(index)}] {e['session'][:44]:44s} "
              f"{dff.shape[0]:4d}x{dff.shape[1]:5d}  {MODEL.split('_30Hz')[0]}", flush=True)

        spike_prob = cascade.predict(MODEL, dff, model_folder=model_folder,
                                     verbosity=0)

        # CASCADE's own per-trace noise metric, useful for QC and for spotting
        # ROIs where no ensemble member really matched
        noise = np.nanmedian(np.abs(np.diff(dff, axis=1)), axis=1) / np.sqrt(1 / 30.0)

        discrete = None
        if DISCRETE_BY_GROUP.get(grp, DO_DISCRETE_DEFAULT):
            try:
                # returns (approximations, spikes) -- the SECOND is the spike
                # times; the first is the reconstructed rate, same length as the
                # trace, which silently looks like one spike per frame
                _, discrete = uds.infer_discrete_spikes(
                    spike_prob, MODEL, model_folder=model_folder, verbosity=0)
                discrete = np.array(
                    [np.asarray(s, dtype=np.int64) for s in discrete],
                    dtype=object)
            except Exception as ex:
                print(f"    discrete inference failed: {type(ex).__name__}: {ex}")

        rate = np.nansum(spike_prob, axis=1) / (dff.shape[1] / 30.0) * 60.0
        print(f"           median inferred rate {np.median(rate):6.1f} spikes/min "
              f"| noise nu {np.median(noise):.2f}")

        np.savez_compressed(out_path, spike_prob=spike_prob.astype(np.float32),
                            noise=noise, session=e["session"], fps=e["fps"],
                            discrete=discrete if discrete is not None else np.array([]),
                            allow_pickle=True)

    print(f"\ndone in {time.time()-t0:.0f}s -> {PRED_DIR}")


if __name__ == "__main__":
    main()
