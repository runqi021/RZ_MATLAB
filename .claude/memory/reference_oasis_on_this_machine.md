---
name: reference_oasis_on_this_machine
description: "OASIS deconvolution on this machine: pyenv is impossible in R2021b and two OASIS entry points crash — use oasisAR1 via the ca_oasis_run.py subprocess"
metadata: 
  node_type: memory
  type: reference
  originSessionId: a1713ea0-2ee8-4943-8d90-ba971a8be969
  modified: 2026-07-25T02:27:21.563Z
---

**`helper.oasis_deconv_and_dff_AR1` cannot run on this machine at all.** Anything calling it — including the OASIS button in `calcium_spike_gui.m` — is silently unusable. Two independent blockers:

1. **pyenv version ceiling.** It bridges through MATLAB's `pyenv`, which in R2021b supports CPython 3.7-3.9 only. Both conda envs holding `oasis-deconv 0.2.1` (`cellpose-gpu`, `oasis`) are **Python 3.10.19**. Not fixable without building a py39 env. Same class of problem as [[project_thermal_ats_breathing]]'s R2021b/py3.13 issue.
2. **Broken LAPACK in those envs.** `numpy.linalg.lstsq` and `np.corrcoef` hard-crash with Windows fatal exception **`0xc06d007f`** (surfaces in MATLAB as `system` status **-1066598273**, and the captured stderr shows only the unrelated cvxpy warning, so it looks like nothing went wrong).

**What crashes vs what works** (numpy 2.2.6, scipy 1.14/1.15, both envs):
- `oasis.functions.deconvolve` — CRASHES. Routes to `estimate_time_constant` -> `numpy.linalg.lstsq`.
- `constrained_oasisAR1` — CRASHES on real dF/F traces, regardless of offsetting. Fine on synthetic AR(1), so a synthetic smoke test will pass and mislead you.
- `oasisAR1(y, g, lam)` — **WORKS**. The plain Cython core.
- `GetSn(y)` — WORKS (Welch PSD, no LAPACK).

**The breakage is wider than LAPACK — plain BLAS is broken too** (found 2026-07-24). In `cellpose-gpu`, even `A.T @ A` on a trivial 900x2 float64 array dies, as does `np.linalg.svd` on a 200x50. The process exits silently with no traceback and no stderr, so a piped command just returns empty output and looks like it produced nothing. If a Python script in this env goes quiet for no reason, suspect a matmul before anything else.

**Use `C:\Users\Admin\.conda\envs\dlc310\python.exe` for any numerical work.** Verified healthy: matmul, `svd`, `eigh`, plus cv2 4.11. `torch113` and `flir` are also healthy. Only `cellpose-gpu` is poisoned — keep it for Cellpose/h5py reads and nothing else.

**The working route:** `ca_oasis_run.py` (repo root) as a **subprocess**, not pyenv. MATLAB writes `-v7` (NOT `-v7.3` — `scipy.io.loadmat` cannot read HDF5 mat), shells out, reads the result back. Call it with the `oasis` env python. Plain `system()` from MATLAB is fine — no PATH scrubbing needed, MATLAB's DLLs are not the problem.

```
[~,~] = system(sprintf('"%s" "%s" --in "%s" --out "%s" --g %.6f --lam 0', ...
                       pyExe, script, fin, fout, g));
```

Pass `g` explicitly (from [[project_ca_spike_detector_260722]]'s kernel fit) so nothing needs estimating — which is what dodges the crashing code path. `lam = 0` avoids pre-censoring small events when a threshold sweep follows.

**Two conventions that bite when consuming the output:**
- OASIS marks the transient **ONSET**; peak-clicked labels sit ~**4 frames** later at 30 fps. Measured, not modelled — the biexponential kernel's argmax predicts 5-9 frames and over-corrects.
- The spike train is mostly zeros, so **MAD is degenerate** (median = 0). Scale by the returned per-trace `sn`, not by MAD.
