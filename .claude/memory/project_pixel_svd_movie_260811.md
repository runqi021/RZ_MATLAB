---
name: project-pixel-svd-movie-260811
description: "Pixel-SVD of the 2P movie — scripts, Gram-matrix method, and the finding that top modes are 10 Hz cardiac that MC does not remove"
metadata: 
  node_type: memory
  type: project
  originSessionId: f145e055-fef7-4364-bc1d-f95f7ce506ec
  modified: 2026-08-11T09:52:35.053Z
---

Pixel-SVD (decompose the calcium MOVIE into spatial+temporal modes), as opposed to
all the `*_svd_*` breath scripts which SVD a behavior AVI — see [[project-breath-svd-pc1]].

**Scripts** (in `SVD\`):
- `Run_SVD_recon_260811.m` — THE current one. Runs a LIST of movie variants
  (raw / mcmc) in one go into `<folderPath>\svd_check_260811\<variant>\`.
  Outputs `SVD_result.mat`, `singular_spectrum.png`, `mode_gallery.png`,
  `recon_k020.mp4` / `recon_k100.mp4` (panels: original | rank-K recon | residual),
  plus `variant_comparison.png` at the root.
- `Run_SVD_mode_psd_260811.m` — companion; PSD of v_k(t) from the saved MATs.
  This is what identifies a mode (cardiac vs breath vs drift). No SVD re-run needed.
- Superseded: `Run_SVD_working_251230.m` (one variant per run, overwrote
  `all_modes\SVD_result.mat`) + `SVD_denoise_recon.m` (snippet, needs
  `SVD_result` already in the workspace, hardcoded path to another dataset).

**Method note — never use `svd(Y,'econ')` on a pixel movie.** Go through the
temporal Gram matrix: `C = Y'*Y` ([T x T], small) → `eig` → V and S=sqrt(lam),
then `U = Y*V./S'`. Exact for the retained modes. 440x440x2970 takes ~30 s
this way; the full econ SVD of a [193600 x 2970] matrix is enormous.

**Trap: `contains(name,'_MC')` is not a safe way to find the MC movie.** In a
processed FOV folder it also matches `*_MC_MC_AVG_for_CP.tif`, a SINGLE-FRAME
projection. Use explicit suffixes + assert the frame count is > 100.

**Frame alignment**: raw is 3000 frames, `_preproc_MC_MC` is 2970 (MC already
tossed 30). Front-drop 30 from raw so both cover the identical frames, otherwise
the spectra are not comparable.

## Finding — `260810_vglut2-soma-g8s\phys\roi1_5x_13lp_3000f_z12_00001`
- Spectrum is shot-noise dominated: tail sits at ~0.034% = 1/T, the white-noise
  floor. Only ~10 modes rise above it. cum@100 modes = 7.5% of variance.
- **Modes 1-2 (~0.9%, 0.8%) are 10.06 Hz CARDIAC pulsation** — red/blue dipole
  eigenimages on cell edges (sub-pixel edge modulation). Modes 5-6 are 9.84 Hz,
  also cardiac. Mode 4 is 5.95 Hz breathing.
- **Rigid+non-rigid MC does NOT touch the cardiac modes.** raw vs mcmc PSDs of
  v1/v2 lie exactly on top of each other; variance is unchanged (0.91→1.00%).
  MC only removed the <2 Hz drift (raw mode 3 = 0.99 Hz wander, gone in mcmc).
  So a 10 Hz within-frame pulsation is out of reach of frame-wise registration —
  suppress it in the trace domain or by rank-removal, not by more MC.
- Rank-20 reconstruction recovers clean somata from an original so photon-starved
  it looks like sparse speckle.
