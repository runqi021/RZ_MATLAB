---
name: project-breath-sign-is-analysis-immune
description: Breath trace polarity flips between sessions but every downstream analysis is sign-immune; SST 260806/260807 use fixedmetric not per-video SVD PC1
metadata: 
  node_type: memory
  type: project
  originSessionId: e48a8451-7533-46a9-bac3-1424231c2bc0
  modified: 2026-08-15T20:41:51.364Z
---

Audited 2026-08-15 because the repr export video looked misaligned with the breath trace.

**Breath trace polarity is NOT stable across sessions.** SST 260806 has temporal
skewness ~-1.7 on all 12 FOVs; 260807 has ~+1.6 on all 12. Onset-triggered averages
are mirror images (corr -0.835, +0.835 when flipped). Within a session the sign IS
consistent, because the fixedmetric pipeline uses one session-level baseline basis.

**Everything downstream is immune to this**, so it has never corrupted a result:
- all analyses read `breath_peak_pc1.mat` -> `insp_onset_idx` / `insp_onsets_train`
  (event times), never the signed `breathTrace`
- coherence phase uses `ref = cos(piecewise_phase_local(peak_idx, foot_idx))` built
  from event indices — `coherence_polar_general_260727.m:137-143`
- the raw waveform `bw` only feeds the PSD for band selection and coherence
  MAGNITUDE, both sign-blind
- the peak detector anchors the same physiological deflection regardless of sign
  (260806 trough at t=0, 260807 peak at t=0)

**Only risk:** overlaying/pooling the breath WAVEFORM itself across sessions (e.g. an
average breath-trace figure) needs a per-session sign normalisation at plot time.

**SST 260806 + 260807 use the fixedmetric pipeline, not per-video SVD PC1.**
Tell-tale in `breath_pc1.mat`: `src_metric='fb'`, `is_baseline=1`, and `varExp`
IDENTICAL across every FOV in the session (89.0% / 79.6%) — the signature of one
shared basis. See [[project_breath_fixedmetric]] and [[project_breath_svd_pc1]].
Breath rates 1.2-2.0 Hz across all 24 FOVs (sane; cf. 0.4 Hz Vgat in §4 of
`analysis_260727\PROGRESS_260731.md`).

`repr_export_RZ.m` computes its OWN ad-hoc SVD and caches it as
`breath_pc1_crop_*.mat` / `breath_pc1_fullframe.mat`. **No analysis ever reads those**,
and it can never clobber `breath_pc1.mat`. See [[project_repr_export_video]].
