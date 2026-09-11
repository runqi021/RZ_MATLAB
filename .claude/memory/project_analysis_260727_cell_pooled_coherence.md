---
name: project-analysis-260727-cell-pooled-coherence
description: analysis_260727/ — PRIMARY is the absolute-time inspiration-triggered PETH with a circular-shift null; coherence and phase are demoted to secondary
metadata: 
  node_type: memory
  type: project
  originSessionId: 1f25dc5a-41b5-40bf-b366-9d77745c5a99
  modified: 2026-07-29T18:11:39.934Z
---

`RZ_MATLAB\analysis_260727\` — breath x calcium, rebuilt around cell identity.
Driver: `run_analysis_260727()` or `run_analysis_260727(rootPath, genotype)`.
Output always `<dataset>\analysis_260727\`. Currently targeting
`D:\260721_Sert_soma_G8s\phys\baseline` (phys was split into baseline/roomAir/reoxygenated).

**PRIMARY = `breath_time\breath_time_peth_260727.m`.** Coherence and phase/Rayleigh
are secondary descriptive measures only. The statistic:
`T_exc = max over search window of (smoothed PETH - its own mean)`, z-scored against
circular-shift shuffles that recompute the SAME max — so the search over latencies is
paid for by the null. No respiratory baseline; each cell vs shuffles of itself.
Significance on the empirical p (rank-based); `M_exc` is only an effect size and is
unreliable below ~40 events.

**"Zero" — the upstream names are a trap.** `foot_idx` comes from
`breath_insp_start_pc1.mat->insp_start_idx` = inspiration ONSET. `peak_idx` comes from
`breath_peak_pc1.mat->insp_onset_idx` = the PEAK (field says "onset", it is not).
267 ms apart. `TRIGGER='onset'|'peak'` → `breath_time\<TRIGGER>\`. Which alignment is
TIGHTER (latency MAD, not M_exc) tells you which landmark a cell locks to. Baseline:
**24 peak / 4 onset / 6 tied**. MAD is quantised so ties are common — report them,
never let a `<` silently assign them.

**Tail plot α lines belong on the NULL curve, not the observed one.** On the observed
survival, y=0.05 means "top 5% of cells I recorded" (a sample quantile) and with N=34
the curve floors at 1/34=0.029 so 0.01/0.005 are unreachable. On the pooled shuffle
null (46k values) they are real critical scores: 1.77 / 3.05 / 4.01 — far above the
normal-theory 1.645 / 2.33 / 2.58, because the null is strongly right-skewed.

**Latency must be SIGNED and to the NEAREST trigger, never the preceding one.**
Preceding-assignment charges an event occurring before a trigger to the previous
breath at nearly a whole cycle, which destroys the precision of pre-inspiratory cells
AND systematically penalises the peak alignment (an event during inspiration has no
preceding *peak* in its own breath). This is not cosmetic: it **reversed** the
onset-vs-peak verdict from 25/34 onset to 28/34 peak. Search range is
[-cycle/2, +cycle/2]; a preferred latency past half a cycle wraps to a negative
latency to the next breath.

**Non-obvious things that were bugs and will be again:**
- exposure must be counted in FRAMES not triggers (33.3 ms frames in 50 ms bins alias)
- only boundary-safe triggers in the trigger vector → the circular FFT correlogram
  becomes exactly the linear one, no wraparound, shuffle machinery unchanged
- zero-event observations MUST be retained (315 of 389 on baseline) or rates and
  recruitment denominators are inflated
- nShuffle must scale with nTests: p floor 1/(n+1) must sit below q/nTests or BH
  cannot reject anything
- fixed bin width + smoothing in SECONDS; an adaptive bin made latency incomparable
  between datasets
- suppression is degenerate at low rates (T_sup <= mean by construction) — gated
- `corr(M_exc, precision_z) = 0.93`: those are ONE axis. Recruitment vs 1/latency-MAD
  is the plot with independent biology.

`cell_metadata_260727.m` writes the canonical join: `roi_metadata.csv` (per ROI:
identity + both coordinate frames + its cell's scores repeated) and
`cell_metadata.csv` (per cell). Read from these, not the scattered CSVs.

**Run `analysis_preflight_260727(rootPath)` FIRST on any new dataset.** It reports
per-recording what exists (meta/dFF/cpSAM/ca_spike_data/breath_pc1/breath_peak/
breath_onset/avi) and names the blocking upstream step. `ensure_pool_260727` builds
cell_link+cell_pool on demand, so no script needs running in order, and
**cell grouping is optional** — cell_link auto-falls back to identity grouping.
Required per recording to contribute: `ca_spike_data.mat` AND both
`breath_peak_pc1.mat` + `breath_insp_start_pc1.mat`. Having the Basler `cam1*.avi`
is NOT enough — it still needs `breath_svd_pc1.m` then the peak/onset GUIs.

See [[feedback_roi_curation_port]], [[project_roi_crossfov_matching]],
[[feedback_roi_review4_save_bugs]].
