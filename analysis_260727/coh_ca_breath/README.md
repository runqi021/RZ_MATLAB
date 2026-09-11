# coh_ca_breath — breath × Ca coherence, now cell-aware

Copy of `..\..\coh_ca_breath_260722\` plus three new scripts that make the analysis
run on **cells** (a neuron matched across recordings) instead of on **observations**
(one ROI in one recording). The original per-ROI scripts are all here and unchanged
in behaviour.

`README_original_260722.md` is the original README, kept for reference.

## The only file you edit: `coh_cfg_260727.m`

```matlab
cfg.rootPath       = 'D:\260721_Sert_soma_G8s\phys';
cfg.genotype       = 'Sert';
cfg.cellOutDirName = 'roi_match_out_260724';   % must match cell_cfg_260727.m
```

Everything else (`cohDir`, `cohData`, `linkFile`, `poolFile`, `curatedFile`) is derived.

## Run order

**1. `coherence_polar_general_260727.m`** — the foundation, run first. Scans `rootPath`
for every recording with `ca_spike_data.mat` + the two breath PC1 files, computes
per-ROI breath × spike coherence, writes `coherence_polar_data.mat`. Now also saves
`recNames` / `roiIdx` / `recPaths` so the join below needs no string parsing.

**2. `cell_link_260727.m`** — joins the curation to the recordings. Computes no
coherence; builds the map and audits it. Read this output before trusting anything
downstream. It asserts that every observation is accounted for exactly once, and warns
loudly about:

- recordings present in one pipeline and not the other,
- ROI counts that disagree between the matcher and `ca_spike_data.mat`,
- non-contiguous `maskL` labels, which would make `roi.roi` a *rank* rather than a
  column index and address the wrong neuron.

**3. `cell_pool_260727.m`** — the lossless container. One file holding, per recording,
the breath waveform / phase / band, and per observation, the exact sparse spike train
and its dF/F. Records how many spikes fell outside the alignment window rather than
just dropping them. Refuses to pool any recording the audit flagged.

**4. `cell_coherence_pooled_260727.m`** — one coherence estimate per cell, pooled over
all its recordings, plus the per-recording estimates and their agreement.

## How pooling avoids distorting anything

Recordings of one cell differ in length, frame rate, and the animal's breathing rate.
Three obvious approaches each break something:

| approach | what it destroys |
|---|---|
| concatenate the traces | invents continuity at the seam, which lands in the breathing band |
| resample to a common fps | distorts spike timing — the thing being measured |
| force one shared frequency band | breathing rate drifts between recordings; the band sits off-peak |

None are used. Each recording is transformed on **its own** sampling grid over **its
own** breath band, and only the cross-spectral quantities are summed:

```
A  = SUM over recordings, tapers, in-band frequencies of  Jx .* conj(Jy)
Bx = SUM |Jx|^2        By = SUM |Jy|^2
coherence = |A| / sqrt(Bx*By)      phase = angle(A)
```

This is the standard multi-segment coherence estimator — what Chronux's `trialave`
does — generalized to unequal-length segments. Nothing is interpolated, resampled or
concatenated.

**Weighting.** Summing weights a recording by how much data it contributes (the
standard choice, and the default). The equal-per-recording alternative is computed too,
as `r_pooled_equalweight`, so the two can be compared rather than argued about.

**Significance is per cell.** `confC = sqrt(1 - alpha^(1/(dim-1)))` with
`dim = K tapers × number of recordings`. A cell pooled from four recordings has a lower
threshold than a cell seen once — that is the entire benefit of pooling. So there is no
single significance circle on the polar plot; each cell carries its own, and the plot
shows the median plus the range. Frequencies within the band contribute further degrees
of freedom that this deliberately does not claim, keeping the threshold conservative.
Cells seen in ≥3 recordings also get a leave-one-recording-out jackknife CI on both
magnitude and phase.

**Comparability.** The foundation reports `mean(coherence)` over the band — a mean of
ratios. The pooled script reports a ratio of sums. These are different estimators and
need not agree exactly, so the per-recording value is computed **both** ways
(`r_rec`, `r_rec_bandmean`), letting any pooled number be traced back to the
foundation's number for the same observation.

## Output figure

`cell_coherence_pooled.png/.pdf`, five panels:

1. polar — one marker per cell, each observation drawn as a small dot tethered to it,
   so a cell whose recordings disagree is visible as a wide fan;
2. pooled vs mean per-recording coherence — does pooling change the answer;
3. phase circular SD vs number of recordings — which cells lock stably;
4. pooled coherence vs number of recordings, with the moving `confC` threshold;
5. pooled coherence vs number of spikes — estimate quality vs evidence.

## Downstream per-ROI scripts (unchanged)

`temporal_phase_perROI.m`, `temporal_phase_hist_population.m`,
`temporal_phase_hist_group_means.m`, `phase_polar_allROI.m`,
`phase_hist_sigROI_overlay.m`, `dFF_heatmaps.m`, `spike_trigger_dFF.m`,
`temporal_phase_single.m`, `temporal_phase_window_gui.m`.

All still operate per ROI-observation. They write into their own output folders under
`rootPath`, and several of those folder names are shared with the 260606/260722 runs —
re-running one overwrites the older figures of the same name.
