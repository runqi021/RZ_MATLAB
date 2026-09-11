---
name: Statistical testing at scale — theoretical options
description: Reference on multiple-comparison strategies as dataset scales from N=30 to ~20k neurons; options and field standards, not committed choices
type: reference
originSessionId: c425b17d-7b64-450c-b474-d258884494d6
---
# Significance testing at scale — reference

**Note:** Reference material on field-standard options and their rationale. Not a recorded decision — no specific pipeline commitment made. Kept for future planning when the dataset grows.

Context: user's dataset is scaling from single-FOV (~3–30 ROIs) toward whole-animal (~20k neurons across many FOVs). Strict Bonferroni doesn't scale, so these are the alternatives the field uses.

## Two complementary per-ROI tests
- **Coherence** `|C| > confC` — tests linear frequency-domain covariance
- **Rayleigh** on spike/phase samples — tests phase-cluster non-uniformity
- They answer different questions, so an AND rule (both pass) is stronger than either alone and lets you loosen α on each

## Options by scale

### Small (N ≤ ~50 pooled tests): Bonferroni is tractable
- α_per_test = 0.05 / N
- At N=30: α ≈ 0.0017; a conservative 0.001 is defensible

### Medium (N ~100–1000): FDR preferred
- Benjamini-Hochberg at q=0.05 controls expected false discovery fraction
- MATLAB: `mafdr(pvals, 'BHFDR', true)`
- Less conservative than Bonferroni, still principled

### Large (N ~20k): three tools typically combined
1. **FDR (Benjamini-Hochberg)** instead of Bonferroni
2. **Shuffle / permutation nulls** — circularly shift breath vs neural signal 500–1000× per ROI, threshold at 99th pct of shuffled statistic. Respects autocorrelation (analytic confC assumes white noise)
3. **Population-level framing** — per-ROI significance becomes a summary statistic, not the endpoint:
   - Fraction sig above shuffle-null floor
   - Rayleigh on population of preferred phases `θ_pref`
   - Effect-size histograms (vector strength / |C|)
   - Stratify by region/cell-type rather than naming individual neurons

## Example recipe one *could* use at 20k scale
```
per-ROI:
  α_coherence = 0.01 (gate, uncorrected)
  Rayleigh p_i on spike phases (intrinsic DOF from n_events)
  shuffle null: ≥500 circshifts → p_shuffle_i

cross-ROI:
  FDR at q=0.05 on Rayleigh and/or shuffle p-values
  AND rule for "phase-locked" label

report:
  - fraction sig by region / cell-type
  - Rayleigh on population θ_pref distribution
  - effect-size (vector strength / |C|) histograms
  - confC kept as diagnostic only, not gate
```

## Key conceptual point
With very large N, huge power detects biologically trivial effects. The scientific claim at scale is about population distributions, not individual neurons. Always report effect sizes alongside p-values; per-neuron pass/fail collapses into a fraction that itself needs a null.
