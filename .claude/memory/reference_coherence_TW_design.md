---
name: Coherence TW design theory (breath vs neural signals)
description: Theoretical framework for choosing Chronux multitaper TW based on record length and breathing peak FWHM — reference only, not a committed decision
type: reference
originSessionId: c425b17d-7b64-450c-b474-d258884494d6
---
# Chronux coherence TW — theoretical rationale

**Note:** Reference material on what *could* be done and the standard rule of thumb. Not a recorded decision — no change committed to the script.

Script discussed: `Breath_summary_timeNphase_260331/breath_combined_summary_260329.m` (currently uses TW=5, f_breath_search=[1.5 2.5])

## Core formulas
- `params.tapers = [TW, 2*TW-1]` — TW is time-bandwidth product (unitless)
- Frequency smoothing half-bandwidth: **W = TW / T_record** (seconds)
- Tapers K = 2·TW − 1 → degrees of freedom
- Standard guideline: **2W should be < FWHM of peak of interest** (avoid smearing)

## Empirical breathing statistics (reference snapshot, scan of 38 sessions on D:/batch_dffQC_test_260325, 2026-04-16)
- Duration T: median **99 s** (36/38 ≤110 s)
- Breathing peak: median **2.03 Hz**, range **1.20–2.72 Hz**
- Breathing FWHM: median **0.161 Hz**, IQR [0.146, 0.198], 90th pct 0.49 Hz
- Minimum observed FWHM 0.055 Hz = artifact of TW=5 smoothing floor
- Implication: `f_breath_search = [1.5 2.5]` in the current script could miss low breathers (peaks span 1.20–2.72 Hz) — widening to `[1.2 3.0]` is an option worth considering

## What different TW values imply at T ≈ 100 s
| TW | K | 2W at T=100s | 2W / median FWHM | Interpretation |
|---|---|---|---|---|
| 3 | 5 | 0.060 Hz | 0.37 | sharp, fewer tapers |
| 4 | 7 | 0.080 Hz | 0.50 | balanced |
| 5 (current) | 9 | 0.100 Hz | 0.63 | starting to broaden peak |
| 7 | 13 | 0.140 Hz | 0.87 | smears peak |

These are trade-offs, not prescriptions. Lower TW = sharper frequency resolution; higher TW = more DOF / smoother estimate.

## How to verify for new datasets
Before changing TW, scan FWHM distribution:
```matlab
bp_files = dir(fullfile(rootPath, '**', '*breath_peak_data.mat'));
% For each: load S.breath, detrend, mtspectrumc with TW=5, find peak in
% f_search=[1.2 3.0], measure half-max width. Then check whether
% 2W = 2*TW/median_T is below median_FWHM.
```
