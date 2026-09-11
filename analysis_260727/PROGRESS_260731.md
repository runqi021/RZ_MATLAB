# 2026-07-31 — Vgat 0730 into the archive, archive made drive-independent, matcher retune audited

Continues `PROGRESS_260729.md`. Written mid-session: the Vgat chain (§4) was still
running when this was recorded, and §5 has not started.

---

## 0. Where everything lives now

The acquisition sessions moved **D: → E:** on 2026-07-30, and E: is an external
drive that gets unplugged. Two things were repointed and one was copied so nothing
downstream depends on a removable drive:

| what | change |
|---|---|
| `archive_into_summary_260729.m` `SRC` | D: → E:, plus the new Vgat source on D: |
| `site_mapping_260729.csv` (31 rows) | `source_path` D: → E: |
| `Ventral_surface_polar_coh_vs_rayleigh_260729.m` `cell_link_sources` | now `fullfile(rootPath, ...)` — reads the archive's own copies |

**The archive is now self-contained.** Grouping artefacts were copied in beside each
session (robocopy, sources untouched):

```
Sert\0721\roi_match_out_260724     138.6 MB   Sert\0721\cell_pooled     22.8 MB
Vglut2\0728\roi_match_out_260727    67.6 MB   Vglut2\0728\cell_pooled    4.0 MB
Vgat\0730\roi_match_out_260727      93.1 MB
```

Verified `curated=1, cell_link=1` for Sert/0721 and Vglut2/0728. **E: was released
after this.** If a script ever needs an acquisition drive again, that is a regression —
the archive should be the only input.

## 1. Archive contents — 13 blocks, 52 recordings

```
ChAT   0124(1) 0521(2) 0522(3) 0523(1) 0524(1)
Sert   0721   6 sites, 10 recordings
Sst    0322(2) 0330(1)
Vgat   0323   3 sites,  3 recordings
Vgat   0730   3 sites,  9 recordings          <- NEW
Vglut2 0224(8)  0728(3)  1124(8)
```

## 2. Vgat 0730 — DONE (9 of 24 recordings archived)

`archive_into_summary_260729('copy')`: **9 copied, 45 skipped**, sources untouched.

**The template's site guess was wrong and was overridden.** It infers site from the
filename prefix, so all ten spike-bearing `roi1_*` recordings landed in `cell1`. Site
must come from `motorPosition`, clustered at the archive's `mergeSiteUm = 100`:

| site | n | stage centre | spread |
|---|---|---|---|
| cell1 | 5 | (−1009, 1141) | 75 µm |
| cell2 | 2 | (−887, 1335) | 38 µm |
| cell3 | 2 | (−1185, 975) | 0 µm |

Sites are ~230 µm apart, so this is three sites, not one. **Judgement call to check:**
`cell1` was formed by single-link chaining and its two extremes are **145 µm** apart —
above the 100 µm merge radius. Defensible because those five span zooms 1.4×–3.2× and a
1.4× FOV comfortably contains all of them, but if you want them split the fix is one
edit to `site_mapping_260729.csv` and a re-run of `'copy'`.

**15 recordings excluded, and the reason is now written into the csv `note` column:**
- 7 × `roi2_*` — **no breath at all** (no `breath_pc1`, no peak, no onset). `cpSAM` is
  done, so these need `breath_svd_pc1` + both GUIs before they can ever enter.
- 6 × `roi1_*` — no `ca_spike_data`
- 1 × `roi1_3x_x950y1000_z70` — has spikes but **0 active ROIs**
- (Vglut2/Sert exclusions unchanged from 260729)

## 3. Matcher retune — audited, one false merge found

`nccMergeThr` 0.60 → **0.40**, `maxPairDist_um` 300 → 400, `maxPairDz_um` 50 → 100,
`patchHalf_um` 30, `targetUmPerPx` 0.5.

Triggered by `Vglut2/1124` cell 13, which was visibly the same neuron as `recB#12` but
scored **0.573** against a 0.60 threshold. Root cause: patches are cut on the Cellpose
**centroid** and scored at **zero shift**. When Cellpose segments the same cell with a
different extent (257 px vs 324 px — the recB mask swallows more of the process) the
centroid moves, the patches misregister by several px, and the correlation collapses.
**It is segmentation disagreement, not appearance change.** The principled fix is a
±3–5 px translation search in the NCC; the threshold drop treats the symptom.

**Position is the check NCC cannot provide.** Estimate the inter-recording frame offset
from the *confirmed* pairs only, then measure each merge's residual after removing it:

| session | result |
|---|---|
| Vglut2/1124 | 12 true pairs ≤ 3.4 µm; **1 false merge at 236 µm** (recA#23↔recB#3, NCC 0.449) |
| Vglut2/0224 | 12 testable cells all ≤ 1.8 µm, **0 false** — including one at NCC 0.439 |

So 0.40 rescues real pairs on both sessions *and* admitted one false merge on 1124.
**Recommendation: gate on offset-corrected residual (<~10 µm) in addition to NCC**,
rather than lowering the threshold further.

**Caveat that matters at scale:** the audit is only as strong as its anchors. On 0224
only **3 of 6 FOVs** were co-registered by a high-NCC cell, so 13 of 25 multi-recording
cells could not be tested at all — position is *silent* there, not clean. My first pass
fitted the registration to all merges including the suspect ones and reported "12 of 25
impossible"; that was an artefact and was withdrawn. Fit on `minNCC ≥ 0.80` only.

Audit script: `scratchpad\check_0224_robust.m` — takes a `roi_match_out_*` path, needs
no matcher re-run.

## 4. Vgat 0730 chain — DONE (4.3 min, all 7 steps OK)

`run_analysis_260727('D:\260730_vgat-g8m_shiverer\phys','Vgat')`. `cell_pooled` copied
into `Ventral_surface_summary\Vgat\0730\cell_pooled`; its row added to
`cell_link_sources`. **`cell_link_260727` reports JOIN IS CLEAN** — no ROI-count
mismatch, no index overflow, no non-contiguous maskL.

```
871 observations -> 709 cells (62 multi-recording, 647 singleton, 93 tossed)
pooled: 268 usable obs, 239 cells with >=1, 29 with >=2, 1429 spikes
510 obs dropped: matcher FOV with no ca_spike_data (13 of 23 recordings)
coherence: 4 of 14 ROIs significant, all at +1.41..+1.75 rad
PETH (12 cells testable): 0 excitation significant, 6 precision (onset) / 7 (peak)
```

### THIS ANIMAL BREATHES AT ~0.4 Hz — five times slower than the rest of the archive

The chain reported breath bands of **0.34–0.74 Hz** on every Vgat/0730 recording. That
is not a pipeline fault. Checked independently against the hand-curated peaks
(`fps / median inter-peak interval`): **0.37, 0.38, 0.38, 0.68 Hz** — it agrees with the
band exactly. Other sessions in the archive run 2–3 Hz (Vglut2 2.2 Hz). Median cycle here
is **2.567 s** versus ~0.45 s for a 2.2 Hz breather.

Consequences to respect:
- **Phase pooling across genotypes stays valid** (phase is normalised per cycle), but any
  ABSOLUTE-TIME comparison against the other genotypes is comparing 2.6 s cycles with
  0.45 s cycles. The PETH window guidance for this session came out −0.5..2.0 s.
- The coherence band is ~0.35 Hz wide at TW=4; check the resolution is adequate at this
  low a frequency before quoting Vgat coherence against the others.
- `roi1_3x_x1350y850_z265` sits at **0.68 Hz**, roughly double the rest, and contributes
  **3 of the 4 significant coherence ROIs**. It is also the deepest (z 265). Worth asking
  whether that recording caught a different physiological state.

### fps 42 and 47 are REAL, not a detection error

Two recordings ran at 42.084 and 47.147 Hz (`scanFrameRate_raw`) against 30 Hz for the
other eight, and TIF page counts match the filename (5000f, 6000f). `detect_session_fps`
is reading them correctly — do not "fix" these.

## 5. Polar plot — DONE

`Ventral_surface_polar_coh_vs_rayleigh_260729.m` →
`Ventral_surface_summary\polar_coh_vs_rayleigh_260729\` (png/pdf/mat +
`polar_scores_percell.csv` 378 rows, `polar_scores_perobs.csv` 389 rows).

```
389 ROI-observations -> 378 cells across 6 groups
  IO 272 | Sert 58 | Vglut2 25 | Vgat 16 | Sst 4 | ChAT 3
cell identity applied: Sert 447->338 masks (16 tossed)
                       Vglut2/0728 155->155 (0 tossed)
                       Vgat/0730 871->709 (93 tossed)
logZ crit(alpha=0.05) = 1.10 -> 80 significant, 51 also pass the shuffle
BH-FDR q=0.05: 56 / 378
agreement: BOTH 18 | coherence only 2 | Rayleigh only 62 | neither 296
```

**Vgat supplies the two strongest cells in the whole archive** — logZ 3.79 and 3.75
(`roi1_3x_x1350y850_z265`/29 and `roi1_2.4x_x1300y900_z250`/32), at +128° and +98°.
Both come from the two deep, faster-breathing recordings (0.64–0.74 and 0.45–0.50 Hz),
i.e. the same pair that supplied 3 of the session's 4 significant coherence ROIs. Vgat
sits in the same +98..+168° sector as Vglut2, while Sert clusters on the opposite side
(−119..−180°).

**Stale file to be aware of:** `polar_scores.csv` in that folder is from 2026-07-29 and
was NOT regenerated — the current run writes `_percell` and `_perobs` instead. Delete or
ignore the old one.

### 5b. Per-cell temporal histograms (Vgat) — DONE

`breath_time_percell_hist_260727` after `coh_cfg_260727('set', <vgat source>, 'Vgat')`.
Window ±2.5 IBI = −6.42..6.38 s at 50 ms bins (median IBI 2.567 s).

```
28 (recording x trigger) tests, 7 significant at BH q<0.05
12 per-cell figures
```

Copied into the archive at
`Ventral_surface_summary\Vgat\0730\analysis_260727\breath_time\`, mirroring how
`Vglut2\1124` carries its own — `percell_hist\` (12 figures), plus `onset\`, `peak\`,
`trigger_compare` and `per_recording_tests.csv`.

Target: `Ventral_surface_polar_coh_vs_rayleigh_260729.m` — the linear-phase / log-Z
version. Panel A radius = band-averaged |C|; panel B radius = log Z with the
occupancy-normalised Rayleigh. `scan_dirs` already includes `Vgat`, so the new session
is picked up automatically once the archive is in place.

`cell_link_sources` now lists Sert/0721, Vglut2/0728 and Vgat/0730, all resolving inside
the archive. **0224 and 1124 still have no `cell_link`** — they have curated matchers but
their chains were never run, so pooling is a no-op for them and their ROIs enter as one
cell each. Running `run_analysis_260727` on those two archive folders would close that gap.

## 5c. Cartoon map — DONE

`Ventral_surface_ROI_on_cartoon_260729.m` → `ROI_on_cartoon.png/.pdf/.mat` in
`C:\Users\Admin\Desktop` (output goes to `fileparts(cartoonMat)`, NOT the summary folder).

```
1414 cells (1536 masks) / 48 FOVs, 35 sites
80 cells significant at logZ >= 1.10 (alpha=0.05)
masks -> cells: Sert 389->296 | Vgat 248->219 | Vglut2 169 (14 matched)
                IO 679, ChAT 39, Sst 12 -- no matcher, 1 cell per mask
```

**Two things had to be fixed to get this right.**

1. **The map has its OWN `cell_link_sources`, and I had only fixed the polar script's.**
   It still pointed at the dead D: paths, the lookup failed with a *warning not an error*,
   and every group silently fell back to "1 cell per mask" — Sert reported 402→402 instead
   of 389→296. Now resolves with `fullfile(sumRoot, ...)` like the polar. **Lesson: grep for
   every `cell_link` registry, there is more than one.**
2. **Vgat/0730 IS a rotated-mount session** (Runqi, 2026-07-31). I had inferred it was not,
   from its stage footprint overlapping the old-mount cloud — that inference was wrong and
   the rule is temporal: **everything from 260721 onward is on the new axis.** Added to
   `rotated_sessions`; rotated FOVs went 13 → 22 and Vgat moved into the same ventral strip
   as Sert/Vglut2. See [[project_rotated_mount_sessions]].

## 6. Also fixed today

- **`pair_behavior_to_phys.m` needs dffQC first.** On a fresh session the TIFFs sit flat
  in `phys\` and the `tifBase == parentName` test matches nothing → "Found 0 FOV
  folders". `Batch_dffQC_260325.m:263-267` is what creates the per-recording folders.
  My earlier claim that the script is "correct in either order" was only about the
  *frame count*, not the folder structure.
- **`Vglut2\0224\pFN\pFN_roi2_z5_25x_...`** has no breath because its interactive crop
  never completed — no `breath_crop.mat`, so `breath_svd_pc1` printed one red
  "crop/extract failed -- skip" and moved on. Video is fine (2001 frames, decodes). The
  `breath_pc1_crop_x*_y*_s*.mat` files there are from an older crop-position scan tool
  and are **not** the pipeline's output.

## 7. Next

1. **Generate `cell_link` for Vglut2 0224 and 1124.** Both have curated matchers but their
   chains were never run, so the polar treated their ROIs as one cell each — 0224 and 1124
   between them contributed 8 of the 18 BOTH-significant cells, all unpooled. Running
   `run_analysis_260727` on those two archive folders is the single change that would most
   affect the numbers above.
2. Decide on the position-residual gate (§3) — the one change that would make the 0.40
   threshold safe rather than lucky.
3. Breath detection for the 7 `roi2_*` Vgat recordings, then re-run the archive template.
4. Cartoon map: **check whether Vgat/0730 shares the single stage zero** the map assumes
   before its ROIs go on it, or it lands in the wrong place the way
   `ro1_2x_x1018y1130` did for Sert. Vgat is now in `scan_dirs`, so the map will pick it up
   automatically the moment it is re-run — this needs checking BEFORE that.
5. §4 flags a real physiological question: does the 0.4 Hz breathing make Vgat comparable
   to the other genotypes in absolute time at all, or should Vgat only be reported in phase?
