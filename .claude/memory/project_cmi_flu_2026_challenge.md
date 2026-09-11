---
name: project-cmi-flu-2026-challenge
description: "2026 CMI-Flu Prediction Challenge — private repo runqi021/cmi-flu-2026 at D:\\RUNQI\\cmi-flu-2026; Runqi competes, research lead is an external assistant, Claude implements"
metadata: 
  node_type: memory
  type: project
  originSessionId: 55c14e42-2167-4625-b80c-027c02c2b116
  modified: 2026-07-28T07:23:00.847Z
---

Runqi is competing in the **2026 CMI-Flu Prediction Challenge** (NIH CMI-Flu
consortium / La Jolla Institute, Bjoern Peters; hosted on Kaggle). Predict human
influenza vaccine responses — antibody breadth, durability, and initial systemic
responses (cytokines, cell composition) — from pre-vaccination multi-omics state.
Runqi also participated in the predecessor CMI-PB challenge in 2024.

**Dates:** Part I opened 2026-07-23, closes 2026-10-06. Part II opens 2026-10-20
(early responses revealed; update breadth/durability), closes 2027-01-19. Results
2027-02-02. Contestants invited as manuscript co-authors. Webinar 2026-07-31.

**Division of labour (set by Runqi):** an external research assistant leads
research, biology, literature and framework; **Claude implements** — repo, data
pipeline, models, tests, submissions. Do not silently redefine the research
question or the validation design.

**Repo:** private `github.com/runqi021/cmi-flu-2026`, local clone
`D:\RUNQI\cmi-flu-2026`. Separate from RZ_MATLAB — it is Python, not MATLAB.
Scaffolded 2026-07-23; `CLAUDE.md` there is the operating contract and there is a
dedicated `cmi-flu-modeler` agent in `.claude/agents/`.

**Why:** this is a sustained multi-month research program with parallel model
branches, not a one-off Kaggle entry.

**How to apply:** work in that repo, follow its `CLAUDE.md`, and read
`docs/claude_briefing.md` for the session-by-session kickoff prompts. The design
rule that overrides convenience: leave-one-cohort-out validation, a mandatory
baseline ladder, and no leaderboard tuning — because in CMI-PB elaborate models
lost to age and baseline titer. See [[reference_cmi_flu_validation_rules]].

**AUTHORITATIVE spec (organizers' "Prediction tasks" sheet + downloaded data,
2026-07-27) → `docs/challenge_tasks_spec.md`.** Kaggle slug is
**`cmi-flu-first-prediction-challenge`** ("CMI-Flu Invited Prediction Challenge"),
NOT `cmi-flu-internal-*` (that's the separate 2025 comp Runqi can't access). Metric =
**column-wise Spearman**. Submission cols: `Task1.1..1.4, Task2.1..2.3`; **40 test
subjects** `2025LJI.SUBxxxx`, sentinel `-99`.
- **1.1** IP10/CXCL10 `analyte=IP10` **fold-change D1/D0**; **1.2** `name=Classical_monocytes`
  abs D1; **1.3** `name=Antibody-secreting_cells_(ASC)` abs D7; **1.4** AIM `stimulation=CON`
  (→ data value **`"Conserved"`**) abs D7.
- **2.1** HAI magnitude = **geomean of 3 vaccine strains** (Victoria/4897/2022,
  DC/27/2023_MDCK, B/Austria/1359417/2021) D28; **2.2** breadth = **geomean of 12
  challenge strains** D28; **2.3** durability = **geomean of the 12 strains at RAW D365**.
- **CRITICAL CORRECTION:** the scored durability target is **raw D365 geomean, NOT the
  peak-controlled residual** — the residual is only an internal Part-2 decomposition.
  Earlier strategy said "use the residual"; that was wrong for the metric. See
  [[project-cmi-flu-strategy-and-immune-models]].
- **Test cohort has BASELINE-ONLY timepoints** (`0/-14/Pre-vacc`) → all Part-1 outcomes
  predicted from pre-vax state; baseline HAI titer (present, all 12 strains) is the
  anchor feature. **Part 2 (Tasks 3.x, opens ~Oct-20)** = same targets but training
  timepoints `0,1,7` → the D1/D7 systems readouts become features then. Task 2.3/3.3 =
  `N*` (D365 collection ongoing, off public LB, scored at end — still submit).

**Full ~12.7GB dataset downloaded + MD5-verified to `D:\RUNQI\cmi-flu-2026\data\raw\`**
(gitignored). Modalities per cohort: aim, bulkBCR, cytokine, ex_vivo_flow, rnaseq
(+BATCH-CORRECTED), serology, vdj; public RNA-seq for 2019/2020/2024_UGA, SDY224_Bcells,
SDY2867, SDY2941; metadata participants/hla/vaccineHistory/investigations. So the earlier
"API missing cytokine/T-cell/HLA/BCR" note is obsolete — Kaggle has everything.
**Kaggle auth mechanism (so next session doesn't refight it):** new `KGAT_` token works
ONLY as an HTTP **Bearer** token to `https://www.kaggle.com/api/v1/...` (the installed
kaggle CLI + kagglehub are version-broken for it); token at `~/.kaggle/access_token`
(37 chars, `KGAT_`+32). **VERIFIED WORKING 2026-07-27** — do NOT "fix" it by creating a
`~/.kaggle/kaggle.json`; that is the old username+key scheme and is the wrong mechanism
here. API confirms slug `cmi-flu-first-prediction-challenge`, deadline **2026-10-07**,
`userHasEntered: true`.

**Python env is `dlc310`, NOT `cmiflu`.** `environment.yml` declares a `cmiflu` env that
was never created; the repo is pip-installed editable into `C:\Users\Admin\.conda\envs\
dlc310\python.exe` (3.10.19). Run tests/scripts with that interpreter.

**Competition rules → `docs/competition_rules.md`** (written 2026-07-27). Binding limits:
**5 submissions/day, 2 final**, team ≤5, one account. **§3.5.d/§3.6.a forbid private
code/data sharing outside a platform-merged Kaggle Team** — the external research assistant
must be a registered teammate. No hand-labeling of the 40 test records. Deps must be
OSI-licensed without commercial-use limits. The rules page's own §1.1/§1.3 title+URL are
**stale** (say `cmi-flu-internal-*`, the 2025 comp) — same lapse as its two unresolved
license checkboxes; the rules body still applies.

**DATA AUDIT DONE 2026-07-27 → `docs/dataset_dictionary.md` § "Authoritative audit"**
(rerun: `scripts/profile_challenge_data.py --full`). Three traps that break naive code:
(1) **`Pre-vacc` = duplicate of day 0 in publicData (100% identical, n=62,285) but the
MEAN of pre-vax draws in 2025LJI** → stacking without dedup double-counts every baseline;
(2) **training support per equally-weighted task**: 2.2=2689 subj/33 cohorts, 2.3=914/8,
2.1=903/5, 1.1=127/4, 1.2=81/3, **1.3=29/1 cohort (LOCO impossible)**, **1.4=0 labels
(no publicData_aim.tsv at all → baseline carry-forward)**; (3) **the 12-strain geomean is
unreproducible in training** — 3 strains absent from publicData (incl. a Task-2.1 vaccine
strain) and no training subject has >8 of 12. Also: test cohort = exactly the 40
submission subjects, age 24-59 (public 18-90 → covariate shift), batch-corrected RNA-seq
only 35/40, bulkBCR 33/40 and public BCR only n=47. Target-definition calls are Tier-0
Q11-Q16 in `research/open_questions.md` — **research lead's call, not mine**.

**STEPS 1-3 BUILT + QC'd 2026-07-27** (89 tests, ruff clean, 30/30 QC). Run:
`build_interim.py` → `build_task_matrices.py` → `run_task_baselines.py` →
`qc_steps_1_3.py`. New code: `ingest/challenge.py`, `outcomes/tasks.py`,
`validation.splits.leave_one_cohort_out_person_safe`, `configs/data.yaml` filled.
**4th trap found while building — THE RESOURCE RE-ENROLS THE SAME PEOPLE ACROSS
COHORTS** (411 UGA + 163 ImmPort labels span >1 study; identity confirmed by
sex-constancy + age tracking the year gap at 99.3%/100%; **49.5% of Task-2.2 rows**).
So `subject_id` = the PERSON, `enrolment_id` = person×study, and LOCO must use the
person-safe splitter (plain `leave_one_cohort_out` now RAISES). Cost: 42.1%/31.4%/5.3%
of naive training volume on 2.1/2.3/2.2. Test cohort is clean (2025LJI SUB1829-9562 vs
ImmPort SUB112829-441405, disjoint).
**Baseline ladder to beat (Spearman, person-safe LOCO):** 2.3=0.673, 2.1=0.648,
2.2=0.572, 1.2=0.581, 1.3=0.311 (OPTIMISTIC single-cohort P2), 1.1=0.181 (SD>mean),
1.4=unscoreable. Baseline titer alone is the bar on HAI tasks. **Task1.1 INVERTS** —
carry-forward = −0.121 there because y=log(post/pre) is anti-correlated with pre
(CLAUDE.md rule 7); never reuse carry-forward on a fold-change target. B0 scores NaN
by design (constant prediction → undefined Spearman). Per-fold ≫ pooled (0.648 vs
0.324) because folds refit scaling; **per-fold LOCO is the right analogue** since the
leaderboard ranks 40 subjects within one cohort.

**🔴 BIGGEST FINDING (2026-07-28) — TASK 1.2's TRAINING TARGET IS A DIFFERENT QUANTITY
FROM THE SCORED TARGET.** `Classical_monocytes`: test 2025LJI = **percentage of CD45+**;
the ONLY three cohorts with D1 data are SDY296 `cells/ul`, SDY301 `cells/ul`, SDY416 `ul`.
**Zero subjects anywhere have D1 monocytes in percentage.** No conversion possible —
SDY296/301 record no parent_population and their 7 populations don't span CD45+. So the
ladder's 0.544-0.581 was validated against COUNTS while the scored quantity is a
FREQUENCY (differ by a per-person denominator varying ~2.5x). **FIX: Task1.2 ships BCF
(carry-forward), not the fitted rung** — carry-forward needs no cross-unit map and is
unit-INVARIANT (measured baseline→post rho: percentage .689/.739, cells/ul .603/.753,
percentile .878, ul .420). Independent confirmation: old fitted preds ranged 24.9–107
for a %-of-CD45+ target (impossible; monocytes are ~2-10%); carry-forward gives 0.94–11.2.
→ `submissions/part1_unitfix_part1.csv`. **ALWAYS CHECK `unit` COLUMN train-vs-test.**
Task1.3 units DO match (both percentage). Task1.1 both pg/ml but platform differs
(Luminex xMAP vs Legendplex) = calibration not quantity.
**Task1.x ladder A–E**: Task1.1 rung **D_offset** (predict absolute log2(D1), subtract
the exactly-known log2(pre); coefficient on baseline fixed at −1) is the ONLY arm
consistent — **4/4 folds, +0.077** — because baseline predicts D1 LEVEL consistently
(+.41/+.86/+.76/+.79) but the fold CHANGE with flipping sign (−.40/+.61/−.19/−.50).
Task1.1's biological hypothesis is DEAD (C_bio_set −0.031, rung B picks a different
feature every fold; **Task1.1 has 0% RNA-seq coverage** so the ISG hypothesis was never
testable). Task1.2: baseline wins, bio −0.115 (0/3), nonlinear −0.248 (0/3).
**Webinar Qs mostly retired by measurement**: Q2 moot — geomean-over-available vs
restrict-to-shared-9 correlate **1.0000** (the 3 missing strains are absent from ALL
training, so the sets are identical); rank-avg alternative 0.96. Q3 non-question (rules
impose no feature restriction). Q4/Q5 already handled conservatively.

**VALIDATION REPAIRS 2026-07-27 (research lead ordered).** (1) **FOLD-GLOBAL LEAK FIXED**:
`rnaseq_baseline(top_variable=)` selected top-k variable genes across ALL subjects incl.
held-out cohorts, then cached it. Default now `None`+warns; selection moved to
`features.blocks.TopVarianceSelector` refit inside the Pipeline; wired via
`_prefix(top_variance=)`. (2) **LOSO added**: `validation.splits.leave_one_season_out`,
`vaccine_season` flows investigations.tsv→subjects.parquet→task matrices as a GROUPING
col; `run_task_baselines.py --protocol {loco,loso,both}` default both. (3) **11 leakage
regression tests** `tests/test_leakage_guards.py` (person/cohort/season/fold-global).
**SURPRISE: LOSO is HIGHER not lower — Task2.2 0.607±0.172 (15 season folds) vs LOCO
0.572±0.283 (33 cohort folds).** LOCO wasn't optimistic, it was NOISY: several LOCO folds
have test n as low as 6 → per-fold Spearman is noise dragging the unweighted mean down.
Tiny-fold problem dominated the shared-season problem. 2.1/2.3 unchanged (each cohort IS
its own season). LOSO still primary (matches predicting a future season) but as a
better-conditioned estimate, not an optimism correction.
**DEADLINE DISCREPANCY UNRESOLVED**: Kaggle API says `2026-10-07T08:00:00Z` (both
`deadline` and `mergerDeadline`); the invitation says **2026-10-06**. NOT a timezone
artifact. **Plan to Oct 6** (earlier = safe); both recorded in `docs/competition_rules.md`.
**RETRACTED OVERCLAIM**: "signal is strain-specific not landscape-level" — failure of the
tested aggregates doesn't establish biology (wrong antigenic coordinate, non-linear
landscape, birth-cohort-varying imprinting reference, subject-varying panels, same-data
interpretation). Permitted: "certain historical-strain titers show cohort-consistent
associations not captured by the simple aggregates tested here."

**CONVERGENT CONCLUSION (2026-07-27): the binding constraint is COHORT HETEROGENEITY,
not method choice.** Two independent designs, same failure: (a) Task2.3 convex blend —
the *oracle* w between baseline and boost FLIPS across cohorts (2017_UGA wants w=1.00
pure-baseline, 2021_UGA w=0.00 pure-boost); inner CV can't learn a cohort-dependent
weight; blend +0.022 [-0.015,+0.060] inconclusive, and ALL arms lose to B3=0.673.
(b) Late fusion (per-modality models, OOF preds, cross-fitted meta) loses: Task2.3
**-0.176 [-0.333,-0.036]** 1/8 folds, Task2.1 -0.052 [-0.125,+0.002] 2/5. I predicted
its 2 real defects (collinear absence-fill, train/test calibration mismatch) were the
cause — **fixing them changed almost nothing** (-0.196→-0.176), so the degradation is
real. Kept `late_fusion_run_FLAWED.log` for audit. **Baseline titer transfers across
cohorts; learned combinations of it with anything else do not.** Untested alternative:
combination rule keyed to OBSERVABLE cohort traits (platform/season/age), not a global w.
**Derived features + interactions: also null.** 2 year-slope probes survive FDR
(q=0.023/0.029) but fail sign-stability; 0/6 hierarchical interactions (all p>0.45).
Aggregate antibody-landscape summaries are null while INDIVIDUAL historical strains
reach q≈0.04 with consistency 1.00 → signal is **strain-specific, not a landscape
gradient**.

**METHODOLOGICAL CHECKPOINT 2026-07-27 (research lead ordered; 3 of my claims retracted).**
(a) **"omics actively hurts" RETRACTED** — `run_model_families.py` compared `full` vs
`omics` against a bar filtered on TASK ONLY, not config, so an 80-subject/2-fold model
was scored against a 903-subject/5-cohort bar. Proof it's fatal: the SAME baseline-only
model scores 0.659 / 0.673 / 0.798 on three subsets of Task2.3. All 12 omics
`beats_bar` values are void. (b) **"we are competitive" RETRACTED** — compared our CV
Spearman to CMI-PB's *count of significant models*, not a score distribution.
(c) **Unit inflation**: rows are ENROLMENTS, persons repeat — Task2.3 = 414 people not
914; Task2.2 = 1,619 not 2,689.
**Correct tool is `scripts/paired_modality_comparison.py`** (same rows, same folds,
paired per-fold delta + cluster bootstrap). Result: 10/11 deltas negative, 5 conclusive
— BUT the decisive 33-fold case is **Task2.2/hai adding 60 more HAI STRAIN titers (not
omics), also harmful (−0.047)**. Degradation tracks FEATURE COUNT not modality
(14-25 feat → −0.05; 47-67 → −0.05..−0.18; 1,507 → −0.40..−0.79). So it's **feature-block
dilution under elastic net**, NOT omics being uninformative. Fix = late fusion / nested
selection / baseline-anchored residuals, NOT imputation or latent factors.
**Single-feature screen** (`screen_single_features.py`, MUST use `--n-perm 1999`):
I had 2 opposing bugs — null mis-aggregated (nanmean vs sum-over-all-cohorts → 1,334/1,521
false positives) and p-floor 0.005 making q<0.10 unreachable. Both fixed. **Corrected
answer: NO single gene/cytokine/cell population reaches q<0.10 on any task** (Task2.2 =
1,521 features × 27 cohorts → 0). Only FDR-survivors are **pre-vax HAI titers to
historical strains on Task2.3**: 4 H3N2 strains NEGATIVE (−0.157..−0.232, sign
consistency 1.00), B/Victoria + recent H1N1 POSITIVE — imprinting-shaped, 3-4 cohorts,
hypothesis not result.
**Imputation GATED**: cytokine/neutralization/rnaseq/tcell_reactive are
**cohort-determined** (>90% of cohorts all-or-nothing) → missingness IS study identity,
nothing to reconstruct from; 19/28 modality pairs have <50 shared subjects or <3 shared
cohorts. Also: 22 of Task2.2's 33 cohorts share a vaccine season → LOCO there is NOT
season-holdout. Docs: `validation_unit_audit.md`, `omics_missingness_review.md`,
`paired_modality_findings.md`, `single_feature_findings.md`, `outputs/audits/`.

**FINAL: 0 of 24 model runs beat baseline titer** (full sweep 2026-07-27). Our mean
over the 6 scoreable tasks = **0.494** (0.424 if Task1.4=0). Calibration: CMI-PB round
2, of 49 models the per-task significant-Spearman counts were 19/12/11/7/2/**0** of 49
— nobody beat chance on the hardest task, and 22/22 literature models failed to beat
baseline. So **not beating baseline is the norm, not our failure**; the round-2
winner's decisive habit (R7) was baseline-fallback, which `make_part1_submission.py`
already does. Position is competitive.
**WHY OMICS FAILED = complete-case join, NOT lack of interactions.** Requiring an
omics block drops Task2.1 from 903 subj/5 cohorts to 80/2 folds at 3,068 features →
indistinguishable from chance (perm p=0.40-0.80 vs 0.00 in `full`). The **round-3
CMI-PB winner used tensor-decomposition imputation (tensorMiss) and never
complete-cased**. So fix order is: (1) block-wise/tensor imputation to keep full n,
(2) module scores (BTM ~350 modules, not 3,000 genes; cell ratios), (3) true
multi-block MCIA/MOFA = the cross-omics "interaction score". At n=80 no feature
representation helps — imputation must come first. NOTE gseapy/scanpy/anndata are in
environment.yml but NOT installed in dlc310.

**FAMILIES 4/6/15 RUN + SUBMITTABLE FILE EXISTS** (2026-07-27):
`submissions/part1_baseline_ladder_part1.csv` — 40 rows, all 7 columns, 0 sentinels,
registry-traced; made by `make_part1_submission.py` (rung per task from recorded CV,
never the leaderboard). **No model family beats baseline titer** — 2.1: enet 0.485 /
boost 0.440 vs bar 0.648; 2.3: boost 0.642 / pcr 0.639 vs bar 0.673. Permutation nulls
all ≈0.00 (p=0.00) so models are far above chance, just not above baseline. **This
reproduces the CMI-PB finding on new data.** Note boost on 2.3 has fold SD 0.142 vs
baseline's 0.223 — lower mean but much more stable; may matter since the leaderboard is
one unseen cohort.
**MODALITY COVERAGE KILLS NAIVE MULTI-OMICS**: within training rows, HAI=100% but
cytokine 0-3.4%, cell_freq 0-11%, RNA-seq 8.9-25%. The test 40 have 100% of everything
— asymmetry is training-side only. Hence two configs in `run_model_families.py`:
`full` (HAI+elisa+tabular, full n) vs `omics` (restricted to rows carrying omics).

**Single-cell objects**: Drive folder has **5** (`lji2025_tcr, sdy2903, sdy3099, sdy3412,
sdy3427`). Have `lji2025_tcr.qs2` (531MB) + `sdy2903.qs2` (9.18GB, resumed 2026-07-27);
**sdy3099/sdy3412/sdy3427 not downloaded**. Not needed for Part-1 tabular targets. Resume
with `gdown.download(id=..., resume=True)` from `base` env (gdown installed there, not
dlc310); folder IDs `1D1kwYPO0...` (superset) / `1Y22S5dNXKf...` (single-cell only).

**CMI-X lineage (5 challenges):** CMI-PB internal 22-23, invited 23-24, public 2024
(Runqi = `runqi021/CMI-PB3`); CMI-Flu internal **2025** (influenza, not pertussis —
`pramodsshinde/CMI-Flu-Internal-Challenge-models`, the direct predecessor), public
2026 (now). Runqi was right that a 2025 round existed and wasn't pertussis.

**Prior work is all open-source** — repo has full archaeology. **The living record
of my understanding is `research/understanding_ledger.md`** (per-attempt, tagged by
verification depth: code-read / verbatim / README) — keep it updated as more is
read. Supporting: `docs/prior_challenges.md` (measured results, verified verbatim),
`research/methods_catalog.md` (techniques by stage), `research/repo_inventory.md`
(~50 repos), `research/prior_model_matrix.csv`.

Verified findings (code + primary-text): 0/22 literature models beat baseline; age
significant on 0/6 round-2 tasks; per-task significant counts 7/19/11/2/12/0.
Round-2 winner = Saonli Basu/UMN (MCIAplus = MCIA factors + baseline + demographics,
LASSO/LMM, **baseline-fallback: ship day-0 baseline when the model loses in CV**);
round-3 winner = Cheng-Chang Wu (tensor-completion imputation → PCA → gamlasso/RF).
Field consensus: impute→normalize/ComBat→MCIA/MOFA reduce→augment-with-baseline→
regularize/RF w/ CV-compare→baseline-anchor→honest validation. Runqi's own 2024
entry validated in-sample (its main flaw) — assess priors by measured Spearman, not
by whose entry it is.
