---
name: project-cmi-flu-strategy-and-immune-models
description: "CMI-Flu 2026 approach strategy + learnings from reproducing leading immune-prediction models (durability, LORIS, MCIA winner, landscape)"
metadata: 
  node_type: memory
  type: project
  originSessionId: 55c14e42-2167-4625-b80c-027c02c2b116
  modified: 2026-07-27T07:12:48.911Z
---

For [[project-cmi-flu-2026-challenge]]. Two linked bodies of work now in the repo
`D:\RUNQI\cmi-flu-2026`:

**The CMI-Flu approach strategy** — `docs/cmi_flu_strategy.md` (the "special
document"). Per-task plan for the 7 Part-I tasks; the winning architecture
(MCIAplus: reduce→augment-with-baseline→regularize + baseline fallback); durability
recipe; feature engineering (BTM modules not raw genes); LOCO validation; what NOT
to do. **The one principle: validation design + target definition + strong baselines
> architecture.**

**Immune-prediction model learning** (the user had me reproduce leading models per
regime, "learn them all"):
- `research/immune_prediction_landscape.md` — the field (Fourati, Cortese/Hagan,
  Immunaut, SCORPIO, LORIS, HR-VILAGE, EVA, single-cell/TCR negative results) + the
  7 high-standard properties.
- `research/durability_model_learnings.md` — **Cortese/Hagan durability model
  reproduced** (Nat Immunol 2025). KEY: durability target = **residual of late titer
  after regressing out peak** (verified orthogonal to magnitude, r=−0.00); Day-7
  BTM-module signature reproduced (M85 platelet activation, M1.0 integrin, M2.0 ECM
  positively predict durability = the megakaryocyte mechanism). Direct template for
  CMI-Flu Part II.
- `reimplementations/immune_reproduction_matrix.md` — per-model status. **LORIS
  reproduced** (external AUC 0.64–0.74 from its 6-feature logistic score). Repos
  cloned to `D:\RUNQI\immune_model_reproductions`.

**Why:** understanding the broader immune-prediction field (not just CMI-PB) grounds
the CMI-Flu modeling in what actually validates externally.

**How to apply:** the durability target to SUBMIT is the **raw D365 geomean** (Task 2.3
ranked by Spearman), NOT the peak-controlled residual — the residual is only an internal
Part-2 decomposition device (see [[project-cmi-flu-2026-challenge]] and
`docs/challenge_tasks_spec.md`). Part-1 features are baseline-only, so anchor
magnitude/breadth/durability to **baseline HAI titer + age + prior-vax + D0 modules**; the
Day-7 signature is a Part-2 feature. Use module scores not raw genes; validate LOCO;
foundation models only if they beat the baseline ladder under leave-one-study-out. A prior
review caught a real leakage bug in
`reimpl_multiomics.R` (in-sample 0.639 → held-out 0.081) — the standing reminder that
all preprocessing must be fit in-fold. See [[reference_cmi_flu_validation_rules]].
