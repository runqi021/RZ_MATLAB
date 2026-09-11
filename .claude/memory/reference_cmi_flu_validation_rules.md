---
name: reference-cmi-flu-validation-rules
description: "The CMI-PB precedent that dictates how CMI-Flu models must be validated — baselines beat elaborate models, random splits overstate transfer"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 55c14e42-2167-4625-b80c-027c02c2b116
  modified: 2026-07-24T02:17:44.773Z
---

Design constraints for [[project-cmi-flu-2026-challenge]], taken from the
published CMI-PB (B. pertussis) predecessor challenges. Verify against the primary
papers before any claim drives a decision — they are recorded as working priors in
`docs/prior_challenges.md`, not as established fact.

1. **Literature signatures did not transfer.** In the CMI-PB resource study >30
   published predictive models were tested and essentially only an **age-based
   model** was predictive.
2. **No single architecture won.** Round 2 evaluated 49 algorithms from 53
   scientists; CatBoost, sparse regression, ensembles and multi-block latent models
   each won different endpoints.
3. **Preprocessing separated the top entries**, not model novelty — imputation,
   batch correction by timepoint, feature selection then dimension reduction,
   MCIA integration, mixed-effects model.
4. **Evaluation was on a newly generated cohort**, not a random split.

**How to apply:** in `cmi-flu-2026`, leave-one-cohort-out is the primary protocol;
repeated subject-grouped K-fold is reported only as an explicitly optimistic
secondary; ungrouped splits are diagnostic. Every model is scored against the
baseline ladder B0–B7, and only counts as beating a baseline when the improvement
exceeds its own fold-to-fold SD. All preprocessing fits inside the fold. Never tune
on the Kaggle leaderboard. These are enforced in `src/cmi_flu/validation/` — the
splitter raises rather than producing a leaky split.

Primary sources to read in full: CMI-PB resource paper (2023); "Putting
computational models of immunity to the test", PLOS Comput Biol (2025), plus the
contestant GitHub repos it required to be deposited.
