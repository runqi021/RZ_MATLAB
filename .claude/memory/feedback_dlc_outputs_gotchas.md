---
name: DLC output gotchas — retrain cleanup, grep patterns, snapshot suffixes
description: Three lessons from a DLC retrain: old model/training-dataset folders must be deleted, grep -v "_R" is buggy because of _Resnet50, and snapshot suffix (best-10 vs best-20) distinguishes old from new outputs
type: feedback
originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---
## 1. Retrain requires cleaning old model + training-dataset folders

`deeplabcut.create_training_dataset(config)` issues an **interactive "overwrite? (yes/no)" prompt** when the iteration folder already exists. Background Python has no stdin → `EOFError: EOF when reading a line` → pipeline dies.

Fix: before `create_training_dataset`, delete:
- `<project>/dlc-models-pytorch/iteration-0/<shuffle>/`
- `<project>/training-datasets/iteration-0/`

Do NOT try to pass `force=True` or similar — `create_training_dataset` doesn't have that kwarg in the PyTorch engine. Matches user's `rerun = delete + redo` preference in `feedback_rerun_logic.md`.

## 2. Never `grep -v "_R"` on DLC filenames

All DLC outputs contain `_Resnet50`, so `grep -v "_R"` strips everything, not just R-side files. Use **`grep -v "dorsal_R"`** to exclude R-side files while keeping L-side.

| Pattern | Matches |
|---|---|
| `*dorsalApr*` | L only (L has `dorsalApr8`; R has `dorsal_RApr10`) |
| `*dorsal_RApr*` | R only |
| `grep -v "dorsal_R"` | correctly excludes R |
| `grep -v "_R"` | wrong — also strips `_Resnet50` |

## 3. Old vs new outputs: use snapshot suffix, not project date

The scorer string encodes the DLC project creation date (e.g., `dorsalApr8shuffle1`), which stays constant across retrains. The **snapshot number changes**: old outputs had `snapshot_best-10`, the 2026-04-12 retrain produced `snapshot_best-20`. When cleaning up old outputs to avoid pattern collisions (e.g., `whisker_L_3w_batch.m`'s `*dorsalApr*.csv` would match both), delete by snapshot suffix:

```
rm "<folder>/*dorsalApr8shuffle1_snapshot_best-10*"
```

**Why:** User caught me about to run `grep -v "_R"` (would have mis-identified L-side files). Separately, the retrain initially failed because I didn't clear old folders. All three lessons came from the same DLC-retrain session.

**How to apply:** For any DLC retrain workflow, always (a) wipe `dlc-models-pytorch/iteration-X` and `training-datasets/iteration-X` before calling `create_training_dataset`, (b) use `dorsal_R` (not `_R`) when L/R filtering, (c) trust snapshot suffix over Apr-date prefix when tracking output vintage.
