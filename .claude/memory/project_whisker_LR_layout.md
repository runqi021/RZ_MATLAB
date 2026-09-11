---
name: Whisker L/R DLC projects, bodyparts, and CSV layout
description: Two-CSV per-video DLC layout, current bodypart sets (vL1/vL3 L-side retrained, vR1/vR2 R-side), DLC project locations and PyTorch training details
type: project
originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---
## Two-CSV per-video layout

Each Basler whisker video is tracked by two DLC networks → two CSVs in the same folder.

| Side | CSV filename pattern | Bodyparts | DLC project folder |
|---|---|---|---|
| **L** | `*DLC_Resnet50_whisker_dorsalApr*.csv` (excludes `_R`) | `vL1, vL3` (retrained 2026-04-12 — vL2 dropped, vL4 dropped) | `C:\Users\Admin\Desktop\whisker_dorsal-RZ-2026-04-08\` |
| **R** | `*DLC_Resnet50_whisker_dorsal_RApr*.csv` | `vR1, vR2` | `C:\260407_KA_electro_dorsal_whisking\whisker_dorsal_R-RZ-2026-04-10\` |

Run stem: strip `DLC_Resnet50.*$` from CSV filename. Same stem on both sides → paired run. Data folder `C:\260407_KA_electro_dorsal_whisking\` has 30 runs (run001..run030) at fps=90 Hz, 5400 frames each, 680×682 px.

## DLC scorer naming caveat

The scorer string in filenames encodes the **project creation date**, not the training date. The L project was created 2026-04-08 so every L CSV keeps `dorsalApr8shuffle1` in the scorer string even after retraining. The snapshot number (e.g., `best-10` vs `best-20`) is what changes between training runs — use that to tell old outputs from new ones.

## PyTorch engine training

DLC 3.x PyTorch engine. Backbone **ResNet-50 GN** (group-norm), batch size 8, AdamW lr=5e-4 with LRListScheduler (drops at epochs 90/120), 200 epochs, snapshots every 25. `pytorch_config.yaml` sits under `dlc-models-pytorch/iteration-0/<shuffle>/train/`. The L-side retrain on 2026-04-12 converged early (test RMSE ≈ 1.06 px at epoch 130) and was stopped manually — `snapshot-best-020.pt` was the model used to analyze.

**Why:** Tracked by two independently-trained DLC networks because the R side was added later (`_RApr10`). The L network was retrained on 2026-04-12 with a 3-bodypart label set (vL1, vL2, vL3) but vL2 is currently excluded from analysis.

**How to apply:** For whisker analysis over this dataset, use the patterns above to pair L and R CSVs. When reading config, always look in the project folder (on Desktop for L, inside data folder for R) — not inside the data folder for both. When DLC's `analyze_videos` outputs look confusingly similar to old ones, check the snapshot suffix, not the "Apr8" prefix.
