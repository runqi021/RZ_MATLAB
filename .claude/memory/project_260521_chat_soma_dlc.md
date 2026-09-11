---
name: project-260521-chat-soma-dlc
description: "260521_ChAT_soma_G8s DLC breathing project — created & trained 2026-05-26, 3 dots, tracking CSVs done"
metadata: 
  node_type: memory
  type: project
  originSessionId: fd1aadc9-ba37-4459-9e34-5423605f4d93
---

DLC breathing-dot tracking for experiment **260521_ChAT_soma_G8s** (created 2026-05-26).

- Data: `Z:\Runqi\behavior\260521_ChAT_soma_G8s` — 10 Basler AVIs, 512×512, **98.6 fps, 3001 frames (~30.4 s)** each.
- DLC project: `...\260521_ChAT_soma_G8s-RZ-2026-05-26\config.yaml`
- 3 bodyparts: `dot1, dot2, dot3`; single animal; 114 train / 6 test labeled images.
- Model: `DLC_Resnet50_260521_ChAT_soma_G8sMay26shuffle1`, snapshot `best-10` (PyTorch engine, resnet50_gn default).
- Eval: test RMSE 1.39 px, mAP/mAR 100. Tracking CSV+h5+pickle written alongside videos.

Ran via `DLC_tracking/dlc_1_setup_and_label.py` then `dlc_2_train_and_analyze.py` (USER SETTINGS edited to this path).

Next: pair CSVs to calcium FOVs ([[reference_pair_behavior_to_phys]]) and compute breathing trace.
