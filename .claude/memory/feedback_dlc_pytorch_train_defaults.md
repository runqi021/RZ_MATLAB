---
name: DLC train_network — call with defaults (PyTorch engine)
description: On this project DLC uses the PyTorch engine; do not pass maxiters/saveiters/displayiters kwargs to train_network
type: feedback
originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---
Call `deeplabcut.train_network(config)` with NO kwargs. Let the training schedule come from `pytorch_config.yaml` defaults.

**Why:** DLC 3.x PyTorch engine does not accept the TF-era kwargs `maxiters`, `saveiters`, `displayiters`. Passing them makes training fail to start. User said: "it will never run like this, so we returned back to default."

**How to apply:** In `DLC_tracking/dlc_2_train_and_analyze.py` (and any similar wrapper), keep the call bare: `deeplabcut.train_network(config)`. If the user wants a custom schedule, edit the project's `pytorch_config.yaml` — not the Python wrapper. Same logic applies to other analogous calls (e.g., don't pass TF-only kwargs to PyTorch-engine APIs).
