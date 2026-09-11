---
name: DLC environment, GPU, and performance notes
description: Conda env path, GPU model, and what does (and doesn't) speed up DLC inference on this machine
type: reference
originSessionId: 8ad2ee8d-8d96-4c84-b84e-fab65cbda69e
---
## Environment

- Conda env: `C:\Users\Admin\.conda\envs\dlc310\python.exe` (DLC 3.0.0rc13, PyTorch engine, cv2 available).
- Run batch files: `DLC_tracking\run_train.bat` calls `dlc_2_train_and_analyze.py`.

## GPU

- **NVIDIA RTX A1000 6 GB Laptop GPU** (driver 581.60).
- Modest Ampere mobile workstation card. Tensor cores present but limited throughput.
- DLC inference: ~14 it/s @ batch size 8 on 680×682 frames → ~6 min per 5400-frame video.

## Speedups that don't help meaningfully on this card

- **FP16 autocast** (`inference.autocast.enabled: true` in `pytorch_config.yaml`): 0% speedup measured. GPU already compute-bound; tensor-core gain eaten by small batch / memory bandwidth.
- **Bigger batch size**: limited by 6 GB VRAM (only ~4 GB free at batch 8).
- **torch.compile**: not tested; compile overhead may not pay off over 30 short videos.

## What to do instead

- Accept ~30 min for full analyze of 30 videos. Don't over-tune — the real bottleneck is the card itself.
- For the custom label-video overlay (read frame + draw crosses + write mp4), OpenCV is single-threaded; ~1 min per video is normal for the A1000 running on CPU side.

**How to apply:** When the user complains DLC analyze is slow, check the GPU first with `nvidia-smi`. If util is 100 %, it's compute-bound and there's no cheap win. Report this rather than suggesting ineffective tweaks like autocast.
