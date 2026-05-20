"""
Run inference with trained ROI QC CNN.

Usage:
    python predict.py --data new_data.h5 --model best_model.pt --out predictions.mat

Reads HDF5 with /traces and /patches (same format as training, no /labels needed).
Saves predictions.mat with:
    prob      [N x 1]  keep probability (0-1)
    pred_keep [N x 1]  binary prediction (1=keep, 0=reject)
    threshold          decision threshold (default 0.5)
"""

import argparse
import h5py
import numpy as np
import torch
from scipy.io import savemat
from model import ROI_QC_CNN


def predict(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Load model
    checkpoint = torch.load(args.model, map_location=device, weights_only=False)
    trace_len  = checkpoint["trace_len"]
    patch_size = checkpoint["patch_size"]
    val_info = checkpoint.get('val_acc', checkpoint.get('metrics', {}).get('acc', '?'))
    print(f"Model trained: epoch {checkpoint['epoch']}, val_acc={val_info}")

    model = ROI_QC_CNN(trace_len=trace_len, patch_size=patch_size).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    # Load data (MATLAB column-major → transpose)
    with h5py.File(args.data, "r") as f:
        traces  = np.array(f["traces"],  dtype=np.float32).T                      # [N, T]
        patches = np.array(f["patches"], dtype=np.float32).transpose(3, 2, 1, 0)  # [N, 2, H, W]
        # Optional: pass through fov_id and roi_label if present
        fov_id    = np.array(f["fov_id"]).squeeze()    if "fov_id"    in f else None
        roi_label = np.array(f["roi_label"]).squeeze()  if "roi_label" in f else None

    N = len(traces)
    print(f"Predicting on {N} ROIs...")

    # Batch inference
    all_probs = []
    batch_size = args.batch_size

    with torch.no_grad():
        for i in range(0, N, batch_size):
            j = min(i + batch_size, N)
            t = torch.from_numpy(traces[i:j]).unsqueeze(1).to(device)   # [B, 1, T]
            p = torch.from_numpy(patches[i:j]).to(device)               # [B, 2, H, W]
            logits = model(t, p).squeeze(-1)
            probs = torch.sigmoid(logits).cpu().numpy()
            all_probs.append(probs)

    prob = np.concatenate(all_probs)
    pred_keep = (prob > args.threshold).astype(np.int8)

    n_keep = pred_keep.sum()
    n_reject = N - n_keep
    print(f"Predictions: {n_keep} keep, {n_reject} reject (threshold={args.threshold})")

    # Save
    out = {
        "prob": prob,
        "pred_keep": pred_keep,
        "threshold": args.threshold,
    }
    if fov_id is not None:
        out["fov_id"] = fov_id
    if roi_label is not None:
        out["roi_label"] = roi_label

    savemat(args.out, out)
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ROI QC CNN inference")
    parser.add_argument("--data", type=str, required=True,
                        help="HDF5 file with /traces and /patches")
    parser.add_argument("--model", type=str, required=True,
                        help="Path to best_model.pt")
    parser.add_argument("--out", type=str, required=True,
                        help="Output .mat file for predictions")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Decision threshold (default 0.5)")
    parser.add_argument("--batch_size", type=int, default=64)
    args = parser.parse_args()
    predict(args)
