"""
Train the ROI QC CNN on exported HDF5 data.

Usage (single split):
    python train.py --data_dir D:/roiQC_cnn_data

Usage (k-fold cross-validation, recommended):
    python train.py --data_dir D:/roiQC_cnn_data --kfold 5

Reads train.h5 and val.h5 from data_dir.
Reports per-class precision/recall/F1 and confusion matrix.
Saves best_model.pt (trained on ALL data after CV).
"""

import argparse
import os
import h5py
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from model import ROI_QC_CNN


# ── Dataset ──────────────────────────────────────────────────────────

class ROI_QC_Dataset(Dataset):
    """Load ROI data from HDF5 (handles MATLAB column-major transpose)."""

    def __init__(self, h5_path=None, augment=False, traces=None, patches=None, labels=None, fov_ids=None):
        if h5_path is not None:
            with h5py.File(h5_path, "r") as f:
                self.traces  = np.array(f["traces"],  dtype=np.float32).T
                self.patches = np.array(f["patches"], dtype=np.float32).transpose(3, 2, 1, 0)
                self.labels  = np.array(f["labels"],  dtype=np.float32).squeeze()
                self.fov_ids = np.array(f["fov_id"]).squeeze() if "fov_id" in f else np.zeros(len(self.labels))
        else:
            self.traces  = traces
            self.patches = patches
            self.labels  = labels
            self.fov_ids = fov_ids if fov_ids is not None else np.zeros(len(labels))
        self.augment = augment

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        trace = self.traces[idx]
        patch = self.patches[idx]
        label = self.labels[idx]

        if self.augment:
            T = len(trace)
            shift = np.random.randint(-T // 4, T // 4)
            trace = np.roll(trace, shift)
            if np.random.rand() > 0.5:
                patch = patch[:, :, ::-1].copy()
            if np.random.rand() > 0.5:
                patch = patch[:, ::-1, :].copy()

        trace = torch.from_numpy(trace).unsqueeze(0)
        patch = torch.from_numpy(patch)
        label = torch.tensor(label, dtype=torch.float32)
        return trace, patch, label


def merge_datasets(ds1, ds2):
    """Merge two datasets into one."""
    return ROI_QC_Dataset(
        traces  = np.concatenate([ds1.traces,  ds2.traces],  axis=0),
        patches = np.concatenate([ds1.patches, ds2.patches], axis=0),
        labels  = np.concatenate([ds1.labels,  ds2.labels],  axis=0),
        fov_ids = np.concatenate([ds1.fov_ids, ds2.fov_ids], axis=0),
    )


def subset_dataset(ds, indices, augment=False):
    """Create a subset of a dataset."""
    return ROI_QC_Dataset(
        traces  = ds.traces[indices],
        patches = ds.patches[indices],
        labels  = ds.labels[indices],
        fov_ids = ds.fov_ids[indices],
        augment = augment,
    )


# ── Metrics ──────────────────────────────────────────────────────────

def compute_metrics(labels, preds):
    """Confusion matrix + per-class precision/recall/F1."""
    tp = int(((preds == 1) & (labels == 1)).sum())
    tn = int(((preds == 0) & (labels == 0)).sum())
    fp = int(((preds == 1) & (labels == 0)).sum())
    fn = int(((preds == 0) & (labels == 1)).sum())

    prec_keep = tp / max(tp + fp, 1)
    rec_keep  = tp / max(tp + fn, 1)
    prec_rej  = tn / max(tn + fn, 1)
    rec_rej   = tn / max(tn + fp, 1)
    f1_keep = 2 * prec_keep * rec_keep / max(prec_keep + rec_keep, 1e-8)
    f1_rej  = 2 * prec_rej  * rec_rej  / max(prec_rej  + rec_rej,  1e-8)
    n = tp + tn + fp + fn
    acc = (tp + tn) / max(n, 1)
    baseline = max(labels.sum(), n - labels.sum()) / max(n, 1)

    return {"tp": tp, "tn": tn, "fp": fp, "fn": fn,
            "prec_keep": prec_keep, "rec_keep": rec_keep, "f1_keep": f1_keep,
            "prec_rej": prec_rej, "rec_rej": rec_rej, "f1_rej": f1_rej,
            "acc": acc, "baseline": baseline}


def print_metrics(m, prefix=""):
    print(f"{prefix}Confusion:  TP={m['tp']}  FP={m['fp']}  FN={m['fn']}  TN={m['tn']}")
    print(f"{prefix}  Keep:   prec={m['prec_keep']:.3f}  recall={m['rec_keep']:.3f}  F1={m['f1_keep']:.3f}")
    print(f"{prefix}  Reject: prec={m['prec_rej']:.3f}  recall={m['rec_rej']:.3f}  F1={m['f1_rej']:.3f}")
    print(f"{prefix}  Accuracy={m['acc']:.3f}  (majority baseline={m['baseline']:.3f})")


# ── Training loop ────────────────────────────────────────────────────

def train_one(train_ds, val_ds, args, device, save_path=None, tag=""):
    """Train one model, return best val metrics."""
    train_ds.augment = True
    val_ds.augment = False

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                              num_workers=0, pin_memory=True)
    val_loader   = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False,
                              num_workers=0, pin_memory=True)

    n_k = int((train_ds.labels == 1).sum())
    n_r = int((train_ds.labels == 0).sum())
    n_kv = int((val_ds.labels == 1).sum())
    n_rv = int((val_ds.labels == 0).sum())
    print(f"{tag}Train: {len(train_ds)} ({n_k} kept, {n_r} rej)  "
          f"Val: {len(val_ds)} ({n_kv} kept, {n_rv} rej)")

    pos_weight = torch.tensor([n_r / max(n_k, 1)], device=device)
    trace_len  = train_ds.traces.shape[1]
    patch_size = train_ds.patches.shape[2]
    model = ROI_QC_CNN(trace_len=trace_len, patch_size=patch_size).to(device)
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.wd)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5, factor=0.5)

    best_val_loss = float("inf")
    best_metrics = None
    patience_ctr = 0

    for epoch in range(1, args.epochs + 1):
        # Train
        model.train()
        t_loss, t_n = 0, 0
        for tr, pa, la in train_loader:
            tr, pa, la = tr.to(device), pa.to(device), la.to(device)
            logits = model(tr, pa).squeeze(-1)
            loss = criterion(logits, la)
            optimizer.zero_grad(); loss.backward(); optimizer.step()
            t_loss += loss.item() * len(la); t_n += len(la)

        # Validate
        model.eval()
        v_loss, v_n = 0, 0
        all_lab, all_pred = [], []
        with torch.no_grad():
            for tr, pa, la in val_loader:
                tr, pa, la = tr.to(device), pa.to(device), la.to(device)
                logits = model(tr, pa).squeeze(-1)
                v_loss += criterion(logits, la).item() * len(la); v_n += len(la)
                all_lab.append(la.cpu().numpy())
                all_pred.append((torch.sigmoid(logits) > 0.5).float().cpu().numpy())

        v_loss /= v_n
        m = compute_metrics(np.concatenate(all_lab), np.concatenate(all_pred))
        scheduler.step(v_loss)

        print(f"{tag}Ep {epoch:3d}  loss={t_loss/t_n:.4f}/{v_loss:.4f}  "
              f"acc={m['acc']:.3f}  keep_F1={m['f1_keep']:.3f}  rej_F1={m['f1_rej']:.3f}  "
              f"lr={optimizer.param_groups[0]['lr']:.1e}")

        if v_loss < best_val_loss:
            best_val_loss = v_loss
            best_metrics = m
            patience_ctr = 0
            if save_path:
                torch.save({"epoch": epoch, "model_state_dict": model.state_dict(),
                             "val_loss": v_loss, "metrics": m,
                             "trace_len": trace_len, "patch_size": patch_size}, save_path)
        else:
            patience_ctr += 1
            if patience_ctr >= args.patience:
                print(f"{tag}Early stop at epoch {epoch}")
                break

    print(f"\n{tag}Best val:")
    print_metrics(best_metrics, prefix=f"{tag}  ")
    return best_metrics


# ── K-fold CV ────────────────────────────────────────────────────────

def kfold_cv(all_ds, args, device, k=5):
    """K-fold CV split by FOV. Trains final model on all data."""
    unique_fovs = np.unique(all_ds.fov_ids)
    np.random.seed(42)
    np.random.shuffle(unique_fovs)

    fold_metrics = []
    for fold in range(k):
        print(f"\n{'='*60}\n  FOLD {fold+1}/{k}\n{'='*60}")
        i0 = fold * len(unique_fovs) // k
        i1 = (fold + 1) * len(unique_fovs) // k
        val_fovs = unique_fovs[i0:i1]

        val_idx   = np.where(np.isin(all_ds.fov_ids, val_fovs))[0]
        train_idx = np.where(~np.isin(all_ds.fov_ids, val_fovs))[0]

        train_sub = subset_dataset(all_ds, train_idx, augment=True)
        val_sub   = subset_dataset(all_ds, val_idx,   augment=False)

        m = train_one(train_sub, val_sub, args, device, tag=f"[F{fold+1}] ")
        fold_metrics.append(m)

    # Summary
    print(f"\n{'='*60}\n  {k}-FOLD CV SUMMARY\n{'='*60}")
    for key in ["acc", "f1_keep", "f1_rej", "prec_keep", "rec_keep", "prec_rej", "rec_rej"]:
        vals = [m[key] for m in fold_metrics]
        print(f"  {key:15s}: {np.mean(vals):.3f} +/- {np.std(vals):.3f}  "
              f"(range {np.min(vals):.3f}-{np.max(vals):.3f})")
    print(f"  {'baseline':15s}: {np.mean([m['baseline'] for m in fold_metrics]):.3f}")

    # Train final model on ALL data (use first fold's val as placeholder)
    print(f"\n{'='*60}\n  FINAL MODEL (trained on all {len(all_ds)} ROIs)\n{'='*60}")
    # Use 90/10 split of all data for final training with early stopping
    n_all = len(all_ds)
    perm = np.random.permutation(n_all)
    n_tr = int(n_all * 0.9)
    final_train = subset_dataset(all_ds, perm[:n_tr], augment=True)
    final_val   = subset_dataset(all_ds, perm[n_tr:], augment=False)
    save_path = os.path.join(args.out_dir, "best_model.pt")
    train_one(final_train, final_val, args, device, save_path=save_path, tag="[Final] ")
    print(f"\nModel saved: {save_path}")

    return fold_metrics


# ── Main ─────────────────────────────────────────────────────────────

def main(args):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    os.makedirs(args.out_dir, exist_ok=True)

    ds_train = ROI_QC_Dataset(os.path.join(args.data_dir, "train.h5"))
    ds_val   = ROI_QC_Dataset(os.path.join(args.data_dir, "val.h5"))

    if args.kfold > 1:
        all_ds = merge_datasets(ds_train, ds_val)
        n_k = int((all_ds.labels == 1).sum())
        n_r = int((all_ds.labels == 0).sum())
        print(f"Total: {len(all_ds)} ROIs ({n_k} kept, {n_r} rejected, "
              f"{len(np.unique(all_ds.fov_ids))} FOVs)")
        kfold_cv(all_ds, args, device, k=args.kfold)
    else:
        save_path = os.path.join(args.out_dir, "best_model.pt")
        train_one(ds_train, ds_val, args, device, save_path=save_path)
        print(f"\nModel saved: {save_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train ROI QC CNN")
    parser.add_argument("--data_dir", type=str, required=True)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--wd", type=float, default=1e-4)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--kfold", type=int, default=1, help="K-fold CV (default 1 = single split)")
    args = parser.parse_args()
    if args.out_dir is None:
        args.out_dir = os.path.join(args.data_dir, "models")
    main(args)
