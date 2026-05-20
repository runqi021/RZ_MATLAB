"""
Two-branch CNN for ROI QC classification.

Branch A: 1D CNN on raw dFF trace → temporal features
Branch B: 2D CNN on avg patch + mask → spatial features
Merge:    concat → FC → sigmoid → keep probability

Usage:
    from model import ROI_QC_CNN
    net = ROI_QC_CNN(trace_len=2000, patch_size=64)
"""

import torch
import torch.nn as nn


class TraceBranch(nn.Module):
    """1D CNN on dFF trace [B, 1, T]."""

    def __init__(self, trace_len=2000, out_dim=128):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv1d(1, 32, kernel_size=15, stride=2, padding=7),
            nn.BatchNorm1d(32),
            nn.ReLU(inplace=True),

            nn.Conv1d(32, 64, kernel_size=11, stride=2, padding=5),
            nn.BatchNorm1d(64),
            nn.ReLU(inplace=True),

            nn.Conv1d(64, 128, kernel_size=7, stride=2, padding=3),
            nn.BatchNorm1d(128),
            nn.ReLU(inplace=True),

            nn.AdaptiveAvgPool1d(1),  # [B, 128, 1]
        )
        self.fc = nn.Linear(128, out_dim)

    def forward(self, x):
        # x: [B, 1, T]
        x = self.conv(x)
        x = x.squeeze(-1)  # [B, 128]
        return self.fc(x)


class PatchBranch(nn.Module):
    """2D CNN on avg patch + mask [B, 2, H, W]."""

    def __init__(self, patch_size=64, out_dim=128):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(2, 32, kernel_size=3, padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),

            nn.Conv2d(32, 64, kernel_size=3, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),

            nn.Conv2d(64, 128, kernel_size=3, padding=1),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),

            nn.AdaptiveAvgPool2d(1),  # [B, 128, 1, 1]
        )
        self.fc = nn.Linear(128, out_dim)

    def forward(self, x):
        # x: [B, 2, H, W]
        x = self.conv(x)
        x = x.view(x.size(0), -1)  # [B, 128]
        return self.fc(x)


class ROI_QC_CNN(nn.Module):
    """Two-branch CNN: trace + patch → keep probability."""

    def __init__(self, trace_len=2000, patch_size=64, branch_dim=128):
        super().__init__()
        self.trace_branch = TraceBranch(trace_len, branch_dim)
        self.patch_branch = PatchBranch(patch_size, branch_dim)

        self.classifier = nn.Sequential(
            nn.Linear(branch_dim * 2, 128),
            nn.ReLU(inplace=True),
            nn.Dropout(0.3),
            nn.Linear(128, 1),
        )

    def forward(self, trace, patch):
        """
        Args:
            trace: [B, 1, T]  raw dFF trace
            patch: [B, 2, H, W]  ch0=avg projection, ch1=binary mask
        Returns:
            logits: [B, 1]  (apply sigmoid for probability)
        """
        t_feat = self.trace_branch(trace)   # [B, branch_dim]
        p_feat = self.patch_branch(patch)   # [B, branch_dim]
        merged = torch.cat([t_feat, p_feat], dim=1)  # [B, 2*branch_dim]
        return self.classifier(merged)
