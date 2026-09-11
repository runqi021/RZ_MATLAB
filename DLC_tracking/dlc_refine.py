"""DLC refinement loop — refine outliers -> merge -> retrain -> re-analyze.

Active-learning loop to fix poorly tracked frames (e.g. a low-confidence
whisker tip). Run AFTER you have already extracted outlier frames (the GUI's
"extract outlier frames" step, or DO_EXTRACT below).

Order of operations:
    (extract outliers) -> refine_labels -> merge_datasets
                       -> create_training_dataset -> train_network -> analyze_videos

Refinement vs labeling:
    refine_labels() loads the network's PREDICTED positions (machinelabels-iter*.h5)
    so you DRAG the wrong dots to fix them. The plain "Label frames" tab shows
    blank frames and does NOT load predictions -- do not use it for refinement.
    merge_datasets() then folds the corrected outliers into CollectedData_<scorer>.

Run in the dlc310 conda env:
    conda activate dlc310
    python dlc_refine.py
"""

import os
import glob
import deeplabcut

# ============================ USER-EDITABLE ============================
# Project config.yaml (root of the DLC project).
CONFIG = r"D:\260615_thermalNbasler\260615_thermalNbasler-RZ-2026-06-15\config.yaml"

# Folder of videos to (optionally) extract outliers from and to re-analyze.
VIDEO_DIR = r"D:\260615_thermalNbasler\260615_thermalNbasler-RZ-2026-06-15\videos"

# Which steps to run. TWO-PASS by default so training never auto-starts:
#   PASS 1 (current): refine ONLY. Correct dots in napari, SAVE (Ctrl+S), close.
#   PASS 2: set DO_REFINE=False and DO_MERGE/DO_RETRAIN/DO_ANALYZE=True, re-run.
DO_EXTRACT  = False   # extract_outlier_frames (skip if GUI already did it)
DO_REFINE   = False   # refine_labels  -> opens napari editor (BLOCKS until closed)
DO_LABEL    = False   # label_frames -> STANDARD labeling GUI (use if refine_labels
                      #   won't save). Label the outlier frames fresh; writes
                      #   CollectedData directly so NO merge_datasets is needed.
DO_MERGE    = False   # merge_datasets -> ONLY needed after refine_labels, NOT after label_frames
# Split so hflip survives: create_training_dataset REGENERATES pytorch_config.yaml
# (wiping any hflip edit). Correct order for custom augmentation:
#   1) DO_CREATE_DATASET=True (run)  2) add hflip to pytorch_config.yaml
#   3) DO_TRAIN=True (run)  -- do NOT re-create in between.
DO_CREATE_DATASET = False  # create_training_dataset (packs labels; wipes hflip)
DO_TRAIN          = False  # train_network only (PyTorch defaults; keeps your hflip)
DO_ANALYZE        = False  # analyze_videos -> new CSVs alongside videos

OUTLIER_ALGO = 'jump'  # only used if DO_EXTRACT: 'jump' | 'uncertain' | 'fitting'
# ======================================================================


def find_videos(video_dir, exts=('avi', 'mp4', 'mov')):
    vids = []
    for ext in exts:
        vids.extend(glob.glob(os.path.join(video_dir, f'*.{ext}')))
    vids.sort()
    return vids


def main():
    assert os.path.isfile(CONFIG), f"config not found: {CONFIG}"
    vids = find_videos(VIDEO_DIR)
    print(f"config : {CONFIG}")
    print(f"videos : {len(vids)} found in {VIDEO_DIR}")

    if DO_EXTRACT:
        print(f"\n=== extract_outlier_frames (algo={OUTLIER_ALGO}) ===", flush=True)
        deeplabcut.extract_outlier_frames(CONFIG, vids, outlieralgorithm=OUTLIER_ALGO)

    if DO_REFINE:
        print("\n=== refine_labels (drag the wrong dots, SAVE Ctrl+S, then close) ===", flush=True)
        deeplabcut.refine_labels(CONFIG)

    if DO_LABEL:
        print("\n=== label_frames (label the outlier frames fresh, SAVE Ctrl+S) ===", flush=True)
        deeplabcut.label_frames(CONFIG)

    if DO_MERGE:
        print("\n=== merge_datasets (fold corrected outliers into training set) ===", flush=True)
        deeplabcut.merge_datasets(CONFIG)

    if DO_CREATE_DATASET:
        print("\n=== create_training_dataset (regenerates pytorch_config.yaml) ===", flush=True)
        deeplabcut.create_training_dataset(CONFIG)
        print(">>> If using hflip/custom augmentation, ADD it to pytorch_config.yaml NOW,", flush=True)
        print(">>> then run again with DO_TRAIN=True (and this step False).", flush=True)

    if DO_TRAIN:
        print("\n=== train_network (PyTorch defaults; uses current pytorch_config.yaml) ===", flush=True)
        deeplabcut.train_network(CONFIG)

    if DO_ANALYZE:
        print(f"\n=== analyze_videos: {len(vids)} videos ===", flush=True)
        deeplabcut.analyze_videos(CONFIG, vids, save_as_csv=True)

    print("\nDone.")


if __name__ == '__main__':
    main()
