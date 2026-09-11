#!/usr/bin/env python3
"""
whisker_roi_check.py  --  evidence for "does ONE per-scope ROI fit ALL runs?".

For each scope it builds a max projection of every run (first --nframes frames),
overlays that scope's cached whisker ROIs (whisker_rois.json), and montages the
runs side by side -> <scope>\whisker_roi_check.png. If the boxes enclose the
whisker marker in every run panel, one ROI per scope is fine; if a run's whisker
sits outside the box, that run (or scope) needs its own ROI.

  C:\\Users\\Admin\\.conda\\envs\\flir\\python.exe whisker_roi_check.py C:\\260613_breathing_thermalNbasler [--nframes 2000]
"""
import os, sys, glob, json, argparse
import numpy as np
import cv2


def avgproj(avi, nframes):
    cap = cv2.VideoCapture(avi); acc = None; c = 0
    while c < nframes:
        ok, fr = cap.read()
        if not ok:
            break
        g = fr[:, :, 0].astype(np.float64)
        acc = g if acc is None else acc + g
        c += 1
    cap.release()
    if acc is None:
        return None
    m = acc / c
    lo, hi = np.percentile(m, 1), np.percentile(m, 99)
    return np.clip((m - lo) / max(hi - lo, 1e-6), 0, 1).astype(np.float32).__mul__(255).astype(np.uint8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--nframes", type=int, default=2000)
    ap.add_argument("--panelw", type=int, default=320)
    args = ap.parse_args()

    scopes = {}
    for avi in glob.glob(os.path.join(args.root, "**", "*.avi"), recursive=True):
        folder = os.path.dirname(avi)
        if not os.path.isfile(os.path.join(folder, "timestamps.csv")):
            continue
        scope = os.path.relpath(folder, args.root).split(os.sep)[0]
        scopes.setdefault(scope, []).append((folder, avi))

    for scope, runs in sorted(scopes.items()):
        cache = os.path.join(args.root, scope, "whisker_rois.json")
        rois = None
        if os.path.isfile(cache):
            rois = json.load(open(cache))["rois"]
        panels = []
        for folder, avi in sorted(runs):
            mx = avgproj(avi, args.nframes)
            if mx is None:
                continue
            vis = cv2.cvtColor(mx, cv2.COLOR_GRAY2BGR)
            if rois:
                for nm, (x, y, w, h) in rois.items():
                    col = (0, 255, 0) if nm.endswith("L") else (0, 128, 255)
                    cv2.rectangle(vis, (x, y), (x + w, y + h), col, 3)
            sc = args.panelw / vis.shape[1]
            vis = cv2.resize(vis, (args.panelw, int(vis.shape[0] * sc)))
            cv2.putText(vis, os.path.basename(folder)[5:], (4, 16),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.4, (255, 255, 0), 1, cv2.LINE_AA)
            panels.append(vis)
            print(f"  {scope}/{os.path.basename(folder)}", flush=True)
        if not panels:
            continue
        h = max(p.shape[0] for p in panels)
        panels = [cv2.copyMakeBorder(p, 0, h - p.shape[0], 0, 6, cv2.BORDER_CONSTANT) for p in panels]
        ncol = min(4, len(panels))
        rows = []
        for i in range(0, len(panels), ncol):
            row = panels[i:i + ncol]
            while len(row) < ncol:
                row.append(np.zeros_like(panels[0]))
            rows.append(np.hstack(row))
        mont = np.vstack(rows)
        out = os.path.join(args.root, scope, "whisker_roi_check.png")
        cv2.imwrite(out, mont)
        print(f"[{scope}] saved {out}  ({len(runs)} runs)", flush=True)


if __name__ == "__main__":
    main()
