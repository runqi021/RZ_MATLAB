---
name: project_make_viewable_video_copies
description: make_viewable_video_copies.py — FFV1 AVI → H.264 mp4 viewing copies at the true fps; handles datasets with NO timestamps.csv by using 2P-triggered imaging fps
metadata: 
  node_type: memory
  type: project
  originSessionId: 3d6a9b12-9acd-41a0-b81b-4030e0e68644
  modified: 2026-07-23T04:44:09.157Z
---

`make_viewable_video_copies.py` (repo root) — makes **viewable** copies of Basler
behaviour videos at the **correct playback rate**. Complements, does not replace,
[[project_fix_avi_timing]].

**Two problems it fixes at once:**
1. Codec — raw runs are FFV1-in-AVI, which VLC / ImageJ / MATLAB `VideoReader`
   cannot decode (see [[feedback_ffv1_video_codec]]). Output is H.264/mp4.
2. Frame rate — the AVI header carries a nominal rate that is not the real one.

**Why a new script instead of fix_avi_timing.py:** that one requires a sibling
`timestamps.csv` and only re-times (lossless stream copy, stays FFV1, still
unviewable). Some experiments have **no timestamps.csv at all** — e.g.
`D:\260721_Sert_soma_G8s\phys` (27 runs, none). Also honours
[[feedback_no_inplace_edits]].

**fps priority chain:** `--fps` → sibling `timestamps.csv` → sibling `*_meta.mat`
`scanFrameRate_raw` (unrounded; NOT the rounded `fps` field) → skip.

**Key insight — the `_meta.mat` branch:** the Basler is 2P-triggered
frame-for-frame, so video frames == ScanImage `framesPerSlice` and the true
playback fps IS the imaging fps. The script *verifies* this per run and warns on
mismatch rather than silently trusting it. On 260721_Sert: exact match in 23/27;
fps was 30.0029 / 42.0883 / 47.0673 depending on run — the AVI header said 30 for
every one. Same triggering assumption as [[project_chat_analysis_breath_alignment]].

**Gotchas:**
- `-r F` BEFORE `-i` re-times correctly ONLY because this re-encodes. On a stream
  copy it just relabels the header and players ignore it — that is exactly why
  fix_avi_timing.py needs the `setts` bitstream filter instead.
- Filename frame counts (`_6000f_`) are UNRELIABLE — names get carried over
  between runs. Trust `_meta.mat` / the actual video, never the name.
- `CAP_PROP_FRAME_COUNT` is a header estimate. One 260721 source claimed 3045
  frames but only 3044 decode (truncated last frame). Count by decoding.

**Usage:** `python make_viewable_video_copies.py <parent|run|file.avi>`
`[--fps F] [--crf N] [--imagej] [--outdir D] [--overwrite] [--dry-run]`
Writes `<stem>_view.mp4` beside each source; originals untouched. `--imagej` also
emits MJPEG AVI for ImageJ installs lacking an mp4 reader.

Result on 260721_Sert_soma_G8s/phys: 27/27, **17.2 GB → 0.34 GB**, CRF 17 gives
~42 dB PSNR (visually indistinguishable), frame counts preserved, MATLAB
`VideoReader` reads them with exact fps/duration.
