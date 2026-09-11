---
name: feedback-matlab-offscreen-figure-reuse
description: "clf(hf) does not make hf the current figure, so subplot() draws into the last visible figure and batch exports come out blank"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 1f25dc5a-41b5-40bf-b366-9d77745c5a99
  modified: 2026-07-28T07:20:14.145Z
---

The common batch-figure pattern in this repo is

```matlab
hOff = figure('Color','w','Visible','off');
for i = 1:N
    if i <= showFirst, hf = figure('Visible','on'); else, hf = hOff; clf(hf); end
    subplot(1,5,1); ...            % <-- BUG
    exportgraphics(hf, ...);
end
```

`clf(hf)` clears `hf` but **does not make it the current figure**. `subplot()` and
bare `plot()` target `gcf`, which is still the last *visible* figure — so every
iteration after `showFirst` piles its content onto that one window and exports a
blank `hf`. The tell is a folder of PNGs with an identical small byte size
(4778 B in `breath_time\percell`, 28 of 34 files).

Fix, either:
- `clf(hf); set(0,'CurrentFigure',hf);` — sets current without raising the window; or
- pass the parent explicitly: `axes('Parent',hf,'Position',...)`.

`cell_activity_summary_260727.m` does NOT have this bug because it uses
`axes('Parent',hf,...)` throughout. `phase_rayleigh_percell_260727.m` DOES have it
(retired pipeline, left unfixed).

**How to apply:** whenever a batch-figure loop reuses one off-screen figure, check
the output file SIZES, not just the count — the script will happily report
"Saved 34 panels" while 28 of them are empty.
