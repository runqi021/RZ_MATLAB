% popsel_gui_260816.m
% -----------------------------------------------------------------------
%  Include/exclude curation with a LIVE population readout.
%
%    left   : the cell's already-rendered per-cell summary
%    right  : population mean dF/F and event histogram over the INCLUDED set,
%             both on a normalised -1..+1 IBI axis, updating on every click
%    stats  : pooled occupancy-weighted Rayleigh over the included cells
%
%  WHY NORMALISED. Cells are averaged on tau/IBI, not seconds. Sert's IBI runs
%  1.74-2.63 s and the archive as a whole 0.43-2.83 s, so a seconds axis would let
%  slow breathers stretch the window and blur the transient. Same rule as the
%  260806 heatmaps: windows in breath cycles.
%
%  WHY THE POOLED STATISTIC IS A SUM. Every event keeps the occupancy weight of
%  ITS OWN recording (1/frames in that phase bin), so pooling across cells is
%  concatenation -- which on binned data is addition. Kish's n_eff follows from
%  the same sums. That makes the statistic exact rather than approximate when the
%  selection changes, and is why nothing is recomputed on a click.
%
%  Decisions are saved to CSV keyed on cell id, so the included set is what the
%  rest of the analysis consumes.
%
%  KEYS: i = include, e = exclude, u = undecided, left/right arrows = prev/next
%
%  Runqi Zhang / 2026-08-16
% -----------------------------------------------------------------------

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));

%% ===================== USER-EDITABLE =====================
sumRoot = 'D:\Ventral_surface_summary';
outDir  = fullfile(sumRoot,'popsel_260816');
% Which cached group to curate; one cache and one decisions CSV per group.
% 'Sert' is the genotype-wide group and now pools 0721 (INTACT) with 0828+0829
% (vagotomised) -- 114 cells across the manipulation. Its CSV holds the 48-cell
% intact 0721 curation. Use 'Sert_vagotomized' (built by
% analysis_260829\popsel_precompute_sessions_260831.m) for the vagotomised set.
GROUP   = 'Sert_vagotomized';
TRIG    = 'peak';        % 'peak' | 'onset' -- which trigger the population uses
altFigDirs = { fullfile(sumRoot,'per-cell-summary_active_260812','_below_active_gate_2evmin')
               fullfile(sumRoot,'per-cell-summary_active_260812','_superseded_by_merge_260814') };
% =========================================================

cacheF = fullfile(outDir, sprintf('popsel_cache_%s.mat', GROUP));
assert(isfile(cacheF), 'cache not found:\n  %s\nrun popsel_precompute_260816.m first', cacheF);
K = load(cacheF);
popsel_launch_260816(K.C, K.tauN, K.ctrN, GROUP, TRIG, outDir, altFigDirs);
