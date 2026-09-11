% render_picked_cells_260817.m
% -----------------------------------------------------------------------
%  Re-render the hand-picked cells of ppt_260815_picked_cells.csv with the
%  CURRENT figure: two triggered averages (onset above peak), no histograms, one
%  trigger line per panel, chronological heatmap rows.
%
%  WHY A SCRIPT AND NOT THE GUI. The picks already exist as a table -- cell,
%  recording, roi and the exact trace window each figure used -- so re-rendering
%  them by hand through the GUI would be re-entering numbers that are already
%  written down, once per cell, with a chance of a typo each time. This reads the
%  windows from the csv, so the new figures differ from the old ones ONLY in the
%  panel changes, not in what part of the trace they show.
%
%  IT IS THE GUI'S PARAMETER BLOCK, COPIED. Every value below is taken from
%  temporal_phase_cell_gui_260812.m, because these renders go into the same
%  folder as the GUI's and have to be indistinguishable from them. If the GUI
%  block changes, this one has to be re-copied -- there is no shared source for
%  it, and pretending otherwise would be worse than saying so here.
%
%  POOLING IS FROM THE REGISTRY, NOT THE CSV. The csv names ONE recording per
%  pick (the one whose trace is drawn); a cell seen in several recordings must
%  still pool all of them, so the observation list is rebuilt from
%  event_latency_data.mat with the same merge overrides per_cell_summary applies.
%  P.traceObs then selects the csv's recording out of that list.
%
%  Runqi Zhang / 2026-08-17
% -----------------------------------------------------------------------

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));
addpath(fullfile(repoRoot,'analysis_260806'));
addpath(fullfile(repoRoot,'2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));
addpath(fullfile(repoRoot,'falloff-analysis-260805'));

%% ===================== USER-EDITABLE =====================
sumRoot = 'D:\Ventral_surface_summary';
regFile = fullfile(sumRoot,'event_latency_260811','event_latency_data.mat');
pickCsv = fullfile(sumRoot,'ppt_260815_picked_cells.csv');
outDir  = fullfile(sumRoot,'gui_renders_260815');
onlyCells = [];        % [] = every pick; otherwise a list of cell numbers
pngDpi  = 600;
savePDF = true;
% =========================================================

% Software OpenGL for the whole batch: the NVIDIA path throws inside
% exportgraphics ~1% of the time and leaves a 9 KB blank behind.
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if ~isfolder(outDir), mkdir(outDir); end

%% ---- GUI parameter block, copied verbatim ----
P = struct();
P.doCoh        = false;
P.nDrop        = 30;          P.fallback_fps = 30;
P.TW_spec      = 6;           P.alpha_sig    = 0.01;
P.minSpikes    = 2;           P.ca_lag_sec   = 0;
P.f_breath_search = [0.2 4];  P.fwhm_factor  = 0.6;   P.min_bw = 0.05;
P.fmin         = 0.05;        P.fmax         = 15;
P.trigWin_sec  = [];          P.trigWinIBI   = 2;
P.ylim_dff     = [];          P.ylim_epc     = [];
P.histBinFrames= 2;
P.nShuffle     = 1200;        P.shiftMinCyc  = 3;
P.pad_um       = 20;          P.clip_pct     = [1 99.5];
P.crop_um      = 70;          P.cropPxUm     = 0.5;
P.scalebar_um  = 10;
P.showRayleigh = false;       P.showAcqTime  = true;
P.trialOverlay = true;        P.dffSpecDeriv = false;
P.guiLayout    = true;        P.boldGenoOnly = true;
P.spkWinIBI    = 1;
P.scalebarThickFrac = 0.04;   P.outlineFrac  = 0.008;
P.outlineAlpha = 0.6;         P.traceBreathCol = [0 0 0];
P.outlineImgPx = 1.2;
P.rayPhaseBins = 36;          P.nShuffleRayleigh = 0;
P.gamma_val    = 1;           P.PixelSizeBase = 1.7778;  P.outlineLW = 0.8;
P.sortMode     = 'none';      P.trace_xlim_sp = [];
P.dffColor     = [0.2 0.7 0.2];
P.onsetCol     = [0.90 0.10 0.10];
P.peakCol      = [0.35 0.75 1.00];
P.TW_coh       = 4;           P.alpha_coh    = 0.001;

%% ---- registry + merges, exactly as per_cell_summary_260812 does it ----
R = load(regFile);
CELL = R.CELL; OBS = R.OBS; REC = R.REC;
[obsMerged, ~, origCells] = ...
    apply_cell_merges_260814(CELL, OBS, cell_merge_overrides_260814());
cellNum = zeros(numel(obsMerged),1);
CELL2 = struct('key',{},'obs',{});
for i = 1:numel(obsMerged)
    cellNum(i)   = min(origCells{i});
    CELL2(i).obs = obsMerged{i};
    CELL2(i).key = CELL(cellNum(i)).key;
end
CELL = CELL2;
fprintf('registry: %d cells after merges\n', numel(CELL));

%% ---- the picks ----
T = readtable(pickCsv,'TextType','string');
if ~isempty(onlyCells), T = T(ismember(T.cell, onlyCells), :); end
% One figure per CELL. The csv has one row per slide image, and a cell can appear
% on more than one slide; the first row wins and the rest are reported, because
% two renders of one cell would differ only by trace window and overwrite each
% other under the same filename.
[~, keep] = unique(T.cell, 'stable');
if numel(keep) < height(T)
    fprintf(2,'%d duplicate cell row(s) in the csv, keeping the first of each\n', ...
            height(T)-numel(keep));
end
T = T(keep,:);
fprintf('%d picked cells to render\n\n', height(T));

ok = 0; failed = strings(0,1);
for k = 1:height(T)
    cid = T.cell(k);
    ci  = find(cellNum == cid, 1);
    if isempty(ci)
        fprintf(2,'cell %d: not in the registry -- skipped\n', cid);
        failed(end+1) = sprintf('cell %d (not in registry)', cid); %#ok<SAGROW>
        continue;
    end

    % every recording this cell appears in
    O = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
    for o = CELL(ci).obs(:)'
        parts = regexp(OBS(o).label, '/', 'split');
        if numel(parts) < 3, continue; end
        roi = str2double(parts{end});
        fp  = REC(OBS(o).rec).folder;
        if ~isfinite(roi) || roi < 1 || ~isfolder(fp), continue; end
        O(end+1) = struct('folder',fp, 'roi',roi, ...
                          'recName',strjoin(parts(3:end-1),'/'), ...
                          'group',parts{1}, 'recDate',parts{2}); %#ok<SAGROW>
    end
    if isempty(O)
        fprintf(2,'cell %d: no usable recordings -- skipped\n', cid);
        failed(end+1) = sprintf('cell %d (no recordings)', cid); %#ok<SAGROW>
        continue;
    end

    % which of them the csv drew the trace from
    it = find(strcmp({O.recName}, char(T.recording(k))), 1);
    if isempty(it)
        fprintf(2,'cell %d: csv recording "%s" not among its %d observation(s) -- using the first\n', ...
                cid, T.recording(k), numel(O));
        it = 1;
    end

    Pr = P;
    Pr.cellId   = cid;
    Pr.traceObs = it;
    Pr.trace_xlim_sp = [T.trace_start_s(k), T.trace_end_s(k)];

    stem = regexprep(sprintf('%s_cell%03d_ROI%02d', O(it).recName, cid, O(it).roi), ...
                     '[\\/:*?"<>|]', '_');
    base = fullfile(outDir, stem);
    fprintf('[%2d/%2d] cell %3d  %-55s roi %2d  %.0f-%.0f s  (%d rec)\n', ...
            k, height(T), cid, O(it).recName, O(it).roi, ...
            T.trace_start_s(k), T.trace_end_s(k), numel(O));
    try
        fg = temporal_phase_cell_fig_260812(O, Pr);
        exportgraphics(fg, [base '.png'], 'Resolution',pngDpi, 'BackgroundColor','white');
        if savePDF
            exportgraphics(fg, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
        end
        close(fg);
        d = dir([base '.png']);
        % A 9,223-byte png is the known OpenGL export stub, not a figure. Check
        % SIZE, not just existence -- a stub file exists perfectly happily.
        if isempty(d) || d.bytes < 20000
            fprintf(2,'   SUSPECT EXPORT: %d bytes\n', max([d.bytes 0]));
            failed(end+1) = sprintf('cell %d (export %d bytes)', cid, max([d.bytes 0])); %#ok<SAGROW>
        else
            ok = ok + 1;
        end
    catch ME
        fprintf(2,'   FAILED: %s\n', ME.message);
        failed(end+1) = sprintf('cell %d (%s)', cid, ME.message); %#ok<SAGROW>
        close(findall(0,'Type','figure'));
    end
end

fprintf('\n%d of %d rendered into %s\n', ok, height(T), outDir);
if ~isempty(failed)
    fprintf(2,'%d problem(s):\n', numel(failed));
    fprintf(2,'   %s\n', failed{:});
end
