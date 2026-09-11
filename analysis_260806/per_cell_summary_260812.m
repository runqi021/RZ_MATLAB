% per_cell_summary_260812.m
% -----------------------------------------------------------------------
%  One temporal-phase summary figure PER CELL, for every archive cell whose
%  pooled Rayleigh log Z clears a threshold.
%
%  THE CELL IS THE UNIT. A cell imaged in several recordings gets ONE figure,
%  with the histograms, triggered averages, per-cycle heatmap and log Z pooled
%  over all of its recordings -- events and cycles concatenated, each recording
%  keeping its own occupancy weights and its own trigger set. The wide dF/F
%  trace and the avg-projection cannot be pooled, so they come from one
%  recording (the first), showing 30 s from the middle of it.
%
%  THE CUT IS APPLIED TO THE NUMBER THE FIGURE PRINTS. A cheap stats-only pass
%  computes each cell's pooled log Z first, using exactly the code that draws the
%  figure, and only cells above logZ_min are rendered. Selecting on the archive's
%  ray_logZ while printing a locally computed one is what previously put logZ<1
%  pages in a logZ>1 folder -- the two are different quantities (the archive
%  value is pooled with a 30 ms ca-lag, and the two ventral CSVs disagree with
%  each other). Here there is NO ca-lag anywhere.
%
%  Renders through analysis_260727\coh_ca_breath\temporal_phase_cell_fig_260812.m,
%  the same function the GUI uses, so the batch and the hand-tuned GUI figure are
%  identical apart from the window/limits the GUI lets you set.
%
%  INPUTS
%    event_latency_data.mat        CELL(.key,.obs) / OBS(.rec,.label) / REC(.folder)
%    breath_time_summary_cells.csv ray_logZ per cell_key -- prefilter only
%
%  Runqi Zhang / 2026-08-12
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot, 'analysis_260727', 'coh_ca_breath'));
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
addpath(fullfile(repoRoot, 'falloff-analysis-260805'));

%% ===================== USER-EDITABLE =====================
rootPath  = 'D:\Ventral_surface_summary';
inFile    = fullfile(rootPath, 'event_latency_260811', 'event_latency_data.mat');
rayCsv    = fullfile(rootPath, 'breath_time_summary_260808', 'breath_time_summary_cells.csv');

% WHICH CELLS. Both criteria are applied to the cell POOLED over its recordings,
% and both are measured by the same stats pass that the figures come from.
%   'logz'   figure_pooled_logZ >= logZ_min
%   'active' pooled nnz(spike_train>0) > activeMinEv -- the archive's ACTIVE
%            definition, straight out of Ventral_surface_breath_time_summary_260808
selMode     = 'active';
logZ_min    = 1;

% ACTIVE = pooled firing rate >= 2 events/min.
%   Was "> 5 events" (activeMinEv), an absolute count that ignored how long the
%   recording ran: durations here span 49 s to 995 s, so a flat 6-event floor
%   meant anywhere from 7.3 to 0.36 events/min purely as a function of length.
%   Rate is pooled (nSpikes/durTot), like every other pooled quantity.
%   MUST match activeMinRateHz in pop_features_260813.m, or T and this folder
%   describe different populations and the IBI split silently mismatches.
activeMinRateHz = 2/60;

% Render only these cell NUMBERS (the ..._cell%03d in the filename). [] = all.
% Use for repairing individual cells without repeating the ~2 h measure pass over
% all 269 -- it restricts the measure pass too, so the parameters are identical
% to the full run by construction rather than by a copied parameter block.
onlyCells   = [];

% Optional breathing-rate gate, applied to the cell's POOLED MEDIAN IBI -- the
% same number printed on the figure. Median not mean: the cycle distribution has
% a tail to 32 s from missed onsets, which a mean would follow. [] = no gate.
ibiMin      = [];
if strcmpi(selMode,'active')
    outDir = fullfile(rootPath, 'per-cell-summary_active_260812');
else
    outDir = fullfile(rootPath, 'per-cell-summary_260812');
end
if ~isempty(ibiMin)
    outDir = fullfile(outDir, sprintf('IBI_ge%g', ibiMin));
end

logZ_pool = -Inf;           % cheap prefilter on the archive's ray_logZ; -Inf = every cell
traceWin  = [];             % [] = 30 s from the middle of the trace recording
overwrite = true;           % false = skip cells already rendered (resumable)
savePDF   = true;
saveProj  = true;

% Figure parameters -- KEEP IN STEP with temporal_phase_cell_gui_260812.m.
P = struct();
P.doCoh        = false;
P.nDrop        = 30;          P.fallback_fps = 30;
P.TW_spec      = 6;           P.alpha_sig    = 0.01;
P.minSpikes    = 2;           P.ca_lag_sec   = 0;     % NO ca-lag, anywhere
P.f_breath_search = [0.2 4];  P.fwhm_factor  = 0.6;   P.min_bw = 0.05;
P.fmin         = 0.05;        P.fmax         = 15;
P.trigWin_sec  = [];          P.trigWinIBI   = 2;
P.ylim_dff     = [];          P.ylim_epc     = [];
P.histBinFrames= 2;
P.nShuffle     = 1200;
P.shiftMinCyc  = 3;
P.pad_um       = 20;          P.clip_pct     = [0.5 99.9];
P.scalebar_um  = 50;
P.rayPhaseBins = 36;
% gamma 1 = LINEAR, no tone curve. The crops are small (some 197x379 px) and the
% contrast stretch is already set by the percentiles of that little window, so a
% gamma on top of it distorts more than it reveals.
P.gamma_val    = 1;           P.PixelSizeBase = 1.7778;  P.outlineLW = 0.8;
P.sortMode     = 'none';
P.dffColor     = [0.2 0.7 0.2];
P.onsetCol     = [0.90 0.10 0.10];
P.peakCol      = [0.35 0.75 1.00];
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if ~isfolder(outDir), mkdir(outDir); end

% Software OpenGL for the whole batch. The hardware path (NVIDIA WGL) throws
%   com.jogamp.opengl.GLException: Error making context ... current
% intermittently inside exportgraphics, roughly once per 100 figures on this
% machine. It is not deterministic and not data-dependent: the 260814 run lost
% cell053 outright and wrote cell054 and cell158 as identical 9,223-byte stubs,
% and took the whole MATLAB session down once. A stub is the dangerous case --
% the file exists, so a folder listing looks complete. Software rendering is
% slower per figure but does not touch the driver.
try
    opengl('software');
catch ME
    fprintf(2,'could not force software OpenGL: %s\n', ME.message);
end

%% ===================== BUILD THE CELL LIST =====================
assert(isfile(inFile), 'Not found:\n  %s\nRun Ventral_surface_event_latency_260811.m first.', inFile);
Din  = load(inFile, 'CELL', 'OBS', 'REC');
CELL = Din.CELL;  OBS = Din.OBS;  REC = Din.REC;
assert(isfield(REC,'folder'), 'REC has no .folder -- rerun Ventral_surface_event_latency_260811.m');

rayLogZ = nan(numel(CELL),1);
if isfile(rayCsv)
    Bt = readtable(rayCsv);
    [tf, loc] = ismember(string({CELL.key}), string(Bt.cell_key));
    rayLogZ(tf) = Bt.ray_logZ(loc(tf));
    fprintf('archive ray_logZ joined for %d of %d cells\n', nnz(tf), numel(CELL));
end

%% ---- cell identity corrections -------------------------------------------
% cell_link left some genuinely identical cells unmerged: it never ran at all on
% Vglut2/0224 and Vglut2/1124 (their curations exist but were never joined), and
% on Vglut2/0810 and Sst/0807 it ran and still missed a pair. Fold those in here.
% The archive on disk is untouched -- the correction lives in a versioned .m file.
%
% A merged cell inherits the LOWEST of its constituent cell numbers, so every
% unmerged cell keeps the number it already had and its filename does not move.
% Its archive ray_logZ is set NaN: the CSV holds one value per OLD cell, so no
% single archive number describes the merged cell. figure_pooled_logZ is measured
% fresh from the pooled data below and is what the figure and the selection use.
[obsMerged, mergeInfo, origCells] = ...
        apply_cell_merges_260814(CELL, OBS, cell_merge_overrides_260814());
nMergedGroups = nnz(strcmp({mergeInfo.status}, 'merged'));
for q = 1:numel(mergeInfo)
    if ~strcmp(mergeInfo(q).status,'merged')
        fprintf(2,'  merge group %d: %s\n', mergeInfo(q).group, mergeInfo(q).status);
    end
end
CELL2 = struct('key',{},'obs',{});
cellNum = zeros(numel(obsMerged),1);
ray2    = nan(numel(obsMerged),1);
for i = 1:numel(obsMerged)
    cellNum(i)  = min(origCells{i});
    CELL2(i).obs = obsMerged{i};
    CELL2(i).key = CELL(cellNum(i)).key;
    if numel(origCells{i}) == 1, ray2(i) = rayLogZ(origCells{i}); end
end
fprintf('cell merges: %d groups applied, %d cells -> %d\n', ...
        nMergedGroups, numel(CELL), numel(CELL2));
CELL = CELL2;  rayLogZ = ray2;
if isfinite(logZ_pool), selCell = find(rayLogZ > logZ_pool);
else,                   selCell = (1:numel(CELL))'; end
fprintf('%d cells enter the measurement pass\n', numel(selCell));

% One entry per CELL, carrying every recording it appears in.
CJ = struct('cell',{},'key',{},'group',{},'obs',{},'pooled',{});
for ii = 1:numel(selCell)
    c = selCell(ii);
    O = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
    for o = CELL(c).obs(:)'
        parts = regexp(OBS(o).label, '/', 'split');
        if numel(parts) < 3, continue; end
        roi = str2double(parts{end});
        fp  = REC(OBS(o).rec).folder;
        if ~isfinite(roi) || roi < 1 || ~isfolder(fp), continue; end
        O(end+1) = struct('folder',fp, 'roi',roi, ...
                          'recName',strjoin(parts(3:end-1),'/'), ...
                          'group',parts{1}, 'recDate',parts{2}); %#ok<SAGROW>
    end
    if isempty(O), continue; end
    % cellNum, not the loop index: a merged cell keeps the lowest of its original
    % numbers so unmerged cells' filenames stay exactly as they were.
    CJ(end+1) = struct('cell',cellNum(c), 'key',CELL(c).key, 'group',O(1).group, ...
                       'obs',O, 'pooled',rayLogZ(c)); %#ok<SAGROW>
end
fprintf('%d cells with usable recordings (%d observations total)\n\n', ...
        numel(CJ), sum(arrayfun(@(x) numel(x.obs), CJ)));

if ~isempty(onlyCells)
    wanted = ismember([CJ.cell], onlyCells);
    missing = setdiff(onlyCells, [CJ.cell]);
    if ~isempty(missing)
        fprintf(2,'onlyCells: %s not in the cell list (merged away, or dropped)\n', ...
                mat2str(missing));
    end
    CJ = CJ(wanted);
    fprintf('onlyCells: restricted to %d cell(s): %s\n\n', numel(CJ), mat2str([CJ.cell]));
    assert(~isempty(CJ), 'onlyCells matched no cells');
end

%% ============== MEASURE FIRST, THEN SELECT ON WHAT WE PRINT ==============
figLogZ = nan(numel(CJ),1);
figNspk = nan(numel(CJ),1);
figRate = nan(numel(CJ),1);   % pooled events/s, for the ACTIVE gate
nRecs   = arrayfun(@(x) numel(x.obs), CJ)';
tM = tic;
for j = 1:numel(CJ)
    try
        Pm = P; Pm.statsOnly = true;
        [~,~,st] = temporal_phase_cell_fig_260812(CJ(j).obs, Pm);
        figLogZ(j) = st.logZ;
        figNspk(j) = st.nSpikes;      % POOLED nnz(spike_train>0) across recordings
        figRate(j) = st.rateHz;       % POOLED nSpikes / durTot
    catch ME
        fprintf(2,'  stats failed cell %d (%s): %s\n', CJ(j).cell, CJ(j).key, ME.message);
    end
    if mod(j,25)==0, fprintf('  measured %d/%d\n', j, numel(CJ)); end
end
fprintf('measured %d cells in %.1f min\n', numel(CJ), toc(tM)/60);

switch lower(selMode)
    case 'active'
        keep = figRate >= activeMinRateHz;
        fprintf('ACTIVE (pooled rate >= %.3g ev/min) : %d of %d cells kept, %d dropped\n', ...
                activeMinRateHz*60, nnz(keep), numel(CJ), nnz(~keep));
        for q = find(~keep)'
            fprintf('   dropped cell %3d  %-10s %5.2f ev/min (%d ev)\n', ...
                    CJ(q).cell, CJ(q).group, figRate(q)*60, figNspk(q));
        end
        fprintf('\n');
    otherwise
        keep = figLogZ >= logZ_min;
        fprintf('pooled logZ >= %g : %d of %d cells kept, %d dropped\n\n', ...
                logZ_min, nnz(keep), numel(CJ), nnz(~keep));
end

%% ===================== RENDER =====================
status = repmat("", numel(CJ), 1);
outPng = repmat("", numel(CJ), 1);
idx = find(keep)';
tAll = tic;
for n = 1:numel(idx)
    j = idx(n);  J = CJ(j);
    stem = regexprep(sprintf('%s_%s_cell%03d', J.group, J.obs(1).recDate, J.cell), ...
                     '[\\/:*?"<>|]', '_');
    base = fullfile(outDir, stem);
    outPng(j) = string([base '.png']);
    if ~overwrite && isfile([base '.png'])
        status(j) = "skipped (exists)";  continue;
    end
    tJ = tic;  figR = [];
    try
        Pr = P; Pr.cellId = J.cell; Pr.traceObs = 1; Pr.trace_xlim_sp = traceWin;
        [figR, projS] = temporal_phase_cell_fig_260812(J.obs, Pr);
        exportgraphics(figR, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
        if savePDF
            exportgraphics(figR, [base '.pdf'], 'ContentType','vector','BackgroundColor','white');
        end
        if saveProj && projS.have_proj
            save_avgproj_png_260812(projS, P, fullfile(outDir, [stem '_avgproj.png']));
        end
        status(j) = "ok";
        fprintf('[%3d/%3d] %5.1fs  %s  (%d rec, %d ev, logZ %.2f)\n', ...
                n, numel(idx), toc(tJ), stem, nRecs(j), figNspk(j), figLogZ(j));
    catch ME
        status(j) = "ERROR: " + string(ME.message);
        fprintf(2,'[%3d/%3d] FAILED %s\n           %s\n', n, numel(idx), stem, ME.message);
    end
    if ~isempty(figR) && isgraphics(figR), close(figR); end
end
status(~keep) = "not rendered (below the " + string(lower(selMode)) + " cut)";

%% ===================== INDEX =====================
T = table([CJ.cell]', string({CJ.key})', string({CJ.group})', nRecs, ...
          figNspk, figLogZ, [CJ.pooled]', keep, status, outPng, ...
    'VariableNames', {'cell_idx','cell_key','group','n_recordings', ...
                      'pooled_n_events','figure_pooled_logZ','archive_ray_logZ', ...
                      'rendered','status','png'});
% A restricted run describes only the cells it touched, so it must NOT overwrite
% the full index -- doing so silently replaced a 269-row record with 3 rows once.
if isempty(onlyCells)
    idxCsv = fullfile(outDir, 'per_cell_summary_index.csv');
else
    idxCsv = fullfile(outDir, 'per_cell_summary_index_subset.csv');
end
writetable(T, idxCsv);
fprintf('  index written    : %s\n', idxCsv);

nOK = nnz(status == "ok");  nBad = nnz(startsWith(status,"ERROR"));
fprintf('\n---- done in %.1f min ----\n', toc(tAll)/60);
fprintf('  rendered %d | dropped %d | failed %d\n', nOK, nnz(~keep), nBad);
fprintf('  %s\n', outDir);
