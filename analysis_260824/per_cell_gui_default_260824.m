% per_cell_gui_default_260824.m
% -----------------------------------------------------------------------
%  One per-CELL temporal-phase summary figure for EVERY cell of the
%  260824 vagotomised Vglut2 session, rendered in the GUI's DEFAULT MODE.
%
%  "GUI default mode" is meant literally: the P block below is a verbatim copy
%  of the USER-EDITABLE block of
%      analysis_260727\coh_ca_breath\temporal_phase_cell_gui_260812.m
%  and the render call is the same call its Render button makes, with every
%  GUI control left at the value it launches with:
%      mode      pooled cell   (all of the cell's recordings)
%      rec       1             (first recording supplies trace + avg-proj)
%      window    0 30          (first 30 s of that recording)
%      trig win  blank         (auto = 2 x median IBI)
%      ylims     blank         (auto)
%      crop um   70   clip %   1 99.5   bar  lower left
%  So a figure here is byte-for-byte the figure the GUI would produce if you
%  loaded that cell and clicked Render without touching anything -- the same
%  look as D:\Ventral_surface_summary\gui_renders_260815.
%
%  NO PERMUTATION PANELS. That is not a switch in this file: P.guiLayout = true
%  is the GUI default and it drops the spike-triggered average and BOTH spike
%  histograms, which are the panels the circular-shift permutation feeds.
%  P.showRayleigh = false likewise removes the log Z line. The statistics are
%  still computed and returned in `stats` -- they are simply not drawn, exactly
%  as in the GUI.
%
%  CELL IDENTITY comes from this session's own curation,
%      <rootPath>\roi_match_out_260824\roi_match_curated.csv
%  which is the numbering cell_link_260727 would produce (its doc says the ids
%  are "identical to the numbering in roi_match_curated.csv"). Reading the CSV
%  directly means this script needs neither cell_link nor cell_pool, and cannot
%  pick up the archive registry, whose cell ids belong to other sessions.
%  Tossed ROIs (grpOf < 0) are already absent from the CSV.
%
%  ROTATED MOUNT. 2026-08-24 is after the 2026-07-21 cutover, so the session was
%  added to the ROT list inside temporal_phase_cell_fig_260812.m. Without that
%  entry the figure's annotation would print raw stage x/y with lateral and
%  rostral swapped.
%
%  OUTPUT  <rootPath>\gui_renders_260824\
%      <recName>_cell%03d_pooled%d.png / .pdf / _avgproj.png
%  the same stem the GUI writes, so a hand-tuned re-render from the GUI lands on
%  the same filename and replaces the default one.
%
%  Runqi Zhang / 2026-08-24
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir); addpath(repoRoot);
addpath(fullfile(repoRoot, 'analysis_260727', 'coh_ca_breath'));
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
addpath(fullfile(repoRoot, 'falloff-analysis-260805'));   % laser_power_calibration
addpath(fullfile(repoRoot, 'analysis_260806'));

%% ===================== USER-EDITABLE =====================
rootPath   = 'C:\260824_Vglut2-soma-g8s_vagotomized\phys';
curatedCsv = fullfile(rootPath, 'roi_match_out_260824', 'roi_match_curated.csv');
genotype   = 'Vglut2';
recDate    = '0824';        % MMDD; keys the rotated-mount lookup in the figure
outDir     = fullfile(rootPath, 'gui_renders_260824');

onlyCells  = [];            % [] = every cell; else a list of cell ids
overwrite  = false;          % false = skip cells whose .png already exists
savePDF    = true;
saveProj   = true;
pngDpi     = 600;           % the GUI's setting

% GUI defaults for the controls (see the header). traceObs = 1 is the "rec"
% dropdown's launch position; the pooled statistics use every recording either
% way, so this only chooses the wide trace and the avg-projection crop.
traceWin   = [0 30];        % the GUI's window box, in seconds. [] = middle 30 s
traceObs   = 1;
barCorner  = 'lower left';

%% ===== FIGURE PARAMETERS -- verbatim from temporal_phase_cell_gui_260812 =====
P = struct();
P.doCoh        = false;
P.nDrop        = 30;          P.fallback_fps = 30;
P.TW_spec      = 6;           P.alpha_sig    = 0.01;
P.minSpikes    = 2;           P.ca_lag_sec   = 0;    % 0 = NO ca-lag compensation
P.f_breath_search = [0.2 4];  P.fwhm_factor  = 0.6;   P.min_bw = 0.05;
P.fmin         = 0.05;        P.fmax         = 15;
P.trigWin_sec  = [];          % [] = auto = P.trigWinIBI x median IBI
P.trigWinIBI   = 2;
P.ylim_dff     = [];          P.ylim_epc     = [];
P.histBinFrames= 2;
P.nShuffle     = 1200;
P.shiftMinCyc  = 3;
P.pad_um       = 20;
P.clip_pct     = [1 99.5];
P.crop_um      = 70;
P.cropPxUm     = 0.5;
P.scalebar_um  = 10;
P.showRayleigh = false;
P.showAcqTime  = true;
P.trialOverlay = true;
P.dffSpecDeriv = false;
P.guiLayout    = true;
P.boldGenoOnly = true;
P.spkWinIBI    = 1;
P.scalebarThickFrac = 0.04;
P.outlineFrac       = 0.008;
P.outlineAlpha      = 0.6;
P.traceBreathCol    = [0 0 0];
P.outlineImgPx = 1.2;
P.rayPhaseBins    = 36;
P.nShuffleRayleigh= 0;
P.gamma_val    = 1;           P.PixelSizeBase = 1.7778;  P.outlineLW = 0.8;
P.sortMode     = 'none';      P.trace_xlim_sp = [];
P.dffColor     = [0.2 0.7 0.2];
P.onsetRowStyle = 'raster';
P.onsetRowTickFrac = 0.5;
P.onsetCol     = [1.00 0.40 0.75];
P.peakCol      = [0.35 0.55 1.00];
P.breathObsByCell = [];
P.TW_coh       = 4;           P.alpha_coh    = 0.001;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

% Software OpenGL. The NVIDIA WGL path throws inside exportgraphics ~1% of the
% time on this machine, and in a batch that silently writes a 9 KB blank over a
% good figure. See reference_matlab_opengl_export_stubs.
try, opengl('software'); catch, end

if ~isfolder(outDir), mkdir(outDir); end

%% ===================== BUILD THE CELL LIST =====================
assert(isfile(curatedCsv), 'curation not found: %s', curatedCsv);
T = readtable(curatedCsv);
req = {'cell_id','fov_name','roi_index'};
for k = 1:numel(req)
    assert(ismember(req{k}, T.Properties.VariableNames), ...
        '%s has no column "%s"', curatedCsv, req{k});
end
cellIds = unique(T.cell_id(:), 'stable');
cellIds = sort(cellIds);
if ~isempty(onlyCells), cellIds = cellIds(ismember(cellIds, onlyCells(:))); end

fprintf('\n=========== per_cell_gui_default_260824 ===========\n');
fprintf('root      : %s\n', rootPath);
fprintf('curation  : %s\n', curatedCsv);
fprintf('cells     : %d   (from %d ROI observations)\n', numel(cellIds), height(T));
fprintf('out       : %s\n\n', outDir);

%% ---- audit the join BEFORE rendering anything ----
% Every (folder, roi) must exist and address a real column of that recording's
% dFF. A silent index error here would draw the wrong neuron, which no downstream
% inspection would catch.
bad = {};
fovList = unique(T.fov_name);
for k = 1:numel(fovList)
    fp = fullfile(rootPath, fovList{k});
    if ~isfolder(fp), bad{end+1} = sprintf('missing folder: %s', fovList{k}); continue; end %#ok<SAGROW>
    dd = dir(fullfile(fp, '*_ch1_dFF.mat'));
    if isempty(dd), bad{end+1} = sprintf('no *_ch1_dFF.mat in %s', fovList{k}); continue; end %#ok<SAGROW>
    if ~isfile(fullfile(fp, 'breath_peak_pc1.mat'))
        bad{end+1} = sprintf('no breath_peak_pc1.mat in %s', fovList{k}); continue; %#ok<SAGROW>
    end
    S = load(fullfile(dd(1).folder, dd(1).name), 'dFF');
    nRoi = size(S.dFF, 2);
    r = T.roi_index(strcmp(T.fov_name, fovList{k}));
    if any(r < 1 | r > nRoi)
        bad{end+1} = sprintf('%s: roi index out of range (dFF has %d ROIs, saw %s)', ...
            fovList{k}, nRoi, mat2str(unique(r(r<1 | r>nRoi))')); %#ok<SAGROW>
    end
end
if ~isempty(bad)
    fprintf(2, 'JOIN AUDIT FAILED:\n');  fprintf(2, '  %s\n', bad{:});
    error('per_cell_gui_default_260824:join', 'fix the above before rendering');
end
fprintf('join audit: OK -- %d FOVs, every ROI index inside its dFF\n\n', numel(fovList));

%% ===================== RENDER =====================
nOK = 0; nSkip = 0; failed = {};
tAll = tic;
for c = 1:numel(cellIds)
    cid = cellIds(c);
    rows = find(T.cell_id == cid);

    OBS = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
    for q = 1:numel(rows)
        rn = T.fov_name{rows(q)};
        OBS(end+1) = struct('folder',  fullfile(rootPath, rn), ...
                            'roi',     T.roi_index(rows(q)), ...
                            'recName', rn, ...
                            'group',   genotype, ...
                            'recDate', recDate); %#ok<SAGROW>
    end

    iT   = min(traceObs, numel(OBS));
    stem = regexprep(sprintf('%s_cell%03d_pooled%d', OBS(iT).recName, cid, numel(OBS)), ...
                     '[\\/:*?"<>|]', '_');
    base = fullfile(outDir, stem);
    if ~overwrite && isfile([base '.png'])
        nSkip = nSkip + 1; continue;
    end

    fprintf('[%3d/%3d] cell %03d  (%d rec)  %s ... ', c, numel(cellIds), cid, numel(OBS), OBS(iT).recName);
    Pr = P;
    Pr.cellId          = cid;
    Pr.traceObs        = iT;
    Pr.trace_xlim_sp   = traceWin;
    Pr.scalebarCorner  = barCorner;

    try
        [figR, proj] = temporal_phase_cell_fig_260812(OBS, Pr);
    catch ME
        fprintf(2, 'RENDER FAILED: %s\n', ME.message);
        failed{end+1} = sprintf('cell %03d: %s', cid, ME.message); %#ok<SAGROW>
        continue;
    end

    try
        exportgraphics(figR, [base '.png'], 'Resolution', pngDpi, 'BackgroundColor','white');
        if savePDF
            exportgraphics(figR, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
        end
        if saveProj && proj.have_proj
            save_avgproj_png_260812(proj, Pr, [base '_avgproj.png']);
        end
    catch ME
        fprintf(2, 'EXPORT FAILED: %s\n', ME.message);
        failed{end+1} = sprintf('cell %03d export: %s', cid, ME.message); %#ok<SAGROW>
        close(figR); continue;
    end
    close(figR);

    % The OpenGL stub writes a ~9 KB blank instead of erroring, so the size is
    % the only evidence the export actually happened. Retry once.
    dpng = dir([base '.png']);
    if isempty(dpng) || dpng.bytes < 50e3
        fprintf(2, 'suspect PNG (%d bytes) -- re-rendering once ... ', ...
                ternary_bytes(dpng));
        try
            [figR, proj] = temporal_phase_cell_fig_260812(OBS, Pr);
            exportgraphics(figR, [base '.png'], 'Resolution', pngDpi, 'BackgroundColor','white');
            close(figR);
            dpng = dir([base '.png']);
        catch ME
            failed{end+1} = sprintf('cell %03d retry: %s', cid, ME.message); %#ok<SAGROW>
        end
        if isempty(dpng) || dpng.bytes < 50e3
            failed{end+1} = sprintf('cell %03d: PNG still %s bytes', cid, ternary_bytes(dpng)); %#ok<SAGROW>
            fprintf(2, 'STILL BAD\n'); continue;
        end
    end

    nOK = nOK + 1;
    fprintf('ok (%.1f MB)\n', dpng.bytes/1e6);
end

fprintf('\n--------------------------------------------------\n');
fprintf('rendered %d / %d cells  (skipped %d, failed %d)  in %.1f min\n', ...
        nOK, numel(cellIds), nSkip, numel(failed), toc(tAll)/60);
if ~isempty(failed)
    fprintf(2, 'failures:\n');  fprintf(2, '  %s\n', failed{:});
end
fprintf('out: %s\n', outDir);

% -----------------------------------------------------------------------
function s = ternary_bytes(d)
if isempty(d), s = '0'; else, s = sprintf('%d', d.bytes); end
end
