% cell_cycle_gui_260831.m
% -----------------------------------------------------------------------
%  IDENTICAL to temporal_phase_cell_gui_260812.m -- same controls, same
%  inputs, same pooled/single-ROI modes, same save behaviour -- with ONE line
%  changed: Render calls cell_cycle_fig_260831 instead of
%  temporal_phase_cell_fig_260812, so you get the PER-CYCLE figure (spaghetti
%  panels, square aligned boxes) rather than the full summary.
%
%  Kept as a copy rather than a flag on the original, so the summary GUI that
%  every existing figure came from is untouched.
%
%  Copied 2026-08-31. The original header follows verbatim.
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
%  Per-CELL temporal-phase summary, driven entirely from the GUI.
%
%  Point it at a recording FOLDER and give a cell (or ROI) number -- that is
%  the whole setup. Nothing has to be edited in this file, and there is NO
%  dependence on coherence_polar_data.mat / coh_cfg any more.
%
%      [ folder ][ Browse ]   [ id ] (o)cell ( )roi   [ Load ]
%      [ pooled cell | single ROI ]
%      [ trig win ][ dFF ylim ][ epc ylim ]
%      [ window ][Apply][Full][ Render figure ][<<Prev][Next>>]
%
%  SAME FIGURE AS THE BATCH.  Rendering calls temporal_phase_cell_fig_260812.m,
%  the very function analysis_260806\per_cell_summary_260812.m calls, so a GUI
%  figure and a batch figure of the same cell are identical apart from the window
%  and y-limits you set here. This file no longer contains any plotting code.
%
%  TWO MODES
%    pooled cell : every recording of the cell goes into the histograms, the
%                  triggered averages, the heatmap, the PSD and the log Z --
%                  exactly what the batch produces. Prev/Next then chooses which
%                  recording supplies the wide trace and the avg-projection.
%    single ROI  : only the loaded (recording, ROI). Pooling one observation is a
%                  no-op, so this is the per-ROI view through the same code path,
%                  not a second implementation that could drift.
%
%  TRIGGER WINDOW.  The three dF/F averages and the per-cycle heatmap share the
%  "trig win" box, which is the TOTAL width in seconds. Blank means auto:
%  P.trigWinIBI (=2) x the median inter-breath interval, IBI being the median gap
%  between inspiration onsets.
%
%  The two SPIKE HISTOGRAMS deliberately do NOT follow that box. Their
%  permutation test is always run on a 1-IBI window -- the width at which each
%  breath contributes exactly once, so "spikes per cycle" means what it says --
%  while the plot is doubled to +/-1 IBI so neighbouring cycles are visible and
%  rhythmicity can be read off. Those outer halves are real measured data, not a
%  copy of the centre, so a spike near the edge is counted in two trigger windows
%  and the bars sum to ~200% over the full 2 IBI. Dotted lines mark the +/-0.5 IBI
%  region the p-value comes from.
%
%  "dFF ylim" and "epc ylim" take "lo hi" and pin the y-axes of the dF/F averages
%  / the spikes-per-cycle histograms so cells can be compared by eye; blank
%  leaves that axis auto-scaled.
%
%  PER CELL, WHEN AVAILABLE.  If a cell_link.mat can be found for the session
%  (walking up from the FOV to <session>\phys\analysis_260727\cell_pooled\),
%  the id is a CELL id and Prev/Next step through that cell's observations --
%  i.e. the same cell as it appears in each FOV, one FOV at a time. With no
%  cell_link the id falls back to a plain ROI index in the loaded folder and
%  Prev/Next step through ROIs.
%
%  NOTE cell_link.mat's own `fov_name` field is truncated at the first dot
%  ("roi1_3x_12" for "roi1_3x_12.5lp_..."), a fileparts-on-dotted-name bug in
%  cell_link_260727.m. `obsT.rec_name` / `obsT.rec_path` are intact, so this
%  script keys off those and never reads link.fov_name.
%
%  NO COHERENCE.  The coherence spectrum and coherence polar panels are gone,
%  along with the coherencyc calls that fed them (P.doCoh = false). The dF/F
%  and breath POWER SPECTRA panel is kept -- that is a spectrum, not coherence.
%
%  NO PHASE PANELS.  The linear breath-phase histogram and the event-phase polar
%  are gone too. Breath locking is shown in SECONDS instead, by two spike
%  histograms side by side: one triggered on inspiration ONSET, one on the
%  breath PEAK, each with its own circular-shift null band. The per-cycle
%  heatmap is likewise PEAK-triggered, with the cycle-average breath waveform
%  drawn across the top of it in place of the old event raster.
%
%  NOTE the dF/F power spectrum is the spectrum of diff(dff)*fps, i.e. of the
%  DERIVATIVE -- that flattens the GCaMP decay but tilts the curve +20 dB/decade,
%  so its peak sits above the true dF/F peak. Display only; no number uses it.
%
%  Window: type "start end" in seconds, Apply to preview, Render to build the
%  full figure. Nothing renders until Render is clicked.
%
%  Requires in each FOV folder: *_ch1_dFF.mat and breath_peak_pc1.mat
%  (breath_insp_start_pc1.mat optional but needed for onset-locked panels).
%
%  Dependencies: Chronux (mtspectrumc), detect_session_fps.m,
%                Image Processing TB (bwboundaries/imdilate).
%
%  Runqi Zhang / 2026-08-12
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
% ONE fileparts, not two: the original lives in analysis_260727\coh_ca_breath and
% needed two to reach the repo root; this copy sits one level down in
% analysis_260831, so two would resolve to the Desktop and every shared
% dependency would silently fail to resolve.
repoRoot  = fileparts(scriptDir); addpath(repoRoot);
addpath(scriptDir);
% the figure function this GUI renders through, plus the canonical summary
% function it calls internally for the crop and the stats
addpath(fullfile(repoRoot, 'analysis_260727', 'coh_ca_breath'));
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));
addpath(fullfile(repoRoot, 'falloff-analysis-260805'));   % laser_power_calibration
addpath(fullfile(repoRoot, 'analysis_260806'));           % registry + cell merges
% the figure itself lives in temporal_phase_cell_fig_260812.m, next to this file,
% and is shared with analysis_260806\per_cell_summary_260812.m

%% ===================== USER-EDITABLE =====================
% Only used to pre-fill the GUI boxes -- everything is changeable in the GUI.
  startFolder = 'D:\Ventral_surface_summary\Vglut2\0824\cell2\roi1_z-10_3x_3000f_16lp_00001';
  startId     = 291;        % pooled id, 282..297
  idMode      = 'cell';
% cell 48 is one of the 260814 merges (was 48 + 52), so the default launch shows
% a genuinely pooled cell with two recordings in the dropdown.

% Cell identity source, in priority order:
%   1. this registry + cell_merge_overrides_260814  -- what the BATCH uses, so a
%      cell id here is the same number as in per-cell-summary_..._cell%03d.png
%   2. the session's cell_link.mat                  -- older per-session grouping
%   3. plain ROI mode
% Priority matters: cell_link was never run for Vglut2/0224 or Vglut2/1124, and
% for Vglut2/0810 and Sst/0807 it ran but missed a pair, so cell_link alone shows
% cells the batch treats as merged. Set regFile = '' to force the old behaviour.
regFile     = 'D:\Ventral_surface_summary\event_latency_260811\event_latency_data.mat';

% One collected output folder for everything the GUI renders. '' would instead
% write into each recording's own <folder>\temporal_phase_cell, which scatters
% renders through the data archive and puts the SAME pooled cell in a different
% place depending on which recording the dropdown happens to be on. The filename
% already carries the recording name and the cell id, so nothing collides here.
outDir      = 'D:\Ventral_surface_summary\figure_260831';
doSave      = true;
pngDpi      = 600;          % PNG export resolution. The PDF is vector and is not
                            % affected. 600 makes the 1.1 px ROI outline and the
                            % single-trial overlays resolve properly; the file
                            % grows roughly (600/200)^2 = 9x.

P = struct();
P.doCoh        = false;       % coherence panels removed; true restores them
P.nDrop        = 30;          P.fallback_fps = 30;
P.TW_spec      = 6;           P.alpha_sig    = 0.01;
P.minSpikes    = 2;           P.ca_lag_sec   = 0;    % 0 = NO ca-lag compensation
P.f_breath_search = [0.2 4];  P.fwhm_factor  = 0.6;   P.min_bw = 0.05;
P.fmin         = 0.05;        P.fmax         = 15;
% ONE window for every triggered panel: averages, heatmap and both histograms.
P.trigWin_sec  = [];          % TOTAL width (s) of every triggered panel; [] = auto
P.trigWinIBI   = 2;           % auto width, in median inter-breath intervals
P.ylim_dff     = [];          % [lo hi] for the dF/F averages;   [] = auto
P.ylim_epc     = [];          % [lo hi] for the spk/cyc % hists; [] = auto
P.histBinFrames= 2;
P.nShuffle     = 1200;        % drives BOTH the null band and the permutation p
P.shiftMinCyc  = 3;
P.pad_um       = 20;
% Contrast stretch percentiles for the avg-proj crop, [lo hi]. GUI uses [1 99.5];
% the batch keeps [0.5 99.9]. There is NO gamma on top -- gamma_val = 1 -- so
% these two numbers are the whole tone mapping. Overridable in the GUI box.
P.clip_pct     = [1 99.5];
% GUI uses crop_um: the TOTAL crop width in microns, centred on the ROI. pad_um
% above is the legacy rule (a margin added to the ROI radius, giving a crop of
% 2*(radius+pad)) and is ignored whenever crop_um is set. The batch still uses
% pad_um = 20. The GUI box overrides crop_um per render.
% 70 um total, which is about what the legacy pad_um = 20 rule produced (crops
% came out 47-247 px at px_um 0.4-1.8, i.e. ~80 um of field). NOTE 15 um would be
% a 13x13 px crop at these zooms -- the number here is a TOTAL WIDTH, not a margin.
P.crop_um      = 70;
% Every crop is resampled to this pixel size (bicubic), so a 70 um crop is always
% 140x140 px whatever zoom it was acquired at. That is what lets the scale bar and
% the ROI outline be specified once, in image pixels, and come out the same width
% on every cell -- native px_um spans ~0.4-1.8 across this archive. Batch keeps
% native pixels ([] = no resample).
P.cropPxUm     = 0.5;
P.scalebar_um  = 10;          % white bar on the avg-proj crop; drawn UNLABELLED.
                              % GUI only -- the batch keeps 50. A 50 um bar would
                              % span most of a 70 um crop.
% Annotation stats, ventral convention (Ventral_surface_polar_coh_vs_rayleigh_260808)
% ---- GUI-ONLY DISPLAY FLAGS -------------------------------------------------
% Every one of these defaults the OTHER way inside temporal_phase_cell_fig_260812,
% so the batch and the 265 figures already on disk are untouched. They exist only
% so the interactive view can differ from the archived figures on purpose.
P.showRayleigh = false;   % no Rayleigh log Z line. logZ is still computed and
                          % returned in stats -- cell selection needs it -- it is
                          % just not a claim this view makes.
P.showAcqTime  = true;    % annotate "Vglut2 260810 15:02": genotype, session and
                          % acquisition clock time, read from the ScanImage epoch
                          % in the RAW tif. Falls back to MMDD, then to genotype
                          % alone, when the acquisition drive is not mounted.
P.trialOverlay = true;    % triggered averages show every single trial faintly
                          % with the mean on top, instead of a +/- SD ribbon.
P.dffSpecDeriv = false;   % plain dF/F spectrum, NOT the derivative.
P.guiLayout    = true;    % two-column layout: drops the spike-triggered average
                          % and BOTH spike histograms; column 2 is the two
                          % triggered averages, onset above peak. Each carries a
                          % single trigger line in its own colour (red onset /
                          % blue peak) -- never both on one panel. The breath
                          % strip shows only the first onset. Everything dropped
                          % is still computed and returned in stats.
P.boldGenoOnly = true;    % the genotype is the ONLY bold text on the figure;
                          % all panel titles, labels and the suptitle go normal.
P.spkWinIBI    = 1;       % spike-triggered panel spans 1 IBI total (+/-0.5 IBI),
                          % narrower than the shared trigger window. It is the one
                          % panel not locked to the breath cycle, so the wider span
                          % mostly showed neighbouring cycles. [] = share the
                          % common window, which is what the batch does.
% Bar and outline as FRACTIONS of the crop, so they hold their proportions if
% crop_um or cropPxUm change. With the crop pinned at 140x140 px these are 5.6 px
% and 1.1 px, identical on every cell.
P.scalebarThickFrac = 0.04;    % bar thickness = 4%   of crop height
P.outlineFrac       = 0.008;   % ROI outline   = 0.8% of crop width
P.outlineAlpha      = 0.6;     % semi-transparent outline, so the cell shows through
P.traceBreathCol    = [0 0 0]; % breath in BLACK in the wide top trace. It is still
                               % plotted BEFORE the dF/F, so it stays underneath.
P.outlineImgPx = 1.2;     % fallback if outlineFrac is cleared. Thickness in IMAGE pixels rather than
                          % points, so it looks the same on every crop size AND
                          % matches between the figure panel and the standalone
                          % avgproj PNG (those measured 3.0 px vs 8.0 px for the
                          % same crop at the fixed-points setting). [] = old way.
% -----------------------------------------------------------------------------
P.rayPhaseBins    = 36;       % 10 deg occupancy bins for the event weights
P.nShuffleRayleigh= 0;        % 0 = do NOT run a second, Rayleigh-specific null.
                              % The only p in the figure is the PSTH permutation
                              % p (P.nShuffle), shown in both places.
% gamma 1 = LINEAR, no tone curve. The crops are small and their contrast stretch
% already comes from the percentiles of that small window; gamma on top distorts
% more than it reveals.
P.gamma_val    = 1;           P.PixelSizeBase = 1.7778;  P.outlineLW = 0.8;
P.sortMode     = 'none';        P.trace_xlim_sp = [];
P.dffColor     = [0.2 0.7 0.2];
P.onsetRowStyle = 'raster';  % per-cycle onset on the heatmap: 'raster' = one
                             % short tick per breath, nothing joining them.
                             % 'line' is the old joined line and is still what
                             % the batch renders.
P.onsetRowTickFrac = 0.5;    % raster tick height, in heatmap rows. 1 = ticks
                             % touch (reads solid again), 0.5 = half a row.
P.onsetCol     = [1.00 0.40 0.75];   % onset marks are PINK (2026-08-17, RZ):
                                     % the same colour as the per-cycle onset
                                     % line on the heatmap, so every onset mark
                                     % on the figure reads as one thing.
P.peakCol      = [0.35 0.55 1.00];   % cornflower, paired to the pink onset:
                                     % same lightness, hue ~219 vs the pink's 325,
                                     % and violet-leaning rather than cyan so the
                                     % two read as one colour pair.
% Per-cell breath-strip restriction: [cellId, whichRecording; ...]. The named
% cell's strip is averaged from that recording ALONE, while the cell stays pooled
% for the dF/F averages, the heatmap, the histograms and the Rayleigh -- only the
% strip narrows. Recording index is the cell's own observation order, the same
% order the "rec" menu lists, so 1 is the first entry in that menu.
%
% 178: two recordings at different breathing rates. Averaging one animal's cycle
% against the other's washed the strip out either side of t=0. Nothing else in
% the figure has that problem -- every other panel pools by trigger and by
% occupancy weight, which is rate-invariant. RZ, 2026-08-17.
%
% The strip title prints "(rec 1 of 2)" whenever this is in force, so a
% restricted strip can never be mistaken for a pooled one.
P.breathObsByCell = [];   % 178 was tried at [178 1] and reverted 2026-08-17:
                          % all cycles, both recordings. The mechanism stays --
                          % add a [cellId whichRec] row to restrict any cell.
% kept only so the render code can still be switched back on
P.TW_coh       = 4;           P.alpha_coh    = 0.001;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

% Software OpenGL: the hardware (NVIDIA WGL) path throws
%   com.jogamp.opengl.GLException: Error making context ... current
% inside exportgraphics, ~1% of the time on this machine. In the batch that cost
% three figures; here it would silently save a 9 KB blank over a good one.
try
    opengl('software');
catch
    % older MATLAB, or already software -- nothing to do
end

%% ===================== LAUNCH GUI =====================
launch_cell_gui(startFolder, startId, idMode, P, doSave, outDir, regFile, pngDpi);

% =======================================================================
% =========================== LOCAL FUNCTIONS ===========================
% =======================================================================

function REG = get_registry(regFile)
%GET_REGISTRY  Flatten event_latency_data.mat into one row per (cell, recording).
%  This is the SAME cell identity the batch renders with: the registry's cells,
%  with cell_merge_overrides_260814 folded in, and a merged cell taking the lowest
%  of its constituent numbers. So a cell id typed into the GUI is the same number
%  as in per-cell-summary_..._cell%03d.png.
%  Returns [] if the registry is unavailable, which drops the GUI back to
%  cell_link and then to ROI mode.
    persistent CACHE CACHEKEY
    REG = [];
    if isempty(regFile) || ~isfile(regFile), return; end
    % Cache key includes the file's SIZE and MODIFIED TIME, not just its path.
    % Keying on the path alone meant that once the registry was read, a rebuilt
    % or appended registry was never picked up for the life of the MATLAB
    % session -- newly registered cells simply did not exist as far as the GUI
    % was concerned, with no error to say so.
    fi  = dir(regFile);
    key = sprintf('%s|%d|%.6f', regFile, fi.bytes, fi.datenum);
    if ~isempty(CACHE) && strcmp(CACHEKEY, key), REG = CACHE; return; end
    try
        D = load(regFile, 'CELL','OBS','REC');
        obsOf = pooled_obs_260814(D.CELL, D.OBS);
    catch ME
        warning('cellgui:noRegistry', ...
                'registry unavailable (%s) -- falling back to cell_link', ME.message);
        return;
    end
    rows = struct('cellNum',{},'fp',{},'recName',{},'roi',{},'group',{},'recDate',{});
    for c = 1:numel(obsOf)
        if isempty(obsOf{c}), continue; end       % merged away into a lower number
        for o = obsOf{c}(:)'
            p = regexp(D.OBS(o).label,'/','split');
            if numel(p) < 3, continue; end
            rows(end+1) = struct('cellNum',c, ...
                'fp',D.REC(D.OBS(o).rec).folder, ...
                'recName',strjoin(p(3:end-1),'/'), 'roi',str2double(p{end}), ...
                'group',p{1}, 'recDate',p{2}); %#ok<AGROW>
        end
    end
    REG = rows;  CACHE = REG;  CACHEKEY = regFile;
end

% -----------------------------------------------------------------------
function k = reg_rows_for_folder(REG, folderPath)
%REG_ROWS_FOR_FOLDER  Registry rows for a FOV folder, by path or by folder name.
%  The registry stores ARCHIVE paths (D:\Ventral_surface_summary\...), but the
%  same recording is often opened from where it was acquired (C:\260810_...\phys\).
%  Exact path first; if that misses, fall back to the leaf folder name, and only
%  when that name is unambiguous across the whole registry -- recording names do
%  repeat between sessions, and quietly matching the wrong session would pool two
%  different animals into one "cell".
    k = [];
    if isempty(REG), return; end
    k = find(strcmpi({REG.fp}, folderPath));
    if ~isempty(k), return; end

    leaf = folderPath;
    leaf = regexprep(leaf, '[\\/]+$', '');
    parts = strsplit(leaf, {'\','/'});
    leaf  = parts{end};
    if isempty(leaf), return; end

    regLeaf = cell(1,numel(REG));
    for i = 1:numel(REG)
        q = strsplit(regexprep(REG(i).fp,'[\\/]+$',''), {'\','/'});
        regLeaf{i} = q{end};
    end
    hit = find(strcmpi(regLeaf, leaf));
    if isempty(hit), return; end
    if numel(unique({REG(hit).fp})) > 1
        warning('cellgui:ambiguousFolder', ...
            ['folder name "%s" appears in %d registry locations -- not matching ' ...
             'by name. Open the archive copy to use registry cell ids.'], ...
            leaf, numel(unique({REG(hit).fp})));
        return;
    end
    k = hit;
end

% -----------------------------------------------------------------------
function items = items_from_registry(REG, folderPath, idVal)
%ITEMS_FROM_REGISTRY  Every recording of registry cell idVal, loaded folder first.
    items = struct('folderPath',{},'roi',{},'group',{},'recDate',{}, ...
                   'recName',{},'cellId',{},'label',{});
    if isempty(REG), return; end
    for k = find([REG.cellNum] == idVal)
        items(end+1) = struct('folderPath',REG(k).fp, 'roi',REG(k).roi, ...
            'group',REG(k).group, 'recDate',REG(k).recDate, ...
            'recName',REG(k).recName, 'cellId',idVal, ...
            'label',sprintf('cell %d | %s ROI%d', idVal, REG(k).recName, REG(k).roi)); %#ok<AGROW>
    end
    if isempty(items), return; end
    same = strcmpi({items.folderPath}, folderPath);
    if any(same), items = [items(same), items(~same)]; end
end

% -----------------------------------------------------------------------
function lk = find_cell_link(folderPath)
%FIND_CELL_LINK  Walk up from a FOV folder looking for the session cell_link.
%  Returns [] when there is none, which is the signal to fall back to ROI mode.
    lk = [];
    d = folderPath;
    for up = 1:4
        d = fileparts(d);
        if isempty(d); return; end
        cand = { fullfile(d,'analysis_260727','cell_pooled','cell_link.mat'), ...
                 fullfile(d,'cell_pooled','cell_link.mat'), ...
                 fullfile(d,'cell_link.mat') };
        for c = cand
            if isfile(c{1})
                S = load(c{1},'link');
                if isfield(S,'link') && isfield(S.link,'obsT')
                    lk = S.link; lk.srcFile = c{1};
                    return;
                end
            end
        end
    end
end

% -----------------------------------------------------------------------
function items = items_for_id(REG, lk, folderPath, idVal, idMode)
%ITEMS_FOR_ID  Build the Prev/Next list.
%  cell mode : every observation of cell idVal, one entry per FOV, ordered so
%              the currently loaded folder comes first. The REGISTRY is tried
%              first so the GUI pools exactly what the batch pools, including the
%              260814 merges; cell_link is only the fallback for recordings the
%              registry does not cover.
%  roi  mode : the single ROI idVal in folderPath.
    items = struct('folderPath',{},'roi',{},'group',{},'recDate',{}, ...
                   'recName',{},'cellId',{},'label',{});

    if strcmpi(idMode,'cell')
        items = items_from_registry(REG, folderPath, idVal);
        if ~isempty(items), return; end
        % The registry is AUTHORITATIVE when it exists. Do not fall through to
        % cell_link on a miss: the two number cells completely differently (the
        % registry is global, cell_link restarts per session), so a fallback
        % silently returns a REAL BUT DIFFERENT cell with no error -- which is
        % exactly what happened to newly registered ids that a stale cache had
        % not picked up yet.
        if ~isempty(REG)
            warning('cellgui:idNotInRegistry', ...
                ['cell %d is not in the registry. NOT falling back to cell_link, ' ...
                 'whose ids mean something different.'], idVal);
            return;
        end
    end

    if strcmpi(idMode,'cell') && ~isempty(lk)
        T = lk.obsT;
        hit = T(T.cell_id == idVal, :);
        if isempty(hit); return; end
        for k = 1:height(hit)
            fp = char(hit.rec_path(k));
            [g, rd] = group_date_from_path(fp);
            items(end+1) = struct('folderPath',fp, 'roi',hit.roi_index(k), ...
                'group',g, 'recDate',rd, 'recName',char(hit.rec_name(k)), ...
                'cellId',idVal, ...
                'label',sprintf('cell %d | %s ROI%d', idVal, ...
                                char(hit.rec_name(k)), hit.roi_index(k))); %#ok<AGROW>
        end
        % put the loaded folder first so Load does not jump you elsewhere
        same = strcmpi({items.folderPath}, folderPath);
        if any(same), items = [items(same), items(~same)]; end
    else
        [g, rd] = group_date_from_path(folderPath);
        pp = strsplit(regexprep(folderPath,'[\\/]+$',''), {'\','/'});
        items(1) = struct('folderPath',folderPath, 'roi',idVal, 'group',g, ...
            'recDate',rd, 'recName',pp{end}, 'cellId',NaN, ...
            'label',sprintf('ROI %d | %s', idVal, pp{end}));
    end
end

function [group, recDate] = group_date_from_path(fp)
%GROUP_DATE_FROM_PATH  Recover Genotype/MMDD when the path is in the archive.
%  Only consumer is the Vglut2/1124 one-frame breath fix in load_trace, so a
%  blank pair is harmless for any other recording.
    group = ''; recDate = '';
    parts = strsplit(regexprep(fp,'[\\/]+$',''), {'\','/'});
    known = {'Vglut2','Vgat','Sst','ChAT','Sert'};
    for i = 1:numel(parts)
        if any(strcmpi(parts{i}, known))
            group = parts{i};
            if i < numel(parts) && ~isempty(regexp(parts{i+1},'^\d{4}$','once'))
                recDate = parts{i+1};
            end
            return;
        end
    end
end

% -----------------------------------------------------------------------
function [t, dff, bw, on_t] = load_trace(folderPath, roi, group, recDate, P)
% Minimal load/align for the preview trace (matches panel-1 alignment).
% on_t = inspiration-onset times (s) for counting breaths in a window.
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
bp = dir(fullfile(folderPath,'breath_peak_pc1.mat'));
ip = dir(fullfile(folderPath,'breath_insp_start_pc1.mat'));
assert(~isempty(df),'No *_ch1_dFF.mat in %s', folderPath);
assert(~isempty(bp),'No breath_peak_pc1.mat in %s', folderPath);
fps = detect_session_fps(folderPath, P.fallback_fps);
D  = load(fullfile(df(1).folder, df(1).name),'dFF'); dff_all = double(D.dFF);
BP = load(fullfile(bp(1).folder, bp(1).name));
assert(roi>=1 && roi<=size(dff_all,2),'ROI %d out of range (1..%d)',roi,size(dff_all,2));
bw = detrend(double(BP.breath(:))); bw(1:min(P.nDrop,numel(bw))) = []; bw = bw - mean(bw);
isV1124 = strcmpi(group,'Vglut2') && strcmp(recDate,'1124');
if isV1124, bw = [bw(1); bw(1:end-1)]; end
T = min(size(dff_all,1), numel(bw));
dff = dff_all(1:T, roi); bw = bw(1:T); t = (0:T-1)'/fps;

on_t = [];
if ~isempty(ip)
    IP = load(fullfile(ip(1).folder, ip(1).name));
    oi = round(IP.insp_start_idx(:)) - P.nDrop;     % to post-toss base
    if isV1124, oi = oi + 1; end                    % rising-edge: breath leads 1 frame
    oi = oi(oi>=1 & oi<=T);
    on_t = (oi-1)/fps;     % frame k -> t=(k-1)/fps, matches the preview trace
end
end

function launch_cell_gui(startFolder, startId, idMode, P, doSave, outDirIn, regFile, pngDpi)
%LAUNCH_CELL_GUI  Folder + cell/ROI entry, window pick, render. No cfg files.
if nargin < 7, regFile = ''; end
if nargin < 8 || isempty(pngDpi), pngDpi = 200; end
REG = get_registry(regFile);      % [] = no registry, fall back to cell_link
curIdx = 1;
items  = struct('folderPath',{},'roi',{},'group',{},'recDate',{}, ...
                'recName',{},'cellId',{},'label',{});
lk = [];                                        % cell_link for the loaded session
t = []; dff = []; bw = []; on_t = []; curFolder = '';

f = figure('Color','w','Name','Per-cell temporal-phase GUI', ...
           'Units','normalized','Position',[0.05 0.30 0.90 0.60]);
ax = axes('Parent',f,'Units','normalized','Position',[0.06 0.40 0.90 0.50]);

% ---- row 1: folder / id / load -----------------------------------------
uicontrol(f,'Style','text','Units','normalized','Position',[0.06 0.30 0.06 0.05], ...
    'String','folder:','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hFold = uicontrol(f,'Style','edit','Units','normalized','Position',[0.115 0.305 0.50 0.055], ...
    'String',startFolder,'FontSize',9,'HorizontalAlignment','left','Callback',@load_folder);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.625 0.305 0.075 0.055], ...
    'String','Browse...','Callback',@browse);
uicontrol(f,'Style','text','Units','normalized','Position',[0.715 0.30 0.03 0.05], ...
    'String','id:','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hId = uicontrol(f,'Style','edit','Units','normalized','Position',[0.745 0.305 0.05 0.055], ...
    'String',num2str(startId),'FontSize',11,'Callback',@load_folder);
hMode = uicontrol(f,'Style','popupmenu','Units','normalized','Position',[0.805 0.305 0.075 0.055], ...
    'String',{'cell','roi'},'FontSize',10, ...
    'Value',1+strcmpi(idMode,'roi'),'Callback',@load_folder);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.89 0.305 0.07 0.055], ...
    'String','Load','FontWeight','bold','Callback',@load_folder);

% ---- row 2: trigger window + manual y-limits ----------------------------
% Window for the dF/F averages and the heatmap. Blank = auto = P.trigWinIBI x the
% median inter-breath interval, TOTAL width. The spike histograms are fixed at
% 1 IBI for the test / 2 IBI for the plot and ignore this box.
% Label box ends before the next label starts: at 0.40 wide this ran to x 0.46
% and sat on top of "trig win", which is a collision that predates the 260815
% controls.
uicontrol(f,'Style','text','Units','normalized','Position',[0.06 0.20 0.32 0.05], ...
    'String','window (s):  start end   e.g.  0 30', ...
    'HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
uicontrol(f,'Style','text','Units','normalized','Position',[0.39 0.20 0.135 0.05], ...
    'String','trig win (s, total):','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hWin = uicontrol(f,'Style','edit','Units','normalized','Position',[0.532 0.205 0.06 0.05], ...
    'String','','FontSize',11,'TooltipString', ...
    ['Total width of the dF/F averages and the heatmap, in seconds. Blank = 2 x median IBI. ' ...
     'The spike histograms ignore this: they are always 1 IBI for the test, 2 IBI for the plot.']);
uicontrol(f,'Style','text','Units','normalized','Position',[0.605 0.20 0.11 0.05], ...
    'String','dFF ylim (lo hi):','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hYdff = uicontrol(f,'Style','edit','Units','normalized','Position',[0.715 0.205 0.075 0.05], ...
    'String','','FontSize',11,'TooltipString', ...
    'y-limits for the three dF/F triggered averages. Blank = auto.');
uicontrol(f,'Style','text','Units','normalized','Position',[0.805 0.20 0.10 0.05], ...
    'String','epc ylim (lo hi):','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hYepc = uicontrol(f,'Style','edit','Units','normalized','Position',[0.905 0.205 0.075 0.05], ...
    'String','','FontSize',11,'TooltipString', ...
    'y-limits for both spike histograms, in spikes/cycle %. Blank = auto.');

% POOLED CELL vs SINGLE ROI. Pooled is what the batch does: every recording of
% the cell, with the recording dropdown choosing which one supplies the trace and
% the crop.
hPool = uicontrol(f,'Style','popupmenu','Units','normalized','Position',[0.06 0.255 0.13 0.05], ...
    'String',{'pooled cell','single ROI'},'FontSize',10, ...
    'TooltipString', ['pooled cell = identical to the batch (all recordings of ' ...
                      'this cell); single ROI = only the selected recording']);

% Which of the cell's recordings to plot. In POOLED mode the statistics still use
% every recording -- this picks the one whose trace and average projection are
% shown (P.traceObs). In SINGLE ROI mode it picks the only recording used at all.
uicontrol(f,'Style','text','Units','normalized','Position',[0.195 0.25 0.035 0.05], ...
    'String','rec:','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hRec = uicontrol(f,'Style','popupmenu','Units','normalized','Position',[0.232 0.255 0.34 0.05], ...
    'String',{'(none)'},'Value',1,'FontSize',9,'Enable','off', ...
    'TooltipString', ['Which recording supplies the trace and the avg-projection crop. ' ...
                      'Pooled statistics still use every recording in the list.'], ...
    'Callback',@pick_rec);

% Crop margin around the ROI, in microns. Applies on the next Render.
uicontrol(f,'Style','text','Units','normalized','Position',[0.578 0.25 0.055 0.05], ...
    'String','crop um:','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hPad = uicontrol(f,'Style','edit','Units','normalized','Position',[0.635 0.255 0.04 0.05], ...
    'String',num2str(P.crop_um),'FontSize',11,'TooltipString', ...
    ['TOTAL width of the avg-projection crop in microns, centred on the ROI. ' ...
     'A fixed physical field of view, so every cell is shown at the same scale ' ...
     'regardless of its size. Applies on the next Render.']);

% Contrast stretch percentiles for the crop. No gamma is applied, so these two
% numbers ARE the tone mapping.
uicontrol(f,'Style','text','Units','normalized','Position',[0.681 0.25 0.042 0.05], ...
    'String','clip %:','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hClip = uicontrol(f,'Style','edit','Units','normalized','Position',[0.725 0.255 0.07 0.05], ...
    'String',sprintf('%g %g', P.clip_pct(1), P.clip_pct(2)),'FontSize',10, ...
    'TooltipString', ['Low and high percentiles of the crop used as black and ' ...
     'white points, e.g. "1 99.5". Contrast is LINEAR between them -- there is ' ...
     'no gamma. Widen (0 100) for the raw range, narrow for more contrast.']);

% Which corner of the crop the scale bar sits in. On the SAME row as the crop
% controls: the row below is fully occupied by the window and the two ylim boxes,
% and putting it there overlapped all three.
uicontrol(f,'Style','text','Units','normalized','Position',[0.801 0.25 0.03 0.05], ...
    'String','bar:','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);
hBarPos = uicontrol(f,'Style','popupmenu','Units','normalized','Position',[0.833 0.255 0.15 0.05], ...
    'String',{'lower left','lower right','upper left','upper right'},'FontSize',9, ...
    'TooltipString','Corner of the avg-projection crop holding the scale bar.');

% ---- row 3: window controls / render / nav ------------------------------
hEdit = uicontrol(f,'Style','edit','Units','normalized','Position',[0.06 0.09 0.14 0.075], ...
    'String','0 30','FontSize',11,'Callback',@apply);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.21 0.09 0.065 0.075], ...
    'String','Apply','Callback',@apply);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.285 0.09 0.055 0.075], ...
    'String','Full','Callback',@(s,e) reset_full());
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.36 0.09 0.13 0.075], ...
    'String','Render figure','FontWeight','bold','Callback',@render);
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.51 0.09 0.085 0.075], ...
    'String','<< Prev','Callback',@(s,e) step(-1));
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.60 0.09 0.085 0.075], ...
    'String','Next >>','Callback',@(s,e) step(1));

% ---- row 4: raw-F trace export -----------------------------------------
% A SEPARATE figure, not a mode of the summary: one panel, the same wide trace
% and the same window, with raw F in counts instead of dF/F. It uses the
% recording chosen in the 'rec' dropdown, because raw F cannot be pooled --
% counts from two recordings at different zoom, laser power and PMT gain are not
% on the same scale, which is the whole reason dF/F exists.
% Renders nothing else and touches nothing else: the summary figure, the batch
% and every other export are unaffected.
uicontrol(f,'Style','pushbutton','Units','normalized','Position',[0.36 0.005 0.13 0.05], ...
    'String','Raw F trace','Callback',@render_rawF, ...
    'TooltipString', ['Single-panel figure: raw F (counts) + breath, for the ' ...
    'recording selected in "rec", over the current window. Not baseline ' ...
    'corrected -- bleaching and focus drift stay in.']);
hRawStat = uicontrol(f,'Style','text','Units','normalized','Position',[0.50 0.0 0.19 0.055], ...
    'String','','HorizontalAlignment','left','BackgroundColor','w','FontSize',8);
hStat = uicontrol(f,'Style','text','Units','normalized','Position',[0.70 0.06 0.28 0.105], ...
    'String','','HorizontalAlignment','left','BackgroundColor','w','FontSize',9);

load_folder();

% ========================================================================
    function browse(~,~)
        d0 = get(hFold,'String');
        if ~isfolder(d0), d0 = pwd; end
        d = uigetdir(d0, 'Pick a recording folder');
        if ischar(d) && isfolder(d)
            set(hFold,'String',d); load_folder();
        end
    end

    function load_folder(~,~)
        fp = strtrim(get(hFold,'String'));
        if ~isfolder(fp)
            set(hStat,'String',sprintf('not a folder:\n%s', fp)); return;
        end
        idVal = str2double(get(hId,'String'));
        if ~isfinite(idVal) || idVal < 1
            set(hStat,'String','id must be a positive integer'); return;
        end
        modes = get(hMode,'String'); mode = modes{get(hMode,'Value')};

        lk = find_cell_link(fp);
        % Cell mode needs SOME identity source. The registry covers the archive;
        % cell_link covers recordings it does not. Only drop to ROI if neither
        % knows this folder.
        regHere = [];
        if ~isempty(REG)
            kHere = reg_rows_for_folder(REG, fp);
            if ~isempty(kHere), regHere = sort(unique([REG(kHere).cellNum])); end
        end
        if strcmpi(mode,'cell') && isempty(lk) && isempty(regHere)
            set(hMode,'Value',2); mode = 'roi';       % automatic fallback
        end

        % which cells does the TYPED folder actually contain?
        % Registry first, so the numbers quoted here match the batch filenames.
        here = [];
        if ~isempty(regHere)
            here = regHere;
        elseif ~isempty(lk)
            mm   = strcmpi(string(lk.obsT.rec_path), fp);
            here = lk.obsT.cell_id(mm);
            here = sort(here(~isnan(here))).';
        end

        items = items_for_id(REG, lk, fp, idVal, mode);
        if isempty(items)
            if strcmpi(mode,'cell')
                set(hStat,'String',sprintf('cell %d not found.\ncells in this folder: %s', ...
                    idVal, mat2str(here)));
            else
                set(hStat,'String','no ROI built');
            end
            return;
        end
        curIdx = 1;
        refresh_rec_list();
        load_current();

        % Loading a cell that is not in the typed folder jumps you to whichever
        % recording does hold it -- say so instead of doing it silently.
        if strcmpi(mode,'cell') && ~any(here == idVal)
            set(hStat,'String', sprintf(['cell %d is NOT in the folder you typed.\n' ...
                'showing %s instead.\ncells in that folder: %s'], ...
                idVal, items(1).recName, mat2str(here)));
        elseif strcmpi(mode,'cell')
            % NOTE: do not read hStat back and re-format it -- uicontrol turns a
            % multi-line String into a 2-D char array, and %s would flatten it
            % column-major into garbage. Always rebuild the message.
            it = items(curIdx);
            set(hStat,'String', sprintf('cell %d in %d FOV(s)\n%s ROI%d\ncells here: %s', ...
                idVal, numel(items), it.recName, it.roi, mat2str(here)));
        end
    end

    function load_current()
        it = items(curIdx);
        curFolder = it.folderPath;
        try
            [t, dff, bw, on_t] = load_trace(curFolder, it.roi, it.group, it.recDate, P);
        catch ME
            set(hStat,'String', sprintf('load error: %s', ME.message)); return;
        end
        yyaxis(ax,'left');  cla(ax);
        yyaxis(ax,'right'); cla(ax);
        yyaxis(ax,'right'); plot(ax, t, bw, '-','Color',[0.6 0.6 0.6],'LineWidth',0.6);
        set(ax,'YColor',[0.6 0.6 0.6],'YTick',[]); ylabel(ax,'breath');
        yyaxis(ax,'left');  plot(ax, t, dff, '-','Color',P.dffColor,'LineWidth',0.8);
        set(ax,'YColor','k'); ylabel(ax,'\DeltaF/F');
        xlim(ax,[t(1) t(end)]); xlabel(ax,'Time (s)'); box(ax,'off');
        if isnan(it.cellId)
            ttl = sprintf('[%d/%d]  ROI%d   %s', curIdx, numel(items), it.roi, it.recName);
        else
            ttl = sprintf('[FOV %d/%d]  cell %d  (ROI%d here)   %s', ...
                          curIdx, numel(items), it.cellId, it.roi, it.recName);
        end
        title(ax, ttl, 'Interpreter','none');
        apply([],[]);
    end

    function step(d)
        if isempty(items); return; end
        curIdx = curIdx + d;
        if curIdx < 1,            curIdx = numel(items); end
        if curIdx > numel(items), curIdx = 1;            end
        refresh_rec_list();      % keep the dropdown showing where we are
        load_current();
    end

    function pick_rec(~,~)
        %PICK_REC  Recording chosen from the dropdown.
        if isempty(items); return; end
        curIdx = get(hRec,'Value');
        load_current();
    end

    function refresh_rec_list()
        %REFRESH_REC_LIST  Repopulate the recording dropdown from items.
        %  Called whenever items change (Load) or curIdx moves (Prev/Next), so the
        %  dropdown and the plotted recording can never disagree.
        if isempty(items)
            set(hRec,'String',{'(none)'},'Value',1,'Enable','off'); return;
        end
        curIdx = min(max(curIdx,1), numel(items));
        lbl = cell(1,numel(items));
        for k = 1:numel(items)
            lbl{k} = sprintf('%d/%d   %s   ROI%d', k, numel(items), ...
                             items(k).recName, items(k).roi);
        end
        set(hRec,'String',lbl,'Value',curIdx,'Enable','on');
    end

    function apply(~,~)
        if isempty(t), return; end
        v = sscanf(strrep(get(hEdit,'String'), ',', ' '), '%f');
        if numel(v) >= 2
            a = max(min(v(1),v(2)), t(1));
            b = min(max(v(1),v(2)), t(end));
            if b > a
                xlim(ax, [a b]);
                mwsp = t>=a & t<=b;
                if nnz(mwsp) > 2
                    yyaxis(ax,'left');  ylim(ax, padlim(dff(mwsp)));
                    yyaxis(ax,'right'); ylim(ax, padlim(bw(mwsp)));
                    yyaxis(ax,'left');
                end
                nB = nnz(on_t>=a & on_t<=b);
                it = items(curIdx);
                if isnan(it.cellId), idstr = sprintf('ROI %d', it.roi);
                else,                idstr = sprintf('cell %d', it.cellId); end
                set(hStat,'String', sprintf('%s | %s\nwindow %.2f-%.2f s | %d breaths', ...
                    idstr, it.recName, a, b, nB));
            else
                set(hStat,'String','end must be > start');
            end
        else
            set(hStat,'String','enter two numbers, e.g. 0 30');
        end
    end

    function reset_full()
        if isempty(t), return; end
        set(hEdit,'String', sprintf('%.0f %.0f', t(1), t(end)));
        xlim(ax,[t(1) t(end)]);
        yyaxis(ax,'left');  ylim(ax, padlim(dff));
        yyaxis(ax,'right'); ylim(ax, padlim(bw));
        yyaxis(ax,'left'); set(hStat,'String','full trace');
    end

    function render(~,~)
        % Renders through temporal_phase_cell_fig_260812 -- the SAME function the
        % batch calls -- so a GUI figure and a batch figure of the same cell are
        % identical apart from the window and the y-limits you set here.
        %
        %   POOLED CELL : every recording of the cell goes in, exactly as the
        %                 batch does it. Prev/Next chooses which recording
        %                 supplies the trace and the projection (P.traceObs).
        %   SINGLE ROI  : only the loaded (recording, ROI). Pooling one thing is
        %                 a no-op, so this is the per-ROI view, same code path.
        if isempty(items), set(hStat,'String','nothing loaded'); return; end
        v = sscanf(strrep(get(hEdit,'String'), ',', ' '), '%f');
        if numel(v) >= 2 && ~isempty(t)
            wsp = [max(min(v(1),v(2)),t(1)), min(max(v(1),v(2)),t(end))];
            if wsp(2) <= wsp(1), wsp = []; end
        else
            wsp = [];                    % [] = 30 s from the middle
        end

        modes  = get(hPool,'String');
        pooled = strcmpi(modes{get(hPool,'Value')}, 'pooled cell');
        it     = items(curIdx);

        if pooled
            use = 1:numel(items);
            Pr_traceObs = curIdx;        % Prev/Next picks the trace recording
        else
            use = curIdx;
            Pr_traceObs = 1;
        end
        OBSs = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
        for q = use
            OBSs(end+1) = struct('folder',items(q).folderPath, 'roi',items(q).roi, ...
                'recName',items(q).recName, 'group',items(q).group, ...
                'recDate',items(q).recDate); %#ok<AGROW>
        end

        Pr = P; Pr.trace_xlim_sp = wsp; Pr.cellId = it.cellId;
        Pr.traceObs = Pr_traceObs;
        % Per-cell breath-strip restriction (P.breathObsByCell). Applied only when
        % the cell is POOLED -- with one recording on screen the strip already
        % comes from that recording and an index would be meaningless. The cell is
        % still pooled everywhere else in the figure; this narrows the strip alone.
        if pooled && isfield(P,'breathObsByCell') && ~isempty(P.breathObsByCell) ...
                  && ~isnan(it.cellId)
            hit = P.breathObsByCell(:,1) == it.cellId;
            if any(hit)
                Pr.breathObs = P.breathObsByCell(find(hit,1), 2);
            end
        end
        tv = sscanf(strrep(get(hWin,'String'), ',', ' '), '%f');
        if ~isempty(tv) && isfinite(tv(1)) && tv(1) > 0, Pr.trigWin_sec = tv(1);
        else,                                            Pr.trigWin_sec = []; end
        Pr.ylim_dff = parse_ylim(get(hYdff,'String'));
        Pr.ylim_epc = parse_ylim(get(hYepc,'String'));
        pv = sscanf(strtrim(get(hPad,'String')), '%f');
        if ~isempty(pv) && isfinite(pv(1)) && pv(1) > 0, Pr.crop_um = pv(1); end
        bp = get(hBarPos,'String');
        Pr.scalebarCorner = bp{get(hBarPos,'Value')};
        cv = sscanf(strrep(get(hClip,'String'), ',', ' '), '%f');
        % Reject anything that is not a valid ordered pair inside [0 100] rather
        % than passing it to prctile, which would error mid-render.
        if numel(cv) >= 2 && all(isfinite(cv(1:2))) && cv(1) >= 0 && ...
           cv(2) <= 100 && cv(2) > cv(1)
            Pr.clip_pct = [cv(1) cv(2)];
        end

        set(hStat,'String','rendering...'); drawnow;
        try
            [figR, proj] = cell_cycle_fig_260831(OBSs, Pr);
        catch ME
            set(hStat,'String', sprintf('render failed:\n%s', ME.message)); return;
        end

        if doSave
            outDir = outDirIn;
            if isempty(outDir), outDir = fullfile(items(curIdx).folderPath,'temporal_phase_cell'); end
            if ~isfolder(outDir), mkdir(outDir); end
            if pooled
                if isnan(it.cellId), idp = sprintf('ROI%02d_pooled', it.roi);
                else,                idp = sprintf('cell%03d_pooled%d', it.cellId, numel(use)); end
                stem = regexprep(sprintf('%s_%s', it.recName, idp), '[\\/:*?"<>|]', '_');
            else
                if isnan(it.cellId), idp = sprintf('ROI%02d', it.roi);
                else,                idp = sprintf('cell%03d_ROI%02d', it.cellId, it.roi); end
                stem = regexprep(sprintf('%s_%s', it.recName, idp), '[\\/:*?"<>|]', '_');
            end
            base = fullfile(outDir, stem);
            exportgraphics(figR, [base '.png'], 'Resolution',pngDpi, 'BackgroundColor','white');
            exportgraphics(figR, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
            if proj.have_proj
                save_avgproj_png_260812(proj, P, fullfile(outDir, sprintf('%s_avgproj.png', stem)));
            end
            set(hStat,'String', sprintf('saved ->\n%s.png', stem));
        else
            set(hStat,'String','rendered');
        end
    end

    function render_rawF(~,~)
    %RENDER_RAWF  Single-panel raw-F trace for the recording in the 'rec' menu.
    %  Deliberately does NOT reuse render()'s OBS list: that list is the POOLED
    %  set, and pooling raw counts across recordings is meaningless. One
    %  recording in, one panel out.
        if isempty(items)
            set(hRawStat,'String','load a cell first'); return;
        end
        it = items(curIdx);
        % Same window parse and same clamp as render(), against the loaded time
        % base t. [] means "30 s from the middle", which the panel handles.
        v = sscanf(strrep(get(hEdit,'String'), ',', ' '), '%f');
        if numel(v) >= 2 && ~isempty(t)
            wsp = [max(min(v(1),v(2)),t(1)), min(max(v(1),v(2)),t(end))];
            if wsp(2) <= wsp(1), wsp = []; end
        else
            wsp = [];
        end

        Pr = P;
        Pr.trace_xlim_sp = wsp;
        Pr.cellId   = it.cellId;
        Pr.traceObs = 1;                    % the one recording handed over below
        OB = struct('folder',it.folderPath, 'roi',it.roi, 'recName',it.recName, ...
                    'group',it.group, 'recDate',it.recDate);

        set(hRawStat,'String','rendering raw F...'); drawnow;
        try
            [figF, infoF] = trace_rawF_panel_260817(OB, Pr);
        catch ME
            set(hRawStat,'String', sprintf('raw F failed: %s', ME.message)); return;
        end

        if doSave
            outDir = outDirIn;
            if isempty(outDir), outDir = fullfile(it.folderPath,'temporal_phase_cell'); end
            if ~isfolder(outDir), mkdir(outDir); end
            if isnan(it.cellId), idp = sprintf('ROI%02d', it.roi);
            else,                idp = sprintf('cell%03d_ROI%02d', it.cellId, it.roi); end
            % _rawF in the stem so it can never overwrite a summary render.
            stem = regexprep(sprintf('%s_%s_rawF', it.recName, idp), '[\\/:*?"<>|]', '_');
            base = fullfile(outDir, stem);
            exportgraphics(figF, [base '.png'], 'Resolution',pngDpi, 'BackgroundColor','white');
            exportgraphics(figF, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
            set(hRawStat,'String', sprintf('saved %s.png\nF %.0f-%.0f counts', ...
                stem, infoF.F_min, infoF.F_max));
        else
            set(hRawStat,'String', sprintf('rendered\nF %.0f-%.0f counts', ...
                infoF.F_min, infoF.F_max));
        end
    end

    function yl = padlim(y)
        lo = min(y); hi = max(y); pad = 0.05*max(hi-lo, eps); yl = [lo-pad hi+pad];
    end
end

% -----------------------------------------------------------------------
function yl = parse_ylim(str)
%PARSE_YLIM  "lo hi" -> [lo hi]; anything else (blank, one number, lo>=hi) -> [].
%  [] is the signal to auto-scale, so a half-typed box never freezes an axis.
yl = [];
v = sscanf(strrep(strtrim(str), ',', ' '), '%f');
if numel(v) >= 2 && all(isfinite(v(1:2))) && v(2) > v(1), yl = [v(1) v(2)]; end
end
