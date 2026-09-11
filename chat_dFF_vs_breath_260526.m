% chat_dFF_vs_breath_260526.m
% -----------------------------------------------------------------------
%  Overview per recording: dF/F (stacked) with breath OVERLAID on top.
%    - breath waveform (blue) at top of the stack, red dots = inspiration onsets
%    - all ROI dF/F traces stacked below, shared/zoomable time axis
%
%  3-line title:
%    1) <date> / <IO or cell#>            (folder structure)
%    2) acquisition time = raw TIFF last-modified time
%    3) coord | zoom | laser(lp) | frames | breath peak f +/- FWHM (lo–hi Hz)
%
%  Alignment matches the coherence script: toss first nDrop breath frames,
%  truncate both to common length. Breath is 2P-frame-triggered (fps = imaging).
%
%  Inputs per folder: ca_spike_data.mat, *_ch1_dFF.mat, *_ch1_meta.mat,
%                     *DLC_*breath_peak_data.mat, raw <recName>.tif
%  Dependencies: detect_session_fps.m, Chronux (mtspectrumc)
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(genpath(fullfile(scriptDir, 'chronux_2_12')));

%% ===================== USER-EDITABLE PARAMETERS ======================
rootPath        = 'C:\Users\Admin\Desktop\ChAT_analysis';
nDrop           = 30;          % breath frames to toss (match calcium)
fallback_fps    = 30;          % used only if meta fps cannot be detected
TW              = 5;           % multitaper time-bandwidth for breath PSD
f_breath_search = [0.2 4];     % Hz, search band for breath peak
fmin            = 0.05;        % Hz, PSD lower bound
doSave          = true;        % save dFF_vs_breath.png into each folder
% =====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

M = dir(fullfile(rootPath, '**', 'ca_spike_data.mat'));
fprintf('Found %d recording(s).\n', numel(M));

nplot = 0;
for i = 1:numel(M)
    fp = M(i).folder;
    [~, nm] = fileparts(fp);

    bp = dir(fullfile(fp, '*DLC*breath_peak_data.mat'));
    df = dir(fullfile(fp, '*_ch1_dFF.mat'));
    if isempty(bp) || isempty(df)
        fprintf('skip (missing breath/dFF): %s\n', fp); continue;
    end

    %% ---- folder structure: <date>/<group>/<recording> ----
    grpDir  = fileparts(fp);              % .../<date>/<group>
    dateDir = fileparts(grpDir);          % .../<date>
    [~, grpName]  = fileparts(grpDir);    % cell1 / cell2 / IO
    [~, dateName] = fileparts(dateDir);   % 0521 ...

    %% ---- acquisition time = raw TIFF last-modified time ----
    tiff = dir(fullfile(fp, [nm '.tif']));
    if isempty(tiff), tiff = dir(fullfile(fp, '*.tif')); end
    if ~isempty(tiff)
        acqTime = datestr(tiff(1).datenum, 'yyyy-mm-dd HH:MM:SS');
    else
        acqTime = '(no tif)';
    end

    %% ---- meta: coordinate / zoom / laser / frames ----
    mt = dir(fullfile(fp, '*_ch1_meta.mat'));
    coordStr = ''; zoomStr = ''; lpStr = ''; frStr = '';
    if ~isempty(mt)
        MM = load(fullfile(mt(1).folder, mt(1).name));
        if isfield(MM,'motorPosition') && numel(MM.motorPosition)>=3
            coordStr = sprintf('x%.0f y%.0f z%.0f', MM.motorPosition(1:3));
        end
        if isfield(MM,'zoomFactor'),     zoomStr = sprintf('%gx',    MM.zoomFactor); end
        if isfield(MM,'laserPower_pct'), lpStr   = sprintf('%glp',   MM.laserPower_pct); end
        if isfield(MM,'framesPerSlice'), frStr   = sprintf('%df',    MM.framesPerSlice); end
    end

    %% ---- load dFF + breath ----
    D   = load(fullfile(df(1).folder, df(1).name), 'dFF');
    dFF = double(D.dFF);
    BP  = load(fullfile(bp(1).folder, bp(1).name));
    fps = detect_session_fps(fp, fallback_fps);

    bw = detrend(double(BP.breath(:))); bw(1:nDrop) = []; bw = bw - mean(bw);
    if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train) == numel(BP.breath)
        ev = double(BP.insp_onsets_train(:) ~= 0);
    else
        ev = zeros(numel(BP.breath),1); oi = round(BP.insp_onset_idx(:));
        ev(oi(oi>=1 & oi<=numel(ev))) = 1;
    end
    ev(1:nDrop) = [];

    T = min([size(dFF,1), numel(bw), numel(ev)]);
    dFF = dFF(1:T,:); bw = bw(1:T); ev = ev(1:T);
    t = (0:T-1)'/fps; ot = t(ev > 0); N = size(dFF,2);

    %% ---- breath peak frequency +/- FWHM (event-train PSD) ----
    be = ev - mean(ev);
    pB.Fs = fps; pB.tapers = [TW, 2*TW-1]; pB.pad = 0;
    pB.fpass = [fmin, fps/2]; pB.err = 0;
    [Sb, fb] = mtspectrumc(be, pB); Sb = Sb(:); fb = fb(:);
    m = fb >= f_breath_search(1) & fb <= f_breath_search(2);
    [~, rel] = max(Sb(m)); ip = find(m,1)+rel-1; f_pk = fb(ip);
    half = Sb(ip)/2;
    lo = ip; while lo>1 && Sb(lo)>half, lo=lo-1; end
    hi = ip; while hi<numel(fb) && Sb(hi)>half, hi=hi+1; end
    f_fwhm = [max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];

    %% ---- figure: dFF stacked + breath overlaid on top ----
    rng = max(dFF,[],1) - min(dFF,[],1);
    sp  = max(prctile(rng, 80), 0.3);

    f = figure('Color','w','Name',nm,'Units','normalized','Position',[0.05 0.08 0.9 0.82]);

    % --- breath panel (own space, ~1/4 height) ---
    ax1 = subplot(4,1,1); hold(ax1,'on');
    plot(ax1, t, bw, 'Color', [0 0.35 1], 'LineWidth', 0.7);
    plot(ax1, ot, interp1(t, bw, ot, 'nearest'), 'r.', 'MarkerSize', 8);
    ylabel(ax1,'breath');
    xlim(ax1, [0 t(end)]);
    ttl = { sprintf('%s / %s', dateName, grpName), ...
            sprintf('acq: %s', acqTime), ...
            sprintf('%s | %s | %s | %s | breath %.2f Hz (%.2f–%.2f Hz)', ...
                    coordStr, zoomStr, lpStr, frStr, f_pk, f_fwhm(1), f_fwhm(2)) };
    title(ax1, ttl, 'Interpreter','none');

    % --- dF/F stacked (bottom 3/4) ---
    ax2 = subplot(4,1,2:4); hold(ax2,'on');
    for r = 1:N, plot(ax2, t, dFF(:,r) + (r-1)*sp, 'LineWidth', 0.5); end
    ylim(ax2, [-sp, N*sp]);
    set(ax2, 'YTick', (0:N-1)*sp, 'YTickLabel', 1:N);
    xlabel(ax2,'Time (s)'); ylabel(ax2,'ROI dF/F (stacked)');
    xlim(ax2, [0 t(end)]);

    linkaxes([ax1 ax2], 'x');
    ax1.Toolbar.Visible = 'off'; ax2.Toolbar.Visible = 'off';

    if doSave
        exportgraphics(f, fullfile(fp,'dFF_vs_breath.png'), ...
            'Resolution',150, 'BackgroundColor','white');
    end
    nplot = nplot + 1;
    fprintf('plotted: %s / %s / %s  (N=%d, breath %.2f Hz)\n', dateName, grpName, nm, N, f_pk);
end
fprintf('\nDone. %d figure(s).\n', nplot);
