% examine_intervening_whisk_epochs_RZ.m
%
% Viewer for the INTERVENING-whisk epochs used in
% whisk_LR_coherence_phase_intervening_RZ.m.
%
% Same definition / same detection pipeline:
%   Intervening whisk = whisk onset (trough of BP whisker angle, inside the
%   whisk-epoch mask) whose nearest inspiration is FARTHER than
%   INSP_EXCLUDE_WIN_S away (i.e. it does NOT fall within a breath).
%
% This script does NOT compute coherence. It opens ONE FIGURE PER intervening
% whisk onset so you can inspect each epoch:
%   top  : L (blue) and R (red) BP whisker angle around the onset
%   bottom: breath (BP) around the onset
% with the whisk onset at t=0, inspiration events marked, the whisk-epoch
% mask shaded, and the coherence window (+/- WIN_SEC/2) marked.
%
% Controls below: cap the number of figures, save them to disk, or restrict
% to a single animal/run.

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [5 50];    % whisker BP
BREATH_BP = [1 20];    % breath BP for inspiration detection

ENV_THR    = 5;        % whisk epoch threshold, deg
ENV_SMOOTH = 0.05;     % s
MIN_DUR    = 1;        % s
MERGE_GAP  = 0.2;      % s

% inspiration detection
MIN_INSP_DIST_S = 0.05;
PROM_MAD        = 0.5;

% intervening rule (must match the coherence script)
INSP_EXCLUDE_WIN_S = 0.100;   % reject whisks within this of inspiration
                              % NOTE: breathing is ~6 Hz (IBI ~0.16 s), so a
                              % whisk can be at most ~0.18 s from a breath.
                              % Values >~0.2 s yield ZERO real intervening
                              % whisks (only flat-tail artifacts). Keep <= 0.15.

% whisk onset detection
MIN_WHISK_DIST_S   = 0.025;   % minimum distance between whisk onsets
WHISK_TROUGH_PROM  = 1.0;     % deg, trough prominence for BP whisker

% coherence window (drawn for reference)
WIN_SEC = 1.0;
MIN_WHISK_FRAC_IN_WINDOW = 0.40;

% ---------- viewer-specific ----------
PLOT_PRE_S  = 0.75;    % seconds shown before onset
PLOT_POST_S = 0.75;    % seconds shown after onset
MAX_FIGS    = 16;      % stop after this many figures (Inf = no cap)
SHOW_FIGS   = true;    % open figures on screen
SAVE_FIGS   = false;   % also save each figure as PNG
SEL_ANIMAL  = "";      % "" = all animals, else e.g. "5916300"
SEL_RUN     = 0;       % 0  = all runs,    else e.g. 1

EXCLUDE = "5840027";

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'intervening_whisk_epochs_view');
if SAVE_FIGS && ~isfolder(outDir), mkdir(outDir); end

ord = 4;
[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

winSamp   = round(WIN_SEC * fpsW);
prePlotN  = round(PLOT_PRE_S  * fpsW);
postPlotN = round(PLOT_POST_S * fpsW);

colL = [0.10 0.30 0.85];   % L blue
colR = [0.85 0.10 0.10];   % R red
colB = [0.20 0.55 0.25];   % breath green

nFigs  = 0;
nInter = 0;
stop_all = false;

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

% ============================ SESSION LOOP ============================

for e = 1:numel(sess)

    if stop_all, break; end

    animal = sess{e}{1};
    kRun   = sess{e}{2};

    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE)), continue; end
    if strlength(SEL_ANIMAL) > 0 && ~strcmp(animal, char(SEL_ANIMAL)), continue; end
    if SEL_RUN > 0 && kRun ~= SEL_RUN, continue; end

    % ============================ WHISKER ============================
    try
        whiskCsv = pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kRun));
        M = dlc_gate_interp(whiskCsv, 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message);
        continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));

    tW = (0:numel(La)-1)' / fpsW;
    N  = numel(tW);

    La0 = fillmissing(La(:), 'linear', 'EndValues', 'nearest');
    Ra0 = fillmissing(Ra(:), 'linear', 'EndValues', 'nearest');
    LaRaw = La0;            % raw filled angle (no BP, no demean) - for display
    RaRaw = Ra0;
    La0 = La0 - mean(La0, 'omitnan');
    Ra0 = Ra0 - mean(Ra0, 'omitnan');

    xL = filtfilt(b_wbp, a_wbp, La0);   % BP only used for detection
    xR = filtfilt(b_wbp, a_wbp, Ra0);

    % ======================= WHISK EPOCH MASK ========================
    env_raw = abs(hilbert(xL));
    env = movmean(env_raw, max(1, round(ENV_SMOOTH * fpsW)));
    is_whisk = env > ENV_THR;
    ep = bool_to_epochs(is_whisk, tW, MIN_DUR, MERGE_GAP);

    m = false(N,1);
    for q = 1:size(ep,1)
        m = m | (tW >= ep(q,1) & tW <= ep(q,2));
    end
    if nnz(m) < 50
        warning('no valid whisk epoch %s n%d', animal, kRun); continue;
    end

    % ============================ BREATH =============================
    try
        noseCsv = pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kRun));
        Pn = thermal_resolve_paths(noseCsv, dataRoot);
        if ~isfile(Pn.breath)
            warning('canonical breath file not found %s n%d', animal, kRun); continue;
        end

        Bs = load(Pn.breath);
        br = Bs.breath(:);
        fb = double(Bs.fps);
        tB = (0:numel(br)-1)' / fb;

        brw = interp1(tB, br, tW, 'linear', NaN);
        breath_finite = isfinite(brw);   % TRUE only where breath was actually recorded
        if nnz(breath_finite) < 50
            warning('too few finite breath samples %s n%d', animal, kRun); continue;
        end
        brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');
        brw_bp = filtfilt(b_br, a_br, brw);

        event_trace = brw_bp;
        prom = PROM_MAD * mad(event_trace, 1);
        if ~isfinite(prom) || prom <= 0
            prom = 0.5 * std(event_trace, 'omitnan');
        end

        [~, locs] = findpeaks(-event_trace, ...
            'MinPeakDistance', round(MIN_INSP_DIST_S * fpsW), ...
            'MinPeakProminence', prom);
        insp_t = tW(locs);
    catch ME
        warning('breath load/detect failed %s n%d: %s', animal, kRun, ME.message);
        continue;
    end
    if isempty(insp_t)
        warning('no inspiration events detected %s n%d', animal, kRun); continue;
    end

    % ======================= DETECT WHISK ONSETS =======================
    tOnL = detect_whisk_troughs(xL, tW, m, fpsW, MIN_WHISK_DIST_S, WHISK_TROUGH_PROM);
    tOnR = detect_whisk_troughs(xR, tW, m, fpsW, MIN_WHISK_DIST_S, WHISK_TROUGH_PROM);
    tOnAll = merge_close_events(sort([tOnL(:); tOnR(:)]), MIN_WHISK_DIST_S);
    if isempty(tOnAll)
        warning('no whisk onsets %s n%d', animal, kRun); continue;
    end

    % ===================== PLOT EACH INTERVENING WHISK =====================
    for ii = 1:numel(tOnAll)

        t0 = tOnAll(ii);
        dInsp = min(abs(insp_t - t0));

        % keep only intervening (not within a breath)
        if dInsp <= INSP_EXCLUDE_WIN_S, continue; end

        c = round(t0 * fpsW) + 1;

        % coherence window must be valid + mostly inside whisk epoch
        cwIdx = c + ((0:winSamp-1)' - floor(winSamp/2));
        if cwIdx(1) < 1 || cwIdx(end) > N, continue; end
        if mean(m(cwIdx)) < MIN_WHISK_FRAC_IN_WINDOW, continue; end
        if ~all(isfinite(xL(cwIdx))) || ~all(isfinite(xR(cwIdx))), continue; end
        % REQUIRE the window to lie inside the actually-recorded breath, so the
        % "intervening" flag cannot come from the flat extrapolated breath tail.
        if ~all(breath_finite(cwIdx)), continue; end

        nInter = nInter + 1;

        if ~SHOW_FIGS && ~SAVE_FIGS, continue; end
        if nFigs >= MAX_FIGS
            stop_all = true;
            fprintf('Reached MAX_FIGS=%d, stopping (more intervening epochs exist).\n', MAX_FIGS);
            break;
        end

        % display window (clamped)
        d0 = max(1, c - prePlotN);
        d1 = min(N, c + postPlotN);
        dispIdx = (d0:d1)';
        tRel = (dispIdx - c) / fpsW;

        segL = LaRaw(dispIdx);    % RAW filled whisker angle
        segR = RaRaw(dispIdx);
        segB = brw(dispIdx);      % RAW (LP<1Hz-subtracted) breath

        tWin = [tRel(1) tRel(end)];

        if SHOW_FIGS, vis = 'on'; else, vis = 'off'; end
        fig = figure('Color','w','Position',[150 120 760 560], 'Visible', vis);

        % ---- top: L/R RAW whisker angle ----
        ax1 = subplot(2,1,1); hold(ax1,'on'); grid(ax1,'on');
        hL = plot(ax1, tRel, segL, '-', 'Color', colL, 'LineWidth', 1.6);
        hR = plot(ax1, tRel, segR, '-', 'Color', colR, 'LineWidth', 1.6);
        xlim(ax1, tWin);
        yl1 = ylim(ax1);                                           % lock to data
        shade_mask_epochs(ax1, ep, t0, tWin, yl1, [0.90 0.90 0.90]);  % whisk epochs
        shade_span(ax1, [-WIN_SEC/2 WIN_SEC/2], yl1, [1.0 0.93 0.75]); % coh window
        xline(ax1, 0, 'k-', 'LineWidth', 1.2);                     % whisk onset
        mark_insp(ax1, insp_t, t0, tWin);                         % inspirations
        ylim(ax1, yl1);
        ylabel(ax1, 'raw whisker angle (deg)');
        legend(ax1, [hL hR], {'L','R'}, 'Location','northeastoutside', 'FontSize',8);
        title(ax1, sprintf(['%s n%d  -  intervening whisk #%d (onset t=%.3f s)\n' ...
              'nearest inspiration %.0f ms away  (rule > %.0f ms)'], ...
              animal, kRun, ii, t0, 1000*dInsp, 1000*INSP_EXCLUDE_WIN_S), ...
              'Interpreter','none');
        box(ax1,'off');

        % ---- bottom: RAW breath ----
        ax2 = subplot(2,1,2); hold(ax2,'on'); grid(ax2,'on');
        plot(ax2, tRel, segB, '-', 'Color', colB, 'LineWidth', 1.6);
        xlim(ax2, tWin);
        yl2 = ylim(ax2);
        shade_mask_epochs(ax2, ep, t0, tWin, yl2, [0.90 0.90 0.90]);
        shade_span(ax2, [-WIN_SEC/2 WIN_SEC/2], yl2, [1.0 0.93 0.75]);
        xline(ax2, 0, 'k-', 'LineWidth', 1.2);
        mark_insp(ax2, insp_t, t0, tWin);
        ylim(ax2, yl2);
        xlabel(ax2, 'time from whisk onset (s)');
        ylabel(ax2, 'raw breath (a.u.)');
        title(ax2, 'raw breath (dashed magenta = inspiration events; grey = whisk epoch; tan = coh window)', ...
              'Interpreter','none');
        box(ax2,'off');

        linkaxes([ax1 ax2], 'x');

        nFigs = nFigs + 1;

        if SAVE_FIGS
            fn = sprintf('intervening_%s_n%d_onset%03d_t%07.3f.png', animal, kRun, ii, t0);
            exportgraphics(fig, fullfile(outDir, fn), 'Resolution',150, 'BackgroundColor','white');
        end
        if ~SHOW_FIGS, close(fig); end
    end
end

fprintf('\n%d intervening-whisk epochs found, %d figures produced.\n', nInter, nFigs);
if SAVE_FIGS, fprintf('Saved PNGs to %s\n', outDir); end

% ============================= HELPERS =============================

function shade_span(ax, xspan, yl, col)
% Vertical shaded band over [xspan(1) xspan(2)], spanning the given y-limits.
    p = patch(ax, [xspan(1) xspan(2) xspan(2) xspan(1)], ...
              [yl(1) yl(1) yl(2) yl(2)], col, ...
              'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
    uistack(p, 'bottom');
end

function shade_mask_epochs(ax, ep, t0, tWin, yl, col)
% Shade whisk-epoch intervals (onset-relative) within the view window.
    for q = 1:size(ep,1)
        a = max(ep(q,1) - t0, tWin(1));
        b = min(ep(q,2) - t0, tWin(2));
        if b <= a, continue; end
        p = patch(ax, [a b b a], [yl(1) yl(1) yl(2) yl(2)], col, ...
                  'EdgeColor','none', 'FaceAlpha',0.7, 'HandleVisibility','off');
        uistack(p, 'bottom');
    end
end

function mark_insp(ax, insp_t, t0, tWin)
% Dashed magenta vertical lines at inspiration events within the window.
    rel = insp_t - t0;
    rel = rel(rel >= tWin(1) & rel <= tWin(2));
    for k = 1:numel(rel)
        xline(ax, rel(k), '--', 'Color',[0.8 0.2 0.8], 'LineWidth',1.0, ...
              'HandleVisibility','off');
    end
end

function tOn = detect_whisk_troughs(x, t, mask, fs, minDistS, prom)
    x = x(:); t = t(:); mask = logical(mask(:));
    ep = bool_to_epochs(mask, t, 0, 0);
    tOn = zeros(0,1);
    minDistSamples = max(1, round(minDistS * fs));
    for q = 1:size(ep,1)
        idx = find(t >= ep(q,1) & t <= ep(q,2));
        if numel(idx) < 5, continue; end
        sig = -x(idx);
        if isfinite(prom) && prom > 0
            [~, locs] = findpeaks(sig, 'MinPeakDistance', minDistSamples, ...
                                       'MinPeakProminence', prom);
        else
            [~, locs] = findpeaks(sig, 'MinPeakDistance', minDistSamples);
        end
        if ~isempty(locs), tOn = [tOn; t(idx(locs))]; end %#ok<AGROW>
    end
    tOn = sort(tOn);
end

function tOut = merge_close_events(tIn, mergeWinS)
    tIn = sort(tIn(:));
    if isempty(tIn), tOut = zeros(0,1); return; end
    clusters = {}; cur = tIn(1);
    for i = 2:numel(tIn)
        if tIn(i) - cur(end) <= mergeWinS
            cur(end+1,1) = tIn(i); %#ok<AGROW>
        else
            clusters{end+1} = cur; %#ok<AGROW>
            cur = tIn(i);
        end
    end
    clusters{end+1} = cur;
    tOut = zeros(numel(clusters),1);
    for i = 1:numel(clusters), tOut(i) = mean(clusters{i}); end
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a = logical(a(:));
    d = diff([false; a; false]);
    s = find(d == 1);
    e = find(d == -1) - 1;
    ep = [t(s), t(e)];
    if isempty(ep), return; end
    ep2 = ep(1,:);
    for i = 2:size(ep,1)
        if ep(i,1) - ep2(end,2) <= mergeGap
            ep2(end,2) = ep(i,2);
        else
            ep2(end+1,:) = ep(i,:); %#ok<AGROW>
        end
    end
    ep = ep2;
    ep = ep(ep(:,2) - ep(:,1) >= minDur, :);
end

function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv'));
    S = {}; key = {};
    for i = 1:numel(d)
        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        k = sprintf('%s_%s', tok{1}, tok{2});
        if any(strcmp(key,k)), continue; end
        key{end+1} = k; %#ok<AGROW>
        S{end+1} = {tok{1}, str2double(tok{2})}; %#ok<AGROW>
    end
    if ~isempty(S)
        an = cellfun(@(c) str2double(c{1}), S);
        rn = cellfun(@(c) c{2}, S);
        [~,o] = sortrows([an(:) rn(:)]);
        S = S(o);
    end
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    [~,ix] = max(arrayfun(@(x) bestnum(x.name), d));
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name, 'best-(\d+)', 'tokens');
    if isempty(tok), n = 0; else, n = str2double(tok{1}{1}); end
end
