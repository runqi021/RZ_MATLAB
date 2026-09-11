% whisk_LR_coherence_phase_intervening_RZ.m
%
% Correct intervening-whisk coherence.
%
% Key definition:
%   Inspiratory whisk = whisk onset within +/-100 ms of inspiration.
%   Intervening whisk = whisk onset whose nearest inspiration is farther than 100 ms away.
%
% Therefore:
%   coherence windows are centered on WHISK ONSETS,
%   not inspiration events.
%
% Whisker:
%   La/Ra -> fill NaNs -> demean
%   BP 5-50 Hz
%   whisk epoch mask from Hilbert envelope of LEFT BP whisker
%
% Breathing:
%   load canonical _breath.mat
%   BP 1-20 Hz
%   detect inspiration events
%
% Coherence:
%   stack fixed-length L/R whisker windows around intervening whisk onsets
%   run Chronux coherencyc across trials
%
% Phase convention:
%   data1 = L, data2 = R
%   phi > 0 means R leads L
%   phi < 0 means L leads R
%
% Plotting:
%   f, C, and phi are treated as the master frequency vectors.
%   Cerr/phistd are used only if their size matches f.
%   The script NEVER truncates f/C/phi to match broken error-bar vectors.

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

% intervening rule
INSP_EXCLUDE_WIN_S = 0.060;   % reject whisks within +/-100 ms of inspiration

% whisk onset detection
MIN_WHISK_DIST_S   = 0.025;   % minimum distance between whisk onsets
WHISK_TROUGH_PROM  = 1.0;     % deg, trough prominence for BP whisker

% coherence windows
WIN_SEC = 1.0;                   % window centered on intervening whisk onset
MIN_WHISK_FRAC_IN_WINDOW = 0.40; % require enough window inside whisk epoch

% Chronux
TW    = 3;
FPASS = [1 60];
ALPHA = 0.001;

EXCLUDE = "5840027";

doSave = false;

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

outDir = fullfile(char(dataRoot), 'whisk_LR_coherence_phase_intervening_CORRECTED');
if doSave && ~isfolder(outDir)
    mkdir(outDir);
end

ord = 4;

[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

winSamp = round(WIN_SEC * fpsW);
relIdx  = (0:winSamp-1)' - floor(winSamp/2);

assert(winSamp >= 2*TW, 'WIN_SEC too short for TW.');

pc.Fs       = fpsW;
pc.tapers   = [TW, 2*TW-1];
pc.pad      = 0;
pc.fpass    = FPASS;
pc.err      = [2 ALPHA];
pc.trialave = 1;

% ============================ ACCUMULATORS ============================

L_win = zeros(winSamp, 0);
R_win = zeros(winSamp, 0);

nSess = 0;
nInter = 0;
nWhiskOnTot = 0;
nWhiskOnNearInsp = 0;
nInspTot = 0;

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

% ============================ SESSION LOOP ============================

for e = 1:numel(sess)

    animal = sess{e}{1};
    kRun   = sess{e}{2};

    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE))
        continue;
    end

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

    La0 = La0 - mean(La0, 'omitnan');
    Ra0 = Ra0 - mean(Ra0, 'omitnan');

    xL = filtfilt(b_wbp, a_wbp, La0);
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
        warning('no valid whisk epoch %s n%d', animal, kRun);
        continue;
    end

    % ============================ BREATH =============================

    try
        noseCsv = pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kRun));
        Pn = thermal_resolve_paths(noseCsv, dataRoot);

        if ~isfile(Pn.breath)
            warning('canonical breath file not found %s n%d', animal, kRun);
            continue;
        end

        Bs = load(Pn.breath);

        br = Bs.breath(:);
        fb = double(Bs.fps);
        tB = (0:numel(br)-1)' / fb;

        brw = interp1(tB, br, tW, 'linear', NaN);
        breath_finite = isfinite(brw);   % TRUE only where breath was actually recorded

        if nnz(breath_finite) < 50
            warning('too few finite breath samples %s n%d', animal, kRun);
            continue;
        end

        brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');
        brw_bp = filtfilt(b_br, a_br, brw);

        event_trace = brw_bp;

        prom = PROM_MAD * mad(event_trace, 1);
        if ~isfinite(prom) || prom <= 0
            prom = 0.5 * std(event_trace, 'omitnan');
        end

        % This follows your existing inspiration detection convention.
        % If your breath sign is flipped, switch -event_trace to event_trace.
        [~, locs] = findpeaks(-event_trace, ...
            'MinPeakDistance', round(MIN_INSP_DIST_S * fpsW), ...
            'MinPeakProminence', prom);

        insp_t = tW(locs);
        nInspTot = nInspTot + numel(insp_t);

    catch ME
        warning('breath load/detect failed %s n%d: %s', animal, kRun, ME.message);
        continue;
    end

    if isempty(insp_t)
        warning('no inspiration events detected %s n%d', animal, kRun);
        continue;
    end

    % ======================= DETECT WHISK ONSETS =======================

    % Protraction onset = trough of BP whisker angle.
    tOnL = detect_whisk_troughs(xL, tW, m, fpsW, MIN_WHISK_DIST_S, WHISK_TROUGH_PROM);
    tOnR = detect_whisk_troughs(xR, tW, m, fpsW, MIN_WHISK_DIST_S, WHISK_TROUGH_PROM);

    % Merge L and R onsets that are basically the same whisk cycle.
    tOnAll = merge_close_events(sort([tOnL(:); tOnR(:)]), MIN_WHISK_DIST_S);

    if isempty(tOnAll)
        warning('no whisk onsets %s n%d', animal, kRun);
        continue;
    end

    nWhiskOnTot = nWhiskOnTot + numel(tOnAll);

    % ====================== KEEP ONLY INTERVENING WHISKS ======================

    sess_added = 0;

    for ii = 1:numel(tOnAll)

        t0 = tOnAll(ii);

        dInsp = min(abs(insp_t - t0));

        % KEY RULE:
        % reject inspiratory whisks, keep true intervening whisks.
        if dInsp <= INSP_EXCLUDE_WIN_S
            nWhiskOnNearInsp = nWhiskOnNearInsp + 1;
            continue;
        end

        c = round(t0 * fpsW) + 1;
        idx = c + relIdx;

        if idx(1) < 1 || idx(end) > N
            continue;
        end

        if mean(m(idx)) < MIN_WHISK_FRAC_IN_WINDOW
            continue;
        end

        % REQUIRE the window inside actually-recorded breath, so "intervening"
        % cannot be an artifact of the flat extrapolated breath tail.
        if ~all(breath_finite(idx))
            continue;
        end

        wl = xL(idx);
        wr = xR(idx);

        if ~all(isfinite(wl)) || ~all(isfinite(wr))
            continue;
        end

        wl = wl - mean(wl);
        wr = wr - mean(wr);

        L_win(:, end+1) = wl; %#ok<SAGROW>
        R_win(:, end+1) = wr; %#ok<SAGROW>

        sess_added = sess_added + 1;
    end

    if sess_added > 0
        nSess = nSess + 1;
        nInter = nInter + sess_added;

        fprintf('  %s n%d: %d intervening whisk windows, %d all whisk onsets, %d insp\n', ...
            animal, kRun, sess_added, numel(tOnAll), numel(insp_t));
    end
end

fprintf('\n%d sessions contributed\n', nSess);
fprintf('%d intervening-whisk windows kept\n', nInter);
fprintf('%d total whisk onsets detected\n', nWhiskOnTot);
fprintf('%d whisk onsets rejected as near inspiration\n', nWhiskOnNearInsp);
fprintf('%d inspiration events total\n', nInspTot);

if nInter < 2
    error('Too few intervening-whisk windows (%d) to estimate coherence.', nInter);
end

% ===================== POPULATION COHERENCE =====================

% Chronux coherencyc (continuous) output order:
% [C12, C, phi, S12, S1, S2, f, confC, phistd, Cerr]
%   C12 = complex coherency, C = magnitude, phi = phase, f = frequencies.
% NOTE: there is NO zerosp output for the continuous version (only in the
% point-process coherencycpt/coherencycpb). Mapping it in shifts f onto S2
% and breaks the frequency axis, so keep this exact order.
[~, C, phi, ~, ~, ~, f, confC, phistd, Cerr] = coherencyc(L_win, R_win, pc);

f   = f(:);
C   = C(:);
phi = phi(:);

% f, C, and phi are the master vectors.
% Do NOT trim them to match Cerr or phistd.
if numel(f) ~= numel(C) || numel(f) ~= numel(phi)
    error('Chronux output mismatch: numel(f)=%d, numel(C)=%d, numel(phi)=%d', ...
        numel(f), numel(C), numel(phi));
end

fprintf('Chronux returned %d frequency bins: %.2f to %.2f Hz\n', ...
    numel(f), min(f), max(f));

% ---------- phase std handling ----------
phistd = squeeze(phistd);

if isempty(phistd)
    phistd = nan(size(f));

elseif isscalar(phistd)
    phistd = repmat(phistd, size(f));

else
    phistd = phistd(:);

    if numel(phistd) ~= numel(f)
        warning('phistd size mismatch: numel(phistd)=%d, numel(f)=%d. Phase band skipped.', ...
            numel(phistd), numel(f));
        phistd = nan(size(f));
    end
end

% ---------- Cerr handling ----------
% Only use Cerr if it exactly matches f.
% If not, skip the coherence error band.
% Never trim f/C/phi to match Cerr.

Clo = nan(size(f));
Chi = nan(size(f));

Cerr = squeeze(Cerr);

if isempty(Cerr)
    warning('Cerr empty. Coherence error band skipped.');

elseif size(Cerr,1) == 2 && size(Cerr,2) == numel(f)
    Clo = Cerr(1,:).';
    Chi = Cerr(2,:).';

elseif size(Cerr,2) == 2 && size(Cerr,1) == numel(f)
    Clo = Cerr(:,1);
    Chi = Cerr(:,2);

else
    warning('Cerr size mismatch: size(Cerr)=[%s], numel(f)=%d. Coherence error band skipped.', ...
        num2str(size(Cerr)), numel(f));
end

% clean / order coherence band if present
tmpLo = min(Clo, Chi);
tmpHi = max(Clo, Chi);
Clo = max(0, min(1, tmpLo));
Chi = max(0, min(1, tmpHi));

phi_deg    = rad2deg(phi);
phistd_deg = rad2deg(phistd);

sig = C >= confC;

[Cpk, ipk] = max(C);

fprintf('Peak coherence C=%.3f at %.2f Hz, phase=%.1f deg\n', ...
    Cpk, f(ipk), phi_deg(ipk));

% ============================== FIGURE ==============================

col  = [0.85 0.10 0.10];
colS = [0.6 0.6 0.6];

fig = figure('Color','w','Position',[180 160 720 740]);

% ---- Panel A: coherence magnitude ----

axA = subplot(2,1,1); hold(axA,'on'); grid(axA,'on');

goodBand = isfinite(f) & isfinite(Clo) & isfinite(Chi);

if any(goodBand)
    fill(axA, [f(goodBand); flipud(f(goodBand))], ...
        [Chi(goodBand); flipud(Clo(goodBand))], ...
        col, ...
        'FaceAlpha',0.20, ...
        'EdgeColor','none');
end

plot(axA, f, C, '-', ...
    'Color', col, ...
    'LineWidth', 2);

yline(axA, confC, 'k--', sprintf('conf %.2f', confC), ...
    'LabelHorizontalAlignment','left');

plot(axA, f(ipk), Cpk, 'o', ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor',col, ...
    'MarkerSize',6);

xlim(axA, FPASS);
ylim(axA, [0 1]);

xlabel(axA, 'frequency (Hz)');
ylabel(axA, 'coherence |C_{LR}|');

title(axA, sprintf(['Bilateral whisking coherence - INTERVENING whisks only\n' ...
    'N=%d windows, %d sessions, exclude +/-%.0f ms from inspiration, win=%.1fs, TW=%d'], ...
    nInter, nSess, 1000*INSP_EXCLUDE_WIN_S, WIN_SEC, TW), ...
    'Interpreter','tex');

box(axA,'off');

% ---- Panel B: phase ----

axB = subplot(2,1,2); hold(axB,'on'); grid(axB,'on');

yline(axB, 0, 'k-');

phi_hi = phi_deg + 1.96 * phistd_deg;
phi_lo = phi_deg - 1.96 * phistd_deg;

goodPhaseBand = isfinite(f) & isfinite(phi_hi) & isfinite(phi_lo);

if any(goodPhaseBand)
    fill(axB, [f(goodPhaseBand); flipud(f(goodPhaseBand))], ...
        [phi_hi(goodPhaseBand); flipud(phi_lo(goodPhaseBand))], ...
        colS, ...
        'FaceAlpha',0.20, ...
        'EdgeColor','none');
end

% all phase, gray
plot(axB, f, phi_deg, '-', ...
    'Color', colS, ...
    'LineWidth', 1);

% significant phase, red
ph_sig = phi_deg;
ph_sig(~sig) = NaN;

plot(axB, f, ph_sig, '-', ...
    'Color', col, ...
    'LineWidth', 2.2);

xlim(axB, FPASS);
ylim(axB, [-180 180]);
yticks(axB, -180:90:180);

xlabel(axB, 'frequency (Hz)');
ylabel(axB, '\Delta phase L,R (deg)');

title(axB, '\Delta phase; positive = R leads L, negative = L leads R; red = significant', ...
    'Interpreter','tex');

box(axB,'off');

% ============================== SAVE ==============================

if doSave

    exportgraphics(fig, fullfile(outDir,'whisk_LR_coherence_phase_intervening_CORRECTED.png'), ...
        'Resolution',200, 'BackgroundColor','white');

    exportgraphics(fig, fullfile(outDir,'whisk_LR_coherence_phase_intervening_CORRECTED.pdf'), ...
        'ContentType','vector', 'BackgroundColor','white');

    save(fullfile(outDir,'whisk_LR_coherence_phase_intervening_CORRECTED_data.mat'), ...
        'f','C','phi','phi_deg','phistd','phistd_deg','Clo','Chi', ...
        'confC','sig','Cpk','ipk', ...
        'nSess','nInter','nWhiskOnTot','nWhiskOnNearInsp','nInspTot', ...
        'WIN_SEC','TW','FPASS','ALPHA','WHISK_BP','BREATH_BP', ...
        'INSP_EXCLUDE_WIN_S','MIN_WHISK_DIST_S','WHISK_TROUGH_PROM', ...
        'ENV_THR','ENV_SMOOTH','MIN_DUR','MERGE_GAP','fpsW','EXCLUDE');

    fprintf('Saved corrected intervening-whisk coherence figure + .mat to %s\n', outDir);
end

% ============================= HELPERS =============================

function tOn = detect_whisk_troughs(x, t, mask, fs, minDistS, prom)
% Detect protraction onsets as troughs of BP whisker angle.
% Runs findpeaks(-x) inside whisking mask epochs only.

    x = x(:);
    t = t(:);
    mask = logical(mask(:));

    ep = bool_to_epochs(mask, t, 0, 0);

    tOn = zeros(0,1);

    minDistSamples = max(1, round(minDistS * fs));

    for q = 1:size(ep,1)

        idx = find(t >= ep(q,1) & t <= ep(q,2));

        if numel(idx) < 5
            continue;
        end

        sig = -x(idx);

        if isfinite(prom) && prom > 0
            [~, locs] = findpeaks(sig, ...
                'MinPeakDistance', minDistSamples, ...
                'MinPeakProminence', prom);
        else
            [~, locs] = findpeaks(sig, ...
                'MinPeakDistance', minDistSamples);
        end

        if ~isempty(locs)
            tOn = [tOn; t(idx(locs))]; %#ok<AGROW>
        end
    end

    tOn = sort(tOn);
end

function tOut = merge_close_events(tIn, mergeWinS)
% Merge events closer than mergeWinS.
% Uses mean time of each cluster.

    tIn = sort(tIn(:));

    if isempty(tIn)
        tOut = zeros(0,1);
        return;
    end

    clusters = {};
    cur = tIn(1);

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

    for i = 1:numel(clusters)
        tOut(i) = mean(clusters{i});
    end
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)

    a = logical(a(:));
    d = diff([false; a; false]);

    s = find(d == 1);
    e = find(d == -1) - 1;

    ep = [t(s), t(e)];

    if isempty(ep)
        return;
    end

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

    S = {};
    key = {};

    for i = 1:numel(d)

        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');

        if isempty(tok)
            continue;
        end

        k = sprintf('%s_%s', tok{1}, tok{2});

        if any(strcmp(key,k))
            continue;
        end

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

    if isempty(tok)
        n = 0;
    else
        n = str2double(tok{1}{1});
    end
end