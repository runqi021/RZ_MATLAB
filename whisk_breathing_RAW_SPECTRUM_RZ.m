% whisk_breath_raw_fft_spectrum_RZ.m
%
% Goal:
%   1) Breathing:
%      Use NEW raw breathing output:
%           <stem>_breath_raw.mat
%      Variables expected:
%           L_raw, R_raw, breath_raw, fps, t
%
%      breathing spectrum = breath_raw, no LP-sub, no invert.
%      Then choose:
%           SPECTRUM_SIGNAL = "raw" : demean only
%           SPECTRUM_SIGNAL = "dt"  : derivative after demean
%
%      For each session, randomly circular-shift the breathing trace,
%      split into 3 chunks, FFT each chunk, then average across chunks.
%
%   2) Whisking:
%      detect whisk epochs using standard BP-envelope method,
%      but compute FFT on raw demeaned whisker angle within each epoch.
%      Then choose raw or derivative with the same SPECTRUM_SIGNAL option.
%
% Plotting style:
%   FFT power stays LINEAR.
%   mean, SEM, mean+SEM, mean-SEM are all computed in LINEAR power.
%   MATLAB plots them on log axes.
%
% Figures:
%   Fig 1: real log-log, power on log y-axis
%   Fig 2: linear frequency x-axis, power on log y-axis
%
% Important:
%   f=0 is excluded from plotting.

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP = [5 50];    % only used for epoch detection
ENV_THR = 5;          % hard envelope threshold, deg
ENV_SMOOTH = 0.05;    % s
MIN_DUR = 1;          % s
MERGE_GAP = 0.2;      % s

BREATH_CHUNKS = 3;    % randomly chop each breathing session into 3 chunks
RNG_SEED = 1;         % reproducible random chopping

FMAX = 80;            % max frequency to plot, Hz
DF = 0.1;             % common frequency grid, Hz

NORMALIZE_EACH = true;   % true = compare spectral shape, false = raw power

SPECTRUM_SIGNAL = "dt";   % "raw" = demeaned signal, "dt" = derivative signal

POWER_FLOOR_FRAC = 1e-12;  % relative floor for log plotting only
YLOG_LIM = [];             % [] = auto, or e.g. [1e-12 1]

EXCLUDE = "5840027";

% ======================================================================

rng(RNG_SEED);

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

fGrid = (0:DF:FMAX)';

[b_wbp, a_wbp] = butter(4, WHISK_BP/(fpsW/2), 'bandpass');

P_breath = zeros(0, numel(fGrid));   % rows = breathing chunks
P_whisk  = zeros(0, numel(fGrid));   % rows = whisk epochs

nBreathSession = 0;
nBreathChunk   = 0;
nEpoch         = 0;

% ============================ SESSION LOOP ============================

ad = dir(char(dataRoot));

for ai = 1:numel(ad)

    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit'))
        continue;
    end

    if any(strcmp(ad(ai).name, EXCLUDE))
        continue;
    end

    animal = ad(ai).name;

    rr = dir(fullfile(char(dataRoot), animal, 'cam1_*'));
    [~,o] = sort({rr.name});
    rr = rr(o);

    for kk = 1:numel(rr)

        % ============================ WHISKER ============================

        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', animal, kk))))
            continue;
        end

        try
            M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk)), 0.6);
        catch ME
            warning('whisk load failed %s n%d: %s', animal, kk, ME.message);
            continue;
        end

        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));

        tW = (0:numel(La)-1)' / fpsW;

        La0 = fillmissing(La(:), 'linear', 'EndValues', 'nearest');
        Ra0 = fillmissing(Ra(:), 'linear', 'EndValues', 'nearest');

        % session-level demean before epoch detection
        La0 = La0 - mean(La0, 'omitnan');
        Ra0 = Ra0 - mean(Ra0, 'omitnan');

        % use BP only for epoch detection
        La_bp5_50 = filtfilt(b_wbp, a_wbp, La0);

        env_raw = abs(hilbert(La_bp5_50));
        env = movmean(env_raw, max(1, round(ENV_SMOOTH * fpsW)));

        is_whisk = env > ENV_THR;
        ep = bool_to_epochs(is_whisk, tW, MIN_DUR, MERGE_GAP);

        % ======================= WHISK EPOCH FFT =========================
        % FFT is computed on raw whisker angle inside each detected epoch.
        % Then prep_spectrum_signal chooses "raw" or "dt".

        for q = 1:size(ep,1)

            idx = tW >= ep(q,1) & tW <= ep(q,2);

            if nnz(idx) < round(MIN_DUR * fpsW)
                continue;
            end

            xL = La0(idx);
            xR = Ra0(idx);

            xL = prep_spectrum_signal(xL, fpsW, SPECTRUM_SIGNAL);
            xR = prep_spectrum_signal(xR, fpsW, SPECTRUM_SIGNAL);

            [fL, pL] = one_sided_fft_power(xL, fpsW, NORMALIZE_EACH);
            [fR, pR] = one_sided_fft_power(xR, fpsW, NORMALIZE_EACH);

            pLi = interp1(fL, pL, fGrid, 'linear', NaN);
            pRi = interp1(fR, pR, fGrid, 'linear', NaN);

            P_whisk(end+1,:) = mean([pLi(:), pRi(:)], 2, 'omitnan')'; %#ok<SAGROW>
            nEpoch = nEpoch + 1;
        end

        % ======================== RAW BREATH FFT =========================
        % Breathing now comes from <stem>_breath_raw.mat:
        %   breath_raw = raw deg C, no LP-sub, no invert.
        % Then prep_spectrum_signal chooses "raw" or "dt".

        try
            noseCsv = pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk));
            Pn = thermal_resolve_paths(noseCsv, dataRoot);

            rawMat = resolve_breath_raw_mat(Pn, noseCsv);

            if strlength(rawMat) == 0 || ~isfile(rawMat)
                warning('raw breath file not found %s n%d', animal, kk);
                continue;
            end

            Bs = load(rawMat);

            if ~isfield(Bs, 'breath_raw')
                warning('breath_raw missing in %s', rawMat);
                continue;
            end

            br = Bs.breath_raw(:);     % raw deg C, no LP-sub, no invert

            if isfield(Bs, 'fps')
                fb = double(Bs.fps);
            else
                warning('fps missing in %s', rawMat);
                continue;
            end

            br = prep_spectrum_signal(br, fb, SPECTRUM_SIGNAL);

            Nbr = numel(br);
            chunkN = floor(Nbr / BREATH_CHUNKS);

            if chunkN < 4
                continue;
            end

            % random circular shift before chopping
            sh = randi(Nbr);
            br = circshift(br, sh);

            nBreathSession = nBreathSession + 1;

            for c = 1:BREATH_CHUNKS

                i1 = (c-1)*chunkN + 1;
                i2 = c*chunkN;

                x = br(i1:i2);
                x = x - mean(x, 'omitnan');

                [fbins, pb] = one_sided_fft_power(x, fb, NORMALIZE_EACH);
                pbi = interp1(fbins, pb, fGrid, 'linear', NaN);

                P_breath(end+1,:) = pbi(:)'; %#ok<SAGROW>
                nBreathChunk = nBreathChunk + 1;
            end

        catch ME
            warning('raw breath load failed %s n%d: %s', animal, kk, ME.message);
        end
    end
end

fprintf('Breathing sessions: %d\n', nBreathSession);
fprintf('Breathing chunks:   %d\n', nBreathChunk);
fprintf('Whisking epochs:    %d\n', nEpoch);
fprintf('Spectrum signal mode: %s\n', SPECTRUM_SIGNAL);

% ============================ LINEAR MEAN +/- SEM =======================
% Everything here is in LINEAR power.
% No dB transform.
% No log(mean).
% No log(std).
% MATLAB log axes handle the plotting.

P_breath(P_breath <= 0) = NaN;
P_whisk(P_whisk <= 0)   = NaN;

% ---------------- breathing ----------------

B_mean = mean(P_breath, 1, 'omitnan');
B_n    = sum(isfinite(P_breath), 1);
B_sem  = std(P_breath, 0, 1, 'omitnan') ./ sqrt(B_n);

B_hi = B_mean + B_sem;
B_lo = B_mean - B_sem;

refB = max(P_breath(:), [], 'omitnan');
if isempty(refB) || ~isfinite(refB) || refB <= 0
    refB = 1;
end

floorB = refB * POWER_FLOOR_FRAC;

B_mean(~isfinite(B_mean) | B_mean <= floorB) = floorB;
B_hi(~isfinite(B_hi)     | B_hi   <= floorB) = floorB;
B_lo(~isfinite(B_lo)     | B_lo   <= floorB) = floorB;

% ---------------- whisking ----------------

W_mean = mean(P_whisk, 1, 'omitnan');
W_n    = sum(isfinite(P_whisk), 1);
W_sem  = std(P_whisk, 0, 1, 'omitnan') ./ sqrt(W_n);

W_hi = W_mean + W_sem;
W_lo = W_mean - W_sem;

refW = max(P_whisk(:), [], 'omitnan');
if isempty(refW) || ~isfinite(refW) || refW <= 0
    refW = 1;
end

floorW = refW * POWER_FLOOR_FRAC;

W_mean(~isfinite(W_mean) | W_mean <= floorW) = floorW;
W_hi(~isfinite(W_hi)     | W_hi   <= floorW) = floorW;
W_lo(~isfinite(W_lo)     | W_lo   <= floorW) = floorW;

% ============================== FIGURE 1 ==============================
% Real log-log.
% f=0 is excluded.

figure('Color','w','Position',[120 160 1100 420]);

ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');
plot_linear_power_spectrum(ax1, fGrid, B_mean, B_hi, B_lo, true);
xline(ax1, 1,  'b--', '1 Hz');
xline(ax1, 20, 'b--', '20 Hz');
set(ax1, 'XScale', 'log', 'YScale', 'log');
xlim(ax1, [0 FMAX]);
if ~isempty(YLOG_LIM), ylim(ax1, YLOG_LIM); end
xlabel(ax1, 'Frequency (Hz)');
ylabel(ax1, 'Power');
title(ax1, sprintf('Breathing %s spectrum, n=%d chunks', SPECTRUM_SIGNAL, nBreathChunk));
axis(ax1, 'square');
box(ax1, 'off');

ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');
plot_linear_power_spectrum(ax2, fGrid, W_mean, W_hi, W_lo, true);
xline(ax2, 5,  'b--', '5 Hz');
xline(ax2, 50, 'b--', '50 Hz');
set(ax2, 'XScale', 'log', 'YScale', 'log');
xlim(ax2, [DF FMAX]);
if ~isempty(YLOG_LIM), ylim(ax2, YLOG_LIM); end
xlabel(ax2, 'Frequency (Hz)');
ylabel(ax2, 'Power');
title(ax2, sprintf('Whisking %s spectrum, n=%d epochs', SPECTRUM_SIGNAL, nEpoch));
axis(ax2, 'square');
box(ax2, 'off');

sgtitle(sprintf('FFT spectra: real log-log, linear mean \\pm SEM, mode=%s', SPECTRUM_SIGNAL));

% ============================== FIGURE 2 ==============================
% Linear frequency, log y.
% f=0 is excluded.

figure('Color','w','Position',[160 220 1100 420]);

ax3 = subplot(1,2,1); hold(ax3,'on'); grid(ax3,'on');
plot_linear_power_spectrum(ax3, fGrid, B_mean, B_hi, B_lo, false);
xline(ax3, 1,  'b--', '1 Hz');
xline(ax3, 20, 'b--', '20 Hz');
set(ax3, 'YScale', 'log');
xlim(ax3, [DF 30]);
if ~isempty(YLOG_LIM), ylim(ax3, YLOG_LIM); end
xlabel(ax3, 'Frequency (Hz)');
ylabel(ax3, 'Power');
title(ax3, sprintf('Breathing %s spectrum, n=%d chunks', SPECTRUM_SIGNAL, nBreathChunk));
axis(ax3, 'square');
box(ax3, 'off');

ax4 = subplot(1,2,2); hold(ax4,'on'); grid(ax4,'on');
plot_linear_power_spectrum(ax4, fGrid, W_mean, W_hi, W_lo, false);
xline(ax4, 5,  'b--', '5 Hz');
xline(ax4, 50, 'b--', '50 Hz');
set(ax4, 'YScale', 'log');
xlim(ax4, [DF 30]);
if ~isempty(YLOG_LIM), ylim(ax4, YLOG_LIM); end
xlabel(ax4, 'Frequency (Hz)');
ylabel(ax4, 'Power');
title(ax4, sprintf('Whisking %s spectrum, n=%d epochs', SPECTRUM_SIGNAL, nEpoch));
axis(ax4, 'square');
box(ax4, 'off');

sgtitle(sprintf('FFT spectra: linear x, log y, linear mean \\pm SEM, mode=%s', SPECTRUM_SIGNAL));

% ============================= HELPERS =============================

function x = prep_spectrum_signal(x, fs, mode)

    x = x(:);
    x = fillmissing(x, 'linear', 'EndValues', 'nearest');
    x = x - mean(x, 'omitnan');

    mode = string(mode);

    switch lower(mode)

        case "raw"
            % use demeaned raw signal

        case {"dt", "derivative", "diff"}
            x = gradient(x) * fs;
            x = x - mean(x, 'omitnan');

        otherwise
            error('Unknown SPECTRUM_SIGNAL mode: %s. Use "raw" or "dt".', mode);
    end
end

function [f, P1] = one_sided_fft_power(x, fs, normalizeEach)

    x = x(:);
    x = x(isfinite(x));

    N = numel(x);

    if N < 4
        f = NaN;
        P1 = NaN;
        return;
    end

    x = x - mean(x, 'omitnan');

    Y = fft(x);

    P2 = abs(Y / N).^2;
    P1 = P2(1:floor(N/2)+1);

    if numel(P1) > 2
        P1(2:end-1) = 2 * P1(2:end-1);
    end

    f = fs * (0:floor(N/2))' / N;
    P1 = P1(:);

    if normalizeEach
        s = sum(P1, 'omitnan');
        if isfinite(s) && s > 0
            P1 = P1 / s;
        end
    end
end

function plot_linear_power_spectrum(ax, f, y, yUp, yLo, useLogX)

    f   = f(:);
    y   = y(:);
    yUp = yUp(:);
    yLo = yLo(:);

    % never plot f=0 for spectrum
    good = isfinite(f) & f > 0 & ...
           isfinite(y) & y > 0 & ...
           isfinite(yUp) & yUp > 0 & ...
           isfinite(yLo) & yLo > 0;

    f   = f(good);
    y   = y(good);
    yUp = yUp(good);
    yLo = yLo(good);

    if isempty(f)
        return;
    end

    % shaded SEM band
    fill(ax, [f; flipud(f)], ...
        [yUp; flipud(yLo)], ...
        [0.5 0.5 0.5], ...
        'FaceAlpha', 0.18, ...
        'EdgeColor', 'none');

    % explicit SEM bounds
    if useLogX
        semilogx(ax, f, yUp, 'k--', 'LineWidth', 1); hold(ax,'on');
        semilogx(ax, f, yLo, 'k--', 'LineWidth', 1);
        semilogx(ax, f, y,   'r',   'LineWidth', 2);
    else
        plot(ax, f, yUp, 'k--', 'LineWidth', 1); hold(ax,'on');
        plot(ax, f, yLo, 'k--', 'LineWidth', 1);
        plot(ax, f, y,   'r',   'LineWidth', 2);
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

function rawMat = resolve_breath_raw_mat(Pn, noseCsv)
% Find <stem>_breath_raw.mat next to the thermal .ats.

    rawMat = "";

    candidates = strings(0,1);

    directFields = {'breath_raw', 'breathRaw', 'raw_breath', 'rawBreath'};

    for i = 1:numel(directFields)
        fn = directFields{i};
        if isfield(Pn, fn)
            try
                candidates(end+1,1) = string(Pn.(fn)); %#ok<AGROW>
            catch
            end
        end
    end

    if isfield(Pn, 'ats')
        atsPath = string(Pn.ats);

        if isfile(atsPath)
            [p, stem] = fileparts(atsPath);
            candidates(end+1,1) = string(fullfile(p, stem + "_breath_raw.mat")); %#ok<AGROW>
        end
    end

    if isfield(Pn, 'breath')
        oldBreath = string(Pn.breath);
        [p, ~] = fileparts(oldBreath);

        if isfolder(p)
            d = dir(fullfile(p, '*_breath_raw.mat'));
            for k = 1:numel(d)
                candidates(end+1,1) = string(fullfile(d(k).folder, d(k).name)); %#ok<AGROW>
            end
        end
    end

    if nargin >= 2
        noseCsv = string(noseCsv);
        [p, ~] = fileparts(noseCsv);

        if isfolder(p)
            d = dir(fullfile(p, '*_breath_raw.mat'));
            for k = 1:numel(d)
                candidates(end+1,1) = string(fullfile(d(k).folder, d(k).name)); %#ok<AGROW>
            end
        end
    end

    candidates = unique(candidates, 'stable');

    for i = 1:numel(candidates)
        if strlength(candidates(i)) > 0 && isfile(candidates(i))
            rawMat = candidates(i);
            return;
        end
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