% whisk_LR_xcorr_all_vs_intervening_RZ.m
%
% Population L/R whisker cross-correlation, computed for two sample sets:
%   ALL          : all whisking samples (whisk-epoch mask)
%   INTERVENING  : whisking samples that are NOT near a breath
%                  (per-sample distance to nearest inspiration > INSP_EXCLUDE_WIN_S)
%
% This is the SAMPLE-MASK version of the reference cc_LR in
% whisk_LR_integral_xcorr_filtered_RZ.m:
%   cc_LR = lagged_corr_fast(xL, xR, mask, ...)   % corr[L(t+lag), R(t)]
% so a single short intervening whisk needs no fixed window -- every
% qualifying sample simply enters the correlation.
%
% Whisker:
%   La/Ra -> fill NaNs -> demean -> BP 5-50 Hz, zero-phase filtfilt
%   whisk epoch mask from Hilbert envelope of LEFT BP whisker, thr = 5 deg
%
% Breathing (only used to define inspiration times for the intervening mask):
%   load canonical _breath.mat, BP 1-20 Hz, inspiration = findpeaks(-breath_bp)
%   the flat extrapolated tail (breath not actually recorded) is excluded.
%
% Lag convention (same as reference):
%   cc_LR = corr[L(t+lag), R(t)];  lag > 0 means R leads L.

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [5 50];   % Hz
BREATH_BP = [1 20];   % Hz

ENV_THR    = 5;       % hard envelope threshold, deg
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0.2;     % s

% inspiration detection
MIN_INSP_DIST_S = 0.05;
PROM_MAD        = 0.5;

% intervening rule (per SAMPLE): a whisking sample is "intervening" if the
% nearest inspiration is farther than this away. Breathing is ~6 Hz here, so
% keep this small (e.g. 0.04-0.08 s) or the intervening set is empty.
INSP_EXCLUDE_WIN_S = 0.060;

XC_MAXLAG   = 0.20;   % s, cross-correlation lag range (whisk-cycle scale)
MIN_SAMPLES = 50;     % min masked samples per lag to accept a session

DO_SHUFFLE  = false;  % circular-shift R control (per condition)
N_SHUFF     = 200;
MIN_SHIFT_S = 5.0;

EXCLUDE = "5840027";

doSave = false;

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_LR_xcorr_all_vs_intervening');
if doSave && ~isfolder(outDir), mkdir(outDir); end

ord = 4;
[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(ord, BREATH_BP/(fpsW/2), 'bandpass');

nLag   = round(XC_MAXLAG * fpsW);
lags_t = (-nLag:nLag)' / fpsW;
nL     = numel(lags_t);
exclSamp = round(INSP_EXCLUDE_WIN_S * fpsW);

% ============================ ACCUMULATORS ============================

xc_all   = zeros(0,nL);   % L-R xcorr, all whisking
xc_inter = zeros(0,nL);   % L-R xcorr, intervening whisking

xc_all_shuf   = zeros(0,nL,N_SHUFF);
xc_inter_shuf = zeros(0,nL,N_SHUFF);

nsess = 0;
nAllSampTot   = 0;
nInterSampTot = 0;

% ============================ SESSION LOOP ============================

ad = dir(char(dataRoot));

for ai = 1:numel(ad)

    if ~ad(ai).isdir || ~all(isstrprop(ad(ai).name,'digit')), continue; end
    if any(strcmp(ad(ai).name, EXCLUDE)), continue; end

    animal = ad(ai).name;
    rr = dir(fullfile(char(dataRoot), animal, 'cam1_*'));
    [~,o] = sort({rr.name}); rr = rr(o);

    for kk = 1:numel(rr)

        if isempty(dir(fullfile(char(whiskDir), sprintf('%s_whisk_n%d*DLC*.csv', animal, kk))))
            continue;
        end

        % ============================ WHISKER ============================
        try
            M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kk)), 0.6);
        catch ME
            warning('whisk load failed %s n%d: %s', animal, kk, ME.message);
            continue;
        end

        La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
        Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));

        t = (0:numel(La)-1)' / fpsW;
        N = numel(t);

        La0 = fillmissing(La(:), 'linear', 'EndValues', 'nearest');
        Ra0 = fillmissing(Ra(:), 'linear', 'EndValues', 'nearest');
        La0 = La0 - mean(La0, 'omitnan');
        Ra0 = Ra0 - mean(Ra0, 'omitnan');

        xL = filtfilt(b_wbp, a_wbp, La0);
        xR = filtfilt(b_wbp, a_wbp, Ra0);

        % ======================= WHISK EPOCH MASK ========================
        env_raw = abs(hilbert(xL));
        env = movmean(env_raw, max(1, round(ENV_SMOOTH * fpsW)));
        ep = bool_to_epochs(env > ENV_THR, t, MIN_DUR, MERGE_GAP);
        m = false(N,1);
        for q = 1:size(ep,1)
            m = m | (t >= ep(q,1) & t <= ep(q,2));
        end
        if nnz(m) < MIN_SAMPLES, continue; end

        % ============================ BREATH =============================
        % (inspiration times only; needed for the intervening mask)
        insp_idx = [];
        breath_finite = false(N,1);
        try
            Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kk)), dataRoot);
            if isfile(Pn.breath)
                Bs = load(Pn.breath);
                br = Bs.breath(:); fb = double(Bs.fps); tBr = (0:numel(br)-1)' / fb;
                brw = interp1(tBr, br, t, 'linear', NaN);
                breath_finite = isfinite(brw);          % recorded region only
                if nnz(breath_finite) < MIN_SAMPLES, continue; end
                brw = fillmissing(brw, 'linear', 'EndValues', 'nearest');
                brw_bp = filtfilt(b_br, a_br, brw);
                prom = PROM_MAD * mad(brw_bp, 1);
                if ~isfinite(prom) || prom <= 0, prom = 0.5*std(brw_bp,'omitnan'); end
                [~, locs] = findpeaks(-brw_bp, ...
                    'MinPeakDistance', round(MIN_INSP_DIST_S * fpsW), ...
                    'MinPeakProminence', prom);
                insp_idx = locs;
            end
        catch ME
            warning('breath load failed %s n%d: %s', animal, kk, ME.message);
            continue;
        end
        if isempty(insp_idx), continue; end

        % ===================== BUILD THE TWO MASKS =======================
        % near-inspiration: within exclSamp of any inspiration sample
        insp_train = false(N,1);
        insp_train(insp_idx(insp_idx>=1 & insp_idx<=N)) = true;
        near_insp = movmax(double(insp_train), 2*exclSamp+1) > 0;

        mask_all   = m & breath_finite;                  % all whisking (recorded breath)
        mask_inter = mask_all & ~near_insp;              % intervening whisking

        if nnz(mask_all) < MIN_SAMPLES, continue; end

        % ============================= XCORR =============================
        cc_all = lagged_corr_fast(xL, xR, mask_all, nLag, MIN_SAMPLES);
        if nnz(mask_inter) >= MIN_SAMPLES
            cc_inter = lagged_corr_fast(xL, xR, mask_inter, nLag, MIN_SAMPLES);
        else
            cc_inter = nan(nL,1);
        end

        % ============================= SHUFFLE ============================
        if DO_SHUFFLE
            cc_all_sh   = nan(nL, N_SHUFF);
            cc_inter_sh = nan(nL, N_SHUFF);
            minShift = round(MIN_SHIFT_S * fpsW); maxShift = N - minShift;
            if maxShift > minShift
                for si = 1:N_SHUFF
                    sh = randi([minShift maxShift]);
                    xR_sh = circshift(xR, sh);
                    cc_all_sh(:,si) = lagged_corr_fast(xL, xR_sh, mask_all, nLag, MIN_SAMPLES);
                    if nnz(mask_inter) >= MIN_SAMPLES
                        cc_inter_sh(:,si) = lagged_corr_fast(xL, xR_sh, mask_inter, nLag, MIN_SAMPLES);
                    end
                end
            end
        end

        % ============================= SAVE ===============================
        if any(isfinite(cc_all))
            xc_all(end+1,:)   = cc_all';   %#ok<SAGROW>
            xc_inter(end+1,:) = cc_inter'; %#ok<SAGROW>
            if DO_SHUFFLE
                xc_all_shuf(end+1,:,:)   = reshape(cc_all_sh,   [1 nL N_SHUFF]); %#ok<SAGROW>
                xc_inter_shuf(end+1,:,:) = reshape(cc_inter_sh, [1 nL N_SHUFF]); %#ok<SAGROW>
            end
            nsess = nsess + 1;
            nAllSampTot   = nAllSampTot   + nnz(mask_all);
            nInterSampTot = nInterSampTot + nnz(mask_inter);
            fprintf('  %s n%d: all=%d, intervening=%d samples\n', ...
                    animal, kk, nnz(mask_all), nnz(mask_inter));
        end
    end
end

fprintf('\n%d sessions. samples: all=%d (%.1fs), intervening=%d (%.1fs)\n', ...
        nsess, nAllSampTot, nAllSampTot/fpsW, nInterSampTot, nInterSampTot/fpsW);

if isempty(xc_all), error('No sessions contributed.'); end

% ============================== FIGURE ==============================

colA = [0.20 0.20 0.20];   % all  : black
colI = [0.85 0.10 0.10];   % intervening : red

fig = figure('Color','w','Position',[160 200 720 560]);
ax = axes(fig); hold(ax,'on'); grid(ax,'on');

if DO_SHUFFLE
    plot_shuffle_band(ax, lags_t, xc_all_shuf,   [0.80 0.80 0.80]);
    plot_shuffle_band(ax, lags_t, xc_inter_shuf, [0.95 0.80 0.80]);
end

[hA, pkA, lagA, nA] = plot_xc(ax, lags_t, xc_all,   colA);
[hI, pkI, lagI, nI] = plot_xc(ax, lags_t, xc_inter, colI);

xline(ax, 0, 'k--');
xlim(ax, [lags_t(1) lags_t(end)]);
xlabel(ax, 'lag (s),  corr[ L(t+lag), R(t) ]   ( + : R leads L )');
ylabel(ax, 'L-R correlation');
legend(ax, [hA hI], ...
    {sprintf('all whisking (peak r=%.3f @ %.0f ms, n=%d)',   pkA, 1000*lagA, nA), ...
     sprintf('intervening (peak r=%.3f @ %.0f ms, n=%d)',    pkI, 1000*lagI, nI)}, ...
    'Location','southoutside', 'FontSize',9);
title(ax, sprintf(['L-R whisker cross-correlation: all vs intervening\n' ...
    'intervening = whisk samples > %.0f ms from any inspiration; mean \\pm SE across sessions'], ...
    1000*INSP_EXCLUDE_WIN_S), 'Interpreter','tex');
box(ax,'off');

if doSave
    exportgraphics(fig, fullfile(outDir,'whisk_LR_xcorr_all_vs_intervening.png'), ...
                   'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'whisk_LR_xcorr_all_vs_intervening.pdf'), ...
                   'ContentType','vector', 'BackgroundColor','white');
    save(fullfile(outDir,'whisk_LR_xcorr_all_vs_intervening_data.mat'), ...
         'lags_t','xc_all','xc_inter','nsess','nAllSampTot','nInterSampTot', ...
         'XC_MAXLAG','INSP_EXCLUDE_WIN_S','WHISK_BP','BREATH_BP','ENV_THR', ...
         'MIN_DUR','MERGE_GAP','MIN_INSP_DIST_S','PROM_MAD','fpsW','EXCLUDE');
    fprintf('Saved figure + .mat to %s\n', outDir);
end

% ============================= HELPERS =============================

function [h, pkv, pklag, n] = plot_xc(ax, lags_t, mat, col)
    n  = sum(any(isfinite(mat),2));
    mu = mean(mat, 1, 'omitnan');
    se = std(mat, 0, 1, 'omitnan') ./ sqrt(max(1,sum(isfinite(mat),1)));
    fill(ax, [lags_t; flipud(lags_t)], [(mu+se)'; flipud((mu-se)')], ...
         col, 'FaceAlpha',0.22, 'EdgeColor','none', 'HandleVisibility','off');
    h = plot(ax, lags_t, mu, '-', 'Color', col, 'LineWidth', 2);
    tmp = abs(mu); tmp(~isfinite(tmp)) = -Inf;
    [~,pki] = max(tmp);
    pkv = mu(pki); pklag = lags_t(pki);
end

function plot_shuffle_band(ax, lags_t, shuf, col)
    if isempty(shuf), return; end
    mu_sh = squeeze(mean(shuf, 1, 'omitnan'));     % nL x nShuff
    lo = prctile(mu_sh, 2.5, 2); hi = prctile(mu_sh, 97.5, 2);
    fill(ax, [lags_t; flipud(lags_t)], [hi(:); flipud(lo(:))], col, ...
         'FaceAlpha',0.5, 'EdgeColor','none', 'HandleVisibility','off');
end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a = logical(a(:));
    d = diff([false; a; false]);
    s = find(d == 1); e = find(d == -1) - 1;
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

function cc = lagged_corr_fast(x, y, mask, nLag, minSamp)
% Fast Pearson corr( x(t+lag), y(t) ) for lag = -nLag:nLag.
% lag > 0 means y/ref leads x/signal.
    x = x(:); y = y(:); mask = logical(mask(:));
    N = numel(x);
    lags = -nLag:nLag;
    cc = nan(numel(lags),1);
    for li = 1:numel(lags)
        k = lags(li);
        if k >= 0
            ix = (1+k):N; iy = 1:(N-k);
        else
            ix = 1:(N+k); iy = (1-k):N;
        end
        v = mask(iy) & isfinite(x(ix)) & isfinite(y(iy));
        if nnz(v) < minSamp, continue; end
        xx = x(ix(v)); yy = y(iy(v));
        xx = xx - mean(xx); yy = yy - mean(yy);
        den = sqrt(sum(xx.^2) * sum(yy.^2));
        if den > 0, cc(li) = sum(xx .* yy) / den; end
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
