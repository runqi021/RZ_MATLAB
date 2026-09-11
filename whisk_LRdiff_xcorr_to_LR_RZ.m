% whisk_LRdiff_xcorr_to_LR_RZ.m
%
% Standalone population cross-correlation of the bilateral difference (L-R)
% against each individual whisker, TWO panels:
%   (1) corr[ (L-R)(t+lag), L(t) ]
%   (2) corr[ (L-R)(t+lag), R(t) ]
%
% Whisker: La/Ra -> fill -> demean -> BP 5-50 Hz (USE_BP) ; L-R = xL - xR.
% Samples taken inside the ABSOLUTE 5-deg-on-L whisk-epoch mask.
% One cross-correlation per session; plotted mean +/- SE across sessions.
%
% Lag convention (as in whisk_LR_integral_xcorr_filtered_RZ.m):
%   corr[ x(t+lag), y(t) ];  lag > 0 means y (the reference) leads x.

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW     = 400;
BP       = [5 50];   % Hz
USE_BP   = true;     % true: BP signal; false: raw demeaned angle

% ABSOLUTE 5-deg-on-L whisk-epoch mask
ENV_THR    = 5;      % deg
ENV_SMOOTH = 0.05;   % s
MIN_DUR    = 1;      % s
MERGE_GAP  = 0.2;    % s

XC_MAXLAG   = 0.20;  % s, cross-correlation lag range
MIN_SAMPLES = 50;

EXCLUDE = "5840027";
doSave  = false;

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_LRdiff_xcorr_to_LR');
if doSave && ~isfolder(outDir), mkdir(outDir); end

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
nLag   = round(XC_MAXLAG*fpsW);
lags_t = (-nLag:nLag)'/fpsW;
nL     = numel(lags_t);

xc_L = zeros(0,nL);   % (L-R) vs L
xc_R = zeros(0,nL);   % (L-R) vs R
nsess = 0;

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

% ============================ SESSION LOOP ============================

for e = 1:numel(sess)

    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE)), continue; end

    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d',animal,kRun)), 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)), -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    t  = (0:numel(La)-1)'/fpsW; N = numel(t);

    LaRaw = fillmissing(La-mean(La,'omitnan'),'linear','EndValues','nearest');
    RaRaw = fillmissing(Ra-mean(Ra,'omitnan'),'linear','EndValues','nearest');
    xLbp = filtfilt(bw,aw, LaRaw);
    xRbp = filtfilt(bw,aw, RaRaw);
    if USE_BP, xL = xLbp; xR = xRbp; else, xL = LaRaw; xR = RaRaw; end

    % ABSOLUTE 5-deg-on-L mask (from BP envelope)
    env = movmean(abs(hilbert(xLbp)), max(1,round(ENV_SMOOTH*fpsW)));
    ep  = detect_abs(env, t, ENV_THR, MIN_DUR, MERGE_GAP);
    m = false(N,1);
    for q=1:size(ep,1), m = m | (t>=ep(q,1) & t<=ep(q,2)); end
    if nnz(m) < MIN_SAMPLES, continue; end

    d = xL - xR;
    cc_L = lagged_corr_fast(d, xL, m, nLag, MIN_SAMPLES);
    cc_R = lagged_corr_fast(d, xR, m, nLag, MIN_SAMPLES);

    if any(isfinite(cc_L))
        xc_L(end+1,:) = cc_L'; %#ok<SAGROW>
        xc_R(end+1,:) = cc_R'; %#ok<SAGROW>
        nsess = nsess + 1;
        fprintf('  %s n%d: %d masked samples\n', animal, kRun, nnz(m));
    end
end

fprintf('\n%d sessions contributed\n', nsess);
if isempty(xc_L), error('No sessions contributed.'); end

% ============================== FIGURE ==============================

col = [0.2 0.2 0.2];
fig = figure('Color','w','Position',[120 220 1020 440]);

ax1 = subplot(1,2,1); hold(ax1,'on'); grid(ax1,'on');
plot_xcorr_panel(ax1, xc_L, lags_t, col, ...
    '(L-R) vs L', 'lag (s),  corr[ (L-R)(t+lag), L(t) ]', 'corr');

ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');
plot_xcorr_panel(ax2, xc_R, lags_t, col, ...
    '(L-R) vs R', 'lag (s),  corr[ (L-R)(t+lag), R(t) ]', 'corr');

if USE_BP, s='BP 5-50 Hz'; else, s='raw demeaned'; end
sgtitle(sprintf('L-R cross-correlation to each whisker (%s);  + lag = whisker leads (L-R);  mean \\pm SE, n=%d', ...
        s, nsess), 'Interpreter','tex');

if doSave
    exportgraphics(fig, fullfile(outDir,'whisk_LRdiff_xcorr_to_LR.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'whisk_LRdiff_xcorr_to_LR.pdf'),'ContentType','vector','BackgroundColor','white');
    save(fullfile(outDir,'whisk_LRdiff_xcorr_to_LR_data.mat'), ...
         'lags_t','xc_L','xc_R','nsess','BP','USE_BP','ENV_THR','XC_MAXLAG','fpsW','EXCLUDE');
    fprintf('Saved figure + .mat to %s\n', outDir);
end

% ============================= HELPERS =============================

function plot_xcorr_panel(ax, mat, lags_t, col, ttl, xlab, ylab)
    n  = size(mat,1);
    mu = mean(mat,1,'omitnan');
    se = std(mat,0,1,'omitnan') ./ sqrt(max(1,sum(isfinite(mat),1)));
    fill(ax, [lags_t; flipud(lags_t)], [(mu+se)'; flipud((mu-se)')], ...
         col, 'FaceAlpha',0.25, 'EdgeColor','none');
    plot(ax, lags_t, mu, 'Color',col, 'LineWidth',2);
    xline(ax, 0, 'k--');
    tmp = abs(mu); tmp(~isfinite(tmp)) = -Inf; [~,pki] = max(tmp);
    title(ax, sprintf('%s\npeak r=%.3f @ %.0f ms, n=%d', ttl, mu(pki), 1000*lags_t(pki), n), ...
          'Interpreter','none');
    xlabel(ax, xlab); ylabel(ax, ylab);
    xlim(ax, [lags_t(1) lags_t(end)]); axis(ax,'square'); box(ax,'off');
end

function cc = lagged_corr_fast(x, y, mask, nLag, minSamp)
% Pearson corr( x(t+lag), y(t) ) for lag = -nLag:nLag. lag>0 means y leads x.
    x = x(:); y = y(:); mask = logical(mask(:));
    N = numel(x); lags = -nLag:nLag; cc = nan(numel(lags),1);
    for li = 1:numel(lags)
        k = lags(li);
        if k >= 0, ix=(1+k):N; iy=1:(N-k); else, ix=1:(N+k); iy=(1-k):N; end
        v = mask(iy) & isfinite(x(ix)) & isfinite(y(iy));
        if nnz(v) < minSamp, continue; end
        xx = x(ix(v)); yy = y(iy(v));
        xx = xx-mean(xx); yy = yy-mean(yy);
        den = sqrt(sum(xx.^2)*sum(yy.^2));
        if den > 0, cc(li) = sum(xx.*yy)/den; end
    end
end

function ep = detect_abs(env, t, thrDeg, minDur, mergeGap)
    a = env(:) > thrDeg;
    d = diff([0; a; 0]); s = find(d==1); e = find(d==-1)-1; ep = [t(s) t(e)];
    if ~isempty(ep)
        m = ep(1,:);
        for i=2:size(ep,1)
            if ep(i,1)-m(end,2) <= mergeGap, m(end,2)=ep(i,2); else, m(end+1,:)=ep(i,:); end %#ok<AGROW>
        end
        ep = m(m(:,2)-m(:,1) >= minDur, :);
    end
end

function S = list_sessions(dirPath)
    d = dir(fullfile(char(dirPath), '*_whisk_n*DLC*.csv')); S = {}; key = {};
    for i=1:numel(d)
        tok = regexp(d(i).name, '^(\d+)_whisk_n(\d+)', 'tokens', 'once');
        if isempty(tok), continue; end
        k = sprintf('%s_%s', tok{1}, tok{2});
        if any(strcmp(key,k)), continue; end
        key{end+1}=k; S{end+1}={tok{1}, str2double(tok{2})}; %#ok<AGROW>
    end
    if ~isempty(S)
        an = cellfun(@(c) str2double(c{1}), S); rn = cellfun(@(c) c{2}, S);
        [~,o]=sortrows([an(:) rn(:)]); S=S(o);
    end
end

function csv = pick_csv(dirPath, prefix)
    d = dir(fullfile(char(dirPath), [char(prefix) '*DLC*.csv']));
    assert(~isempty(d), 'no DLC csv matching %s* in %s', prefix, dirPath);
    bn = arrayfun(@(x) bestnum(x.name), d); [~,ix]=max(bn);
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
