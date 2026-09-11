% whisk_LR_signed_movsum_distribution_RZ.m
%
% Standalone: distribution of the SIGNED, moving-sum-integrated bilateral
% whisker difference   int(L - R)   (no abs), pooled across all sessions.
%
%   diff_inst = xL - xR                                   (deg)
%   int_LR    = movsum(diff_inst, winSamp) / fpsW         (deg*s)
%   winSamp   = round(INT_WIN_S * fpsW)
%
% Whisker: La/Ra -> fill NaNs -> demean -> BP 5-50 Hz, zero-phase filtfilt.
% Samples taken inside the whisk-epoch mask (Hilbert envelope of L > ENV_THR);
% set ENV_THR = 0 to use every sample.
%
% NOTE: at fpsW = 400, INT_WIN_S = 0.002 s -> winSamp = round(0.8) = 1 sample,
% i.e. essentially no integration (just the per-sample signed L-R scaled by
% 1/fpsW). Raise INT_WIN_S for a real integration window.

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP = [5 50];   % Hz
USE_BP   = true;     % true: (L-R) of BP angle; false: (L-R) of raw demeaned angle

INT_WIN_S = 0.02;   % moving-sum integration window (s)

ENV_THR    = 3;      % deg, whisk-epoch gate (0 = all samples)
ENV_SMOOTH = 0.05;   % s
MIN_DUR    = 1;      % s
MERGE_GAP  = 0.5;    % s

NBINS   = 100;
XLIM    = 0;         % symmetric x-limit (deg*s); 0 = auto (99.5th pct of |value|)

EXCLUDE = "5840027";
doSave  = false;

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_LR_signed_movsum_distribution');
if doSave && ~isfolder(outDir), mkdir(outDir); end

ord = 4;
[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2), 'bandpass');
winSamp = max(1, round(INT_WIN_S * fpsW));

intLR = [];          % pooled signed int(L-R) samples (deg*s)
nSess = 0;

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

% ============================ SESSION LOOP ============================

for e = 1:numel(sess)

    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE)), continue; end

    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kRun)), 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    tW = (0:numel(La)-1)'/fpsW; N = numel(tW);

    La0 = fillmissing(La(:),'linear','EndValues','nearest'); La0 = La0 - mean(La0,'omitnan');
    Ra0 = fillmissing(Ra(:),'linear','EndValues','nearest'); Ra0 = Ra0 - mean(Ra0,'omitnan');

    xLbp = filtfilt(b_wbp,a_wbp,La0);
    xRbp = filtfilt(b_wbp,a_wbp,Ra0);
    if USE_BP, xL = xLbp; xR = xRbp; else, xL = La0; xR = Ra0; end

    % whisk-epoch mask (from L BP envelope)
    if ENV_THR > 0
        env = movmean(abs(hilbert(xLbp)), max(1,round(ENV_SMOOTH*fpsW)));
        ep = bool_to_epochs(env > ENV_THR, tW, MIN_DUR, MERGE_GAP);
        m = false(N,1);
        for q=1:size(ep,1), m = m | (tW>=ep(q,1)&tW<=ep(q,2)); end
    else
        m = true(N,1);
    end
    if nnz(m) < 50, continue; end

    % SIGNED moving-sum integral of (L - R)
    d = movsum(xL - xR, winSamp, 'omitnan') / fpsW;     % deg*s
    d = d(m & isfinite(d));
    intLR = [intLR; d(:)]; %#ok<AGROW>
    nSess = nSess + 1;
    fprintf('  %s n%d: %d samples, median int(L-R)=%.4g, mean=%.4g deg*s\n', ...
        animal, kRun, numel(d), median(d), mean(d));
end

fprintf('\n%d sessions, %d pooled samples. int(L-R): median=%.4g, mean=%.4g, std=%.4g deg*s (win=%d samp=%.4g s)\n', ...
    nSess, numel(intLR), median(intLR), mean(intLR), std(intLR), winSamp, winSamp/fpsW);

% ============================== FIGURE ==============================

if XLIM <= 0, XLIM = prctile(abs(intLR), 100); end
edges = linspace(-XLIM, XLIM, NBINS+1);

fig = figure('Color','w','Position',[200 220 660 470]); hold on; grid on;
histogram(intLR, edges, 'Normalization','pdf', ...
          'FaceColor',[0.25 0.45 0.75], 'EdgeColor','w');
med = median(intLR); mn = mean(intLR);
xline(0,   'k-',  'LineWidth',1.0);
xline(med, 'r-',  sprintf('median %.3g', med), 'LineWidth',1.5);
xline(mn,  'r--', sprintf('mean %.3g',   mn),  'LineWidth',1.2, 'LabelHorizontalAlignment','left');
xlim([-XLIM XLIM]);
if USE_BP, sig = 'BP 5-50 Hz'; else, sig = 'raw demeaned'; end
xlabel('signed int(L - R)   (deg\cdots)');
ylabel('probability density');
title(sprintf('Distribution of signed int(L-R)  (%s, win=%.3g s = %d samp)\n%d samples, %d sessions%s', ...
    sig, winSamp/fpsW, winSamp, numel(intLR), nSess, ...
    ternary(ENV_THR>0, sprintf(', whisk epochs >%.0f deg', ENV_THR), ', all samples')), ...
    'Interpreter','tex');
box off;

if doSave
    exportgraphics(fig, fullfile(outDir,'whisk_LR_signed_movsum_distribution.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'whisk_LR_signed_movsum_distribution.pdf'),'ContentType','vector','BackgroundColor','white');
    save(fullfile(outDir,'whisk_LR_signed_movsum_distribution_data.mat'), ...
         'intLR','nSess','WHISK_BP','USE_BP','INT_WIN_S','winSamp','ENV_THR', ...
         'MIN_DUR','MERGE_GAP','fpsW','EXCLUDE');
    fprintf('Saved figure + .mat to %s\n', outDir);
end

% ============================= HELPERS =============================

function out = ternary(c,a,b), if c, out=a; else, out=b; end, end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a = logical(a(:)); d = diff([false; a; false]);
    s = find(d==1); e = find(d==-1)-1; ep = [t(s), t(e)];
    if isempty(ep), return; end
    ep2 = ep(1,:);
    for i=2:size(ep,1)
        if ep(i,1)-ep2(end,2)<=mergeGap, ep2(end,2)=ep(i,2);
        else, ep2(end+1,:)=ep(i,:); end %#ok<AGROW>
    end
    ep = ep2; ep = ep(ep(:,2)-ep(:,1)>=minDur,:);
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
    [~,ix] = max(arrayfun(@(x) bestnum(x.name), d));
    csv = fullfile(d(ix).folder, d(ix).name);
end

function n = bestnum(name)
    tok = regexp(name, 'best-(\d+)', 'tokens');
    if isempty(tok), n = 0; else, n = str2double(tok{1}{1}); end
end
