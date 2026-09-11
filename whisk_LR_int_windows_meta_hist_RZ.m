% whisk_LR_int_windows_meta_hist_RZ.m
%
% Standalone: ONE-panel meta histogram of int(L-R) = movsum(xL-xR, window),
% overlaid for several integration windows, pooled across ALL sessions.
%
%   int(L-R) = movsum(xL - xR, winSamp)        % plain moving sum, no /fps (deg)
%   windows: 20, 50, 100, 500, 2000 ms
%
% xL/xR = BP 5-50 Hz whisker angle, samples inside the ABSOLUTE 5-deg-on-L
% whisk-epoch mask (ENV_THR = 0 -> all samples).

clear; clc;

% ============================ USER-EDITABLE ============================
dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW     = 400;
BP       = [5 50];   % Hz
USE_BP   = true;

WINS_MS  = [0 5 10 20 50 100 500 1000];   % integration windows (ms); 0 = no integration (1 sample)

ENV_THR    = 5;      % deg, whisk-epoch gate (0 = all samples)
ENV_SMOOTH = 0.05;   % s
MIN_DUR    = 1;      % s
MERGE_GAP  = 0.2;    % s

NBINS   = 400;
XLIM    = 0;         % symmetric x-limit (deg); 0 = auto

EXCLUDE = "";        % "" = all 17 sessions
doSave  = false;
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_LR_int_windows_meta_hist');
if doSave && ~isfolder(outDir), mkdir(outDir); end

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');
nW = numel(WINS_MS);
winSamp = max(1, round(WINS_MS/1000*fpsW));
pool = cell(1,nW); for k=1:nW, pool{k} = []; end
nSess = 0;

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
    tW = (0:numel(La)-1)'/fpsW; N = numel(tW);

    LaRaw = fillmissing(La-mean(La,'omitnan'),'linear','EndValues','nearest');
    RaRaw = fillmissing(Ra-mean(Ra,'omitnan'),'linear','EndValues','nearest');
    xLbp = filtfilt(bw,aw,LaRaw); xRbp = filtfilt(bw,aw,RaRaw);
    if USE_BP, xL = xLbp; xR = xRbp; else, xL = LaRaw; xR = RaRaw; end

    if ENV_THR > 0
        env = movmean(abs(hilbert(xLbp)), max(1,round(ENV_SMOOTH*fpsW)));
        ep = bool_to_epochs(env > ENV_THR, tW, MIN_DUR, MERGE_GAP);
        if isempty(ep), continue; end
        epIdx = zeros(size(ep));
        for q=1:size(ep,1)
            epIdx(q,1) = find(tW>=ep(q,1),1,'first');
            epIdx(q,2) = find(tW<=ep(q,2),1,'last');
        end
    else
        epIdx = [1 N];
    end

    d = xL - xR;
    % Moving average computed WITHIN each epoch only; keep just center samples
    % whose FULL window fits inside the epoch, so the average never leaks into
    % the quiet (non-whisk) regions. Epochs shorter than the window contribute 0.
    for q=1:size(epIdx,1)
        seg = d(epIdx(q,1):epIdx(q,2)); L = numel(seg);
        for k=1:nW
            W = winSamp(k); h = floor(W/2);
            if L < W, continue; end
            v = movmean(seg, W, 'omitnan');
            pool{k} = [pool{k}; v((1+h):(L-h))]; %#ok<AGROW>
        end
    end
    nSess = nSess + 1;
end

fprintf('%d sessions contributed.\n', nSess);
for k=1:nW
    if isempty(pool{k})
        fprintf('  win %4d ms (%4d samp): NO samples (all epochs shorter than window)\n', WINS_MS(k), winSamp(k));
    else
        fprintf('  win %4d ms (%4d samp): %7d samples, std=%.2f deg\n', ...
            WINS_MS(k), winSamp(k), numel(pool{k}), std(pool{k}));
    end
end

% ============================== FIGURE ==============================
ne = find(~cellfun(@isempty, pool));      % windows that actually have samples
if isempty(ne), error('No window produced samples (epochs too short).'); end
if XLIM <= 0
    XLIM = max(cellfun(@(p) prctile(abs(p),99), pool(ne)));
end
edges = linspace(-XLIM, XLIM, NBINS+1);
cols  = turbo(nW);

fig = figure('Color','w','Position',[200 200 760 500]); hold on; grid on;
h = []; lab = {};
for k=1:nW
    if isempty(pool{k}), continue; end
    hk = histogram(pool{k}, edges, 'Normalization','pdf', 'DisplayStyle','stairs', ...
                   'EdgeColor',cols(k,:), 'LineWidth',1.8);
    if WINS_MS(k)==0, nm='no avg (L-R)'; else, nm=sprintf('%g ms',WINS_MS(k)); end
    h(end+1) = hk; %#ok<AGROW>
    lab{end+1} = sprintf('%s  (n=%d, std %.2f)', nm, numel(pool{k}), std(pool{k})); %#ok<AGROW>
end
xline(0,'k-','LineWidth',1.0);
xlim([-XLIM XLIM]);
if USE_BP, s='BP 5-50 Hz'; else, s='raw demeaned'; end
xlabel('movmean(L-R, window)   (deg)  -- averaged WITHIN epoch only'); ylabel('probability density');
legend(h, lab, 'Location','northeast');
title(sprintf('Distribution of movmean(L-R) vs window (%s)\n%d sessions%s  (window must fit inside epoch)', ...
    s, nSess, tern(ENV_THR>0, sprintf(', whisk epochs >%.0f deg',ENV_THR), ', all samples')), ...
    'Interpreter','none');
box off;

if doSave
    exportgraphics(fig, fullfile(outDir,'whisk_LR_int_windows_meta_hist.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'whisk_LR_int_windows_meta_hist.pdf'),'ContentType','vector','BackgroundColor','white');
    save(fullfile(outDir,'whisk_LR_int_windows_meta_hist_data.mat'),'pool','WINS_MS','winSamp','nSess','BP','USE_BP','ENV_THR','fpsW','EXCLUDE');
    fprintf('Saved figure + .mat to %s\n', outDir);
end

% ============================= HELPERS =============================
function out = tern(c,a,b), if c, out=a; else, out=b; end, end

function ep = bool_to_epochs(a, t, minDur, mergeGap)
    a = logical(a(:)); d = diff([false; a; false]);
    s = find(d==1); e = find(d==-1)-1; ep = [t(s), t(e)];
    if isempty(ep), return; end
    ep2 = ep(1,:);
    for i=2:size(ep,1)
        if ep(i,1)-ep2(end,2)<=mergeGap, ep2(end,2)=ep(i,2); else, ep2(end+1,:)=ep(i,:); end %#ok<AGROW>
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
    tok = regexp(name,'best-(\d+)','tokens'); if isempty(tok), n=0; else, n=str2double(tok{1}{1}); end
end
