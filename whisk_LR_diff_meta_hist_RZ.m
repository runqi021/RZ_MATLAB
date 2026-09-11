% whisk_LR_diff_meta_hist_RZ.m
%
% Standalone: ONE-panel meta histogram of L-R (signed bilateral whisker
% difference, deg), pooled across ALL sessions.
%
% L-R = xL - xR of the BP 5-50 Hz whisker angle, taken inside the
% ABSOLUTE 5-deg-on-L whisk-epoch mask.  Set ENV_THR = 0 to use every sample.

clear; clc;

% ============================ USER-EDITABLE ============================
dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";

fpsW     = 400;
BP       = [5 50];   % Hz
USE_BP   = true;     % true: L-R of BP angle; false: raw demeaned angle

ENV_THR    = 5;      % deg, whisk-epoch gate (0 = all samples)
ENV_SMOOTH = 0.05;   % s
MIN_DUR    = 1;      % s
MERGE_GAP  = 0.2;    % s

NBINS   = 120;
XLIM    = 0;         % symmetric x-limit (deg); 0 = auto (99.5th pct of |L-R|)

EXCLUDE = "";        % "" = all 17 sessions; "5840027" = drop that animal
doSave  = false;
% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));

outDir = fullfile(char(dataRoot), 'whisk_LR_diff_meta_hist');
if doSave && ~isfolder(outDir), mkdir(outDir); end

[bw,aw] = butter(3, BP/(fpsW/2),'bandpass');

diffLR = [];   % pooled L-R samples (deg)
nSess  = 0;

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
        m = false(N,1);
        for q=1:size(ep,1), m = m | (tW>=ep(q,1)&tW<=ep(q,2)); end
    else
        m = true(N,1);
    end
    if nnz(m) < 50, continue; end

    d = xL - xR;
    d = d(m & isfinite(d));
    diffLR = [diffLR; d(:)]; %#ok<AGROW>
    nSess = nSess + 1;
    fprintf('  %s n%d: %d samples, median=%.2f deg\n', animal, kRun, numel(d), median(d));
end

fprintf('\n%d sessions, %d pooled samples. L-R: median=%.2f, mean=%.2f, std=%.2f deg\n', ...
    nSess, numel(diffLR), median(diffLR), mean(diffLR), std(diffLR));

% ============================== FIGURE ==============================
if XLIM <= 0, XLIM = prctile(abs(diffLR), 99.5); end
edges = linspace(-XLIM, XLIM, NBINS+1);

fig = figure('Color','w','Position',[200 220 660 470]); hold on; grid on;
histogram(diffLR, edges, 'Normalization','pdf', 'FaceColor',[0.25 0.45 0.75], 'EdgeColor','none');
xline(0, 'k-', 'LineWidth',1.0);
xline(median(diffLR), 'r-', sprintf('median %.2f', median(diffLR)), 'LineWidth',1.5);
xlim([-XLIM XLIM]);
if USE_BP, s='BP 5-50 Hz'; else, s='raw demeaned'; end
xlabel('L - R  (deg)'); ylabel('probability density');
title(sprintf('Meta distribution of L-R whisker angle (%s)\n%d samples, %d sessions%s', ...
    s, numel(diffLR), nSess, tern(ENV_THR>0, sprintf(', whisk epochs >%.0f deg', ENV_THR), ', all samples')), ...
    'Interpreter','none');
box off;

if doSave
    exportgraphics(fig, fullfile(outDir,'whisk_LR_diff_meta_hist.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'whisk_LR_diff_meta_hist.pdf'),'ContentType','vector','BackgroundColor','white');
    save(fullfile(outDir,'whisk_LR_diff_meta_hist_data.mat'),'diffLR','nSess','BP','USE_BP','ENV_THR','fpsW','EXCLUDE');
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
