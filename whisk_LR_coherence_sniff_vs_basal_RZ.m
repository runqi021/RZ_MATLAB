% whisk_LR_coherence_sniff_vs_basal_RZ.m
%
% Bilateral (L/R) whisker coherence spectrum + delta phase, computed
% SEPARATELY for BASAL and SNIFFING breathing epochs -> TWO figures.
%
% Replicates the Moore 2013 / Deschenes 2016 approach: segment the data by
% breathing state BEFORE coherence, stack fixed-length windows as Chronux
% trials, trial-average coherency.
%
% Whisker:  La/Ra -> fill -> demean -> BP 5-50 Hz filtfilt
% Breathing-state: instantaneous breath rate from Hilbert phase of BP breath
%   basal  : rate < BASAL_HZ
%   sniff  : rate > SNIFF_HZ
% Only WHISKING samples (Hilbert-envelope mask) inside recorded breath are used.
%
% Phase convention: data1=L, data2=R; phi>0 => R leads L.

clear; clc;

% ============================ USER-EDITABLE ============================

dataRoot = "D:\260615_thermalNbasler";
whiskDir = "D:\260615_thermalNbasler\whisk\260615_whisk-RZ-2026-06-17\videos";
noseDir  = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos";

fpsW = 400;

WHISK_BP  = [5 50];   % Hz
BREATH_BP = [1 15];   % Hz

ENV_THR    = 0;       % deg, whisk-epoch gate
ENV_SMOOTH = 0.05;    % s
MIN_DUR    = 1;       % s
MERGE_GAP  = 0.2;     % s

% breathing state (mouse)
BASAL_HZ = 6;
SNIFF_HZ = 7;
FREQ_SMOOTH_S = 0.10;  % smoothing of instantaneous breath rate

% coherence windows (per state) -- non-overlapping, stacked as trials
WIN_BASAL_S = 1.5;    % longer window for slow basal rhythm
WIN_SNIFF_S = 1.0;
TW    =1;            % time-bandwidth product (K = 2*TW-1 tapers)
FPASS = [1 60];
ALPHA = 0.05;
MIN_FRAC_IN_STATE = 0.90;   % window must be >=90% in-state & whisking & breath

EXCLUDE = "5840027";
doSave  = false;

% ======================================================================

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot, fullfile(repoRoot,'thermal_breathing'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

outDir = fullfile(char(dataRoot), 'whisk_LR_coherence_sniff_vs_basal');
if doSave && ~isfolder(outDir), mkdir(outDir); end

ord = 4;
[b_wbp, a_wbp] = butter(ord, WHISK_BP/(fpsW/2),  'bandpass');
[b_br,  a_br ] = butter(3,   BREATH_BP/(fpsW/2), 'bandpass');

winB = round(WIN_BASAL_S*fpsW);
winS = round(WIN_SNIFF_S*fpsW);
assert(winS >= 2*TW && winB >= 2*TW, 'window too short for TW');

pcB = chronux_params(fpsW, TW, FPASS, ALPHA);
pcS = pcB;

% accumulators: stacked windows
Lb = zeros(winB,0); Rb = zeros(winB,0);   % basal
Ls = zeros(winS,0); Rs = zeros(winS,0);   % sniff
nSess = 0; nWb = 0; nWs = 0;

sess = list_sessions(whiskDir);
assert(~isempty(sess), 'no *_whisk_n*.csv in %s', whiskDir);

% ============================ SESSION LOOP ============================

for e = 1:numel(sess)

    animal = sess{e}{1}; kRun = sess{e}{2};
    if strlength(EXCLUDE) > 0 && strcmp(animal, char(EXCLUDE)), continue; end

    % ---------------------------- WHISKER ----------------------------
    try
        M = dlc_gate_interp(pick_csv(whiskDir, sprintf('%s_whisk_n%d', animal, kRun)), 0.6);
    catch ME
        warning('whisk load failed %s n%d: %s', animal, kRun, ME.message); continue;
    end
    La = rad2deg(unwrap(atan2(-(M(:,6)-M(:,3)),  -(M(:,5)-M(:,2)))));
    Ra = rad2deg(unwrap(atan2(-(M(:,12)-M(:,9)),  M(:,11)-M(:,8))));
    tW = (0:numel(La)-1)'/fpsW; N = numel(tW);
    La0 = fillmissing(La(:),'linear','EndValues','nearest'); La0 = La0-mean(La0,'omitnan');
    Ra0 = fillmissing(Ra(:),'linear','EndValues','nearest'); Ra0 = Ra0-mean(Ra0,'omitnan');
    xL = filtfilt(b_wbp,a_wbp,La0); xR = filtfilt(b_wbp,a_wbp,Ra0);

    env = movmean(abs(hilbert(xL)), max(1,round(ENV_SMOOTH*fpsW)));
    ep = bool_to_epochs(env>ENV_THR, tW, MIN_DUR, MERGE_GAP);
    m = false(N,1);
    for q=1:size(ep,1), m = m | (tW>=ep(q,1)&tW<=ep(q,2)); end
    if nnz(m) < 50, continue; end

    % ---------------------------- BREATH -----------------------------
    try
        Pn = thermal_resolve_paths(pick_csv(noseDir, sprintf('%s_nose_n%d', animal, kRun)), dataRoot);
        if ~isfile(Pn.breath), continue; end
        Bs = load(Pn.breath); br = Bs.breath(:); fb = double(Bs.fps); tBr = (0:numel(br)-1)'/fb;
        brw = interp1(tBr, br, tW, 'linear', NaN);
        breath_finite = isfinite(brw);
        if nnz(breath_finite) < 50, continue; end
        brw = fillmissing(brw,'linear','EndValues','nearest');
        brw_bp = filtfilt(b_br,a_br,brw);
    catch ME
        warning('breath failed %s n%d: %s', animal, kRun, ME.message); continue;
    end

    % instantaneous breath rate from Hilbert phase
    ph = unwrap(angle(hilbert(brw_bp)));
    fr = [0; diff(ph)] * fpsW/(2*pi);                 % Hz
    fr = movmean(fr, max(1,round(FREQ_SMOOTH_S*fpsW)));
    basal_mask = m & breath_finite & (fr>0) & (fr<BASAL_HZ);
    sniff_mask = m & breath_finite & (fr>SNIFF_HZ);

    % --------------- tile windows inside each state ------------------
    nb = tile_into(xL, xR, basal_mask, winB, MIN_FRAC_IN_STATE);
    ns = tile_into(xL, xR, sniff_mask, winS, MIN_FRAC_IN_STATE);
    if ~isempty(nb), Lb = [Lb nb.L]; Rb = [Rb nb.R]; nWb = nWb + size(nb.L,2); end %#ok<AGROW>
    if ~isempty(ns), Ls = [Ls ns.L]; Rs = [Rs ns.R]; nWs = nWs + size(ns.L,2); end %#ok<AGROW>

    if (~isempty(nb)) || (~isempty(ns))
        nSess = nSess + 1;
        fprintf('  %s n%d: basal win %d, sniff win %d\n', animal, kRun, ...
            size_or0(nb), size_or0(ns));
    end
end

fprintf('\n%d sessions. basal windows=%d (%.0fs), sniff windows=%d (%.0fs)\n', ...
    nSess, nWb, nWb*WIN_BASAL_S, nWs, nWs*WIN_SNIFF_S);

% ============================ COHERENCE + FIGS ============================

make_fig('BASAL',  Lb, Rb, pcB, WIN_BASAL_S, TW, nWb, FPASS, ALPHA, ...
         fullfile(outDir,'coherence_basal'), doSave);
make_fig('SNIFF',  Ls, Rs, pcS, WIN_SNIFF_S, TW, nWs, FPASS, ALPHA, ...
         fullfile(outDir,'coherence_sniff'), doSave);

if doSave
    save(fullfile(outDir,'coherence_sniff_vs_basal_data.mat'), ...
         'nSess','nWb','nWs','WIN_BASAL_S','WIN_SNIFF_S','TW','FPASS','ALPHA', ...
         'BASAL_HZ','SNIFF_HZ','WHISK_BP','BREATH_BP');
    fprintf('Saved coherence figures + .mat to %s\n', outDir);
end

% ============================= HELPERS =============================

function p = chronux_params(Fs, TW, fpass, alpha)
    p.Fs=Fs; p.tapers=[TW 2*TW-1]; p.pad=0; p.fpass=fpass; p.err=[2 alpha]; p.trialave=1;
end

function out = tile_into(xL, xR, mask, winSamp, minFrac)
% Non-overlapping windows fully (>=minFrac) inside mask; demeaned columns.
    out = [];
    N = numel(xL); s = 1; Lc = []; Rc = [];
    while s+winSamp-1 <= N
        idx = s:(s+winSamp-1);
        if mean(mask(idx)) >= minFrac
            wl = xL(idx); wr = xR(idx);
            if all(isfinite(wl)) && all(isfinite(wr))
                Lc(:,end+1) = wl - mean(wl); %#ok<AGROW>
                Rc(:,end+1) = wr - mean(wr); %#ok<AGROW>
            end
            s = s + winSamp;            % step a full window past in-state block
        else
            s = s + round(winSamp/4);   % slide forward looking for an in-state block
        end
    end
    if ~isempty(Lc), out.L = Lc; out.R = Rc; end
end

function n = size_or0(s)
    if isempty(s), n = 0; else, n = size(s.L,2); end
end

function make_fig(name, Lw, Rw, pc, winSec, TW, nW, FPASS, ALPHA, stem, doSave)
    col = [0.10 0.30 0.85]; colS = [0.6 0.6 0.6];
    fig = figure('Color','w','Position',[160 160 620 720], 'Name',name);
    if size(Lw,2) < 2
        annotation(fig,'textbox',[0.1 0.45 0.8 0.1],'String', ...
            sprintf('%s: only %d windows -- not enough for coherence', name, size(Lw,2)), ...
            'EdgeColor','none','FontSize',12,'HorizontalAlignment','center');
        return;
    end
    [~, C, phi, ~,~,~, f, confC, phistd, Cerr] = coherencyc(Lw, Rw, pc);
    f=f(:); C=C(:); phi=phi(:); phistd=phistd(:);
    Clo=Cerr(1,:)'; Chi=Cerr(2,:)';
    phid=rad2deg(phi); phisd=rad2deg(phistd);
    sig = C>=confC;
    [Cpk,ipk]=max(C);

    axA=subplot(2,1,1); hold(axA,'on'); grid(axA,'on');
    fill(axA,[f;flipud(f)],[Chi;flipud(Clo)],col,'FaceAlpha',0.2,'EdgeColor','none');
    plot(axA,f,C,'-','Color',col,'LineWidth',2);
    yline(axA,confC,'k--',sprintf('conf %.2f',confC),'LabelHorizontalAlignment','left');
    plot(axA,f(ipk),Cpk,'o','MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize',6);
    xlim(axA,FPASS); ylim(axA,[0 1]);
    xlabel(axA,'frequency (Hz)'); ylabel(axA,'coherence |C_{LR}|');
    title(axA,sprintf('%s: L-R coherence  (N=%d win, %.0fs, win=%.1fs, TW=%d)\npeak %.3f @ %.2f Hz', ...
        name, nW, nW*winSec, winSec, TW, Cpk, f(ipk)),'Interpreter','tex');
    box(axA,'off');

    axB=subplot(2,1,2); hold(axB,'on'); grid(axB,'on');
    yline(axB,0,'k-');
    fill(axB,[f;flipud(f)],[phid+1.96*phisd;flipud(phid-1.96*phisd)],colS,'FaceAlpha',0.2,'EdgeColor','none');
    plot(axB,f,phid,'-','Color',colS,'LineWidth',1);
    ps=phid; ps(~sig)=NaN; plot(axB,f,ps,'-','Color',col,'LineWidth',2.2);
    xlim(axB,FPASS); ylim(axB,[-180 180]); yticks(axB,-180:90:180);
    xlabel(axB,'frequency (Hz)'); ylabel(axB,'\Delta phase L,R (deg)');
    title(axB,'\Delta phase ( + R leads L ); bold = significant','Interpreter','tex');
    box(axB,'off');

    if doSave
        exportgraphics(fig,[stem '.png'],'Resolution',200,'BackgroundColor','white');
        exportgraphics(fig,[stem '.pdf'],'ContentType','vector','BackgroundColor','white');
    end
end

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
