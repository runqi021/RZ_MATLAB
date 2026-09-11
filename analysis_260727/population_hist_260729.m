function population_hist_260729()
%% population_hist_260729  Population Ca-event histogram, breath-triggered.
% -----------------------------------------------------------------------
% Pools every ACTIVE cell across every dataset into one histogram per genotype,
% triggered on inspiration ONSET (left) and the inspiratory PEAK (right).
%
%   rows    = genotype, plus an ALL row at the bottom
%   left    = Ca events vs time from inspiration onset
%   right   = Ca events vs time from the inspiratory peak
%
% Three things per panel, same grammar as the per-cell figure:
%   orange bars = observed Ca events
%   grey line   = shuffle mean
%   grey band   = 95% of shuffles
%
% ------------------------- FIXED 3 s WINDOW -------------------------------
% +/-1.5 s, the SAME for every dataset. Deliberately not scaled to each animal's
% breath period: the median IBI is 2.13 s in 260721_Sert but 0.57 s in the older
% ventral recordings, so an IBI-scaled window would put different amounts of
% physiology on the same axis and the genotypes could not be compared. A fixed
% window means 0 is the same landmark and 1 s is the same second everywhere.
% The cost is that a fast-breathing animal fits ~2.6 cycles in the window and a
% slow one only ~0.7, so side-peaks appear for some genotypes and not others --
% that is a real difference in breath rate, not an artefact.
%
% ----------------------- ACTIVE CELL CRITERION ----------------------------
% Matches the existing ventral analysis exactly:
%     active ROI  <=>  nnz(spike_train > 0) > 5      (i.e. >= 6 events)
% from spike_trigger_dFF.m:34 and temporal_phase_perROI.m:56. Nothing else is
% included. Note this is DELIBERATELY looser than the >=20 events the per-cell PETH
% test requires: pooling means a 6-event cell still contributes 6 real events,
% whereas a per-cell shuffle test on 6 events has no power at all.
%
% Output: <ARCHIVE>\analysis_260729\population_hist.png / .pdf / .csv
%
% Runqi Zhang / 2026-07-29
close all;

here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here,'coh_ca_breath')); addpath(fileparts(here));

%% ===================== USER-EDITABLE PARAMETERS ======================
OUTROOT = 'D:\Ventral_surface_summary\analysis_260729';

% {root, genotype-or-empty}. Empty genotype = infer from the first path component
% under the root, which is how the archive is organised.
SRC = { ...
    'D:\Ventral_surface_summary',  ''
    'D:\260721_Sert_soma_G8s\phys\baseline',           'Sert'
    'D:\260728_vglut2_soma-g8s\phys',                  'Vglut2'
    };

win_s        = 1.5;     % FIXED window half-width -> 3 s total, same for all datasets
binWidth_s   = 0.050;
smoothWidth_s= 0.150;
activeMinEv  = 5;       % ventral's criterion: ACTIVE = nnz(spike_train>0) > this
ampFrac      = 0.20;    % trigger QC, amplitude only (long pauses never rejected)
nShuffle     = 500;     % per ROI; the pooled band is far tighter than any single one
nDrop        = 30;      % breath frames tossed up front, matching the rest of the pipeline
rngSeed      = 260729;
doSave       = true;
% =====================================================================

EXCL = coh_cfg_260727().excludeRecordings;   % recordings on a different stage zero
rng(rngSeed);
edges = -win_s : binWidth_s : win_s;
ctrs  = edges(1:end-1) + binWidth_s/2;
nB    = numel(ctrs);
smB   = max(1, round(smoothWidth_s/binWidth_s));

fprintf('\n============ population_hist_260729 ============\n');
fprintf('window +/-%.2f s (%.1f s total), %.0f ms bins (%d), active = >%d events\n', ...
        win_s, 2*win_s, 1000*binWidth_s, nB, activeMinEv);

%% ---- gather every recording ----
recs = struct('path',{},'genotype',{});
for s = 1:size(SRC,1)
    root = SRC{s,1};
    if ~isfolder(root), fprintf('  (missing, skipped) %s\n', root); continue; end
    hits = dir(fullfile(root,'**','ca_spike_data.mat'));
    for h = 1:numel(hits)
        g = SRC{s,2};
        if isempty(g)
            rel = erase(hits(h).folder, [root filesep]);
            parts = split(string(rel), filesep);
            g = char(parts(1));
        end
        if startsWith(g,'_') || startsWith(g,'analysis') || contains(g,'test'), continue; end
        [~,rn] = fileparts(hits(h).folder);
        if any(strcmp(rn, EXCL)), continue; end   % excluded: different zero reference
        recs(end+1) = struct('path',hits(h).folder,'genotype',string(g)); %#ok<AGROW>
    end
end
assert(~isempty(recs), 'No ca_spike_data.mat found in any source.');
fprintf('found %d recordings with spikes across %d genotypes\n', ...
        numel(recs), numel(unique([recs.genotype])));

%% ---- accumulate, per genotype x trigger ----
gens = unique([recs.genotype]);
A = struct();                       % A.(gen).(trig) = accumulators
for g = gens(:)'
    for tg = {'onset','peak'}
        A.(matlab.lang.makeValidName(g)).(tg{1}) = ...
            struct('cnt',zeros(1,nB),'shMu',zeros(1,nB),'shVar',zeros(1,nB), ...
                   'nRoi',0,'nEv',0,'nRec',0,'nTrig',0,'ibi',[]);
    end
end
nSkipNoTrig = 0; nSkipNoActive = 0;

for r = 1:numel(recs)
    p = recs(r).path;  gf = matlab.lang.makeValidName(recs(r).genotype);
    bp = fullfile(p,'breath_peak_pc1.mat');  ip = fullfile(p,'breath_insp_start_pc1.mat');
    if ~isfile(bp) || ~isfile(ip), nSkipNoTrig = nSkipNoTrig + 1; continue; end
    try
        CA = load(fullfile(p,'ca_spike_data.mat'),'roi_spikes');
        ev = arrayfun(@(x) nnz(x.spike_train>0), CA.roi_spikes);
        act = find(ev > activeMinEv);                       % <-- ventral's criterion
        if isempty(act), nSkipNoActive = nSkipNoActive + 1; continue; end

        fs  = detect_session_fps(p, 30);
        nCa = numel(CA.roi_spikes(1).spike_train);
        BP  = load(bp);  IP = load(ip);
        bw  = detrend(double(BP.breath(:)));
        bw(1:min(nDrop,numel(bw))) = [];
        pk  = round(BP.insp_onset_idx(:)) - nDrop;
        ft  = round(IP.insp_start_idx(:)) - nDrop;
        T   = min(numel(bw), nCa);
        pk  = pk(pk>=1 & pk<=T);  ft = ft(ft>=1 & ft<=T);
        if numel(ft) < 3, nSkipNoTrig = nSkipNoTrig + 1; continue; end
        bwz = (bw(1:T) - median(bw(1:T))) / max(mad(bw(1:T),1)*1.4826, eps);

        % trigger QC + boundary-safe, both landmarks
        ft = sort(ft);  amp = nan(numel(ft)-1,1);  pkf = nan(numel(ft)-1,1);
        for i = 1:numel(ft)-1
            q = pk(pk>ft(i) & pk<ft(i+1));
            if ~isempty(q), amp(i) = bwz(q(1)) - bwz(ft(i)); pkf(i) = q(1); end
        end
        good = amp > ampFrac*median(amp,'omitnan');
        lo = ceil(edges(1)*fs); hi = floor(edges(end)*fs);
        TRG = struct('onset', ft(good), 'peak', pkf(good));
        ibi = median(diff(ft))/fs;

        for tg = {'onset','peak'}
            f = TRG.(tg{1});  f = f(~isnan(f));
            f = f(f+lo>=1 & f+hi<=T);
            if numel(f) < 3, continue; end
            trig = zeros(T,1); trig(f) = 1;
            m  = (lo:hi)';  bo = discretize(m/fs, edges);
            kp = ~isnan(bo);  m = m(kp);  bo = bo(kp);
            lagIdx = mod(m,T)+1;
            acc = A.(gf).(tg{1});
            for a = act(:)'
                e = double(CA.roi_spikes(a).spike_train(:));
                if numel(e)<T, e(end+1:T,1)=0; end
                e = e(1:T);
                ccf = real(ifft(conj(fft(trig)).*fft(e)));   % boundary-safe -> linear
                acc.cnt = acc.cnt + accumarray(bo, ccf(lagIdx), [nB 1])';
                SH = zeros(nShuffle, nB);
                for s2 = 1:nShuffle
                    d = randi(T)-1;
                    SH(s2,:) = accumarray(bo, ccf(mod(lagIdx-1+d,T)+1), [nB 1])';
                end
                SH = movmean(SH, smB, 2);
                acc.shMu  = acc.shMu  + mean(SH,1);
                acc.shVar = acc.shVar + var(SH,0,1);         % independent -> variances add
                acc.nRoi  = acc.nRoi + 1;
                acc.nEv   = acc.nEv + ev(a);
            end
            acc.nRec = acc.nRec + 1;  acc.nTrig = acc.nTrig + numel(f);
            acc.ibi  = [acc.ibi; ibi];
            A.(gf).(tg{1}) = acc;
        end
    catch ME
        fprintf(2,'  ERROR %s: %s\n', p, ME.message);
    end
end
fprintf('skipped: %d without breath triggers, %d with no active ROI\n', nSkipNoTrig, nSkipNoActive);

%% ---- figure ----
gf = arrayfun(@(g) matlab.lang.makeValidName(g), gens);
keep = arrayfun(@(f) A.(f).onset.nRoi > 0, gf);
gens = gens(keep);  gf = gf(keep);
nRow = numel(gens) + 1;                       % + ALL row

fig = figure('Color','w','Name','population Ca-event histogram', ...
             'Units','centimeters','Position',[1 1 22 3.4*nRow + 2.2]);
set(fig,'DefaultAxesFontSize',8);
col = [0.90 0.45 0.10];
yTop = 0.92; yBot = 0.09; hgt = (yTop-yBot)/nRow;
rows = {};

for r = 1:nRow
    for t = 1:2
        tg = ternary(t==1,'onset','peak');
        ax = axes('Parent',fig,'Position',[0.11+0.46*(t-1), yTop-r*hgt+0.035, 0.37, hgt*0.66]); %#ok<LAXES>
        hold(ax,'on'); box(ax,'on');
        if r <= numel(gens)
            a = A.(gf(r)).(tg);  ttl = char(gens(r));
        else
            a = struct('cnt',zeros(1,nB),'shMu',zeros(1,nB),'shVar',zeros(1,nB), ...
                       'nRoi',0,'nEv',0,'nRec',0,'nTrig',0,'ibi',[]);
            for q = 1:numel(gf)
                b = A.(gf(q)).(tg);
                a.cnt=a.cnt+b.cnt; a.shMu=a.shMu+b.shMu; a.shVar=a.shVar+b.shVar;
                a.nRoi=a.nRoi+b.nRoi; a.nEv=a.nEv+b.nEv; a.nRec=a.nRec+b.nRec;
                a.ibi=[a.ibi;b.ibi];
            end
            ttl = 'ALL';
        end
        if a.nRoi == 0, axis(ax,'off'); continue; end
        sd = sqrt(a.shVar);
        fill(ax,[ctrs fliplr(ctrs)],[a.shMu+1.96*sd fliplr(max(a.shMu-1.96*sd,0))], ...
             [.80 .80 .80],'EdgeColor','none');
        plot(ax, ctrs, a.shMu, '-','Color',[.40 .40 .40],'LineWidth',1);
        bar(ax, ctrs, movmean(a.cnt,smB), 1, 'FaceColor',col,'EdgeColor','none');
        xline(ax, 0, 'k-','LineWidth',1.2);
        mib = median(a.ibi);
        if isfinite(mib) && mib < win_s
            for k = [-2 -1 1 2]
                if abs(k*mib) < win_s, xline(ax, k*mib, ':','Color',[.35 .35 .35]); end
            end
        end
        xlim(ax,[edges(1) edges(end)]);
        if r == nRow, xlabel(ax, sprintf('time from %s (s)', tg)); else, set(ax,'XTickLabel',[]); end
        if t == 1, ylabel(ax,'# Ca events'); end
        title(ax, sprintf('%s  |  %s  |  %d active ROI, %d rec, %d events, IBI %.2f s', ...
              upper(tg), ttl, a.nRoi, a.nRec, a.nEv, mib), 'Interpreter','none','FontSize',7);
        rows(end+1,:) = {string(ttl), string(tg), a.nRoi, a.nRec, a.nEv, mib, ...
                         max(movmean(a.cnt,smB)), max(a.shMu+1.96*sd)}; %#ok<AGROW>
    end
end
sgtitle({sprintf('population Ca-event histogram   |   fixed %.0f s window   |   active ROI = >%d events', ...
                 2*win_s, activeMinEv), ...
         'orange = observed    grey line = shuffle mean    grey band = 95% of shuffles'}, ...
        'Interpreter','none','FontSize',9);

%% ---- save ----
if ~isfolder(OUTROOT), mkdir(OUTROOT); end
if doSave
    exportgraphics(fig, fullfile(OUTROOT,'population_hist.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(OUTROOT,'population_hist.pdf'),'ContentType','vector','BackgroundColor','white');
    T = cell2table(rows,'VariableNames', ...
        {'group','trigger','n_active_roi','n_recordings','n_events','median_IBI_s', ...
         'peak_observed','peak_of_95_band'});
    writetable(T, fullfile(OUTROOT,'population_hist.csv'));
    fprintf('\n'); disp(T);
    fprintf('Saved population_hist.png/.pdf/.csv to\n  %s\n', OUTROOT);
end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
