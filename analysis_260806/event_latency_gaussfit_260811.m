% event_latency_gaussfit_260811.m
% -----------------------------------------------------------------------
%  TIME-DOMAIN tuning, one session at a time. NO PHASE ANYWHERE.
%
%  Per ROI, per trigger (inspiration ONSET and inspiratory PEAK):
%
%    1. WINDOW = one mean IBI centred on the trigger, and that IBI is THIS
%       RECORDING'S own mean inter-trigger interval -- not a session or archive
%       mean. Breath rate differs between recordings, so a shared window would
%       show one recording more than a cycle and another less.
%    2. every event is assigned to its NEAREST trigger and given a SIGNED
%       latency:  < 0 LEADS the trigger,  > 0 LAGS it. One event, one trigger.
%    3. PERMUTATION TEST FIRST. Nothing is described until the latency
%       distribution is shown to be CONCENTRATED:
%           statistic  IQR of the nearest-trigger latencies (SMALLER = tighter)
%           null       nPerm circular shifts of the BREATH TRIGGER train,
%                      latencies recomputed from scratch against each shift
%           p          (1 + #{IQR_null <= IQR_obs}) / (1 + nPerm)
%       Only cells with p < alpha are described. The rest are reported as
%       "not distinguishable from uniform", which is an answer, not a gap.
%
%       THE STATISTIC IS THE THING WE REPORT. An earlier version tested
%       max(smoothed histogram) - mean on an FFT cross-correlogram and then
%       reported median/IQR of the nearest-trigger latencies -- two different
%       quantities. It passed cells it should not have. Two causes, both real:
%         (a) BIN CLIPPING. ceil(t/binS) was clamped to nb, so the last trigger
%             and the last events piled into bin nb of BOTH trains. Coincident
%             pile-up in both trains contributes pairs at LAG 0 -- exactly the
%             observed window -- inflating the observed statistic alone.
%         (b) ONE LUCKY BIN. With ~20 events over 25 bins the tallest bin is
%             noise, and a max-based statistic scores it as structure.
%       Measured on Vglut2/1124 pFN roi5_1400-1230-0 r4 (21 events, latencies
%       filling 72% of the window): correlogram null p = 0.007, honest
%       trigger-shift null p = 0.282 with the same statistic, p = 0.106 with
%       IQR. The correlogram is gone. Latencies are recomputed against every
%       shifted trigger train, which costs more and cannot drift from what is
%       reported.
%    4. THEN the description: a KERNEL DENSITY ESTIMATE of the latencies plus
%       the MEDIAN and the INTERQUARTILE RANGE. Both are read off the latencies
%       themselves, so nothing is imposed on the data:
%           centre = median(L)          -/+ ms, - leads the trigger, + lags
%           spread = IQR = Q3 - Q1      and Q1, Q3 are reported separately
%       This replaces the Gaussian mu/sigma. A Gaussian returns a centre and a
%       width for ANY histogram, symmetric or not, and the latency
%       distributions here are visibly skewed (a calcium event can lag a long
%       way but cannot lead by more than the window). The median sits where
%       half the events are, the IQR brackets the middle half, and neither
%       assumes symmetry. The KDE is drawn as the smooth curve only -- the
%       numbers do not come from it, so its bandwidth cannot move them.
%
%  WHY SHIFT THE BREATH, NOT THE SPIKES. One trigger train serves every ROI in a
%  recording, so the nPerm shifted trigger trains are built ONCE per recording and
%  every ROI is tested against the same nPerm nulls. Circular shifts preserve the
%  trigger train's own rate and rhythm, so the null keeps everything except the
%  alignment. Shifts within minShiftIBI cycles of 0 (or of the full duration,
%  which wraps to the same place) are rejected: a shift of half a breath is not an
%  independent null, it is the same alignment moved by half a cycle.
%
%  WHY TIME, NOT PHASE. Phase maps inspiration (~15-21% of the cycle) onto half
%  the axis, and a GCaMP lag -- a fixed number of MILLISECONDS -- becomes a phase
%  offset that depends on the breath rate. In time the lag is a constant shift of
%  the centre, identical for every cell. NOTHING HERE IS LAG-COMPENSATED.
%
%  Runqi Zhang / 2026-08-11

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
addpath(fileparts(scriptDir));

%% ===================== USER-EDITABLE =====================
sessionDir = 'D:\Ventral_surface_summary\Vglut2\1124';
outDir     = 'D:\Ventral_surface_summary\breath_trig_heatmap_260806';
nBins      = 25;        % bins across the full window (one IBI), for DISPLAY only
nPerm      = 1000;      % circular-shift permutations
caLagSec   = 0.10;      % GCaMP LEAD COMPENSATION, seconds. ca_spike_data events sit
                        % ON the dF/F peak, so they LAG the true spike by the rise
                        % time; spikes are moved EARLIER by round(caLagSec*fps)
                        % frames. 0.1 s = 3 frames @30 Hz, the value the
                        % svd_breath_motion family settled on. Set 0 to disable --
                        % this shifts every median by the same amount and cannot
                        % create or destroy significance.
fix1124    = true;      % Vglut2/1124 used the breath RISING EDGE, so breath leads
                        % calcium by one frame there. Delay breath 1 frame, exactly
                        % as the svd_breath_motion scripts do. Session-specific.
ctrlSwapBreath = false; % ARTIFACT CONTROL. true = pair every recording with a
                        % DIFFERENT recording's breath train. Real coupling must
                        % vanish; whatever survives is method artifact.
alphaPerm  = 0.005;     % a cell is described only if p < this.  NOTE the floor:
                        % with nPerm shifts the smallest attainable p is
                        % 1/(nPerm+1), so at nPerm = 300 that is 0.00332 and the
                        % next value up is 0.00664. A 0.005 cut therefore admits
                        % ONLY cells where not one null beat the observed IQR.
                        % Raise nPerm to resolve the threshold more finely.
minShiftIBI = 2;        % reject shifts within this many breath cycles of no-shift
nDrop      = 30;        % breath frames tossed to align with calcium
minEvents  = 20;        % minimum in-window events
perPage    = 24;
doSave     = true;
% =========================================================

TRIG = {'onset','peak'};
tCol = [0.85 0.20 0.10; 0.15 0.45 0.85];

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
rng(0);
outTag = sessTag(sessionDir);
if ctrlSwapBreath, outTag = [outTag '_CTRLswapBreath']; end   %% never clobber the real run
hits = dir(fullfile(sessionDir,'**','ca_spike_data.mat'));
assert(~isempty(hits), 'no ca_spike_data under %s', sessionDir);

C = struct('rec',{},'site',{},'roi',{},'ibi',{},'nCyc',{},'binMs',{}, ...
           'y',{},'ctrs',{},'nEv',{},'evPerCyc',{}, ...
           'p',{},'T',{},'med',{},'q1',{},'q3',{},'iqr',{},'kx',{},'kd',{});
for h = 1:numel(hits)
    fp = hits(h).folder;  pp = strsplit(fp, filesep);
    rec = pp{end};  site = pp{end-1};
    % CONTROL: take the breath from a DIFFERENT recording. Everything else -- the
    % events, the window, the test -- is untouched, so any surviving significance
    % is produced by the method rather than by the animal.
    fpB = fp;
    if ctrlSwapBreath, fpB = hits(mod(h, numel(hits)) + 1).folder; end
    bp = fullfile(fpB,'breath_peak_pc1.mat');  ip = fullfile(fpB,'breath_insp_start_pc1.mat');
    if ~isfile(bp) || ~isfile(ip), fprintf('  skip (no breath): %s\n', rec); continue; end
    PKr = load(bp,'insp_onset_idx');  PK = sort(round(PKr.insp_onset_idx(:)) - nDrop);
    ONr = load(ip,'insp_start_idx');  ON = sort(round(ONr.insp_start_idx(:)) - nDrop);
    PK = PK(PK>=1);  ON = ON(ON>=1);
    if numel(PK) < 5 || numel(ON) < 5, continue; end

    fps = detect_session_fps(fp, 30);
    S   = load(fullfile(fp,'ca_spike_data.mat'),'roi_spikes');
    % Vglut2/1124 rising-edge fix: delay the breath by one frame.
    if fix1124 && contains(fp, fullfile('Vglut2','1124')), PK = PK + 1;  ON = ON + 1; end
    nLag = round(caLagSec*fps);                     % GCaMP lead comp, frames
    T   = max([numel(S.roi_spikes(1).spike_train), max(PK), max(ON)]);
    if ctrlSwapBreath                               % wrap a foreign train into range
        PK = unique(mod(PK-1, T) + 1);  ON = unique(mod(ON-1, T) + 1);
    end
    trg = {ON, PK};

    Trec = T/fps;                                   % recording duration, seconds
    for q = 1:2
        tg   = trg{q};
        ibi  = mean(diff(tg))/fps;                  % THIS recording's mean IBI
        binS = ibi/nBins;                           % display bin width, seconds
        half = floor(nBins/2);
        ctrs = (-half:half)*binS*1000;              % ms
        edg  = [ctrs - binS*500, ctrs(end) + binS*500];
        nCyc = numel(tg) - 1;
        tT   = tg(:).'/fps;

        % ---- the nPerm shifted trigger trains, built ONCE for every ROI ----
        % Shifts within minShiftIBI cycles of no-shift are rejected -- a shift of
        % half a breath is the same alignment displaced, not an independent null.
        shMin = minShiftIBI*ibi;
        sh = zeros(nPerm,1); k = 0;
        while k < nPerm
            v = Trec*rand;
            if v > shMin && v < Trec - shMin, k = k+1; sh(k) = v; end
        end
        tShift = arrayfun(@(v) sort(mod(tT + v, Trec)), sh, 'uni', 0);

        for r = 1:numel(S.roi_spikes)
            st = find(S.roi_spikes(r).spike_train > 0) - nLag;   % lead-shift
            st = st(st >= 1);
            if numel(st) < minEvents, continue; end
            tE = st(:)/fps;

            % ---- observed nearest-trigger latencies ----
            L = near_lat(tE, tT, Trec);
            L = 1000*L(abs(L) <= ibi/2);            % ms
            if numel(L) < minEvents, continue; end
            y  = histcounts(L, edg) / nCyc;         % EVENTS PER CYCLE per bin
            Qo = prctile(L,[25 75]);  iqrObs = Qo(2) - Qo(1);

            % ---- permutation: recompute the latencies against each shift ----
            % No correlogram, no binning, no clipping. The statistic tested is the
            % IQR that gets reported, so the test cannot pass a cell whose
            % reported spread is indistinguishable from flat.
            iqrNul = nan(nPerm,1);
            for s = 1:nPerm
                Ls = near_lat(tE, tShift{s}, Trec);
                Ls = 1000*Ls(abs(Ls) <= ibi/2);
                if numel(Ls) >= 4
                    Qn = prctile(Ls,[25 75]);  iqrNul(s) = Qn(2) - Qn(1);
                end
            end
            p = (1 + nnz(iqrNul <= iqrObs)) / (1 + sum(~isnan(iqrNul)));

            c = struct('rec',rec,'site',site,'roi',r,'ibi',ibi,'nCyc',nCyc, ...
                       'binMs',binS*1000,'y',y,'ctrs',ctrs,'nEv',numel(L), ...
                       'evPerCyc',numel(L)/nCyc,'p',p,'T',median(iqrNul,'omitnan'), ...
                       'med',NaN,'q1',NaN,'q3',NaN,'iqr',NaN,'kx',[],'kd',[]);

            % ---- describe ONLY if it passed ----
            if p < alphaPerm
                c.med = median(L);  c.q1 = Qo(1);  c.q3 = Qo(2);  c.iqr = iqrObs;
                % KDE for the CURVE only. Scaled to the histogram: density (per ms)
                % x events x bin width -> the same "events per cycle" the bars use.
                [kd, kx] = ksdensity(L, linspace(ctrs(1), ctrs(end), 300));
                c.kx = kx;  c.kd = kd * numel(L) * c.binMs / nCyc;
            end

            k = find(strcmp({C.rec},rec) & [C.roi]==r, 1);
            if q == 1
                C(end+1) = c;                                       %#ok<SAGROW>
                C(end).p=[p NaN]; C(end).T=[c.T NaN]; C(end).ibi=[ibi NaN];
                C(end).med=[c.med NaN]; C(end).q1=[c.q1 NaN]; C(end).q3=[c.q3 NaN];
                C(end).iqr=[c.iqr NaN]; C(end).kx={c.kx}; C(end).kd={c.kd};
                C(end).y={y}; C(end).ctrs={ctrs};
                C(end).nEv=[c.nEv NaN]; C(end).evPerCyc=[c.evPerCyc NaN];
            elseif ~isempty(k)
                C(k).p(2)=p; C(k).T(2)=c.T; C(k).ibi(2)=ibi;
                C(k).med(2)=c.med; C(k).q1(2)=c.q1; C(k).q3(2)=c.q3; C(k).iqr(2)=c.iqr;
                C(k).kx{2}=c.kx; C(k).kd{2}=c.kd;
                C(k).y{2}=y; C(k).ctrs{2}=ctrs;
                C(k).nEv(2)=c.nEv; C(k).evPerCyc(2)=c.evPerCyc;
            end
        end
        fprintf('%-42s %-5s IBI %.3f s | bin %.0f ms | %d cycles\n', ...
                rec(1:min(42,end)), TRIG{q}, ibi, binS*1000, nCyc);
    end
end
assert(~isempty(C), 'no ROI cleared minEvents = %d', minEvents);

P = vertcat(C.p);  MD = vertcat(C.med);  IQ = vertcat(C.iqr);
Q1 = vertcat(C.q1); Q3 = vertcat(C.q3);
both = all(~isnan(P),2);
fprintf('\n%d ROIs tested | %d permutations, alpha %.3f\n', nnz(both), nPerm, alphaPerm);
for q = 1:2
    sg = both & P(:,q) < alphaPerm & ~isnan(MD(:,q));
    fprintf(['  %-5s significant %3d / %3d  |  median %+.0f..%+.0f (median %+.0f) ms | ' ...
             'IQR %.0f..%.0f (median %.0f) ms | lead %d / lag %d\n'], ...
        TRIG{q}, nnz(sg), nnz(both), min(MD(sg,q)), max(MD(sg,q)), median(MD(sg,q)), ...
        min(IQ(sg,q)), max(IQ(sg,q)), median(IQ(sg,q)), nnz(MD(sg,q)<0), nnz(MD(sg,q)>0));
end

%% ---- per-cell pages ----
sig = find(both & any(P < alphaPerm, 2));
[~,o] = sort(MD(sig,2)); sig = sig(o);
nPage = max(1, ceil(numel(sig)/perPage));
for pg = 1:nPage
    sub = sig((pg-1)*perPage+1 : min(pg*perPage, numel(sig)));
    if isempty(sub), continue; end
    fh = figure('Color','w','Position',[30 30 1500 950]);
    tl = tiledlayout(fh, 4, 6, 'TileSpacing','compact','Padding','compact');
    for k = 1:numel(sub)
        c = C(sub(k)); ax = nexttile(tl); hold(ax,'on');
        yTop = max([cellfun(@max, c.y), eps]) * 1.15;
        for q = 1:2
            stairs(ax, c.ctrs{q}, c.y{q}, 'Color',[tCol(q,:) 0.45], 'LineWidth',0.8);
            if ~isnan(c.med(q))
                % No IQR band on this page -- KDE and median line only. The
                % quartiles are still computed, and live in the summary figure
                % and the .csv.
                plot(ax, c.kx{q}, c.kd{q}, '-', 'Color',tCol(q,:), 'LineWidth',1.4);
                xline(ax, c.med(q), '--', 'Color',tCol(q,:), 'LineWidth',1.0);
            end
        end
        xlim(ax, c.ctrs{2}([1 end])); ylim(ax, [0 yTop]);
        set(ax,'TickDir','out','FontSize',7); box(ax,'off');
        if mod(k-1,6)==0, ylabel(ax,'epc','FontSize',7); end
        if k > numel(sub)-6, xlabel(ax,'latency (ms)','FontSize',7); end
        title(ax, sprintf('%s r%d  n=%d', shortrec(c.rec), c.roi, c.nEv(2)), 'FontSize',7);
        subtitle(ax, sprintf('on %s | pk %s', fmt(c.med(1),c.p(1)), fmt(c.med(2),c.p(2))), ...
                 'FontSize',6.5, 'Color',[0.3 0.3 0.3]);
    end
    title(tl, sprintf(['%s  page %d/%d  |  red ONSET, blue PEAK  |  window = that recording''s ' ...
        'mean IBI, y = epc (events per cycle per bin)  |  curve = KDE, dashed = median  |  ' ...
        '%d-shift permutation, described only if p<%.3f'], ...
        outTag, pg, nPage, nPerm, alphaPerm), 'FontWeight','bold','Interpreter','none');
    set(findall(fh,'Type','axes'),'Toolbar',[]);   %% kill the toolbar exportgraphics warns about
    if doSave
        b = fullfile(outDir, sprintf('event_latency_%s_cells_p%02d', outTag, pg));
        exportgraphics(fh,[b '.png'],'Resolution',200,'BackgroundColor','white');
        exportgraphics(fh,[b '.pdf'],'ContentType','vector','BackgroundColor','white');
    end
end

%% ---- summary ----
fh = figure('Color','w','Position',[60 60 1250 520]);
tl = tiledlayout(fh,1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl,1); hold(ax,'on');
for q=1:2, sg = both & P(:,q)<alphaPerm; histogram(ax, MD(sg,q), 14,'FaceColor',tCol(q,:),'FaceAlpha',0.55); end
xline(ax,0,'k-','LineWidth',1.2); xlabel(ax,'median latency (ms)  [- leads, + lags]'); ylabel(ax,'cells');
box(ax,'off'); set(ax,'TickDir','out'); title(ax,'centre  (median)'); legend(ax,TRIG,'Box','off');
ax = nexttile(tl,2); hold(ax,'on');
for q=1:2, sg = both & P(:,q)<alphaPerm; histogram(ax, IQ(sg,q), 14,'FaceColor',tCol(q,:),'FaceAlpha',0.55); end
xlabel(ax,'IQR (ms)'); ylabel(ax,'cells'); box(ax,'off'); set(ax,'TickDir','out');
title(ax,'spread  (Q3 - Q1)'); legend(ax,TRIG,'Box','off');
ax = nexttile(tl,3); hold(ax,'on');
for q=1:2
    sg = find(both & P(:,q)<alphaPerm);
    ev = arrayfun(@(c) c.evPerCyc(q), C(sg));
    % whisker = the actual quartiles, so the asymmetry the median/IQR captures is visible
    for i = 1:numel(sg)
        plot(ax, [Q1(sg(i),q) Q3(sg(i),q)], IQ(sg(i),q)*[1 1], '-', 'Color',[tCol(q,:) 0.35], 'LineWidth',0.8);
    end
    scatter(ax, MD(sg,q), IQ(sg,q), 20+300*ev/max([ev 1]), tCol(q,:), 'filled', ...
            'MarkerFaceAlpha',0.6,'MarkerEdgeColor','w');
end
xline(ax,0,'k-','LineWidth',1.2); xlabel(ax,'median (ms), whisker = Q1..Q3'); ylabel(ax,'IQR (ms)');
box(ax,'off'); set(ax,'TickDir','out'); title(ax,'centre vs spread  (area \propto events/cycle)');
title(tl, sprintf('%s  |  %d of %d ROIs significant (p<%.3f, %d shifts)  |  median & IQR, KDE curve', ...
      outTag, nnz(both & any(P<alphaPerm,2)), nnz(both), alphaPerm, nPerm), ...
      'FontWeight','bold','Interpreter','none');

set(findall(fh,'Type','axes'),'Toolbar',[]);
if doSave
    b = fullfile(outDir, sprintf('event_latency_%s_summary', outTag));
    exportgraphics(fh,[b '.png'],'Resolution',200,'BackgroundColor','white');
    exportgraphics(fh,[b '.pdf'],'ContentType','vector','BackgroundColor','white');
    writetable(table({C.rec}',{C.site}',[C.roi]', vertcat(C.ibi)*1000, P, MD, Q1, Q3, IQ, ...
        vertcat(C.nEv), vertcat(C.evPerCyc), [C.nCyc]', ...
        'VariableNames',{'recording','site','roi','ibi_ms','p_perm','median_ms', ...
                         'q1_ms','q3_ms','iqr_ms','n_events','ev_per_cycle','n_cycles'}), ...
        fullfile(outDir, sprintf('event_latency_%s.csv', outTag)));
    fprintf('saved event_latency_%s_* to %s\n', outTag, outDir);
end

%% ---- local ----
function L = near_lat(tE, ts, Trec)
% Signed latency from each event to its NEAREST trigger, treating the trigger
% train as circular with period Trec so events before the first trigger and
% after the last one wrap instead of being thrown to a distant neighbour --
% which matters for the shifted nulls, where the wrap point moves every time.
% Nearest is resolved by midpoints, so it is exact, not interpolated.
ts = ts(:).';
te = [ts - Trec, ts, ts + Trec];
k  = discretize(tE(:), [-inf, (te(1:end-1) + te(2:end))/2, inf]);
L  = tE(:) - te(k).';
end

function s = fmt(md, p)
if isnan(md), s = sprintf('n.s. p=%.2f', p);
else,         s = sprintf('%+.0f ms (p=%.4f)', md, p);
end
end

function t = sessTag(p), pp = strsplit(p, filesep); t = sprintf('%s_%s', pp{end-1}, pp{end}); end
function s = shortrec(r), s = regexprep(r,'_\d{5}$',''); if numel(s)>16, s = s(1:16); end, end
