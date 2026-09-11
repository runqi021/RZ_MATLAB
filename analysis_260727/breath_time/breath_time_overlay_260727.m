function breath_time_overlay_260727()
%% breath_time_overlay_260727  Measure the breath cycle before analysing anything.
%  Picks the PETH window, and screens the triggers for junk detections.
% -----------------------------------------------------------------------
% Step 0 of the ABSOLUTE-TIME analysis, which replaces the phase warping.
%
% WHY TIME INSTEAD OF PHASE. Cycle-interpolated phase stretches every breath onto
% a common 0..2pi axis. Measured here, inspiration is only ~15% of the cycle, so
% the phase axis is occupied ~12x more densely in one half than the other and an
% unmodulated cell looks strongly modulated until an ECDF correction is applied.
% In absolute time from inspiration onset there is no stretch, no occupancy
% imbalance, and the null is flat. Nothing has to be explained away.
%
% TWO THINGS THIS SCRIPT EXISTS TO CATCH
%
% 1. HOW LONG ALIGNMENT SURVIVES. A fixed time window is only meaningful while
%    breath duration is reproducible. Panels 1-3 measure that directly, so the
%    window is chosen from data instead of assumed.
%
% 2. JUNK TRIGGERS. Every detected inspiration onset becomes a PETH trigger, and
%    a trigger that is not really an inspiration dilutes every PETH equally. On
%    this dataset ~10% of detected cycles have almost no breath amplitude: during
%    a genuine respiratory pause the PC1 trace is flat, and the peak/foot detector
%    puts spurious low-amplitude events into it. The effect is that a long apnoea
%    does NOT appear as one long interval -- it gets subdivided into several fake
%    short cycles, which is why a duration histogram alone will not show you your
%    pauses. Panel 4 finds them; panel 8 shows what removing them does.
%
% Reads cell_pool.mat only. Writes figures + a cycle table, nothing else.
% Output: <phys>\analysis_260727\breath_time\breath_time_overlay.png/.pdf
%                                            breath_cycles.csv
%
% Runqi Zhang / 2026-07-27
close all;

%% ---- path setup ----
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir), 'coh_ca_breath'));

cfg = coh_cfg_260727();
outDir = fullfile(cfg.outRoot, 'breath_time');
if ~isfolder(outDir), mkdir(outDir); end

%% ===================== USER-EDITABLE PARAMETERS ======================
maxLag_s   = 6.0;     % how far out to LOOK (the window you use will be shorter)
preLag_s   = 0.5;     % how far before inspiration onset to show
binWidth_s = 0.050;   % PETH bin width. Keep >= the slowest frame period
                      %   (33.3 ms at 30 fps).
ampFrac    = 0.20;    % TRIGGER QC: a cycle whose peak-minus-foot breath amplitude
                      %   is below this fraction of its recording's MEDIAN cycle
                      %   amplitude is treated as a failed detection, not a breath.
nWaveShow  = 400;     % breath waveforms drawn in panel 1
doSave     = true;
% =====================================================================

fprintf('\n=========== breath_time_overlay_260727 ===========\n');
pool = ensure_pool_260727();
rec = pool.rec; obs = pool.obs;

edges = -preLag_s : binWidth_s : maxLag_s;
ctrs  = edges(1:end-1) + binWidth_s/2;
tGrid = edges(1):binWidth_s:edges(end);

%% ---- 1. every cycle: durations, amplitude, waveform ----
W=[]; Ti=[]; Te=[]; Tc=[]; Amp=[]; RecI=[]; FootFr={}; GoodC={};
for k = 1:numel(rec)
    if ~rec(k).usable, FootFr{k}=[]; GoodC{k}=[]; continue; end %#ok<AGROW>
    fs = rec(k).fps;  T = rec(k).T;
    f = sort(rec(k).foot_idx(:));  p = sort(rec(k).peak_idx(:));
    bw = rec(k).bw(:);
    bwz = (bw - median(bw)) / max(mad(bw,1)*1.4826, eps);   % robust z
    a = nan(numel(f)-1,1);
    for i = 1:numel(f)-1
        pk = p(p > f(i) & p < f(i+1));
        if isempty(pk), continue; end
        a(i) = bwz(pk(1)) - bwz(f(i));
        Ti(end+1,1)=(pk(1)-f(i))/fs; Te(end+1,1)=(f(i+1)-pk(1))/fs; %#ok<AGROW>
        Tc(end+1,1)=(f(i+1)-f(i))/fs; Amp(end+1,1)=a(i); RecI(end+1,1)=k; %#ok<AGROW>
        idx = f(i) + round(tGrid*fs);
        ok = idx>=1 & idx<=T;  w = nan(1,numel(tGrid));  w(ok) = bwz(idx(ok));
        W(end+1,:) = w; %#ok<AGROW>
    end
    thr = ampFrac * median(a,'omitnan');
    FootFr{k} = f; %#ok<AGROW>
    GoodC{k}  = [a > thr; false];   % per foot: does the cycle STARTING here pass
end
nCyc = numel(Tc);
recMedAmp = accumarray(RecI, Amp, [numel(rec) 1], @(x) median(x,'omitnan'), NaN);
good = Amp > ampFrac * recMedAmp(RecI);
fprintf('cycles: %d | insp %.3f s (CV %.2f) | cycle %.3f s (CV %.2f)\n', ...
    nCyc, median(Ti), std(Ti)/mean(Ti), median(Tc), std(Tc)/mean(Tc));
fprintf('TRIGGER QC: %d of %d cycles (%.1f%%) fail the %.0f%%-of-median amplitude test\n', ...
    nnz(~good), nCyc, 100*mean(~good), 100*ampFrac);
fprintf('  failed cycles: duration med %.2f s, max %.2f s | passed: med %.2f s, max %.2f s\n', ...
    median(Tc(~good)), max(Tc(~good)), median(Tc(good)), max(Tc(good)));

%% ---- 2. pooled PETH, with and without the trigger QC ----
[rateAll, expAll] = pooled_peth(obs, rec, edges, ctrs, FootFr, GoodC, false);
[rateQC , ~     ] = pooled_peth(obs, rec, edges, ctrs, FootFr, GoodC, true);

%% ---- figure ----
fig = figure('Color','w','Name','breath time-domain diagnostics', ...
             'Units','centimeters','Position',[1 1 40 17]);
set(fig,'DefaultAxesFontSize',8);
col = cfg.genotype_color;

a1 = subplot(2,4,1); hold(a1,'on'); box(a1,'on');
sh = randperm(size(W,1), min(nWaveShow,size(W,1)));
plot(a1, tGrid, W(sh,:)', '-','Color',[.78 .78 .78],'LineWidth',0.3);
plot(a1, tGrid, median(W,1,'omitnan'), 'r-','LineWidth',2);
xline(a1,0,'k-'); xline(a1,median(Ti),'b--'); xline(a1,median(Tc),'k--');
xlim(a1,[edges(1) edges(end)]);
xlabel(a1,'time from insp onset (s)'); ylabel(a1,'breath (robust z)');
title(a1,{sprintf('%d cycles overlaid',size(W,1)),'median is flat after the 1st breath'},'FontSize',7);

a2 = subplot(2,4,2); hold(a2,'on'); box(a2,'on'); grid(a2,'on');
eb = 0:0.05:ceil(max(Tc));
histogram(a2, Ti, eb, 'FaceColor',[.20 .45 .80],'EdgeColor','none','DisplayName','inspiration');
histogram(a2, Tc, eb, 'FaceColor',[.35 .35 .35],'EdgeColor','none','FaceAlpha',0.5,'DisplayName','full cycle');
set(a2,'YScale','log'); xlim(a2,[0 ceil(max(Tc))]);
legend(a2,'Location','northeast');
xlabel(a2,'duration (s)'); ylabel(a2,'# cycles (log)');
title(a2,{sprintf('full cycle = foot->foot, med %.2f s',median(Tc)), ...
          sprintf('LOG y and full range: max %.1f s is now visible',max(Tc))},'FontSize',7);

a3 = subplot(2,4,3); hold(a3,'on'); box(a3,'on'); grid(a3,'on');
ts = sort(Tc);  plot(a3, ts, (1:numel(ts))/numel(ts), 'k-','LineWidth',1.5);
for fr = [0.05 0.5 0.95]
    v = prctile(Tc,fr*100);
    plot(a3,[v v],[0 fr],':','Color',[.5 .5 .5]); plot(a3,[0 v],[fr fr],':','Color',[.5 .5 .5]);
    text(a3,v,fr,sprintf('  %.2f s',v),'FontSize',7,'VerticalAlignment','bottom');
end
xlim(a3,[0 maxLag_s]); ylim(a3,[0 1]);
xlabel(a3,'time from insp onset (s)'); ylabel(a3,'fraction of cycles ended');
title(a3,'alignment holds while this is flat','FontSize',7);

a4 = subplot(2,4,4); hold(a4,'on'); box(a4,'on'); grid(a4,'on');
scatter(a4, Tc(good), Amp(good), 8, [.35 .35 .35], 'filled','MarkerFaceAlpha',0.35);
scatter(a4, Tc(~good), Amp(~good), 14, [.85 .10 .10], 'filled');
set(a4,'XScale','log');
xlabel(a4,'cycle duration (s)'); ylabel(a4,'peak - foot amplitude (robust z)');
title(a4,{sprintf('TRIGGER QC: %d of %d fail (%.1f%%)',nnz(~good),nCyc,100*mean(~good)), ...
          'red = flat breath, i.e. a detection inside a pause'},'FontSize',7);

a5 = subplot(2,4,5); hold(a5,'on'); box(a5,'on'); grid(a5,'on');
bar(a5, ctrs, rateAll, 1, 'FaceColor',col,'EdgeColor','none');
xline(a5,0,'k-'); xline(a5,median(Ti),'b--'); xline(a5,median(Tc),'k--');
xlim(a5,[edges(1) edges(end)]);
xlabel(a5,'time from insp onset (s)'); ylabel(a5,'events / s of exposure');
title(a5,'pooled PETH, ALL triggers','FontSize',7);

a6 = subplot(2,4,6); hold(a6,'on'); box(a6,'on'); grid(a6,'on');
m = ctrs <= 2.0;
bar(a6, ctrs(m), rateQC(m), 1, 'FaceColor',col,'EdgeColor','none');
xline(a6,0,'k-'); xline(a6,median(Ti),'b--');
xlabel(a6,'time from insp onset (s)'); ylabel(a6,'events / s of exposure');
title(a6,'first cycle, QC-passed triggers only','FontSize',7);

a7 = subplot(2,4,7); hold(a7,'on'); box(a7,'on'); grid(a7,'on');
plot(a7, ctrs, expAll/max(expAll), 'k-','LineWidth',1.2);
ylim(a7,[0 1.05]); xlim(a7,[edges(1) edges(end)]);
xlabel(a7,'time from insp onset (s)'); ylabel(a7,'relative exposure');
title(a7,{'exposure counted in FRAMES, not triggers','its ripple is the aliasing this divides out'},'FontSize',7);

a8 = subplot(2,4,8); hold(a8,'on'); box(a8,'on'); grid(a8,'on');
plot(a8, ctrs, rateAll, '-','Color',[.6 .6 .6],'LineWidth',1.4);
plot(a8, ctrs, rateQC , '-','Color',col,'LineWidth',1.8);
xline(a8,0,'k-'); xlim(a8,[edges(1) 2.5]);
legend(a8,{'all triggers','QC-passed only'},'Location','northeast','FontSize',7);
xlabel(a8,'time from insp onset (s)'); ylabel(a8,'events / s of exposure');
title(a8,'what removing the junk triggers does','FontSize',7);

sgtitle(sprintf(['breath diagnostics  |  %d cycles  |  insp %.0f ms (CV %.2f), cycle %.2f s (CV %.2f), max %.1f s  ' ...
                 '|  %d junk triggers removed (%.1f%%)'], nCyc, 1000*median(Ti), std(Ti)/mean(Ti), ...
                 median(Tc), std(Tc)/mean(Tc), max(Tc), nnz(~good), 100*mean(~good)),'FontSize',9);

%% ---- guidance + save ----
fprintf('\n---- window guidance ----\n');
fprintf('  next onset: 5%%%% %.2f s | 50%%%% %.2f s | 95%%%% %.2f s\n', prctile(Tc,5), prctile(Tc,50), prctile(Tc,95));
fprintf('  6 x median inspiration = %.2f s = %.2f cycles (NOT 2-3)\n', 6*median(Ti), 6*median(Ti)/median(Tc));
fprintf('  recommended PETH window: -0.5 .. 2.0 s\n');
pk1 = max(rateAll(ctrs>0 & ctrs<1)); pk2 = max(rateQC(ctrs>0 & ctrs<1));
b1 = mean(rateAll(ctrs<0)); b2 = mean(rateQC(ctrs<0));
fprintf('  pooled peak/baseline: ALL triggers %.2f  ->  QC-passed %.2f\n', pk1/b1, pk2/b2);

if doSave
    exportgraphics(fig, fullfile(outDir,'breath_time_overlay.png'),'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'breath_time_overlay.pdf'),'ContentType','vector','BackgroundColor','white');
    Tt = table(RecI, string({rec(RecI).name})', Tc, Ti, Te, Amp, good, ...
        'VariableNames',{'rec_index','rec_name','cycle_s','insp_s','exp_s','amplitude_z','passes_qc'});
    writetable(Tt, fullfile(outDir,'breath_cycles.csv'));
    fprintf('\nSaved breath_time_overlay.png/.pdf + breath_cycles.csv to %s\n', outDir);
end
end

%% ========================= helpers =========================
function [rate, expSec] = pooled_peth(obs, rec, edges, ctrs, FootFr, GoodC, useQC)
cnt = zeros(1,numel(ctrs)); expSec = zeros(1,numel(ctrs));
for i = 1:numel(obs)
    o = obs(i);
    if ~o.usable, continue; end
    k = o.rec; if ~rec(k).usable || isempty(FootFr{k}), continue; end
    fs = rec(k).fps; T = rec(k).T;
    f = FootFr{k};
    if useQC, f = f(GoodC{k}); end
    evFr = find(o.spikes);
    if isempty(evFr) || isempty(f), continue; end
    loFr = ceil(edges(1)*fs); hiFr = floor(edges(end)*fs);
    for j = 1:numel(f)
        d = (evFr - f(j))/fs;
        d = d(d >= edges(1) & d < edges(end));
        if ~isempty(d), cnt = cnt + histcounts(d, edges); end
        fr = (f(j)+loFr):(f(j)+hiFr);  fr = fr(fr>=1 & fr<=T);
        if isempty(fr), continue; end
        expSec = expSec + histcounts((fr-f(j))/fs, edges)/fs;
    end
end
rate = cnt ./ max(expSec, eps);
end
