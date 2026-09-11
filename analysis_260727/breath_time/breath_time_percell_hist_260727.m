function breath_time_percell_hist_260727()
%% breath_time_percell_hist_260727  Per-cell, PER-RECORDING Ca-event histograms.
% -----------------------------------------------------------------------
% One figure per cell. Simple.
%
%     rows    = the recordings that cell appears in (+ a pooled row at the bottom
%               when there is more than one)
%     left    = histogram of Ca events vs time from inspiration ONSET
%     right   = histogram of Ca events vs time from the inspiratory PEAK
%
% That is the whole figure: counts of calcium events per time bin, a line at zero.
% Each panel title carries that recording's event count and its OWN test result
% (M_z and shuffle p), so a cell seen in three recordings is tested three times and
% you can see directly whether the modulation reproduces across them -- rather than
% only seeing the pooled answer.
%
% This replaces breath_time_peth_percell_compare_260727.m for reading; that one
% shows the shuffle envelopes and null distributions and is for auditing the
% statistics, not for looking at data.
%
% NOTE ON WHAT IS PLOTTED. The bars are raw EVENT COUNTS, which is what "histogram"
% means and is the easiest thing to sanity-check. The statistics underneath still
% use the exposure-normalised rate (events per second of imaging actually observed
% in that bin), because raw counts alias against the frame rate. For a single
% recording at one frame rate the two have the same shape, so the bars and the test
% agree; across pooled recordings at different frame rates they can differ slightly.
%
% Input : cell_pool.mat
% Output: <dataset>\analysis_260727\breath_time\percell_hist\cell_###.png
%         <dataset>\analysis_260727\breath_time\per_recording_tests.csv
%
% Runqi Zhang / 2026-07-28
close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir),'coh_ca_breath'));
cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
winIBI       = 2.5;    % window = +/- winIBI x median inter-breath interval, so the
                       %   trigger is CENTRED at 0 and the panel spans 5 IBI. Side peaks
                       %   at +/-1 and +/-2 IBI then show how many cycles the locking
                       %   persists for; how fast they wash out is the IBI jitter.
binWidth_s   = 0.050;
smoothWidth_s= 0.150;
ampFrac      = 0.20;   % trigger QC, amplitude only (long pauses never rejected)
minEvents_test = 10;   % a recording is TESTED at or above this; below it the
                       %   histogram is still drawn, just not tested
nShuffle     = 2000;
fdr_q        = 0.05;
which_cells  = 'tested';   % 'sig' | 'tested' | 'all'
showFirst    = 4;
rngSeed      = 260728;
doSave       = true;
% =====================================================================

rng(rngSeed);
outDir = fullfile(cfg.outRoot,'breath_time');
sub    = fullfile(outDir,'percell_hist');
fprintf('\n======= breath_time_percell_hist_260727 =======\n');
pool = ensure_pool_260727();
rec = pool.rec; obs = pool.obs;
keyOf = strings(max([pool.obsT.cell_id;0]),1);   % stable content-derived cell key
kk = ~isnan(pool.obsT.cell_id);
keyOf(pool.obsT.cell_id(kk)) = pool.obsT.cell_key(kk);

%% ---- window from the measured cycle ----
cyc = [];
for k = 1:numel(rec)
    if ~rec(k).usable, continue; end
    f = sort(rec(k).foot_idx(:));
    if numel(f) > 1, cyc = [cyc; diff(f)/rec(k).fps]; end %#ok<AGROW>
end
medCyc = median(cyc);
edges  = -winIBI*medCyc : binWidth_s : winIBI*medCyc;
ctrs   = edges(1:end-1) + binWidth_s/2;
nB     = numel(ctrs);
srchM  = ctrs >= 0 & ctrs <= medCyc;
smB    = max(1, round(smoothWidth_s/binWidth_s));
fprintf('median IBI %.3f s | window %.2f..%.2f s (+/-%.1f IBI) | %.0f ms bins (%d)\n', ...
        medCyc, edges(1), edges(end), winIBI, 1000*binWidth_s, nB);

%% ---- trigger sets per recording, both landmarks ----
TRG = struct('onset',{},'peak',{});
for k = 1:numel(rec)
    e = struct('onset',[],'peak',[]);
    if rec(k).usable
        fs = rec(k).fps; T = rec(k).T;
        f  = sort(rec(k).foot_idx(:));  p = sort(rec(k).peak_idx(:));
        if numel(f) >= 2
            bwz = (rec(k).bw(:) - median(rec(k).bw)) / max(mad(rec(k).bw,1)*1.4826, eps);
            amp = nan(numel(f)-1,1); pk = nan(numel(f)-1,1);
            for i = 1:numel(f)-1
                q = p(p > f(i) & p < f(i+1));
                if ~isempty(q), amp(i) = bwz(q(1)) - bwz(f(i)); pk(i) = q(1); end
            end
            ok = amp > ampFrac*median(amp,'omitnan');
            lo = ceil(edges(1)*fs); hi = floor(edges(end)*fs);
            on = f(ok);   on = on(on+lo>=1 & on+hi<=T);
            pe = pk(ok);  pe = pe(~isnan(pe));  pe = pe(pe+lo>=1 & pe+hi<=T);
            e.onset = on;  e.peak = pe;
        end
    end
    TRG(k) = e; %#ok<AGROW>
end

%% ---- per (cell, recording, trigger) test ----
fprintf('testing each recording separately...\n');
Rrec = struct('cell_id',{},'rec_name',{},'trigger',{},'n_events',{},'n_trig',{}, ...
              'counts',{},'rate',{},'Mz',{},'p',{},'lat_s',{},'tested',{},'obsIdx',{}, ...
              'sh_mean',{},'sh_var',{},'sh_lo',{},'sh_hi',{},'sh_glob',{});
for i = 1:numel(obs)
    o = obs(i);
    if ~o.usable || isnan(o.cell_id), continue; end
    k = o.rec; if isnan(k) || ~rec(k).usable, continue; end
    T = rec(k).T; fs = rec(k).fps;
    ev = full(double(o.spikes(:))); if numel(ev)<T, ev(end+1:T,1)=0; end; ev = ev(1:T);
    for tg = {'onset','peak'}
        f = TRG(k).(tg{1});
        if numel(f) < 2, continue; end
        [cnt, expSec, ccf, lagIdx, binOf] = one_peth(ev, f, T, fs, edges, nB);
        S = struct('cell_id',o.cell_id,'rec_name',rec(k).name,'trigger',string(tg{1}), ...
                   'n_events',nnz(ev>0),'n_trig',numel(f),'counts',cnt,'rate',cnt./max(expSec,eps), ...
                   'Mz',NaN,'p',NaN,'lat_s',NaN,'tested',false,'obsIdx',i, ...
                   'sh_mean',nan(1,nB),'sh_var',nan(1,nB),'sh_lo',nan(1,nB), ...
                   'sh_hi',nan(1,nB),'sh_glob',NaN);
        if S.n_events > 0
            % ---- shuffle: keep the FULL per-bin count matrix, not just the statistic ----
            % The bars are counts, so the band has to be in counts too. Exposure is
            % identical under a circular shift, so shuffled counts are directly
            % comparable to the observed bars.
            SH = zeros(nShuffle, nB);
            nl = nan(nShuffle,1);
            [Tobs, ~, S.lat_s] = peth_stat(S.rate, ctrs, srchM, smB);
            for s = 1:nShuffle
                d = randi(T)-1;
                idx = mod(lagIdx-1+d, T)+1;
                cs  = accumarray(binOf, ccf(idx), [nB 1])';
                SH(s,:) = cs;
                nl(s) = peth_stat(cs./max(expSec,eps), ctrs, srchM, smB);
            end
            % SMOOTH the shuffles exactly as the test smooths the observed data.
            % Without this the band is built from raw 50 ms counts (very spiky) while
            % the p-value comes from the 150 ms-smoothed rate, so the drawn band and
            % the reported p disagree -- the peak can be significant yet sit well
            % inside a band computed on unsmoothed counts.
            SHs = movmean(SH, smB, 2);
            S.sh_mean = mean(SHs,1);
            S.sh_var  = var(SHs,0,1);
            S.sh_lo   = prctile(SHs,  2.5, 1);     % POINTWISE 95% band
            S.sh_hi   = prctile(SHs, 97.5, 1);
            % GLOBAL band: 95th percentile of each shuffle's LARGEST excursion above
            % its own mean. A pointwise band is exceeded in ~5% of bins by chance;
            % this is the multiple-comparison-corrected version, and it is the band
            % the p-value actually corresponds to.
            S.sh_glob = prctile(max(SHs - S.sh_mean, [], 2), 95);
            nl = nl(isfinite(nl));
            if S.n_events >= minEvents_test && ~isempty(nl) && std(nl) > 0
                S.Mz = (Tobs - mean(nl))/std(nl);
                S.p  = (1 + nnz(nl >= Tobs))/(1 + numel(nl));
                S.tested = true;
            end
        end
        Rrec(end+1) = S; %#ok<AGROW>
    end
end
tf = [Rrec.tested];
q = nan(1,numel(Rrec));
if any(tf), q(tf) = bh_fdr([Rrec(tf).p], fdr_q); end
fprintf('  %d (recording x trigger) tests, %d significant at BH q<%.2f\n', ...
        nnz(tf), nnz(q <= fdr_q), fdr_q);

%% ---- which cells to draw ----
allC = unique([Rrec.cell_id]);
switch lower(which_cells)
    case 'sig', sel = unique([Rrec(tf & q <= fdr_q).cell_id]);
    case 'all', sel = allC;
    otherwise,  sel = unique([Rrec(tf).cell_id]);
end
assert(~isempty(sel),'No cells selected.');
mz = arrayfun(@(c) max([Rrec([Rrec.cell_id]==c).Mz, -Inf]), sel);
[~,o2] = sort(mz,'descend'); sel = sel(o2);

if ~isfolder(sub), mkdir(sub);
else, old = dir(fullfile(sub,'cell_*.png')); for z=1:numel(old), delete(fullfile(old(z).folder,old(z).name)); end
end

%% ---- one figure per cell ----
col = cfg.genotype_color;
hOff = figure('Color','w','Visible','off');
for si = 1:numel(sel)
    c = sel(si);
    idx = find([Rrec.cell_id] == c);
    recNames = unique([Rrec(idx).rec_name],'stable');
    nR = numel(recNames);
    nRow = nR + (nR > 1);                       % extra pooled row when >1 recording

    if si <= showFirst, hf = figure('Color','w','Visible','on');
    else, hf = hOff; clf(hf); set(0,'CurrentFigure',hf); end
    set(hf,'Units','centimeters','Position',[1 1 20 3.2*nRow + 2.2]);

    yTop = 0.90; yBot = 0.10; hgt = (yTop-yBot)/nRow;
    for r = 1:nRow
        for t = 1:2
            tg = ternary(t==1,"onset","peak");
            ax = axes('Parent',hf,'Position',[0.10+0.47*(t-1), yTop-r*hgt+0.035, 0.38, hgt*0.70]); %#ok<LAXES>
            hold(ax,'on'); box(ax,'on');
            if r <= nR
                j = idx([Rrec(idx).rec_name] == recNames(r) & [Rrec(idx).trigger] == tg);
                ttl = recNames(r);
            else
                j = idx([Rrec(idx).trigger] == tg);      % pooled row
                ttl = "POOLED";
            end
            if isempty(j), axis(ax,'off'); continue; end
            cnts = sum(cat(1, Rrec(j).counts), 1);
            shMu = sum(cat(1, Rrec(j).sh_mean), 1);

            % ---- shuffle band, in COUNTS so it overlays the bars directly ----
            if numel(j) == 1
                shLo = Rrec(j).sh_lo;  shHi = Rrec(j).sh_hi;
                bandNote = '';
            else
                % Pooled row: shuffles are independent across recordings, so variances
                % add. Percentiles do NOT add, hence the normal approximation here --
                % flagged in the title rather than passed off as exact.
                sd   = sqrt(sum(cat(1, Rrec(j).sh_var), 1));
                shLo = shMu - 1.96*sd;  shHi = shMu + 1.96*sd;
                bandNote = ' (band approx)';
            end
            % THREE things on the panel, nothing else:
            %   orange bars = observed Ca events
            %   grey line   = shuffle mean
            %   grey band   = 95% of shuffles
            % All three are SMOOTHED the same way (150 ms), so they are directly
            % comparable and match what the test sees. Smoothing the bars but not the
            % band (or vice versa) is what made the earlier version unreadable.
            fill(ax, [ctrs fliplr(ctrs)], [shHi fliplr(max(shLo,0))], [.80 .80 .80], ...
                 'EdgeColor','none');
            plot(ax, ctrs, shMu, '-', 'Color',[.40 .40 .40], 'LineWidth',1);
            bar(ax, ctrs, movmean(cnts, smB), 1, 'FaceColor',col, 'EdgeColor','none');

            % ---- guides: the trigger, and +/-1, +/-2 IBI ----
            xline(ax, 0, 'k-', 'LineWidth',1.2);
            for mIBI = [-2 -1 1 2]
                xline(ax, mIBI*medCyc, ':', 'Color',[.35 .35 .35], 'LineWidth',0.8);
            end
            xlim(ax,[edges(1) edges(end)]);
            if r == nRow
                xlabel(ax, sprintf('time from %s (s)   dotted = +/-1, +/-2 IBI', tg));
            else
                set(ax,'XTickLabel',[]);
            end
            if t == 1, ylabel(ax,'# Ca events'); end
            nev = sum([Rrec(j).n_events]);
            if numel(j) == 1 && Rrec(j).tested
                st = sprintf('  M_z=%.1f  p=%.3f', Rrec(j).Mz, Rrec(j).p);
            else
                st = '';
            end
            title(ax, sprintf('%s  |  %s  |  n=%d%s%s', upper(tg), shorten(ttl), nev, st, bandNote), ...
                  'Interpreter','none','FontSize',7);
        end
    end
    ck = '';  if c <= numel(keyOf), ck = char(keyOf(c)); end
    sgtitle({sprintf('cell %d   |   key %s   |   %d recording(s)', c, ck, nR), ...
             'orange = observed Ca events    grey line = shuffle mean    grey band = 95% of shuffles'}, ...
            'Interpreter','none','FontSize',9);
    if doSave
        exportgraphics(hf, fullfile(sub, sprintf('cell_%03d.png', c)), 'Resolution',150);
    end
end
if ishandle(hOff), close(hOff); end
fprintf('Saved %d per-cell histogram figures to %s\n', numel(sel), sub);

%% ---- per-recording table ----
if doSave
    T = table([Rrec.cell_id]', [Rrec.rec_name]', [Rrec.trigger]', [Rrec.n_events]', ...
              [Rrec.n_trig]', [Rrec.Mz]', [Rrec.p]', q', [Rrec.tested]', 1000*[Rrec.lat_s]', ...
        'VariableNames',{'cell_id','rec_name','trigger','n_events','n_triggers', ...
                         'Mz','p','q','tested','latency_ms'});
    writetable(T, fullfile(outDir,'per_recording_tests.csv'));
    fprintf('Saved per_recording_tests.csv to %s\n', outDir);
end
end

%% ========================= helpers =========================
function [cnt, expSec, ccf, lagIdx, binOf] = one_peth(ev, f, T, fs, edges, nB)
trig = zeros(T,1); trig(f) = 1;
ccf  = real(ifft(conj(fft(trig)) .* fft(ev)));   % boundary-safe triggers -> linear
lo = ceil(edges(1)*fs); hi = floor(edges(end)*fs);
m  = (lo:hi)';  bo = discretize(m/fs, edges);
keep = ~isnan(bo); m = m(keep); bo = bo(keep);
lagIdx = mod(m,T)+1;  binOf = bo;
cnt    = accumarray(bo, ccf(lagIdx), [nB 1])';
expSec = (accumarray(bo, 1, [nB 1]) * numel(f) / fs)';
end

function [Tex, rbar, lat] = peth_stat(rate, ctrs, srchM, smB)
rbar = mean(rate,'omitnan');
rs = rate; if smB > 1, rs = movmean(rate, smB); end
d = rs - rbar; d(~srchM) = -Inf;
[Tex, ix] = max(d);
lat = ctrs(ix);
end

function q = bh_fdr(p, ~)
p = p(:); n = numel(p);
[ps,ord] = sort(p);
qs = min(1, ps.*n./(1:n)');
for i = n-1:-1:1, qs(i) = min(qs(i), qs(i+1)); end
q = nan(n,1); q(ord) = qs; q = q';
end

function s = shorten(x)
x = char(x);
if numel(x) > 30, s = [x(1:14) '..' x(end-13:end)]; else, s = x; end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
