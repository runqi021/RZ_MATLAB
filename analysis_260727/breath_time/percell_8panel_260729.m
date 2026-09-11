function percell_8panel_260729()
%% percell_8panel_260729  One 8-panel figure per active cell.
% -----------------------------------------------------------------------
%   row 1   onset dF/F avg | peak dF/F avg | onset event hist | peak event hist
%   row 2   onset dF/F map | peak dF/F map | power spectrum   | spike-triggered dF/F
%
% Everything breath-triggered uses the SAME fixed 3 s window (+/-1.5 s) as
% population_hist_260729, so panels are comparable across cells and across animals
% whose breath rates differ (median IBI 2.13 s in 260721_Sert vs 0.57 s in the older
% ventral recordings). Zero is the trigger in every triggered panel.
%
% PANEL NOTES
%  dF/F avg      mean across ALL accepted breaths, shaded +/-SEM. The grey band is
%                the circular-shift shuffle 95%, so "is this rise real" is readable
%                without a p-value.
%  event hist    Ca event counts, smoothed identically to the shuffle band -- the
%                same three-element grammar as the per-cell histogram figure.
%  dF/F map      one row per breath, CHRONOLOGICAL order, every accepted breath.
%                Deliberately NOT sorted by event latency and NOT restricted to
%                event-containing cycles: doing either manufactures a diagonal band
%                out of noise, which is what made the older heatmaps unreadable.
%  power spectrum of the dF/F trace, 0.1-14 Hz, LINEAR axes as requested. Chronux
%                multitaper (TW=4), computed per recording then averaged on a common
%                frequency grid, since fps differs between recordings. The breath
%                rate is marked, so a peak there is breath-locked and a peak
%                elsewhere is the cell's own rhythm.
%  spike-trig    dF/F averaged on the cell's own detected events. This is partly
%                circular -- the events were detected FROM this dF/F -- so read it
%                as the event KERNEL / a detection sanity check, not as evidence.
%
% ACTIVE CELL = nnz(spike_train>0) > 5 pooled over the cell's recordings, matching
% the existing ventral criterion (spike_trigger_dFF.m:34, temporal_phase_perROI.m:56).
%
% Output: <dataset>\analysis_260727\breath_time\percell_8panel\cell_###.png
%
% Runqi Zhang / 2026-07-29
close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);
addpath(fullfile(fileparts(scriptDir),'coh_ca_breath'));
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));
cfg = coh_cfg_260727();

%% ===================== USER-EDITABLE PARAMETERS ======================
win_s        = 1.5;      % fixed +/-1.5 s -> 3 s window, all triggered panels
binWidth_s   = 0.050;
smoothWidth_s= 0.150;
activeMinEv  = 5;        % ventral criterion: active = nnz(spike_train>0) > this
ampFrac      = 0.20;     % trigger QC, amplitude only
staWin_s     = 1.0;      % spike-triggered dF/F half-window
fBand        = [0.1 14]; % power spectrum band, linear axes
TW           = 4;        % Chronux time-bandwidth
nShuffle     = 500;
showFirst    = 4;
rngSeed      = 260729;
doSave       = true;
% =====================================================================

rng(rngSeed);
outDir = fullfile(cfg.outRoot,'breath_time','percell_8panel');
fprintf('\n========== percell_8panel_260729 ==========\n');
pool = ensure_pool_260727();
rec = pool.rec; obs = pool.obs; cells = pool.cells;
EXCL = cfg.excludeRecordings;

edges = -win_s : binWidth_s : win_s;
ctrs  = edges(1:end-1) + binWidth_s/2;
nB    = numel(ctrs);
smB   = max(1, round(smoothWidth_s/binWidth_s));
fGrid = linspace(fBand(1), fBand(2), 400);
staTgrid = -staWin_s : binWidth_s/2 : staWin_s;   % fixed 25 ms grid, fps-independent
keyOf = strings(max([pool.obsT.cell_id;0]),1);
kk = ~isnan(pool.obsT.cell_id);
keyOf(pool.obsT.cell_id(kk)) = pool.obsT.cell_key(kk);

%% ---- trigger sets, both landmarks, per recording ----
TRG = struct('onset',{},'peak',{},'ibi',{});
for k = 1:numel(rec)
    e = struct('onset',[],'peak',[],'ibi',NaN);
    if rec(k).usable && ~any(strcmp(rec(k).name, EXCL))
        fs = rec(k).fps; T = rec(k).T;
        f = sort(rec(k).foot_idx(:)); p = sort(rec(k).peak_idx(:));
        if numel(f) >= 3
            bwz = (rec(k).bw(:)-median(rec(k).bw))/max(mad(rec(k).bw,1)*1.4826,eps);
            amp = nan(numel(f)-1,1); pkf = nan(numel(f)-1,1);
            for i = 1:numel(f)-1
                q = p(p>f(i) & p<f(i+1));
                if ~isempty(q), amp(i)=bwz(q(1))-bwz(f(i)); pkf(i)=q(1); end
            end
            g = amp > ampFrac*median(amp,'omitnan');
            lo = ceil(edges(1)*fs); hi = floor(edges(end)*fs);
            on = f(g);            on = on(on+lo>=1 & on+hi<=T);
            pe = pkf(g); pe = pe(~isnan(pe)); pe = pe(pe+lo>=1 & pe+hi<=T);
            e.onset = on; e.peak = pe; e.ibi = median(diff(f))/fs;
        end
    end
    TRG(k) = e; %#ok<AGROW>
end

%% ---- which cells ----
nEvOf = zeros(numel(cells),1);
for c = 1:numel(cells)
    ii = cells{c};
    for q = ii(:)'
        if obs(q).usable && ~any(strcmp(rec(obs(q).rec).name, EXCL))
            nEvOf(c) = nEvOf(c) + obs(q).n_spikes_used;
        end
    end
end
sel = find(nEvOf > activeMinEv);
[~,o] = sort(nEvOf(sel),'descend'); sel = sel(o);
fprintf('%d active cells (>%d events pooled)\n', numel(sel), activeMinEv);
assert(~isempty(sel),'No active cells.');

if ~isfolder(outDir), mkdir(outDir);
else, old = dir(fullfile(outDir,'cell_*.png')); for z=1:numel(old), delete(fullfile(old(z).folder,old(z).name)); end
end

col = cfg.genotype_color;
hOff = figure('Color','w','Visible','off');
for si = 1:numel(sel)
    c = sel(si);  ii = cells{c};
    D = gather_cell(ii, obs, rec, TRG, EXCL, edges, ctrs, nB, smB, ...
                    staTgrid, fGrid, TW, nShuffle);
    if isempty(D.onset.trials) && isempty(D.peak.trials), continue; end

    if si <= showFirst, hf = figure('Color','w','Visible','on');
    else, hf = hOff; clf(hf); set(0,'CurrentFigure',hf); end
    set(hf,'Units','centimeters','Position',[1 1 36 17]);
    W = 0.195; H = 0.33; X0 = 0.065; GX = 0.055; Y1 = 0.55; Y2 = 0.09;
    axp = @(cc,rr) axes('Parent',hf,'Position',[X0+(cc-1)*(W+GX), ternary(rr==1,Y1,Y2), W, H]); %#ok<LAXES>
    AX = gobjects(0);   % every TRIGGERED axes, so xlim can be enforced last

    for t = 1:2
        tg = ternary(t==1,'onset','peak');  S = D.(tg);
        % --- row 1: dF/F average with shuffle band ---
        a = axp(t,1); hold(a,'on'); box(a,'on');
        if ~isempty(S.trials)
            mu = mean(S.trials,1,'omitnan');
            se = std(S.trials,0,1,'omitnan')./sqrt(max(sum(~isnan(S.trials),1),1));
            fill(a,[ctrs fliplr(ctrs)],[S.dffShHi fliplr(S.dffShLo)],[.82 .82 .82],'EdgeColor','none');
            fill(a,[ctrs fliplr(ctrs)],[mu+se fliplr(mu-se)],col,'FaceAlpha',0.25,'EdgeColor','none');
            plot(a, ctrs, mu, '-','Color',col,'LineWidth',1.8);
        end
        xline(a,0,'k-','LineWidth',1.2);
        ylabel(a,'dF/F'); title(a,sprintf('%s  dF/F avg (n=%d breaths)',upper(tg),size(S.trials,1)),'FontSize',7.5);
        AX(end+1) = a;

        % --- row 1 cols 3-4: event histogram ---
        a = axp(t+2,1); hold(a,'on'); box(a,'on');
        fill(a,[ctrs fliplr(ctrs)],[S.evShHi fliplr(max(S.evShLo,0))],[.80 .80 .80],'EdgeColor','none');
        plot(a, ctrs, S.evShMu, '-','Color',[.40 .40 .40],'LineWidth',1);
        bar(a, ctrs, movmean(S.evCnt,smB), 1, 'FaceColor',col,'EdgeColor','none');
        xline(a,0,'k-','LineWidth',1.2);
        ylabel(a,'# Ca events'); title(a,sprintf('%s  event hist (%d events)',upper(tg),round(sum(S.evCnt))),'FontSize',7.5);
        xlabel(a,sprintf('time from %s (s)',tg));
        AX(end+1) = a;

        % --- row 2: dF/F heatmap, chronological ---
        a = axp(t,2); hold(a,'on');
        if ~isempty(S.trials)
            imagesc(a, ctrs, 1:size(S.trials,1), S.trials);
            cl = prctile(S.trials(:),[2 98]); if cl(2)<=cl(1), cl = cl(1)+[0 1]; end
            caxis(a, cl); colormap(a, parula); axis(a,'tight');
            plot(a,[0 0],ylim(a),'w-','LineWidth',1.2);
        end
        set(a,'YDir','normal');
        xlabel(a,sprintf('time from %s (s)',tg)); ylabel(a,'breath # (chronological)');
        AX(end+1) = a;
        title(a,sprintf('%s  dF/F per breath',upper(tg)),'FontSize',7.5);
    end

    % --- row 2 col 3: power spectrum ---
    a = axp(3,2); hold(a,'on'); box(a,'on'); grid(a,'on');
    if ~isempty(D.spec)
        plot(a, fGrid, D.spec, '-','Color',col,'LineWidth',1.4);
        if isfinite(D.fBreath)
            xline(a, D.fBreath, '--','Color',[.15 .15 .15],'LineWidth',1.2);
            text(a, D.fBreath, max(D.spec)*0.95, sprintf('  breath %.2f Hz',D.fBreath), 'FontSize',7);
        end
    end
    xlim(a, fBand); xlabel(a,'frequency (Hz)'); ylabel(a,'power');
    title(a,'dF/F power spectrum (linear)','FontSize',7.5);

    % --- row 2 col 4: spike-triggered dF/F ---
    a = axp(4,2); hold(a,'on'); box(a,'on'); grid(a,'on');
    if ~isempty(D.sta)
        mu = mean(D.sta,1,'omitnan');
        se = std(D.sta,0,1,'omitnan')./sqrt(max(sum(~isnan(D.sta),1),1));
        fill(a,[D.staT fliplr(D.staT)],[mu+se fliplr(mu-se)],col,'FaceAlpha',0.25,'EdgeColor','none');
        plot(a, D.staT, mu, '-','Color',col,'LineWidth',1.8);
    end
    xline(a,0,'k-','LineWidth',1.2);
    xlabel(a,'time from Ca event (s)'); ylabel(a,'dF/F');
    title(a,sprintf('spike-triggered dF/F (n=%d)',size(D.sta,1)),'FontSize',7.5);

    % enforce the triggered-window limits LAST: any later plotting call on a
    % different axes must not be able to leave these at auto-scale.
    for h = AX(:)', xlim(h, [edges(1) edges(end)]); end
    ck = ''; if c <= numel(keyOf), ck = char(keyOf(c)); end
    sgtitle({sprintf('cell %d   |   key %s   |   %d recording(s), %d events', ...
                     c, ck, numel(ii), nEvOf(c)), ...
             'orange = observed   grey band = 95% of circular-shift shuffles   window fixed at 3 s'}, ...
            'Interpreter','none','FontSize',9);
    if doSave
        exportgraphics(hf, fullfile(outDir, sprintf('cell_%03d.png', c)), 'Resolution',140);
    end
end
if ishandle(hOff), close(hOff); end
fprintf('Saved %d 8-panel figures to %s\n', numel(sel), outDir);
end

%% ========================= helpers =========================
function D = gather_cell(ii, obs, rec, TRG, EXCL, edges, ctrs, nB, smB, staTgrid, fGrid, TW, nShuffle)
D = struct('onset',blank(nB),'peak',blank(nB),'spec',[],'fBreath',NaN,'sta',[],'staT',[]);
specAcc = []; ibis = []; staAcc = []; staT = [];
for tg = {'onset','peak'}
    S = D.(tg{1});
    for q = ii(:)'
        o = obs(q);
        if ~o.usable, continue; end
        k = o.rec; if isnan(k) || ~rec(k).usable, continue; end
        if any(strcmp(rec(k).name, EXCL)), continue; end
        f = TRG(k).(tg{1});  if numel(f) < 3, continue; end
        T = rec(k).T; fs = rec(k).fps;
        dff = double(o.dff(:)); if numel(dff) < T, dff(end+1:T,1) = NaN; end
        ev  = full(double(o.spikes(:)));   % stored SPARSE in cell_pool; fft needs full
        if numel(ev) < T, ev(end+1:T,1) = 0; end
        ev  = ev(1:T);

        % dF/F trials, sampled at the bin centres (nearest frame, no interpolation)
        off = round(ctrs*fs);
        idx = f(:) + off;                          % nTrig x nB
        ok  = idx>=1 & idx<=T;
        tr  = nan(numel(f), nB);
        tr(ok) = dff(idx(ok));
        S.trials = [S.trials; tr];

        % events + shuffle band, in counts
        trig = zeros(T,1); trig(f) = 1;
        ccf  = real(ifft(conj(fft(trig)).*fft(ev)));
        lo = ceil(edges(1)*fs); hi = floor(edges(end)*fs);
        m = (lo:hi)'; bo = discretize(m/fs, edges);
        kp = ~isnan(bo); m = m(kp); bo = bo(kp);
        lag = mod(m,T)+1;
        S.evCnt = S.evCnt + accumarray(bo, ccf(lag), [nB 1])';
        SH = zeros(nShuffle, nB);
        for s2 = 1:nShuffle
            d = randi(T)-1;
            SH(s2,:) = accumarray(bo, ccf(mod(lag-1+d,T)+1), [nB 1])';
        end
        SH = movmean(SH, smB, 2);
        S.evShMu = S.evShMu + mean(SH,1);
        S.evShVar = S.evShVar + var(SH,0,1);

        % dF/F shuffle band: shift the dF/F trace, not the triggers
        SHd = zeros(nShuffle, nB);
        for s2 = 1:nShuffle
            d = randi(T)-1;
            dsh = circshift(dff(1:T), d);
            j = idx; j(~ok) = 1;
            v = dsh(j); v(~ok) = NaN;
            SHd(s2,:) = mean(v,1,'omitnan');
        end
        S.dffSh = [S.dffSh; SHd];

        if strcmp(tg{1},'onset')
            % power spectrum of the dF/F, per recording, on a common grid
            try
                pr.Fs = fs; pr.tapers = [TW 2*TW-1]; pr.pad = 0;
                pr.fpass = [fGrid(1) min(fGrid(end), fs/2)]; pr.err = 0;
                x = dff(1:T); x(isnan(x)) = 0; x = x - mean(x);
                [Sp,fp] = mtspectrumc(x, pr);
                specAcc = [specAcc; interp1(fp, Sp, fGrid, 'linear', NaN)]; %#ok<AGROW>
            catch
            end
            ibis = [ibis; TRG(k).ibi]; %#ok<AGROW>
            % spike-triggered dF/F, on a FIXED TIME grid. A frame-count window
            % (-w:w) has a different width at 30 / 42 / 47 fps, so trials from
            % different recordings could not be stacked.
            staT = staTgrid;
            sOff = round(staT*fs);
            e = find(ev > 0);
            e = e(e + sOff(1) >= 1 & e + sOff(end) <= T);
            if ~isempty(e)
                staAcc = [staAcc; dff(e(:) + sOff)]; %#ok<AGROW>
            end
        end
    end
    sd = sqrt(S.evShVar);
    S.evShLo = S.evShMu - 1.96*sd;  S.evShHi = S.evShMu + 1.96*sd;
    if ~isempty(S.dffSh)
        S.dffShLo = prctile(S.dffSh, 2.5, 1);  S.dffShHi = prctile(S.dffSh, 97.5, 1);
    else
        S.dffShLo = nan(1,nB); S.dffShHi = nan(1,nB);
    end
    D.(tg{1}) = S;
end
if ~isempty(specAcc), D.spec = mean(specAcc,1,'omitnan'); end
if ~isempty(ibis),   D.fBreath = 1/median(ibis,'omitnan'); end
D.sta = staAcc;  D.staT = staT;
end

function S = blank(nB)
S = struct('trials',zeros(0,nB),'evCnt',zeros(1,nB),'evShMu',zeros(1,nB), ...
           'evShVar',zeros(1,nB),'evShLo',nan(1,nB),'evShHi',nan(1,nB), ...
           'dffSh',zeros(0,nB),'dffShLo',nan(1,nB),'dffShHi',nan(1,nB));
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
