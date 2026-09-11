% Ventral_surface_dFF_heatmaps_svd_260606.m
% -----------------------------------------------------------------------
%  Population dF/F + spike figures for the Ventral_surface_summary dataset,
%  split SIG vs NON-SIG by the CURRENT (SVD) breath x Ca-spike coherence
%  (coherence_polar_svd_260606).  Matches the current method:
%    - breath from PC1 (breath_peak_pc1 / breath_insp_start_pc1)
%    - phase = onset 0, peak pi  (piecewise_phase_local, NOT uniform-linear)
%    - time figures triggered on inspiration ONSET
%    - spikes lead-shifted by ca_lag_sec (7 frames @30Hz, data-driven)
%    - Vglut2/1124 rising-edge +1-frame fix
%    - PHASE panels 0..2pi (no 4pi tile)
%
%  Six figures (each SIG | NON-SIG):
%    Fig 1  ROI-mean dF/F heatmap, TIME    (+/- t_win s, onset-triggered)
%    Fig 2  ROI-mean dF/F heatmap, PHASE   (0..2pi)
%    Fig 3  pooled spike-count hist, TIME  (each spike once, nearest onset)
%    Fig 4  pooled spike-phase hist + cos, PHASE
%    Fig 5  single-cycle dF/F, spike-cycles only, TIME, sorted by dt-to-nearest
%    Fig 6  single-cycle dF/F, spike-cycles only, PHASE, sorted by earliest spike phase
% -----------------------------------------------------------------------

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath')); repoRoot = fileparts(scriptDir);
addpath(repoRoot); addpath(fullfile(repoRoot,'2p_breathing_coherence'));

%% ===================== USER-EDITABLE =====================
rootPath = 'D:\Ventral_surface_summary';
cohData  = fullfile(rootPath, 'coherence_polar_svd_260606', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'dFF_heatmaps_svd_260606');

doSave = false;
t_win        = 3;            % +/- s for time figures
t_dt         = 0.05;         % s per column on common time grid
nPhaseBins   = 100;          % bins per cycle (dF/F phase grid)
nHistBins    = 24;           % bins per cycle for spike-phase histogram
prcLim       = [0.5 99.5];   % percentile clamp for color limits
climManual   = [];           % [] = auto; [lo hi] forces
ca_lag_sec   = 0.1;        % spike lead (shift earlier 3 frames @30Hz, 0.1 s)
spike_tick   = true;         % overlay spike ticks on single-cycle heatmaps
tick_color   = [0 0.4 1];
tick_lw      = 0.8;
fallback_fps = 30;
nDrop        = 30;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
if ~isfolder(outDir), mkdir(outDir); end

%% ---- coherence results + folder map ----
assert(isfile(cohData), 'coherence_polar_data.mat not found: %s', cohData);
S = load(cohData, 'PP','labels','confC');
PP = S.PP; labels = S.labels; confC = S.confC;
nROI = numel(labels);
fprintf('Loaded %d ROI labels.  confC=%.3f\n', nROI, confC);

allca = dir(fullfile(rootPath,'**','ca_spike_data.mat'));
folderMap = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(allca)
    pp = strsplit(allca(i).folder, filesep); bn = pp{end};
    if ~isKey(folderMap, bn), folderMap(bn) = allca(i).folder; end
end

t_grid     = -t_win:t_dt:t_win;  nW = numel(t_grid);
phase_1cyc = linspace(0, 2*pi, nPhaseBins+1); phase_1cyc = phase_1cyc(1:end-1);

% per-ROI summary
T_mat  = nan(nROI, nW);            % Fig1: ROI-mean dF/F, time
P_mat  = nan(nROI, nPhaseBins);    % Fig2: ROI-mean dF/F, phase
spk_t_pool  = cell(nROI,1);        % Fig3: nearest-onset dt (each spike once)
spk_ph_pool = cell(nROI,1);        % Fig4: spike phase
spk_ph_occ  = cell(nROI,1);        % Fig4: per-ROI phase occupancy (dwell) for compensation
cyc_T = cell(nROI,1); cyc_dt = cell(nROI,1); cyc_spk   = cell(nROI,1);   % Fig5 (time)
cyc_P = cell(nROI,1); cyc_phf= cell(nROI,1); cyc_spkph = cell(nROI,1);   % Fig6 (phase)

%% ---- per ROI ----
for ii = 1:nROI
    parts = regexp(labels{ii}, '/', 'split');
    if numel(parts) < 4, continue; end
    recName = strjoin(parts(3:end-1), '/'); rid = str2double(parts{end});
    if ~isKey(folderMap, recName), continue; end
    fp = folderMap(recName);
    try
        df = dir(fullfile(fp,'*_ch1_dFF.mat'));
        bp = dir(fullfile(fp,'breath_peak_pc1.mat'));
        ip = dir(fullfile(fp,'breath_insp_start_pc1.mat'));
        if isempty(df) || isempty(bp) || isempty(ip), continue; end

        fps = detect_session_fps(fp, fallback_fps);
        D   = load(fullfile(df(1).folder, df(1).name), 'dFF');
        BP  = load(fullfile(bp(1).folder, bp(1).name));
        IP  = load(fullfile(ip(1).folder, ip(1).name));
        CA  = load(fullfile(fp, 'ca_spike_data.mat'), 'roi_spikes');
        dff_all = double(D.dFF);
        if rid<1 || rid>size(dff_all,2) || rid>numel(CA.roi_spikes), continue; end
        tr = dff_all(:, rid);
        spk_train = double(CA.roi_spikes(rid).spike_train(:));
        if numel(spk_train) < numel(tr), spk_train(end+1:numel(tr)) = 0; end
        spk_train = spk_train(1:numel(tr));

        nB = numel(BP.breath);
        ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;  % PEAK
        ev(1:min(nDrop,numel(ev))) = [];
        ef = zeros(nB,1); fi = round(IP.insp_start_idx(:)); ef(fi(fi>=1 & fi<=nB)) = 1;   % ONSET
        ef(1:min(nDrop,numel(ef))) = [];
        % Vglut2/1124: rising-edge 2P trigger -> breath leads calcium by 1 frame; delay events.
        if contains(fp, fullfile('Vglut2','1124'))
            ev = [0; ev(1:end-1)]; ef = [0; ef(1:end-1)];
        end

        Tn = min([numel(ev), numel(ef), numel(tr), numel(spk_train)]);
        ev = ev(1:Tn); ef = ef(1:Tn); tr = tr(1:Tn); spk_train = spk_train(1:Tn);
        peakIdx = find(ev>0); on = find(ef>0);          % on = inspiration ONSET
        phi = piecewise_phase_local(peakIdx, on, Tn);   % onset=0, peak=pi (cumulative)

        lag = round(ca_lag_sec*fps);                    % GCaMP lead-comp
        spike_idx = find(spk_train>0) - lag; spike_idx = spike_idx(spike_idx>=1 & spike_idx<=Tn);

        win = round(t_win*fps); T_local = (-win:win)/fps;
        valid_on = on(on - win >= 1 & on + win <= Tn);

        %% Fig1: ROI-mean dF/F, onset-triggered TIME
        if ~isempty(valid_on)
            E = zeros(numel(valid_on), 2*win+1);
            for k = 1:numel(valid_on), E(k,:) = tr(valid_on(k)-win : valid_on(k)+win); end
            T_mat(ii,:) = interp1(T_local, mean(E,1), t_grid, 'linear', NaN);
        end

        %% Fig2: ROI-mean dF/F, PHASE (onset=0,peak=pi via cumulative phi)
        nCyc = numel(on) - 1;
        if nCyc >= 1
            dff_phase = nan(nCyc, nPhaseBins);
            for k = 1:nCyc
                i1 = on(k); i2 = on(k+1); if i2 <= i1, continue; end
                ps = phi(i1:i2); ds = tr(i1:i2); gd = ~isnan(ps);
                if nnz(gd) < 2, continue; end
                ps = ps(gd); ds = ds(gd); [ps,iu] = unique(ps); ds = ds(iu);
                dff_phase(k,:) = interp1(ps, ds, ps(1)+phase_1cyc, 'linear', NaN);
            end
            P_mat(ii,:) = mean(dff_phase, 1, 'omitnan');
        end

        %% Fig3: nearest-onset spike dt (each spike once)
        if ~isempty(spike_idx) && ~isempty(on)
            dt_all = nan(numel(spike_idx),1);
            for s = 1:numel(spike_idx)
                [~, mi] = min(abs(double(on) - double(spike_idx(s))));
                dt_all(s) = (spike_idx(s) - on(mi)) / fps;
            end
            spk_t_pool{ii} = dt_all(abs(dt_all) <= t_win);
        end

        %% Fig4: spike phase
        sp_ph = mod(phi(spike_idx), 2*pi); spk_ph_pool{ii} = sp_ph(~isnan(sp_ph));
        spk_ph_occ{ii} = histcounts(mod(phi(~isnan(phi)),2*pi), linspace(0,2*pi,nHistBins+1));  % phase dwell

        %% Fig5: single-cycle dF/F (spike-cycles only), TIME, dt to nearest spike
        nv = numel(valid_on);
        Tsnips = nan(nv, nW); Tdt = nan(nv,1); Tspk = cell(nv,1); nc = 0;
        for e = 1:nv
            c = valid_on(e);
            sw = spike_idx(spike_idx >= c-win & spike_idx <= c+win);
            if isempty(sw), continue; end
            bef = spike_idx(spike_idx < c); aft = spike_idx(spike_idx >= c);
            dt_bef = Inf; dt_aft = Inf;
            if ~isempty(bef), dt_bef = (c - bef(end))/fps; end
            if ~isempty(aft), dt_aft = (aft(1) - c)/fps; end
            if dt_bef <= dt_aft, dt_sign = -dt_bef; else, dt_sign = dt_aft; end
            nc = nc + 1;
            Tsnips(nc,:) = interp1(T_local, tr(c-win:c+win), t_grid, 'linear', NaN);
            Tdt(nc)      = dt_sign;
            Tspk{nc}     = (sw - c)/fps;
        end
        cyc_T{ii} = Tsnips(1:nc,:); cyc_dt{ii} = Tdt(1:nc); cyc_spk{ii} = Tspk(1:nc);

        %% Fig6: single-cycle dF/F (spike-cycles only), PHASE, earliest spike phase
        ncMax = max(numel(on)-1, 0);
        Psnips = nan(ncMax, nPhaseBins); Pfirst = nan(ncMax,1); Pspk = cell(ncMax,1); np = 0;
        for k = 1:numel(on)-1
            i1 = on(k); i2 = on(k+1); if i2 <= i1, continue; end
            sc = spike_idx(spike_idx >= i1 & spike_idx < i2);
            if isempty(sc), continue; end
            ps = phi(i1:i2); ds = tr(i1:i2); gd = ~isnan(ps);
            if nnz(gd) < 2, continue; end
            ps = ps(gd); ds = ds(gd); [ps,iu] = unique(ps); ds = ds(iu);
            np = np + 1;
            Psnips(np,:) = interp1(ps, ds, ps(1)+phase_1cyc, 'linear', NaN);
            spk_ph_cyc   = mod(phi(sc) - ps(1), 2*pi);               % 0 = onset
            Pfirst(np)   = min(spk_ph_cyc);
            Pspk{np}     = spk_ph_cyc;
        end
        cyc_P{ii} = Psnips(1:np,:); cyc_phf{ii} = Pfirst(1:np); cyc_spkph{ii} = Pspk(1:np);
    catch ME
        warning('  ERROR %s: %s', labels{ii}, ME.message);
    end
end

%% ---- masks + sort orders ----
sig = PP.r >= confC; nonsig = ~sig;
[~, os]  = sort(PP.th(sig));    ord_sig    = sub2abs(sig,    os);
[~, ons] = sort(PP.th(nonsig)); ord_nonsig = sub2abs(nonsig, ons);
fprintf('SIG n=%d | NON-SIG n=%d\n', sum(sig), sum(nonsig));

gclim = clamp_from(prcLim, climManual, [T_mat(:); P_mat(:)]);
dt_bin = t_grid(2) - t_grid(1);
peth_edges = [t_grid(:) - dt_bin/2; t_grid(end) + dt_bin/2];
peth_ctrs  = t_grid;
edges_hist = linspace(0, 2*pi, nHistBins+1); ctrs_hist = (edges_hist(1:end-1)+edges_hist(2:end))/2;

%%
% %% ============ FIG 1: dF/F TIME ============
% f1 = figure('Color','w','Name','dFF heatmap (time)','Units','centimeters','Position',[2 2 26 12]);
% draw_heatmap(subplot(1,2,1), T_mat(ord_sig,:),    t_grid, sprintf('SIG  n=%d', sum(sig)),     'time from insp onset (s)', gclim, '\DeltaF/F', 'time');
% draw_heatmap(subplot(1,2,2), T_mat(ord_nonsig,:), t_grid, sprintf('NON-SIG  n=%d', sum(nonsig)), 'time from insp onset (s)', gclim, '\DeltaF/F', 'time');
% sgtitle('onset-triggered dF/F (ROI mean)  |  rows sorted by preferred phase');

% %% ============ FIG 2: dF/F PHASE ============
% f2 = figure('Color','w','Name','dFF heatmap (phase)','Units','centimeters','Position',[2 2 26 12]);
% draw_heatmap(subplot(1,2,1), P_mat(ord_sig,:),    phase_1cyc, sprintf('SIG  n=%d', sum(sig)),     'breath phase (rad)', gclim, '\DeltaF/F', 'phase');
% draw_heatmap(subplot(1,2,2), P_mat(ord_nonsig,:), phase_1cyc, sprintf('NON-SIG  n=%d', sum(nonsig)), 'breath phase (rad)', gclim, '\DeltaF/F', 'phase');
% sgtitle('breath-phase dF/F (ROI mean)  |  onset=0, peak=\pi');
% 
% %% ============ FIG 3: spike-count hist TIME ============
% spk_t_sig = vertcat(spk_t_pool{sig}); spk_t_nonsig = vertcat(spk_t_pool{nonsig});
% f3 = figure('Color','w','Name','spike hist (time)','Units','centimeters','Position',[2 2 26 12]);
% draw_spike_hist_time(subplot(1,2,1), spk_t_sig,    peth_edges, peth_ctrs, sprintf('SIG  %d spikes', numel(spk_t_sig)));
% draw_spike_hist_time(subplot(1,2,2), spk_t_nonsig, peth_edges, peth_ctrs, sprintf('NON-SIG  %d spikes', numel(spk_t_nonsig)));
% sgtitle('onset-triggered spike-count histogram (pooled, each spike once)');

%% ============ FIG 4: spike-phase hist + cosine ============
spk_ph_sig = vertcat(spk_ph_pool{sig}); spk_ph_nonsig = vertcat(spk_ph_pool{nonsig});
occ_sig    = sum(vertcat(spk_ph_occ{sig}),    1);   % pooled phase dwell (sig)
occ_nonsig = sum(vertcat(spk_ph_occ{nonsig}), 1);   % pooled phase dwell (non-sig)
f4 = figure('Color','w','Name','spike hist (phase)','Units','centimeters','Position',[2 2 26 12]);
draw_spike_hist_phase(subplot(1,2,1), spk_ph_sig,    occ_sig,    edges_hist, ctrs_hist, sprintf('SIG  %d spikes', numel(spk_ph_sig)));
draw_spike_hist_phase(subplot(1,2,2), spk_ph_nonsig, occ_nonsig, edges_hist, ctrs_hist, sprintf('NON-SIG  %d spikes', numel(spk_ph_nonsig)));
sgtitle('breath-phase spike histogram + cosine fit  |  onset=0, peak=\pi');

%% ============ FIG 5: single-cycle dF/F, TIME (sorted by dt) ============
[T_sig,    ~, Tspk_sig   ] = pool_and_sort(cyc_T, cyc_dt, cyc_spk, sig);
[T_nonsig, ~, Tspk_nonsig] = pool_and_sort(cyc_T, cyc_dt, cyc_spk, nonsig);
gclim_cyc = clamp_from(prcLim, climManual, [T_sig(:); T_nonsig(:)]);
f5 = figure('Color','w','Name','single-cycle dFF (time)','Units','centimeters','Position',[2 2 26 12]);
draw_cycle_heatmap(subplot(1,2,1), T_sig,    t_grid, Tspk_sig,    gclim_cyc, tick_color, tick_lw, spike_tick, ...
    sprintf('SIG  %d spike-cycles', size(T_sig,1)),    'time from insp onset (s)', 'time');
draw_cycle_heatmap(subplot(1,2,2), T_nonsig, t_grid, Tspk_nonsig, gclim_cyc, tick_color, tick_lw, spike_tick, ...
    sprintf('NON-SIG  %d spike-cycles', size(T_nonsig,1)), 'time from insp onset (s)', 'time');
sgtitle('single-cycle dF/F, spike-containing cycles only, sorted by dt to nearest spike');

%% ============ FIG 6: single-cycle dF/F, PHASE (sorted by earliest spike phase) ============
[P_sig,    ~, Pspk_sig   ] = pool_and_sort(cyc_P, cyc_phf, cyc_spkph, sig);
[P_nonsig, ~, Pspk_nonsig] = pool_and_sort(cyc_P, cyc_phf, cyc_spkph, nonsig);
gclim_cycP = clamp_from(prcLim, climManual, [P_sig(:); P_nonsig(:)]);
f6 = figure('Color','w','Name','single-cycle dFF (phase)','Units','centimeters','Position',[2 2 26 12]);
draw_cycle_heatmap(subplot(1,2,1), P_sig,    phase_1cyc, Pspk_sig,    gclim_cycP, tick_color, tick_lw, spike_tick, ...
    sprintf('SIG  %d spike-cycles', size(P_sig,1)),    'breath phase (rad)', 'phase');
draw_cycle_heatmap(subplot(1,2,2), P_nonsig, phase_1cyc, Pspk_nonsig, gclim_cycP, tick_color, tick_lw, spike_tick, ...
    sprintf('NON-SIG  %d spike-cycles', size(P_nonsig,1)), 'breath phase (rad)', 'phase');
sgtitle('single-cycle dF/F, spike-containing cycles only, sorted by earliest spike phase  |  onset=0, peak=\pi');

%% ---- save ----
if doSave
    out_pairs = {f1,'dFF_heatmap_time'; f2,'dFF_heatmap_phase'; ...
                 f3,'spike_hist_time';  f4,'spike_hist_phase'; ...
                 f5,'cycle_dFF_sorted_time'; f6,'cycle_dFF_sorted_phase'};
    
    for k = 1:size(out_pairs,1)
        fh = out_pairs{k,1}; nm = out_pairs{k,2};
        exportgraphics(fh, fullfile(outDir,[nm '.png']), 'Resolution',200,'BackgroundColor','white');
        exportgraphics(fh, fullfile(outDir,[nm '.pdf']), 'ContentType','vector','BackgroundColor','white');
    end
    save(fullfile(outDir,'dFF_heatmap_data.mat'), 'T_mat','P_mat','spk_t_pool','spk_ph_pool', ...
         'T_sig','T_nonsig','P_sig','P_nonsig','Tspk_sig','Tspk_nonsig','Pspk_sig','Pspk_nonsig', ...
         't_grid','phase_1cyc','sig','labels','PP','confC','gclim','gclim_cyc','gclim_cycP');
    fprintf('Saved 6 figures + .mat to %s\n', outDir);

end
%% ========================= LOCAL FUNCTIONS ==========================
function idx = sub2abs(mask, localOrd)
    abs_idx = find(mask); idx = abs_idx(localOrd);
end

function c = clamp_from(prcLim, climManual, pool)
    if ~isempty(climManual), c = climManual; return; end
    f = pool(isfinite(pool));
    if isempty(f), c = [0 1]; return; end
    c = prctile(f, prcLim);
    if ~(c(2) > c(1)), c = [min(f) max(f)]; end
end

function [M_sorted, dt_sorted, spk_sorted] = pool_and_sort(M_cell, dt_cell, spk_cell, mask)
% Pool rows from ROIs selected by mask, sort by dt_cell ascending.
    keep = find(mask);
    M_all = []; dt_all = []; spk_all = {};
    for k = keep(:)'
        if isempty(M_cell{k}), continue; end
        M_all  = [M_all;  M_cell{k}];        %#ok<AGROW>
        dt_all = [dt_all; dt_cell{k}(:)];    %#ok<AGROW>
        spk_all = [spk_all; spk_cell{k}(:)]; %#ok<AGROW>
    end
    if isempty(M_all), M_sorted = M_all; dt_sorted = dt_all; spk_sorted = spk_all; return; end
    dt_sort = dt_all; dt_sort(isnan(dt_sort)) = Inf;
    [~, ord] = sort(dt_sort, 'ascend');
    M_sorted = M_all(ord, :); dt_sorted = dt_all(ord); spk_sorted = spk_all(ord);
end

function draw_heatmap(ax, M, xax, ttl, xlab, clim, cbLabel, mode)
    imagesc(ax, xax, 1:size(M,1), M); hold(ax,'on');
    colormap(ax, flipud(gray(256))); caxis(ax, clim);
    set(ax,'YDir','reverse'); axis(ax,'tight');
    if strcmp(mode,'phase')
        xline(ax, pi, 'k--', 'LineWidth',0.8); xlim(ax,[0 2*pi]);
        set(ax,'XTick',[0 pi 2*pi], 'XTickLabel',{'0','\pi','2\pi'});
    else
        xline(ax, 0, 'k--', 'LineWidth',0.8);
    end
    xlabel(ax, xlab); ylabel(ax, 'ROI (sorted by pref. phase)');
    title(ax, ttl); cb = colorbar(ax); cb.Label.String = cbLabel;
end

function draw_spike_hist_time(ax, dt_pool, edges, ctrs, ttl)
    hold(ax,'on');
    if ~isempty(dt_pool)
        counts = histcounts(dt_pool, edges);
        bar(ax, ctrs, counts, 1, 'FaceColor',[0.3 0.3 0.8], 'FaceAlpha',0.85, 'EdgeColor','none');
    end
    xline(ax, 0, 'k--', 'LineWidth',1); xlim(ax, [ctrs(1) ctrs(end)]);
    xlabel(ax,'time from insp onset (s)'); ylabel(ax,'spike count');
    title(ax, ttl); box(ax,'on');
end

function draw_spike_hist_phase(ax, ph_pool, occ_pool, edges, ctrs, ttl)
    hold(ax,'on');
    if ~isempty(ph_pool)
        cnt = histcounts(ph_pool, edges);
        if numel(occ_pool) == numel(cnt) && any(occ_pool>0)         % occupancy compensation
            Ob = mean(occ_pool(occ_pool>0)); pos = occ_pool>0;
            w = ones(size(cnt)); w(pos) = Ob ./ occ_pool(pos); cnt = cnt .* w;
        end
        prob = cnt / max(sum(cnt),1);
        mu_fit = angle(sum(prob .* exp(1i*ctrs)));                   % circular mean of (reweighted) dist
        a_fit  = mean(prob); b_fit = 2 * mean(prob .* cos(ctrs - mu_fit));
        bar(ax, ctrs, prob, 1, 'FaceColor',[0.3 0.3 0.8], 'FaceAlpha',0.85, 'EdgeColor','none');
        x_fit = linspace(0, 2*pi, 400); plot(ax, x_fit, a_fit + b_fit*cos(x_fit - mu_fit), 'r-', 'LineWidth',1.8);
    end
    xline(ax, pi, 'k--', 'LineWidth',1); xlim(ax, [0 2*pi]);
    set(ax,'XTick',[0 pi 2*pi], 'XTickLabel',{'0','\pi','2\pi'});
    xlabel(ax,'breath phase (rad)'); ylabel(ax,'P(spike)');
    title(ax, ttl); box(ax,'on');
end

function draw_cycle_heatmap(ax, M, xax, spk_cell, clim, tickCol, tickLW, doTicks, ttl, xlab, mode)
% Per-cycle stacked dFF heatmap with optional spike ticks per row (single cycle).
    if isempty(M) || size(M,1) == 0
        text(ax, 0.5, 0.5, 'no cycles', 'Units','normalized','HorizontalAlignment','center');
        axis(ax,'off'); title(ax, ttl); return;
    end
    imagesc(ax, xax, 1:size(M,1), M); hold(ax,'on');
    colormap(ax, flipud(gray(256))); caxis(ax, clim);
    set(ax,'YDir','reverse'); axis(ax,'tight');
    if doTicks
        for r = 1:size(M,1)
            sp = spk_cell{r}; if isempty(sp), continue; end
            for s = 1:numel(sp)
                plot(ax, [sp(s) sp(s)], [r-0.4 r+0.4], '-', 'Color', tickCol, 'LineWidth', tickLW);
            end
        end
    end
    if strcmp(mode,'phase')
        xline(ax, pi, 'k--', 'LineWidth',0.8); xlim(ax,[0 2*pi]);
        set(ax,'XTick',[0 pi 2*pi], 'XTickLabel',{'0','\pi','2\pi'});
    else
        xline(ax, 0, 'k--', 'LineWidth',0.8);
    end
    xlabel(ax, xlab); ylabel(ax, 'cycle (sorted)');
    title(ax, ttl); cb = colorbar(ax); cb.Label.String = '\DeltaF/F';
end

function phi = piecewise_phase_local(peak_idx, foot_idx, T)
phi = nan(T,1);
events = [peak_idx(:); foot_idx(:)];
types  = [ones(numel(peak_idx),1); zeros(numel(foot_idx),1)];
[events, ord] = sort(events); types = types(ord);
keep = true(size(events));
for i = 2:numel(events), if types(i) == types(i-1), keep(i) = false; end, end
events = events(keep); types = types(keep);
if numel(events) < 2, return; end
phases = nan(size(events)); phi_cur = types(1) * pi;
for i = 1:numel(events), phases(i) = phi_cur; phi_cur = phi_cur + pi; end
for i = 1:numel(events)-1
    a = events(i); b = events(i+1);
    if a < 1 || b > T || b <= a, continue; end
    phi(a:b) = linspace(phases(i), phases(i+1), b - a + 1);
end
end
