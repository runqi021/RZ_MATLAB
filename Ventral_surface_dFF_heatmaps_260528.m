% Ventral_surface_dFF_heatmaps_260528.m
% -----------------------------------------------------------------------
%  Eight population figures for the Ventral_surface_summary dataset,
%  split into SIG vs NON-SIG ROIs by breath-peak x Ca-spike coherence
%  (saved by Ventral_surface_coherence_polar_260528.m).
%
%    Fig 1  ROI-mean dF/F heatmap, TIME    (sig | non-sig), +/- t_win s
%    Fig 2  ROI-mean dF/F heatmap, PHASE   (sig | non-sig), [0, 4pi] tile
%    Fig 3  pooled spike-count hist, TIME  (sig | non-sig)
%    Fig 4  pooled spike-phase hist+cos,   (sig | non-sig), [0, 4pi] tile
%    Fig 5  ROI-mean spike rate (Hz), TIME (sig | non-sig)
%    Fig 6  ROI-mean spike rate (Hz), PHASE(sig | non-sig), [0, 4pi] tile
%    Fig 7  SINGLE-CYCLE dF/F, spike-cycles only, sorted by dt-to-nearest,
%           TIME (sig | non-sig)  -- matches breath_combined_summary
%           local_build_stack_sorted (R4,C1 in the 5x4 grid)
%    Fig 8  SINGLE-CYCLE dF/F, spike-cycles only, sorted by min spike phase,
%           PHASE (sig | non-sig) -- matches local_build_stack_sorted_phase
%
%  Assumes coherence_polar_data.mat lives at
%     <rootPath>/coherence_polar_260528/coherence_polar_data.mat
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));

%% ===================== USER-EDITABLE =================================
rootPath = 'D:\Ventral_surface_summary';
dataMat  = fullfile(rootPath, 'coherence_polar_260528', 'coherence_polar_data.mat');
outDir   = fullfile(rootPath, 'coherence_polar_260528');

t_win        = 3;           % +/- s for time heatmap
t_dt         = 0.05;        % s per column on common time grid
nPhaseBins   = 100;         % bins per cycle (dF/F phase grid)
nHistBins    = 24;          % bins per cycle for spike-phase histogram
prcLim       = [0.5 99.5];    % percentile clamp for color limits (raw dF/F)
climManual   = [];          % [] = auto; [lo hi] forces

spike_tick   = true;        % overlay spike ticks on cycle heatmaps (Fig 7/8)
tick_color   = [0 0.4 1];
tick_lw      = 0.8;

fallback_fps = 30;
nDrop        = 30;
% =====================================================================

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');
if ~isfolder(outDir), mkdir(outDir); end

%% ---- load coherence results ----
assert(isfile(dataMat), 'Run Ventral_surface_coherence_polar_260528.m first (missing %s)', dataMat);
S = load(dataMat, 'PP','labels','confC');
PP = S.PP; labels = S.labels; confC = S.confC;

nROI = numel(labels);
fprintf('Loaded %d ROI labels.  confC=%.3f\n', nROI, confC);

t_grid     = -t_win:t_dt:t_win;
nW         = numel(t_grid);
phase_1cyc = linspace(0, 2*pi, nPhaseBins+1); phase_1cyc = phase_1cyc(1:end-1);

% per-ROI summary matrices
T_mat  = nan(nROI, nW);
P_mat  = nan(nROI, nPhaseBins);
ST_mat = nan(nROI, nW);
SP_mat = nan(nROI, nPhaseBins);
spk_t_pool  = cell(nROI,1);
spk_ph_pool = cell(nROI,1);

% per-ROI cycle-level cell arrays (only cycles with >=1 spike)
cyc_T   = cell(nROI,1);   cyc_dt       = cell(nROI,1);   cyc_spk   = cell(nROI,1);
cyc_P   = cell(nROI,1);   cyc_ph_first = cell(nROI,1);   cyc_spkph = cell(nROI,1);

%% ---- per ROI: locate folder, load dFF / breath / ca_spike_data ----
scan_map = struct('IO','ChAT','ChAT','ChAT','Vglut2','Vglut2','Vgat','Vgat','Sst','Sst');

for ii = 1:nROI
    lab = labels{ii};
    parts = regexp(lab, '/', 'split');
    if numel(parts) < 4, continue; end
    gname = parts{1}; dateStr = parts{2}; fov = parts{3}; rid = str2double(parts{4});
    if ~isfield(scan_map, gname), continue; end
    sname = scan_map.(gname);

    hit = dir(fullfile(rootPath, sname, dateStr, '*', fov, 'ca_spike_data.mat'));
    if isempty(hit), continue; end
    fp = hit(1).folder;

    try
        df = dir(fullfile(fp,'*_ch1_dFF.mat'));
        bp = dir(fullfile(fp,'*breath_peak_data.mat'));
        if isempty(df) || isempty(bp), continue; end

        fps = detect_session_fps(fp, fallback_fps);
        D   = load(fullfile(df(1).folder, df(1).name), 'dFF');
        BPc = load(fullfile(bp(1).folder, bp(1).name));
        CA  = load(fullfile(fp, 'ca_spike_data.mat'));
        dff_all = double(D.dFF);
        if rid<1 || rid>size(dff_all,2), continue; end
        tr = dff_all(:, rid);
        spk_train = double(CA.roi_spikes(rid).spike_train(:));
        if numel(spk_train) < numel(tr), spk_train(end+1:numel(tr)) = 0; end
        spk_train = spk_train(1:numel(tr));

        if isfield(BPc,'insp_onsets_train') && numel(BPc.insp_onsets_train)==numel(BPc.breath)
            ev = double(BPc.insp_onsets_train(:) ~= 0);
        else
            ev = zeros(numel(BPc.breath),1); oi = round(BPc.insp_onset_idx(:));
            ev(oi(oi>=1 & oi<=numel(ev))) = 1;
        end
        ev(1:min(nDrop,numel(ev))) = [];
        Tn = min([numel(ev), numel(tr), numel(spk_train)]);
        ev = ev(1:Tn); tr = tr(1:Tn); spk_train = spk_train(1:Tn);
        on = find(ev>0);
        spike_idx = find(spk_train > 0);
        win = round(t_win*fps);
        T_local = (-win:win)/fps;

        %% ===== ROI-mean dF/F + Hz (time, phase) =====
        valid_on = on(on - win >= 1 & on + win <= Tn);
        if ~isempty(valid_on)
            E  = zeros(numel(valid_on), 2*win+1);
            ES = zeros(numel(valid_on), 2*win+1);
            for k = 1:numel(valid_on)
                E (k,:) = tr       (valid_on(k)-win : valid_on(k)+win);
                ES(k,:) = spk_train(valid_on(k)-win : valid_on(k)+win);
            end
            T_mat (ii,:) = interp1(T_local, mean(E, 1),         t_grid, 'linear', NaN);
            ST_mat(ii,:) = interp1(T_local, mean(ES,1) * fps,   t_grid, 'linear', NaN);
        end
        nCyc = numel(on) - 1;
        if nCyc >= 1
            dff_phase = nan(nCyc, nPhaseBins);
            spk_phase = nan(nCyc, nPhaseBins);
            for k = 1:nCyc
                i1 = on(k); i2 = on(k+1);
                if i2 > i1
                    chunk    = tr(i1:i2);
                    schunk   = spk_train(i1:i2);
                    orig_ph  = linspace(0, 2*pi, numel(chunk));
                    dff_phase(k,:) = interp1(orig_ph, chunk,  phase_1cyc, 'linear');
                    cyc_dur_s = (i2 - i1) / fps;
                    spk_phase(k,:) = interp1(orig_ph, schunk, phase_1cyc, 'linear') ...
                                     * numel(schunk) / cyc_dur_s;
                end
            end
            P_mat (ii,:) = mean(dff_phase, 1, 'omitnan');
            SP_mat(ii,:) = mean(spk_phase, 1, 'omitnan');
        end

        %% ===== pooled spike times/phases (for Fig 3/4) =====
        if ~isempty(spike_idx) && ~isempty(on)
            dt_all = nan(numel(spike_idx),1);
            for s = 1:numel(spike_idx)
                [~, mi] = min(abs(double(on) - double(spike_idx(s))));
                dt_all(s) = (spike_idx(s) - on(mi)) / fps;
            end
            spk_t_pool{ii} = dt_all(abs(dt_all) <= t_win);
        end
        if ~isempty(spike_idx) && numel(on) >= 2
            ph_all = nan(numel(spike_idx),1);
            for s = 1:numel(spike_idx)
                c = spike_idx(s);
                bef = on(on <= c); aft = on(on >  c);
                if ~isempty(bef) && ~isempty(aft)
                    ph_all(s) = 2*pi * (c - bef(end)) / (aft(1) - bef(end));
                end
            end
            spk_ph_pool{ii} = ph_all(~isnan(ph_all));
        end

        %% ===== single-cycle dF/F, spike-cycles only (Fig 7) =====
        Tsnips = []; Tdt = []; Tspk = {};
        for e = 1:numel(valid_on)
            c = valid_on(e);
            spk_in_win = spike_idx(spike_idx >= c - win & spike_idx <= c + win);
            if isempty(spk_in_win), continue; end
            snip = tr(c-win : c+win);
            dt_bef = Inf; dt_aft = Inf;
            bef = spike_idx(spike_idx <  c);
            aft = spike_idx(spike_idx >= c);
            if ~isempty(bef), dt_bef = (c - bef(end))/fps; end
            if ~isempty(aft), dt_aft = (aft(1) - c)/fps; end
            if dt_bef <= dt_aft, dt_sign = -dt_bef; else, dt_sign = dt_aft; end
            Tsnips(end+1,:) = interp1(T_local, snip, t_grid, 'linear', NaN);   %#ok<AGROW>
            Tdt(end+1,1)    = dt_sign;                                          %#ok<AGROW>
            Tspk{end+1,1}   = (spk_in_win - c) / fps;                           %#ok<AGROW>
        end
        cyc_T{ii} = Tsnips; cyc_dt{ii} = Tdt; cyc_spk{ii} = Tspk;

        %% ===== single-cycle dF/F, spike-cycles only (Fig 8: phase) =====
        Psnips = []; Pfirst = []; Pspk = {};
        for c = 1:numel(on)-1
            i1 = on(c); i2 = on(c+1);
            if i2 <= i1, continue; end
            spk_in_cyc = spike_idx(spike_idx >= i1 & spike_idx < i2);
            if isempty(spk_in_cyc), continue; end
            chunk    = tr(i1:i2);
            orig_ph  = linspace(0, 2*pi, numel(chunk));
            ph_row   = interp1(orig_ph, chunk, phase_1cyc, 'linear');
            spk_ph_cyc = 2*pi * (double(spk_in_cyc) - i1) / (i2 - i1);
            Psnips(end+1,:) = ph_row;                                           %#ok<AGROW>
            Pfirst(end+1,1) = min(spk_ph_cyc);                                  %#ok<AGROW>
            Pspk{end+1,1}   = spk_ph_cyc;                                       %#ok<AGROW>
        end
        cyc_P{ii} = Psnips; cyc_ph_first{ii} = Pfirst; cyc_spkph{ii} = Pspk;
    catch ME
        warning('  ERROR %s: %s', lab, ME.message);
    end
end

%% ---- masks + sort orders ----
sig    = PP.r >= confC;
nonsig = ~sig;
[~, ord_sig_t]    = sort(PP.th(sig));
[~, ord_nonsig_t] = sort(PP.th(nonsig));
ord_sig    = sub2abs(sig,    ord_sig_t);
ord_nonsig = sub2abs(nonsig, ord_nonsig_t);

% Tile phase to [0, 4pi)
P_tile   = [P_mat,  P_mat ];
SP_tile  = [SP_mat, SP_mat];
phs_tile = [phase_1cyc, phase_1cyc + 2*pi];

% color clamps
gclim       = clamp_from(prcLim, climManual, [T_mat(:);  P_mat(:) ]);
gclim_event = clamp_from(prcLim, climManual, [ST_mat(:); SP_mat(:)]);

%% ============ FIG 1: ROI-mean dFF, TIME ============
f1 = figure('Color','w','Name','dFF heatmap (time, ROI-mean)','Units','centimeters','Position',[2 2 26 12]);
draw_heatmap(subplot(1,2,1), T_mat(ord_sig,:),    t_grid, sprintf('SIG  n=%d', sum(sig)),     'time (s)', gclim, '\DeltaF/F', 'time');
draw_heatmap(subplot(1,2,2), T_mat(ord_nonsig,:), t_grid, sprintf('NON-SIG  n=%d', sum(nonsig)), 'time (s)', gclim, '\DeltaF/F', 'time');
sgtitle('breath-triggered dF/F (ROI mean)  |  rows sorted by preferred phase');

%% ============ FIG 2: ROI-mean dFF, PHASE ============
f2 = figure('Color','w','Name','dFF heatmap (phase, ROI-mean)','Units','centimeters','Position',[2 2 26 12]);
draw_heatmap(subplot(1,2,1), P_tile(ord_sig,:),    phs_tile, sprintf('SIG  n=%d', sum(sig)),     'phase (rad)', gclim, '\DeltaF/F', 'phase');
draw_heatmap(subplot(1,2,2), P_tile(ord_nonsig,:), phs_tile, sprintf('NON-SIG  n=%d', sum(nonsig)), 'phase (rad)', gclim, '\DeltaF/F', 'phase');
sgtitle('breath-phase dF/F (ROI mean, tiled [0, 4\pi])');

%% ============ FIG 3: pooled spike-count histogram, TIME ============
spk_t_sig    = vertcat(spk_t_pool{sig});
spk_t_nonsig = vertcat(spk_t_pool{nonsig});
dt_bin     = t_grid(2) - t_grid(1);
peth_edges = [t_grid(:) - dt_bin/2; t_grid(end) + dt_bin/2];
peth_ctrs  = t_grid;
f3 = figure('Color','w','Name','spike hist (time)','Units','centimeters','Position',[2 2 26 12]);
draw_spike_hist_time(subplot(1,2,1), spk_t_sig,    peth_edges, peth_ctrs, sprintf('SIG  %d spikes', numel(spk_t_sig)));
draw_spike_hist_time(subplot(1,2,2), spk_t_nonsig, peth_edges, peth_ctrs, sprintf('NON-SIG  %d spikes', numel(spk_t_nonsig)));
sgtitle('breath-triggered spike-count histogram (pooled)');

%% ============ FIG 4: pooled spike-phase hist + cosine fit ============
spk_ph_sig    = vertcat(spk_ph_pool{sig});
spk_ph_nonsig = vertcat(spk_ph_pool{nonsig});
edges_hist = linspace(0, 2*pi, nHistBins+1);
ctrs_hist  = (edges_hist(1:end-1) + edges_hist(2:end))/2;
f4 = figure('Color','w','Name','spike hist (phase)','Units','centimeters','Position',[2 2 26 12]);
draw_spike_hist_phase(subplot(1,2,1), spk_ph_sig,    edges_hist, ctrs_hist, sprintf('SIG  %d spikes', numel(spk_ph_sig)));
draw_spike_hist_phase(subplot(1,2,2), spk_ph_nonsig, edges_hist, ctrs_hist, sprintf('NON-SIG  %d spikes', numel(spk_ph_nonsig)));
sgtitle('breath-phase spike histogram + cosine fit (tiled [0, 4\pi])');

%% ============ FIG 5: ROI-mean spike rate (Hz), TIME ============
f5 = figure('Color','w','Name','spike rate (time)','Units','centimeters','Position',[2 2 26 12]);
draw_heatmap(subplot(1,2,1), ST_mat(ord_sig,:),    t_grid, sprintf('SIG  n=%d', sum(sig)),     'time (s)',     gclim_event, 'spike rate (Hz)', 'time');
draw_heatmap(subplot(1,2,2), ST_mat(ord_nonsig,:), t_grid, sprintf('NON-SIG  n=%d', sum(nonsig)), 'time (s)',   gclim_event, 'spike rate (Hz)', 'time');
sgtitle('breath-triggered calcium event rate (Hz, ROI mean)');

%% ============ FIG 6: ROI-mean spike rate (Hz), PHASE ============
f6 = figure('Color','w','Name','spike rate (phase)','Units','centimeters','Position',[2 2 26 12]);
draw_heatmap(subplot(1,2,1), SP_tile(ord_sig,:),    phs_tile, sprintf('SIG  n=%d', sum(sig)),     'phase (rad)', gclim_event, 'spike rate (Hz)', 'phase');
draw_heatmap(subplot(1,2,2), SP_tile(ord_nonsig,:), phs_tile, sprintf('NON-SIG  n=%d', sum(nonsig)), 'phase (rad)', gclim_event, 'spike rate (Hz)', 'phase');
sgtitle('breath-phase calcium event rate (Hz, ROI mean, tiled [0, 4\pi])');

%% ============ FIG 7: SINGLE-CYCLE dFF, spike-cycles only, sorted by dt ============
[T_sig,    ~, Tspk_sig   ] = pool_and_sort(cyc_T, cyc_dt, cyc_spk, sig);
[T_nonsig, ~, Tspk_nonsig] = pool_and_sort(cyc_T, cyc_dt, cyc_spk, nonsig);
gclim_cyc = clamp_from(prcLim, climManual, [T_sig(:); T_nonsig(:)]);
f7 = figure('Color','w','Name','single-cycle dFF (time)','Units','centimeters','Position',[2 2 26 12]);
draw_cycle_heatmap(subplot(1,2,1), T_sig,    t_grid, Tspk_sig,    gclim_cyc, tick_color, tick_lw, spike_tick, ...
    sprintf('SIG  %d spike-cycles', size(T_sig,1)),    'time (s)', 'time');
draw_cycle_heatmap(subplot(1,2,2), T_nonsig, t_grid, Tspk_nonsig, gclim_cyc, tick_color, tick_lw, spike_tick, ...
    sprintf('NON-SIG  %d spike-cycles', size(T_nonsig,1)), 'time (s)', 'time');
sgtitle('single-cycle dF/F, spike-containing cycles only, sorted by dt to nearest spike');

%% ============ FIG 8: SINGLE-CYCLE dFF, spike-cycles only, sorted by min spike phase ============
[P_sig,    ~, Pspk_sig   ] = pool_and_sort(cyc_P, cyc_ph_first, cyc_spkph, sig);
[P_nonsig, ~, Pspk_nonsig] = pool_and_sort(cyc_P, cyc_ph_first, cyc_spkph, nonsig);
gclim_cycP = clamp_from(prcLim, climManual, [P_sig(:); P_nonsig(:)]);
P_sig_tile    = [P_sig,    P_sig];
P_nonsig_tile = [P_nonsig, P_nonsig];
f8 = figure('Color','w','Name','single-cycle dFF (phase)','Units','centimeters','Position',[2 2 26 12]);
draw_cycle_heatmap(subplot(1,2,1), P_sig_tile,    phs_tile, Pspk_sig,    gclim_cycP, tick_color, tick_lw, spike_tick, ...
    sprintf('SIG  %d spike-cycles', size(P_sig,1)),    'phase (rad)', 'phase');
draw_cycle_heatmap(subplot(1,2,2), P_nonsig_tile, phs_tile, Pspk_nonsig, gclim_cycP, tick_color, tick_lw, spike_tick, ...
    sprintf('NON-SIG  %d spike-cycles', size(P_nonsig,1)), 'phase (rad)', 'phase');
sgtitle('single-cycle dF/F, spike-containing cycles only, sorted by nearest-onset spike phase');

%% ============ FIG 9: DUPLICATING histogram (every cycle counts) ============
% Same accounting as Fig 5 (each spike contributes once per cycle window it
% falls in), shown as a 1D pooled bar histogram. Should display the
% recurring breath rhythm across the +/-3 s window.
psth_sig    = [];
psth_nonsig = [];
sig_idx_v    = find(sig);
nonsig_idx_v = find(nonsig);
for k = sig_idx_v(:)'
    if ~isempty(cyc_spk{k}), psth_sig = [psth_sig; vertcat(cyc_spk{k}{:})]; end %#ok<AGROW>
end
for k = nonsig_idx_v(:)'
    if ~isempty(cyc_spk{k}), psth_nonsig = [psth_nonsig; vertcat(cyc_spk{k}{:})]; end %#ok<AGROW>
end
f9 = figure('Color','w','Name','spike hist (time, duplicating PSTH)','Units','centimeters','Position',[2 2 26 12]);
draw_spike_hist_time(subplot(1,2,1), psth_sig,    peth_edges, peth_ctrs, sprintf('SIG  %d (spike,cycle) pairs', numel(psth_sig)));
draw_spike_hist_time(subplot(1,2,2), psth_nonsig, peth_edges, peth_ctrs, sprintf('NON-SIG  %d (spike,cycle) pairs', numel(psth_nonsig)));
sgtitle('breath-triggered PSTH (every cycle counts every nearby spike)');

%% ============ FIG 10: NON-DUPLICATING per-ROI heatmap (nearest-peak only) ============
% Same accounting as Fig 3 histogram (each spike counted once, at its
% nearest breath peak). Per-ROI histcounts of nearest-peak dt, stacked.
T_nearest_mat = zeros(nROI, numel(t_grid));
for ii = 1:nROI
    if ~isempty(spk_t_pool{ii})
        T_nearest_mat(ii,:) = histcounts(spk_t_pool{ii}, peth_edges);
    end
end
gclim_nearest = clamp_from(prcLim, climManual, T_nearest_mat(:));
f10 = figure('Color','w','Name','spike-count heatmap (time, nearest-peak)','Units','centimeters','Position',[2 2 26 12]);
draw_heatmap(subplot(1,2,1), T_nearest_mat(ord_sig,:),    t_grid, sprintf('SIG  n=%d', sum(sig)),     'time (s)', gclim_nearest, 'spike count (nearest-peak)', 'time');
draw_heatmap(subplot(1,2,2), T_nearest_mat(ord_nonsig,:), t_grid, sprintf('NON-SIG  n=%d', sum(nonsig)), 'time (s)', gclim_nearest, 'spike count (nearest-peak)', 'time');
sgtitle('per-ROI nearest-peak spike count heatmap (each spike counted once)');

%% ---- save ----
out_pairs = {f1,'dFF_heatmap_time'; f2,'dFF_heatmap_phase'; ...
             f3,'spike_hist_time';  f4,'spike_hist_phase'; ...
             f5,'spike_rate_heatmap_time'; f6,'spike_rate_heatmap_phase'; ...
             f7,'cycle_dFF_sorted_time';   f8,'cycle_dFF_sorted_phase'; ...
             f9,'spike_hist_time_PSTH';    f10,'spike_count_heatmap_nearest_time'};
for k = 1:size(out_pairs,1)
    fh = out_pairs{k,1}; nm = out_pairs{k,2};
    exportgraphics(fh, fullfile(outDir,[nm '.png']), 'Resolution',200,'BackgroundColor','white');
    exportgraphics(fh, fullfile(outDir,[nm '.pdf']), 'ContentType','vector','BackgroundColor','white');
end
save(fullfile(outDir,'dFF_heatmap_data.mat'), ...
     'T_mat','P_mat','ST_mat','SP_mat','spk_t_pool','spk_ph_pool', ...
     'T_sig','T_nonsig','P_sig','P_nonsig','Tspk_sig','Tspk_nonsig','Pspk_sig','Pspk_nonsig', ...
     't_grid','phase_1cyc','sig','labels','PP','confC','gclim','gclim_event','gclim_cyc','gclim_cycP');
fprintf('Saved 10 figures + .mat to %s\n', outDir);

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
        M_all  = [M_all;  M_cell{k}];           %#ok<AGROW>
        dt_all = [dt_all; dt_cell{k}(:)];       %#ok<AGROW>
        spk_all = [spk_all; spk_cell{k}(:)];    %#ok<AGROW>
    end
    if isempty(M_all)
        M_sorted = M_all; dt_sorted = dt_all; spk_sorted = spk_all; return;
    end
    dt_sort = dt_all; dt_sort(isnan(dt_sort)) = Inf;
    [~, ord] = sort(dt_sort, 'ascend');
    M_sorted   = M_all(ord, :);
    dt_sorted  = dt_all(ord);
    spk_sorted = spk_all(ord);
end

function draw_heatmap(ax, M, xax, ttl, xlab, clim, cbLabel, mode)
    imagesc(ax, xax, 1:size(M,1), M); hold(ax,'on');
    colormap(ax, flipud(gray(256))); caxis(ax, clim);
    set(ax,'YDir','reverse'); axis(ax,'tight');
    if strcmp(mode,'phase')
        xline(ax, 2*pi, 'k--', 'LineWidth',0.8); xlim(ax,[0 4*pi]);
        set(ax,'XTick',[0 pi 2*pi 3*pi 4*pi], 'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
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
    xline(ax, 0, 'k--', 'LineWidth',1);
    xlim(ax, [ctrs(1) ctrs(end)]);
    xlabel(ax,'time from inspiration (s)'); ylabel(ax,'spike count');
    title(ax, ttl); box(ax,'on');
end

function draw_spike_hist_phase(ax, ph_pool, edges, ctrs, ttl)
    hold(ax,'on');
    if ~isempty(ph_pool)
        cnt   = histcounts(ph_pool, edges);
        prob  = cnt / max(sum(cnt),1);
        mu_fit = angle(mean(exp(1i * ph_pool)));
        a_fit  = mean(prob);
        b_fit  = 2 * mean(prob .* cos(ctrs - mu_fit));
        ctrs_dup = [ctrs, ctrs + 2*pi];
        prob_dup = [prob, prob];
        bar(ax, ctrs_dup, prob_dup, 1, 'FaceColor',[0.3 0.3 0.8], 'FaceAlpha',0.85, 'EdgeColor','none');
        x_fit = linspace(0, 4*pi, 500);
        y_fit = a_fit + b_fit * cos(x_fit - mu_fit);
        plot(ax, x_fit, y_fit, 'r-', 'LineWidth',1.8);
    end
    xline(ax, 2*pi, 'k--', 'LineWidth',1);
    xlim(ax, [0 4*pi]);
    set(ax,'XTick',[0 pi 2*pi 3*pi 4*pi], 'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
    xlabel(ax,'breath phase (rad)'); ylabel(ax,'P(spike)');
    title(ax, ttl); box(ax,'on');
end

function draw_cycle_heatmap(ax, M, xax, spk_cell, clim, tickCol, tickLW, doTicks, ttl, xlab, mode)
% Per-cycle stacked dFF heatmap with optional blue spike ticks per row.
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
                if strcmp(mode,'phase')
                    plot(ax, [sp(s)       sp(s)],       [r-0.4 r+0.4], '-', 'Color', tickCol, 'LineWidth', tickLW);
                    plot(ax, [sp(s)+2*pi  sp(s)+2*pi],  [r-0.4 r+0.4], '-', 'Color', tickCol, 'LineWidth', tickLW);
                else
                    plot(ax, [sp(s) sp(s)], [r-0.4 r+0.4], '-', 'Color', tickCol, 'LineWidth', tickLW);
                end
            end
        end
    end

    if strcmp(mode,'phase')
        xline(ax, 2*pi, 'k--', 'LineWidth',0.8); xlim(ax,[0 4*pi]);
        set(ax,'XTick',[0 pi 2*pi 3*pi 4*pi], 'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
    else
        xline(ax, 0, 'k--', 'LineWidth',0.8);
    end
    xlabel(ax, xlab); ylabel(ax, 'cycle (sorted by dt)');
    title(ax, ttl); cb = colorbar(ax); cb.Label.String = '\DeltaF/F';
end
