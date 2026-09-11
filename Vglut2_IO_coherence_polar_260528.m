% Vglut2_IO_coherence_polar_260528.m
% -----------------------------------------------------------------------
%  Single-group preliminary analysis for the Vglut2/1124/IO folder. Same
%  math as Ventral_surface_coherence_polar_260528 + dFF_heatmaps but:
%    - no sig / non-sig split (all ROIs treated as one black cloud)
%    - one panel per figure (not side-by-side)
%    - nothing saved (figures only)
%
%  Walks every ca_spike_data.mat under <root>, includes ROIs with
%  >= minSpikes spikes, then plots:
%    Fig 1   polar coherence (one black cloud)
%    Fig 2   dF/F heatmap, breath-triggered, TIME (ROI mean, sorted by pref. phase)
%    Fig 3   dF/F heatmap, breath-phase (tiled [0, 4pi])
%    Fig 4   pooled spike-count histogram, TIME (nearest-peak)
%    Fig 5   pooled spike-phase histogram + cosine fit
%    Fig 6   spike-rate Hz heatmap, TIME
%    Fig 7   spike-rate Hz heatmap, PHASE
%    Fig 8   single-cycle dF/F TIME, spike-cycles only, sorted by dt-nearest
%    Fig 9   single-cycle dF/F PHASE, spike-cycles only, sorted by min spike phase
%    Fig 10  pooled PSTH (duplicating: every spike counted per cycle window)
%    Fig 11  per-ROI nearest-peak spike-count heatmap (non-duplicating)
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));
addpath(genpath(fullfile(scriptDir, 'chronux_2_12')));

%% ===================== USER-EDITABLE =================================
root = 'D:\Ventral_surface_summary\Vglut2\1124\IO';
dot_color = [0 0 0];

% coherence
nDrop           = 30;
fallback_fps    = 30;
minSpikes       = 2;
TW              = 5;
alpha_sig       = 0.001;
f_breath_search = [0.2 4];
fwhm_factor     = 0.6;
min_bw          = 0.05;
fmin            = 0.05;
fmax            = 15;

% heatmaps / histograms
t_win        = 3;
t_dt         = 0.05;
nPhaseBins   = 100;
nHistBins    = 24;
prcLim       = [1 99.5];
climManual   = [];
spike_tick   = true;
tick_color   = [0 0.4 1];
tick_lw      = 0.8;
% =====================================================================

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

%% ---- walk folder, collect coherence + heatmap matrices per ROI ----
PP = init_coll(); labels = {}; confC = NaN;
T_mat = []; P_mat = []; ST_mat = []; SP_mat = [];
spk_t_pool = {}; spk_ph_pool = {};
cyc_T = {}; cyc_dt = {}; cyc_spk = {};
cyc_P = {}; cyc_ph_first = {}; cyc_spkph = {};

t_grid     = -t_win:t_dt:t_win;
nW         = numel(t_grid);
phase_1cyc = linspace(0, 2*pi, nPhaseBins+1); phase_1cyc = phase_1cyc(1:end-1);

allMat = dir(fullfile(root, '**', 'ca_spike_data.mat'));
fprintf('Found %d recording(s) under %s\n', numel(allMat), root);

for kk = 1:numel(allMat)
    folderPath = allMat(kk).folder;
    [~, recName] = fileparts(folderPath);
    try
        bp = dir(fullfile(folderPath, '*breath_peak_data.mat'));
        df = dir(fullfile(folderPath, '*_ch1_dFF.mat'));
        if isempty(bp) || isempty(df), fprintf('  skip (no breath/dFF): %s\n', recName); continue; end

        fps = detect_session_fps(folderPath, fallback_fps);
        CA  = load(fullfile(folderPath, 'ca_spike_data.mat'));
        D   = load(fullfile(df(1).folder, df(1).name), 'dFF');
        BP  = load(fullfile(bp(1).folder, bp(1).name));
        dff_all = double(D.dFF);
        nROI = numel(CA.roi_spikes);
        nCa  = numel(CA.roi_spikes(1).spike_train);

        nB = numel(BP.breath);
        if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==nB
            evt = double(BP.insp_onsets_train(:) ~= 0);
        else
            evt = zeros(nB,1); oi = round(BP.insp_onset_idx(:));
            evt(oi(oi>=1 & oi<=nB)) = 1;
        end
        evt(1:min(nDrop,numel(evt))) = [];
        bw = detrend(double(BP.breath(:)));
        bw(1:min(nDrop,numel(bw))) = []; bw = bw - mean(bw);
        T = min([numel(evt), numel(bw), nCa, size(dff_all,1)]);
        evt = evt(1:T); bw = bw(1:T); dff_all = dff_all(1:T,:);
        if sum(evt) < 3, continue; end
        on = find(evt>0);
        win = round(t_win*fps); T_local = (-win:win)/fps;

        %% breath PSD -> coherence band
        pB.Fs=fps; pB.tapers=[TW,2*TW-1]; pB.pad=0;
        pB.fpass=[fmin,min(fmax,fps/2)]; pB.err=0;
        [Sb,fb] = mtspectrumc(bw, pB); Sb=Sb(:); fb=fb(:);
        m = fb>=f_breath_search(1) & fb<=f_breath_search(2);
        [~,rl]=max(Sb(m)); ip=find(m,1)+rl-1; f_pk=fb(ip);
        h=Sb(ip)/2; lo=ip; while lo>1&&Sb(lo)>h, lo=lo-1; end
        hi=ip;            while hi<numel(fb)&&Sb(hi)>h, hi=hi+1; end
        f_fwhm=[max(fb(lo),f_breath_search(1)), min(fb(hi),f_breath_search(2))];
        bwd=max(diff(f_fwhm)*fwhm_factor, min_bw);
        band=[max(f_pk-bwd/2,f_breath_search(1)), min(f_pk+bwd/2,f_breath_search(2))];

        pc.Fs=fps; pc.tapers=[TW,2*TW-1]; pc.pad=0;
        pc.fpass=band; pc.err=[2,alpha_sig];
        br_p = evt - mean(evt);

        nInc = 0;
        for rid = 1:nROI
            spk_train = double(CA.roi_spikes(rid).spike_train(:));
            if numel(spk_train)<T, spk_train(end+1:T) = 0; end
            spk_train = spk_train(1:T);
            if sum(spk_train) < minSpikes, continue; end
            nInc = nInc + 1;

            % coherence
            [PP, confC] = add_coh(PP, br_p, spk_train - mean(spk_train), pc, band, confC);
            labels{end+1} = sprintf('%s#%d', recName, rid); %#ok<SAGROW>

            tr = dff_all(:, rid);
            ii = numel(labels);

            % preallocate row in matrices
            T_mat(ii,:)  = nan(1, nW);  %#ok<SAGROW>
            P_mat(ii,:)  = nan(1, nPhaseBins); %#ok<SAGROW>
            ST_mat(ii,:) = nan(1, nW);  %#ok<SAGROW>
            SP_mat(ii,:) = nan(1, nPhaseBins); %#ok<SAGROW>
            spk_t_pool{ii,1}  = []; spk_ph_pool{ii,1} = []; %#ok<SAGROW>
            cyc_T{ii,1} = []; cyc_dt{ii,1} = []; cyc_spk{ii,1} = {}; %#ok<SAGROW>
            cyc_P{ii,1} = []; cyc_ph_first{ii,1} = []; cyc_spkph{ii,1} = {}; %#ok<SAGROW>

            % time-triggered dF/F + spike rate
            spike_idx = find(spk_train > 0);
            valid_on = on(on - win >= 1 & on + win <= T);
            if ~isempty(valid_on)
                E  = zeros(numel(valid_on), 2*win+1);
                ES = zeros(numel(valid_on), 2*win+1);
                for k = 1:numel(valid_on)
                    E (k,:) = tr       (valid_on(k)-win : valid_on(k)+win);
                    ES(k,:) = spk_train(valid_on(k)-win : valid_on(k)+win);
                end
                T_mat (ii,:) = interp1(T_local, mean(E, 1),       t_grid, 'linear', NaN);
                ST_mat(ii,:) = interp1(T_local, mean(ES,1) * fps, t_grid, 'linear', NaN);
            end

            % phase-triggered dF/F + spike rate
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

            % pooled nearest-peak spike times
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

            % single-cycle dFF -- time
            Tsnips = []; Tdt = []; Tspk = {};
            for e = 1:numel(valid_on)
                c = valid_on(e);
                spk_in_win = spike_idx(spike_idx >= c - win & spike_idx <= c + win);
                if isempty(spk_in_win), continue; end
                snip = tr(c-win : c+win);
                dt_bef = Inf; dt_aft = Inf;
                bef = spike_idx(spike_idx <  c); aft = spike_idx(spike_idx >= c);
                if ~isempty(bef), dt_bef = (c - bef(end))/fps; end
                if ~isempty(aft), dt_aft = (aft(1) - c)/fps; end
                if dt_bef <= dt_aft, dt_sign = -dt_bef; else, dt_sign = dt_aft; end
                Tsnips(end+1,:) = interp1(T_local, snip, t_grid, 'linear', NaN); %#ok<AGROW>
                Tdt(end+1,1)    = dt_sign; %#ok<AGROW>
                Tspk{end+1,1}   = (spk_in_win - c) / fps; %#ok<AGROW>
            end
            cyc_T{ii} = Tsnips; cyc_dt{ii} = Tdt; cyc_spk{ii} = Tspk;

            % single-cycle dFF -- phase
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
                Psnips(end+1,:) = ph_row; %#ok<AGROW>
                Pfirst(end+1,1) = min(spk_ph_cyc); %#ok<AGROW>
                Pspk{end+1,1}   = spk_ph_cyc; %#ok<AGROW>
            end
            cyc_P{ii} = Psnips; cyc_ph_first{ii} = Pfirst; cyc_spkph{ii} = Pspk;
        end
        fprintf('[%d] %-50s band [%.2f %.2f] Hz  %d/%d ROI\n', ...
                kk, recName, band(1), band(2), nInc, nROI);
    catch ME
        warning('  ERROR %s: %s', recName, ME.message);
    end
end

if isempty(PP.r), error('No ROIs collected.'); end
nROI = numel(labels);
fprintf('\nTotal: %d ROIs, confC=%.3f, %d >= confC.\n', nROI, confC, sum(PP.r >= confC));

%% ---- sort rows by preferred phase (ascending) ----
[~, ord] = sort(PP.th);

% phase tiling [0, 4pi]
P_tile   = [P_mat,  P_mat ];
SP_tile  = [SP_mat, SP_mat];
phs_tile = [phase_1cyc, phase_1cyc + 2*pi];

% color clamps
gclim       = clamp_from(prcLim, climManual, [T_mat(:);  P_mat(:) ]);
gclim_event = clamp_from(prcLim, climManual, [ST_mat(:); SP_mat(:)]);

%% ============ FIG 1: POLAR COHERENCE (one black cloud) ============
f1 = figure('Color','w','Name','Vglut2/IO coherence polar','Units','centimeters','Position',[2 2 14 13]);
set(f1,'DefaultAxesFontSize',8,'DefaultTextFontSize',8);
ax = polaraxes(f1,'Position',[0.10 0.10 0.78 0.78]); hold(ax,'on');
thc = linspace(0,2*pi,360);
polarplot(ax, thc, repmat(confC,1,360), 'k--','LineWidth',1);
for k = 1:numel(PP.r)
    if isnan(PP.th(k)) || isnan(PP.r(k)), continue; end
    if PP.r(k) >= confC
        polarplot(ax, [PP.th(k) PP.th(k)], [PP.rlo(k) PP.rhi(k)], '-', 'Color', dot_color, 'LineWidth', 1.0);
        if ~isnan(PP.dphi(k))
            arc = linspace(PP.th(k)-PP.dphi(k), PP.th(k)+PP.dphi(k), 30);
            polarplot(ax, arc, PP.r(k)*ones(size(arc)), '-', 'Color', dot_color, 'LineWidth', 1.0);
        end
    end
    polarplot(ax, PP.th(k), PP.r(k), 'o', ...
              'MarkerFaceColor', dot_color, 'MarkerEdgeColor','k', ...
              'MarkerSize', 6, 'LineWidth', 0.4);
end
ax.RLim=[0 1]; ax.ThetaZeroLocation='right'; ax.ThetaDir='counterclockwise';
ax.RAxisLocation=180; ax.FontSize=8;
title(ax, sprintf('Vglut2/1124/IO   N=%d ROI   confC=%.2f   \\alpha=%.3f', ...
                  nROI, confC, alpha_sig), 'Interpreter','none');

%% ============ FIG 2: dFF heatmap TIME ============
f2 = figure('Color','w','Name','dFF heatmap (time)','Units','centimeters','Position',[2 2 14 12]);
draw_heatmap(gca, T_mat(ord,:), t_grid, sprintf('all ROIs (n=%d), sorted by pref. phase', nROI), 'time (s)', gclim, '\DeltaF/F', 'time');

%% ============ FIG 3: dFF heatmap PHASE ============
f3 = figure('Color','w','Name','dFF heatmap (phase)','Units','centimeters','Position',[2 2 14 12]);
draw_heatmap(gca, P_tile(ord,:), phs_tile, sprintf('all ROIs (n=%d), tiled [0, 4\\pi]', nROI), 'phase (rad)', gclim, '\DeltaF/F', 'phase');

%% ============ FIG 4: pooled nearest-peak spike-count histogram TIME ============
spk_t_all = vertcat(spk_t_pool{:});
dt_bin     = t_grid(2) - t_grid(1);
peth_edges = [t_grid(:) - dt_bin/2; t_grid(end) + dt_bin/2];
peth_ctrs  = t_grid;
f4 = figure('Color','w','Name','spike hist (time, nearest)','Units','centimeters','Position',[2 2 14 10]);
draw_spike_hist_time(gca, spk_t_all, peth_edges, peth_ctrs, sprintf('all ROIs (n=%d), %d spikes', nROI, numel(spk_t_all)));

%% ============ FIG 5: pooled spike-phase histogram + cosine fit ============
spk_ph_all = vertcat(spk_ph_pool{:});
edges_hist = linspace(0, 2*pi, nHistBins+1);
ctrs_hist  = (edges_hist(1:end-1) + edges_hist(2:end))/2;
f5 = figure('Color','w','Name','spike hist (phase)','Units','centimeters','Position',[2 2 14 10]);
draw_spike_hist_phase(gca, spk_ph_all, edges_hist, ctrs_hist, sprintf('all ROIs (n=%d), %d spikes', nROI, numel(spk_ph_all)));

%% ============ FIG 6: spike rate (Hz) heatmap TIME ============
f6 = figure('Color','w','Name','spike rate (time)','Units','centimeters','Position',[2 2 14 12]);
draw_heatmap(gca, ST_mat(ord,:), t_grid, sprintf('all ROIs (n=%d)', nROI), 'time (s)', gclim_event, 'spike rate (Hz)', 'time');

%% ============ FIG 7: spike rate (Hz) heatmap PHASE ============
f7 = figure('Color','w','Name','spike rate (phase)','Units','centimeters','Position',[2 2 14 12]);
draw_heatmap(gca, SP_tile(ord,:), phs_tile, sprintf('all ROIs (n=%d), tiled [0, 4\\pi]', nROI), 'phase (rad)', gclim_event, 'spike rate (Hz)', 'phase');

%% ============ FIG 8: single-cycle dFF TIME, spike-cycles only, sorted by dt-nearest ============
[T_all,    ~, Tspk_all] = pool_and_sort(cyc_T, cyc_dt, cyc_spk, true(nROI,1));
gclim_cyc = clamp_from(prcLim, climManual, T_all(:));
f8 = figure('Color','w','Name','single-cycle dFF (time)','Units','centimeters','Position',[2 2 14 12]);
draw_cycle_heatmap(gca, T_all, t_grid, Tspk_all, gclim_cyc, tick_color, tick_lw, spike_tick, ...
    sprintf('%d spike-cycles', size(T_all,1)), 'time (s)', 'time');

%% ============ FIG 9: single-cycle dFF PHASE, sorted by min spike phase ============
[P_all,    ~, Pspk_all] = pool_and_sort(cyc_P, cyc_ph_first, cyc_spkph, true(nROI,1));
gclim_cycP = clamp_from(prcLim, climManual, P_all(:));
P_all_tile = [P_all, P_all];
f9 = figure('Color','w','Name','single-cycle dFF (phase)','Units','centimeters','Position',[2 2 14 12]);
draw_cycle_heatmap(gca, P_all_tile, phs_tile, Pspk_all, gclim_cycP, tick_color, tick_lw, spike_tick, ...
    sprintf('%d spike-cycles', size(P_all,1)), 'phase (rad)', 'phase');

%% ============ FIG 10: PSTH duplicating histogram TIME ============
psth_all = [];
for k = 1:nROI
    if ~isempty(cyc_spk{k}), psth_all = [psth_all; vertcat(cyc_spk{k}{:})]; end %#ok<AGROW>
end
f10 = figure('Color','w','Name','PSTH (time, duplicating)','Units','centimeters','Position',[2 2 14 10]);
draw_spike_hist_time(gca, psth_all, peth_edges, peth_ctrs, ...
    sprintf('all ROIs (n=%d), %d (spike,cycle) pairs', nROI, numel(psth_all)));

%% ============ FIG 11: per-ROI nearest-peak count heatmap TIME (non-duplicating) ============
T_nearest_mat = zeros(nROI, numel(t_grid));
for ii = 1:nROI
    if ~isempty(spk_t_pool{ii})
        T_nearest_mat(ii,:) = histcounts(spk_t_pool{ii}, peth_edges);
    end
end
gclim_nearest = clamp_from(prcLim, climManual, T_nearest_mat(:));
f11 = figure('Color','w','Name','nearest-peak count heatmap (time)','Units','centimeters','Position',[2 2 14 12]);
draw_heatmap(gca, T_nearest_mat(ord,:), t_grid, sprintf('all ROIs (n=%d), nearest-peak', nROI), ...
    'time (s)', gclim_nearest, 'spike count', 'time');

fprintf('11 figures rendered. Nothing saved.\n');

%% ========================= LOCAL FUNCTIONS ==========================
function C = init_coll()
    C = struct('th',[],'r',[],'rlo',[],'rhi',[],'dphi',[]);
end

function [C, confC] = add_coh(C, x, y, pc, band, confC)
    [~, Cxy, phi, ~,~,~, f, cC, phistd, Cerr] = coherencyc(x, y, pc);
    if isnan(confC), confC = cC; end
    f = f(:); mb = f>=band(1) & f<=band(2); if ~any(mb), mb = true(size(f)); end
    C.th(end+1,1)   = angle(mean(exp(1i*(-phi(mb)))));
    C.r(end+1,1)    = mean(Cxy(mb));
    C.rlo(end+1,1)  = max(0, mean(Cerr(1,mb)));
    C.rhi(end+1,1)  = min(1, mean(Cerr(2,mb)));
    C.dphi(end+1,1) = 1.96*mean(phistd(mb));
end

function c = clamp_from(prcLim, climManual, pool)
    if ~isempty(climManual), c = climManual; return; end
    f = pool(isfinite(pool));
    if isempty(f), c = [0 1]; return; end
    c = prctile(f, prcLim);
    if ~(c(2) > c(1)), c = [min(f) max(f)]; end
end

function [M_sorted, dt_sorted, spk_sorted] = pool_and_sort(M_cell, dt_cell, spk_cell, mask)
    keep = find(mask);
    M_all = []; dt_all = []; spk_all = {};
    for k = keep(:)'
        if isempty(M_cell{k}), continue; end
        M_all  = [M_all;  M_cell{k}];       %#ok<AGROW>
        dt_all = [dt_all; dt_cell{k}(:)];   %#ok<AGROW>
        spk_all = [spk_all; spk_cell{k}(:)];%#ok<AGROW>
    end
    if isempty(M_all)
        M_sorted = M_all; dt_sorted = dt_all; spk_sorted = spk_all; return;
    end
    dt_sort = dt_all; dt_sort(isnan(dt_sort)) = Inf;
    [~, ord] = sort(dt_sort, 'ascend');
    M_sorted = M_all(ord, :); dt_sorted = dt_all(ord); spk_sorted = spk_all(ord);
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
