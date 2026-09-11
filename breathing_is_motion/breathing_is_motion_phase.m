%% breathing_is_motion_phase.m
% Breathing-phase modulation for a single ROI, in the manner of
% Batch_breath_trig_combined_260328.m, driven by shift_breath_master.mat
% (output of breathing_is_motion_analysis.m).
%
% Breathing phase = linear 0 -> 2*pi between consecutive inspiration onsets.
% Three plots (rows sorted by calcium peaks; dFF heatmaps in grayscale):
%   1) dFF snippets around inspiration (TIME), sorted by signed dt to nearest spike
%   2) per-cycle dFF resampled to phase (PHASE), sorted by spike phase
%   3) spike probability P(Ca spike | phase) over [0,4*pi] + cosine fit
%        a + b*cos(phi - mu)  via circular mean:  mu = preferred phase,
%        b = 2*B*R (R = resultant length), Rayleigh p for non-uniformity.

close all; clear; clc;

%% ---------------- user params ----------------
folderPath  = "C:\Users\Admin\Desktop\260521_ChAT_soma_G8s\roi5_7x_x-1200y200z-30_3000f_23lp_00001";
roi         = 3;        % ROI to analyze (needs calcium spikes)
win_sec     = 1.5;      % half-window (s) for the time-domain snippets
nPhaseGrid  = 60;       % phase columns for the phase heatmap (resampling)
nBins       = 20;       % phase bins for the spike-probability histogram
prcLim      = [1 100];   % percentile clamp for heatmap color scale
% ----------------------------------------------

%% Load master
mhits = dir(fullfile(folderPath, 'shift_breath_master.mat'));
if isempty(mhits)
    error(['No shift_breath_master.mat in:\n  %s\n' ...
           'Run breathing_is_motion_analysis.m first.'], folderPath);
end
L = load(fullfile(folderPath, mhits(1).name));
fps    = L.fps;
dFF    = L.dFF;
onsets = sort(L.insp_onsets(:));
T = size(dFF,1);

if roi < 1 || roi > size(dFF,2)
    error('roi=%d out of range (1..%d).', roi, size(dFF,2));
end
y = dFF(:, roi);

ca = [];
if isfield(L,'ca_idx') && numel(L.ca_idx) >= roi
    ca = L.ca_idx{roi}(:);
    ca = ca(ca >= 1 & ca <= T);
end
fprintf('ROI %d: %d inspiration onsets, %d calcium spikes\n', roi, numel(onsets), numel(ca));
if isempty(ca)
    warning(['ROI %d has no calcium spikes - spike-based sorting and the ' ...
             'P(spike|phase) fit are not meaningful. Pick a spiking ROI.'], roi);
end

%% Breathing phase (0..2*pi), linear peak-to-peak
phase = nan(T,1);
for k = 1:numel(onsets)-1
    a = onsets(k); b = onsets(k+1);
    idx = a:b-1;
    phase(idx) = 2*pi * (idx - a) / (b - a);
end
nCyc = numel(onsets) - 1;

%% ===== Panel 1 data: TIME snippets around inspiration peaks =====
win   = round(win_sec * fps);
t_snip = (-win:win) / fps;
valid_on = onsets(onsets > win & onsets <= T - win);
nV = numel(valid_on);
snip = nan(nV, 2*win+1);
for b = 1:nV
    c = valid_on(b);
    snip(b,:) = y(c-win : c+win);
end

% sort key: nearest spike = min absolute dt (s) to a spike, regardless of side
dt_near = inf(nV,1);
for b = 1:nV
    c = valid_on(b);
    before = ca(ca <  c);  after = ca(ca >= c);
    dt_bef = inf; dt_aft = inf;
    if ~isempty(before), dt_bef = (c - before(end)) / fps; end
    if ~isempty(after),  dt_aft = (after(1) - c)    / fps; end
    dt_near(b) = min(dt_bef, dt_aft);
end
[~, ord_t] = sort(dt_near, 'ascend');
snip_s = snip(ord_t, :);

% spikes within +/- win of each valid onset (sec rel to onset) for the time raster
ca_rel = cell(nV,1);
for b = 1:nV
    c = valid_on(b);
    ca_rel{b} = (ca(ca >= c-win & ca <= c+win) - c) / fps;
end

%% ===== Panel 2 data: per-cycle dFF resampled to phase (batch-style) =====
phgrid = (0.5:nPhaseGrid-0.5) / nPhaseGrid * 2*pi;
M_ph = nan(nCyc, nPhaseGrid);
for k = 1:nCyc
    a = onsets(k); b = onsets(k+1);
    seg = y(a:b-1);
    if numel(seg) >= 2
        M_ph(k,:) = interp1(linspace(0,2*pi,numel(seg)), seg, phgrid, 'linear');
    end
end
% keep valid (non all-NaN) cycles; per-cycle spike phases + min-phase sort key
valid_ph_idx   = find(~all(isnan(M_ph),2));
spk_ph_per_cyc = cell(numel(valid_ph_idx),1);
min_spk_ph     = inf(numel(valid_ph_idx),1);
for cc = 1:numel(valid_ph_idx)
    k = valid_ph_idx(cc);
    a = onsets(k); b = onsets(k+1);
    spk = ca(ca >= a & ca < b);
    if ~isempty(spk)
        spk_ph_per_cyc{cc} = 2*pi * (spk - a + 0.5) / (b - a);
        min_spk_ph(cc)     = min(spk_ph_per_cyc{cc});
    end
end
[~, sort_ph_ord] = sort(min_spk_ph, 'ascend');
dff_ph_sorted    = M_ph(valid_ph_idx(sort_ph_ord), :);

%% ===== Panel 3 data: spike probability + circular-mean cosine fit =====
phi_sp = phase(ca);  phi_sp = phi_sp(~isnan(phi_sp));
edges  = linspace(0, 2*pi, nBins+1);
ctrs   = (edges(1:end-1) + edges(2:end)) / 2;
counts = histcounts(phi_sp, edges);
counts_pc = counts / max(nCyc,1);          % spikes per cycle per bin
ctrs_ext  = [ctrs, ctrs + 2*pi];
cnts_ext  = [counts_pc, counts_pc];

if ~isempty(phi_sp)
    r_vec = mean(exp(1i*phi_sp));
    mu    = angle(r_vec);                  % preferred (modulated) phase
    Rlen  = abs(r_vec);                    % resultant length (0..1)
    B_cos = numel(phi_sp) / (nBins * max(nCyc,1));
    A_cos = 2 * B_cos * Rlen;              % cosine amplitude
    th    = linspace(0, 4*pi, 300);
    fit_ext = B_cos + A_cos * cos(th - mu);
    n_sp  = numel(phi_sp);
    p_ray = exp(-n_sp * Rlen^2);           % Rayleigh test (non-uniformity)
else
    mu = NaN; Rlen = NaN; A_cos = NaN; B_cos = NaN; th = []; fit_ext = []; p_ray = NaN;
end

fprintf('ROI %d: preferred phase mu=%.1f deg, R=%.3f, p_Rayleigh=%.4f\n', ...
        roi, rad2deg(mod(mu,2*pi)), Rlen, p_ray);

%% ===== Figure: 3 panels =====
fig = figure('Color','w','Position',[120 120 1500 460]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% shared color limits for both heatmaps (same prcLim over combined dFF)
clim_dff = prctile([snip_s(:); dff_ph_sorted(:)], prcLim);

% --- Panel 1: TIME heatmap (grayscale, sorted by nearest spike) ---
ax1 = nexttile;
imagesc(ax1, t_snip, 1:nV, snip_s); set(ax1,'YDir','reverse');
colormap(ax1, flipud(gray(256))); if clim_dff(2) > clim_dff(1), caxis(ax1, clim_dff); end
hold(ax1,'on'); xline(ax1, 0, 'k--', 'LineWidth', 1); hold(ax1,'off');
xlim(ax1,[-win_sec win_sec]); ylim(ax1,[0.5 nV+0.5]);
xlabel(ax1,'Time from inspiration (s)'); ylabel(ax1,'Breath # (sorted by nearest spike)');
title(ax1, sprintf('ROI#%d  \\DeltaF/F (time)', roi));
cb1 = colorbar(ax1); cb1.Label.String = '\DeltaF/F';

% --- Panel 2: PHASE heatmap over [0,4pi] (grayscale, sorted by spike phase) ---
phase_2cyc = [phgrid, phgrid + 2*pi];
dff_2x     = [dff_ph_sorted, dff_ph_sorted];
nR = size(dff_2x,1);
ax2 = nexttile;
imagesc(ax2, phase_2cyc, 1:nR, dff_2x); set(ax2,'YDir','reverse');
colormap(ax2, flipud(gray(256))); if clim_dff(2) > clim_dff(1), caxis(ax2, clim_dff); end
hold(ax2,'on'); xline(ax2, 2*pi, 'k--', 'LineWidth', 0.8); hold(ax2,'off');
xlim(ax2,[0 4*pi]); ylim(ax2,[0.5 nR+0.5]);
xticks(ax2,[0 pi 2*pi 3*pi 4*pi]); xticklabels(ax2,{'0','\pi','2\pi','3\pi','4\pi'});
xlabel(ax2,'Breathing phase (rad)'); ylabel(ax2,'Cycle (sorted by spike phase)');
title(ax2, sprintf('ROI#%d  \\DeltaF/F (phase)', roi));
cb2 = colorbar(ax2); cb2.Label.String = '\DeltaF/F';

% --- Panel 3: spike probability [0,4pi] + cosine fit ---
ax3 = nexttile; hold(ax3,'on');
bar(ax3, ctrs_ext, cnts_ext, 1, 'FaceColor',[0.8 0.8 0.8], 'EdgeColor',[0.5 0.5 0.5]);
if ~isempty(fit_ext)
    plot(ax3, th, fit_ext, 'b-', 'LineWidth', 2);
    %xline(ax3, mod(mu,2*pi),        'r--', 'LineWidth', 1.5);
    %xline(ax3, mod(mu,2*pi) + 2*pi, 'r--', 'LineWidth', 1.5);
end
xline(ax3, 2*pi, 'k--', 'LineWidth', 1);
xlabel(ax3,'Breathing phase (rad)'); ylabel(ax3,'P(Ca spike | phase)  [spikes/cycle/bin]');
xticks(ax3,[0 pi 2*pi 3*pi 4*pi]); xticklabels(ax3,{'0','\pi','2\pi','3\pi','4\pi'});
xlim(ax3,[0 4*pi]);
ylim([0 0.3]);
title(ax3, sprintf('Spike probability (Phase)'));
box(ax3,'on');

outFig = fullfile(folderPath, sprintf('phase_modulation_ROI%02d.png', roi));
exportgraphics(fig, outFig, 'Resolution',150);
fprintf('Saved: %s\n', outFig);

%% Save results
phase_result = struct('roi',roi,'mu_rad',mod(mu,2*pi),'mu_deg',rad2deg(mod(mu,2*pi)), ...
    'Rlen',Rlen,'A_cos',A_cos,'B_cos',B_cos,'p_rayleigh',p_ray, ...
    'ctrs',ctrs,'counts_per_cycle',counts_pc,'n_spikes',numel(phi_sp), ...
    'win_sec',win_sec,'nPhaseGrid',nPhaseGrid,'nBins',nBins);
outMat = fullfile(folderPath, sprintf('phase_modulation_ROI%02d.mat', roi));
save(outMat, '-struct', 'phase_result');
fprintf('Saved: %s\n', outMat);
