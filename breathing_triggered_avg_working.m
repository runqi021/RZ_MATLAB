%%
close all; clear; clc;

folderPath = "D:\RUNQI\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi3_z20_512x512_3x_2000f_00001";

% Auto-find DLC CSV and SAM output
csv_hits = dir(fullfile(folderPath, '*snapshot_best-20.csv'));
sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));

if isempty(csv_hits) || isempty(sam_hits)
    error('Could not find DLC CSV or SAM output in:\n  %s', folderPath);
end

dlc_csv = fullfile(folderPath, csv_hits(1).name);
SAM     = fullfile(folderPath, sam_hits(1).name);

fprintf('DLC CSV : %s\n', csv_hits(1).name);
fprintf('SAM     : %s\n', sam_hits(1).name);

%% Load F
SAMload = load(SAM);
F_raw = SAMload.F;
% toss first second (system steady state)
F = F_raw;
F(1:30, :) = [];

% Inverted scanner signal to trigger bassler acq (Falling edge, 2P acq start, bassler acq start, synched).
% dFF
dFFout = helper.dFF_RZ(F);
dFF = dFFout.dFF;

%% ------------------- BASIC PARAMS ---------------------------
[T, N_roi] = size(dFF);

fps_img  = 30;   % FPS
fps_breath = fps_img; % stim trigger by flyback
nDrop = 30;      % for breathing video

t_img = (0:T-1)' / fps_img;
roi_ids = arrayfun(@(k) sprintf('%02d', k), 1:N_roi, 'UniformOutput', false);

%% load calcium spike
in = load(fullfile(folderPath, 'ca_spike_data.mat'));
ifSpk = in.ifSpike; roiSpk_id = find(ifSpk == 1);
ca_spk_data = in.roi_spikes;

nSpkROI = numel(roiSpk_id);

ca_spk_train = zeros(T, nSpkROI, 'double');
ca_spk_id = {}; t_spk = {};

for k = 1:nSpkROI
    i = roiSpk_id(k);
    ca_spk_train(:, k) = ca_spk_data(i).spike_train;
    ca_spk_id{k} = ca_spk_data(i).spike_idx;
    t_spk{k} = ca_spk_data(i).spike_t;
end

%% load breath peak
breath_peak = load(fullfile(folderPath, 'breath_peak_data.mat'));

insp_onsets = breath_peak.insp_onset_idx;
Sb = breath_peak.insp_onsets_train;

insp_onsets(insp_onsets<30)=[];
Sb(1:30)=[];

%% ------------------- DLC BREATHING INPUT --------------------
% DLC CSV has 3 header lines: scorer / bodyparts / coords
dlc_raw = readmatrix(dlc_csv, 'NumHeaderLines', 3);

data_dlc = dlc_raw;
data_dlc(1:30, :) = [];
frame_idx = data_dlc(:,1);
t_breath  = frame_idx / fps_breath;

% Build breath signal using dot/coord/combine settings saved in breath_peak_data.
% Column layout: dot1=[2,3,4], dot2=[5,6,7], dot3=[8,9,10], dot4=[11,12,13]
dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
dot_idx  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);

fp          = breath_peak.findpeak_params;
nF_dlc      = size(data_dlc, 1);
nActiveDots = numel(fp.dot_selection);
traces      = NaN(nF_dlc, nActiveDots);

fprintf('Breath signal: %d dots, combine=%s, invert=%d, lk_thr=%.2f\n', ...
    nActiveDots, fp.combine_method, fp.inverted, fp.likelihood_thr);

for d = 1:nActiveDots
    di  = dot_idx.(fp.dot_selection{d});
    xc  = data_dlc(:, dot_cols(di,1));
    yc  = data_dlc(:, dot_cols(di,2));
    pc  = data_dlc(:, dot_cols(di,3));
    bad = pc < fp.likelihood_thr;

    switch fp.coord_types{di}
        case 'x',         sig =  xc;
        case 'y',         sig =  yc;
        case '-x',        sig = -xc;
        case '-y',        sig = -yc;
        case 'magnitude', sig = sqrt(xc.^2 + yc.^2);
        otherwise,        sig =  xc;
    end
    sig(bad) = NaN;
    if sum(~isnan(sig)) >= 2
        sig = fillmissing(sig, 'linear', 'EndValues', 'nearest');
    end
    traces(:,d) = sig;
    fprintf('  %s -> %s\n', fp.dot_selection{d}, fp.coord_types{di});
end

switch fp.combine_method
    case 'sum',  breath = sum(traces,  2, 'omitnan');
    otherwise,   breath = mean(traces, 2, 'omitnan');   % 'average'
end

if fp.inverted, breath = -breath; end

breath = detrend(breath);
breath = (breath - mean(breath)) / std(breath);

fprintf('Breathing trace: %d samples after dropping first %d frames\n', ...
        numel(breath), nDrop);


%% Chronux spectrum
TW_spec = 4;
params_spec.Fs     = fps_img;
params_spec.tapers = [TW_spec, 2*TW_spec - 1];
params_spec.pad    = 0;
params_spec.fpass  = [0.1 15];
params_spec.err    = [2 0.05];

[Sk_raw, fk_raw] = mtspectrumc(breath, params_spec);
figure('Color','White');
plot(fk_raw, Sk_raw,                'k',  'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([0 10]);
grid('on');
axis square;

% Save breathing spectrum
outFig = fullfile(folderPath, 'breathing_spectrum_raw.png');
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved breathing spectrum:\n  %s\n', outFig);
%% ------------------- FIGURE 1: BREATH + STACKED dFF (Top 10 ROIs) --------
T_img = numel(breath);

% global amplitude
amp = max(dFF(:)) - min(dFF(:));
if amp == 0, amp = 1; end; gap = 0.8 * amp;

% build stacked matrix (highest at top)
Y = zeros(T_img, N_roi);
for k = 1:N_roi
    Y(:,k) = dFF(:,k) + gap*(N_roi - k);
end

% y-ticks and labels
yt  = 0:gap:gap*(N_roi-1);                    % positions (bottom -> top)
ytl = roi_ids(end:-1:1);                      % labels (bottom = lowest ROI)

figure('Position',[200 200 900 700]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% ---- TOP PANEL: breathing + breathing events (Sb) ----
ax1 = nexttile(1);
plot(ax1, t_img, -breath, 'k', 'LineWidth', 1); hold(ax1,'on');

% --- breathing spike overlay (Sb) ---
% assume Sb is [T_img x 1] or [1 x T_img], binary 0/1
idx_b = find(Sb(:) ~= 0);
if ~isempty(idx_b)
    t_b = t_img(idx_b);
    y_b = breath(idx_b);
    plot(ax1, t_b, -y_b, 'ro', 'MarkerFaceColor','r', ...
         'MarkerSize', 3);   % breathing event markers
end
xlim(ax1, [t_img(1) t_img(end)])
xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Breathing (z)');
title(ax1, 'Breathing trace (DLC) with detected inspirations');
grid(ax1,'off');
hold(ax1,'off');

% BOTTOM PANEL: dFF
ax2 = nexttile(2, [2 1]);
hold(ax2,'on');

for kk = 1:N_roi
    % color spiking ROIs orange
    if ismember(kk, roiSpk_id)
        trace_color = [1 0 0];
    else
        trace_color = [0 0 0];
    end
    plot(ax2, t_img, Y(:,kk), 'Color', trace_color, 'LineWidth', 0.5);
end

xlim(ax2, [t_img(1) t_img(end)]);
ylim(ax2, [-gap, gap*(N_roi-1)+gap]);

yticks(ax2, yt);
yticklabels(ax2, ytl);

xlabel(ax2, 'Time (s)');
ylabel(ax2, 'ROI#');
title(ax2, '\DeltaF/F');
grid(ax2,'off');

% SCALE BAR
dFF_scale = 0.2;     % 0.2 dF/F vertical bar

% position vertical bar near right
x0 = 1*(t_img(end));
y0 = gap*(N_roi-1) + 0.2*gap;

plot(ax2, [x0 x0], [dFF_scale 2*dFF_scale], 'k', 'LineWidth', 2);
text(ax2, x0+0.5, dFF_scale + 0.05*dFF_scale, ...
    sprintf('%.2f \\DeltaF/F', dFF_scale), ...
    'VerticalAlignment','bottom','HorizontalAlignment','center', ...
    'FontSize', 10);
hold(ax2,'off');

% Save breathing spectrum
outFig = fullfile(folderPath, 'breath_N_dFF.png');
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved breathing and dFF traces:\n  %s\n', outFig);

% ---- Figure 1b: middle 10 s zoom ----
mid_fr    = round(T_img / 2);
mid_start = mid_fr - 5 * fps_img;
mid_end   = mid_fr + 5 * fps_img - 1;
mid_start = max(mid_start, 1);
mid_end   = min(mid_end,   T_img);

t_mid = t_img(mid_start:mid_end);

figure('Position',[200 200 900 700]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

ax1b = nexttile(1);
plot(ax1b, t_img, -breath, 'k', 'LineWidth', 1); hold(ax1b,'on');
if ~isempty(idx_b)
    plot(ax1b, t_b, -y_b, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 3);
end
xlim(ax1b, [t_mid(1) t_mid(end)]);
xlabel(ax1b, 'Time (s)');
ylabel(ax1b, 'Breathing (z)');
title(ax1b, 'Breathing trace (DLC) with detected inspirations');
grid(ax1b,'off');
hold(ax1b,'off');

ax2b = nexttile(2, [2 1]);
hold(ax2b,'on');

for kk = 1:N_roi
    if ismember(kk, roiSpk_id)
        trace_color = [1 0 0];
    else
        trace_color = [0 0 0];
    end
    plot(ax2b, t_img, Y(:,kk), 'Color', trace_color, 'LineWidth', 1);
end

xline(ax2b, t_b, 'Color', [1 0.1 0.1], 'LineWidth', 0.3);

xlim(ax2b, [t_mid(1) t_mid(end)]);
ylim(ax2b, [-gap, gap*(N_roi-1)+gap]);
yticks(ax2b, yt);
yticklabels(ax2b, ytl);
xlabel(ax2b, 'Time (s)');
ylabel(ax2b, 'ROI#');
title(ax2b, '\DeltaF/F');
grid(ax2b,'off');

plot(ax2b, [t_mid(end) t_mid(end)], [dFF_scale 2*dFF_scale], 'k', 'LineWidth', 2);
text(ax2b, t_mid(end)+0.1, dFF_scale + 0.05*dFF_scale, ...
    sprintf('%.2f \\DeltaF/F', dFF_scale), ...
    'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize', 10);

hold(ax2b,'off');

outFig = fullfile(folderPath, 'breath_N_dFF_mid10s.png');
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved middle 10 s zoom:\n  %s\n', outFig);

%% Figure 2: Plot breath and detected calcium spk
fs   = fps_img;

% build stacked matrix (highest at top)
Y = zeros(T_img, nSpkROI);
for k = 1:nSpkROI
    i = roiSpk_id(k);
    Y(:,k) = dFF(:,i) + 1.5*gap*(nSpkROI - k);
end

% y-tick positions and labels for spiking ROIs (ascending position, bottom=last ROI)
yt_spk  = fliplr(1.5*gap*(nSpkROI - (1:nSpkROI)));
ytl_spk = arrayfun(@(id) sprintf('ROI#%02d', id), roiSpk_id(end:-1:1), 'UniformOutput', false);

figure('Position',[200 200 900 700]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1);
plot(ax1, t_img, -breath, 'Color', 'k', 'LineWidth',1.2); xlim([0 30]);
xline(ax1, t_b, 'Color', [1 0.1 0.1], 'LineWidth', 0.3);

ax2 = nexttile(2, [2 1]);
plot(ax2, t_img, Y, 'Color', 'k', 'LineWidth',0.8); xlim([0 30]);
xline(ax2, t_b, 'Color', [1 0.1 0.1], 'LineWidth', 0.3);

off = 1.5*gap*(nSpkROI - (1:nSpkROI));

tickCenter = off - 0.50*gap;
tickHalf   = 0.12*gap;

hold on
for k = 1:nSpkROI
    spk = find(ca_spk_train(:,k));
    if isempty(spk), continue; end

    x = spk;
    xx = [x(:) x(:)]';
    yy = [(tickCenter(k)-tickHalf)*ones(numel(spk),1), ...
          (tickCenter(k)+tickHalf)*ones(numel(spk),1)]';

    plot(xx/fs, yy, 'k', 'LineWidth', 1);
end
yticks(ax2, yt_spk);
yticklabels(ax2, ytl_spk);
ylabel(ax2, 'ROI#');
xlabel('Time (s)')

outFig = fullfile(folderPath, 'breath_ca_spk_30s.png');
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved breath+ca spk (30s):\n  %s\n', outFig);

figure('Position',[200 200 900 700]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1);
plot(ax1, t_img, -breath, 'Color', 'k', 'LineWidth',1.2); xlim([0 T_img/fs]);
xline(ax1, t_b, 'Color', [1 0.1 0.1], 'LineWidth', 0.3);

ax2 = nexttile(2, [2 1]);
plot(ax2, t_img, Y, 'Color', 'k', 'LineWidth',0.8); xlim([0 T_img/fs]);
xline(ax2, t_b, 'Color', [1 0.1 0.1], 'LineWidth', 0.3);

off = 1.5*gap*(nSpkROI - (1:nSpkROI));

tickCenter = off - 0.50*gap;
tickHalf   = 0.12*gap;

hold on
for k = 1:nSpkROI
    spk = find(ca_spk_train(:,k));
    if isempty(spk), continue; end

    x = spk;
    xx = [x(:) x(:)]';
    yy = [(tickCenter(k)-tickHalf)*ones(numel(spk),1), ...
          (tickCenter(k)+tickHalf)*ones(numel(spk),1)]';

    plot(xx/fs, yy, 'k', 'LineWidth', 1);
end
yticks(ax2, yt_spk);
yticklabels(ax2, ytl_spk);
ylabel(ax2, 'ROI#');
xlabel('Time (s)')

outFig = fullfile(folderPath, 'breath_ca_spk_full.png');
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved breath+ca spk (full):\n  %s\n', outFig);

%% ---- Figure 3. Inspiration-triggered dFF ----

win  = round(0.5 * fs);          % +/-0.5 s
t_evt = (-win:win) / fs;

T_img = size(dFF,1);
insp_ev = insp_onsets(insp_onsets > win & insp_onsets < T_img - win);  % local only
nEv = numel(insp_ev);

fprintf('Event-triggered : using %d inspiration onsets.\n', nEv);

dff_evt_seg_all = cell(N_roi,1);

for i = 1:N_roi
    roi_trace = dFF(:, i);               % RAW dF/F

    seg = zeros(nEv, 2*win+1);
    for e = 1:nEv
        c = insp_ev(e);
        seg(e,:) = roi_trace(c-win : c+win);
    end

    dff_evt_seg_all{i} = seg;
end

% ---------- Per-ROI event-triggered overlays (all trials + mean, RAW dFF) ----------

figure('Name', (sprintf('Event-triggered dFF (n=%d breaths)', nEv)));

for i = 1:N_roi
    seg = dff_evt_seg_all{i};    % nEv x (2*win+1)
    if isempty(seg)
        continue;
    end

    % average for this ROI (RAW dF/F)
    dff_avg = mean(seg, 1);

    subplot(ceil(N_roi/10), 10, i);
    hold on;

    % all individual trials
    plot(t_evt, seg', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.7);
    % overlay mean
    plot(t_evt, dff_avg, 'Color', 'k', 'LineWidth', 2);

    % reference lines
    xline(0, '--', 'Color', [0.4 0.4 0.4]);  % inspiration onset at 0
    yline(0, '--',  'Color', [0.6 0.6 0.6]);

    xlabel('Time from inspiration onset (s)');
    ylabel('\DeltaF/F');
    title(sprintf('ROI#%d', i));

    hold off;
    axis square;

end

outFig = fullfile(folderPath, 'triggered_avg.png');
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved triggered average:\n  %s\n', outFig);

%% ---- Figure 4. Per-ROI calcium-triggered breathing overlays ----

win_ca = round(0.5 * fs);      % +/-0.5 s window for calcium-triggered breathing
t_ca_evt = (-win_ca:win_ca) / fs;

for k = 1:nSpkROI
    roi_id = roiSpk_id(k);
    ca_frames = ca_spk_id{k};       % calcium spike frame indices for this ROI
    nCaEv = numel(ca_frames);

    if nCaEv == 0, continue; end

    % --- LEFT: breath-triggered dFF overlay ---
    seg_breath = dff_evt_seg_all{roi_id};   % nEv x (2*win+1)
    dff_avg = mean(seg_breath, 1);

    % --- RIGHT: calcium-triggered breathing segments ---
    ca_valid = ca_frames(ca_frames > win_ca & ca_frames < T_img - win_ca);
    nCaValid = numel(ca_valid);
    seg_ca = zeros(nCaValid, 2*win_ca+1);
    for e = 1:nCaValid
        c = ca_valid(e);
        seg_ca(e,:) = breath(c-win_ca : c+win_ca);
    end
    breath_avg = mean(seg_ca, 1);

    % --- Figure ---
    fig4 = figure('Color','w', 'Name', sprintf('ROI#%02d triggered avg', roi_id), ...
                  'Position', [100 100 900 400]);

    % LEFT subplot: breath-triggered dFF
    subplot(1,2,1);
    hold on;
    plot(t_evt, seg_breath', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.7);
    plot(t_evt, dff_avg, 'Color', 'k', 'LineWidth', 2);
    xline(0, '--', 'Color', [0.4 0.4 0.4]);
    yline(0, '--', 'Color', [0.6 0.6 0.6]);
    xlabel('Time from inspiration onset (s)');
    ylabel('\DeltaF/F');
    title(sprintf('ROI#%02d  breath-trig dFF (n=%d)', roi_id, nEv));
    axis square;
    hold off;

    % RIGHT subplot: calcium-triggered breathing
    subplot(1,2,2);
    hold on;
    if nCaValid > 0
        plot(t_ca_evt, seg_ca', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.7);
        plot(t_ca_evt, breath_avg, 'Color', 'k', 'LineWidth', 2);
    end
    xline(0, '--', 'Color', [0.4 0.4 0.4]);
    yline(0, '--', 'Color', [0.6 0.6 0.6]);
    xlabel('Time from Ca^{2+} spike (s)');
    ylabel('Breathing (z-score)');
    title(sprintf('ROI#%02d  Ca-trig breath (n=%d)', roi_id, nCaValid));
    axis square;
    hold off;

    % Save
    outPng = fullfile(folderPath, sprintf('ROI%02d_triggered_avg.png', roi_id));
    exportgraphics(fig4, outPng, 'ContentType','vector');
    outFigFile = fullfile(folderPath, sprintf('ROI%02d_triggered_avg.fig', roi_id));
    savefig(fig4, outFigFile);
    fprintf('Saved triggered avg for ROI#%02d:\n  %s\n  %s\n', roi_id, outPng, outFigFile);
end

%% ---- Figure 5: dt_last analysis (raster + histogram + permutation) ----
% For each spiking ROI, compute time from LAST inspiration (at or before)
% each calcium spike.  dt_last >= 0 always.

win5   = round(2.5 * fs);
w_ms   = 150;
w_sec  = w_ms / 1000;
M_shuf = 3000;

b_frames = sort(insp_onsets(:));
Delta_k  = diff(b_frames) / fs;       % breath interval durations (s)

% Collect per-ROI results
dtlast_results = struct();

for k = 1:nSpkROI
    roi_id   = roiSpk_id(k);
    ca_ev    = ca_spk_id{k}(:);        % calcium spike frames for this ROI
    nCaTotal = numel(ca_ev);

    if nCaTotal == 0, continue; end

    % ---- Assign each Ca event to its breath interval ----
    % bin_idx(j) = i means ca_ev(j) falls in [b_frames(i), b_frames(i+1))
    [~, ~, bin_idx] = histcounts(ca_ev, b_frames);
    keep_ca = (bin_idx >= 1) & (bin_idx <= numel(b_frames) - 1);
    ca_in   = ca_ev(keep_ca);
    bin_in  = bin_idx(keep_ca);
    nEv_in  = numel(ca_in);

    if nEv_in == 0
        fprintf('ROI#%02d : no calcium events within breath intervals, skipping.\n', roi_id);
        continue;
    end

    % ---- dt_last: time from last inspiration (at or before) ----
    dt_last = (ca_in - b_frames(bin_in)) / fs;   % always >= 0

    % ---- Build insp_rel for raster ----
    insp_rel_in = cell(nEv_in, 1);
    for e = 1:nEv_in
        c      = ca_in(e);
        nearby = insp_onsets(insp_onsets >= (c-win5) & insp_onsets <= (c+win5));
        insp_rel_in{e} = (nearby - c) / fs;
    end

    % ---- Histogram bins (one-sided, 0 to max interval) ----
    half_max   = max(Delta_k);
    bin_number = round(half_max / (1/fs));
    if bin_number < 10, bin_number = 10; end
    if mod(bin_number,2), bin_number = bin_number + 1; end

    bin_edges_dt = linspace(0, half_max, bin_number);
    bin_ctrs     = (bin_edges_dt(1:end-1) + bin_edges_dt(2:end)) / 2;
    nBins        = numel(bin_ctrs);
    H_obs        = histcounts(dt_last, bin_edges_dt, 'Normalization','probability') * 100;

    % ---- Observed scan statistic ----
    S_obs = breath_scan_stat(dt_last, w_sec);

    % ---- Shuffle (interval-based) ----
    S_shuf = nan(M_shuf, 1);
    H_shuf = nan(M_shuf, nBins);
    for m = 1:M_shuf
        % Place each calcium event uniformly within its breath interval
        dt_last_sh    = rand(nEv_in, 1) .* Delta_k(bin_in);
        S_shuf(m)     = breath_scan_stat(dt_last_sh, w_sec);
        H_shuf(m,:)   = histcounts(dt_last_sh, bin_edges_dt, 'Normalization','probability') * 100;
    end
    p_val   = (1 + sum(S_shuf >= S_obs)) / (M_shuf + 1);
    n_close = sum(dt_last <= w_sec);
    fprintf('ROI#%02d : S_obs=%d  median(S_null)=%.1f  p=%.4f  (w=%.2fs, M=%d, n=%d)\n', ...
        roi_id, S_obs, median(S_shuf), p_val, w_sec, M_shuf, nEv_in);

    % ---- Color by significance ----
    sig_color = [1 0.1 0.1];
    ns_color  = [0.1 0.1 0.1];
    obs_color = double(p_val < 0.05) * sig_color + double(p_val >= 0.05) * ns_color;

    % ---- Figure: 1x3 layout ----
    fig5 = figure('Color', 'w', ...
        'Name', sprintf('dt_last: ROI#%02d', roi_id), ...
        'Position', [0 0 1100 400]);

    % ---- Subplot 1: Raster sorted by dt_last ascending ----
    ax_rast = subplot(1,3,1);
    hold(ax_rast, 'on');

    [~, sort_idx] = sort(dt_last, 'ascend');

    for row = 1:nEv_in
        e       = sort_idx(row);
        t_ticks = insp_rel_in{e};
        for tt = 1:numel(t_ticks)
            plot(ax_rast, [t_ticks(tt) t_ticks(tt)], [row-0.4  row+0.4], ...
                'k', 'LineWidth', 1);
        end
    end

    xline(ax_rast, 0, 'r--', 'LineWidth', 1);
    xlim(ax_rast, [-2.5 2.5]);
    ylim(ax_rast, [0 nEv_in+1]);
    set(ax_rast, 'YDir', 'reverse');
    xlabel(ax_rast, 'Time from Ca^{2+} event (s)');
    ylabel(ax_rast, sprintf('Ca^{2+} event # \n (sorted by dt_{last})'), 'HorizontalAlignment', 'center');
    title(ax_rast, sprintf('ROI#%02d  (%d events in %d breaths)', roi_id, nEv_in, sum(Sb)));
    hold(ax_rast, 'off');

    % ---- Subplot 2: Histogram vs shuffle envelope ----
    ax_hist = subplot(1,3,2);
    hold(ax_hist, 'on');

    % Shuffle envelope (2.5-97.5 percentile band)
    H_mean = mean(H_shuf, 1);
    H_lo   = prctile(H_shuf,  2.5, 1);
    H_hi   = prctile(H_shuf, 97.5, 1);

    be = bin_edges_dt;
    [xs_hi, ys_hi] = stairs(be, [H_hi  H_hi(end)]);
    [xs_lo, ys_lo] = stairs(be, [H_lo  H_lo(end)]);

    fill(ax_hist, [xs_hi; flipud(xs_lo)], [ys_hi; flipud(ys_lo)], ...
         [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    stairs(ax_hist, be, [H_mean H_mean(end)], 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
    stairs(ax_hist, be, [H_obs  H_obs(end)],  'Color', obs_color, 'LineWidth', 2.5);

    xlabel(ax_hist, 'Time from last inspiration (s)');
    ylabel(ax_hist, 'Proportion of Ca^{2+} events (%)');
    legend(ax_hist, {'95% shuffle envelope', 'Shuffle mean', 'Observed'}, 'Location', 'northeast');
    xlim(ax_hist, [0 half_max]);
    ylim(ax_hist, [0 40]);
    axis(ax_hist, 'square');
    hold(ax_hist, 'off');

    % ---- Subplot 3: Scan statistic null distribution ----
    ax_scan = subplot(1,3,3);
    hold(ax_scan, 'on');

    S_shuf_pct = S_shuf / nEv_in * 100;
    S_obs_pct  = S_obs  / nEv_in * 100;

    histogram(ax_scan, S_shuf_pct, 'Normalization', 'percentage', ...
              'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    xline(ax_scan, S_obs_pct, 'Color', obs_color, 'LineWidth', 2);

    xlabel(ax_scan, sprintf('Max %%events in a %d-ms sliding window', w_ms));
    ylabel(ax_scan, '% Shuffles');

    xlim(ax_scan, [0 100]);
    ylim(ax_scan, [0 30]);
    yl_sc = ylim(ax_scan);  xl_sc = xlim(ax_scan);
    text(ax_scan, mean(xl_sc), yl_sc(2) + 2, ...
         sprintf('%d/%d  (p=%.3f)', S_obs, nEv_in, p_val), ...
         'Color', obs_color, 'FontSize', 11, ...
         'HorizontalAlignment','center', 'VerticalAlignment','middle');
    axis(ax_scan, 'square');
    box(ax_scan, 'on');
    hold(ax_scan, 'off');

    % ---- Save figure ----
    outFig_fig = fullfile(folderPath, sprintf('calcium_breath_dtlast_ROI%02d.fig', roi_id));
    savefig(fig5, outFig_fig);
    outFig_png = fullfile(folderPath, sprintf('calcium_breath_dtlast_ROI%02d.png', roi_id));
    exportgraphics(fig5, outFig_png, 'ContentType', 'vector');
    fprintf('Saved dt_last figure for ROI#%02d:\n  %s\n  %s\n', roi_id, outFig_fig, outFig_png);

    % ---- Store per-ROI results ----
    dtlast_results(k).roi_id    = roi_id;
    dtlast_results(k).dt_last   = dt_last;
    dtlast_results(k).ca_in     = ca_in;
    dtlast_results(k).bin_in    = bin_in;
    dtlast_results(k).nEv_in    = nEv_in;
    dtlast_results(k).S_obs     = S_obs;
    dtlast_results(k).S_shuf    = S_shuf;
    dtlast_results(k).p_val     = p_val;
    dtlast_results(k).n_close   = n_close;
    dtlast_results(k).H_obs     = H_obs;
    dtlast_results(k).H_shuf    = H_shuf;
    dtlast_results(k).bin_edges = bin_edges_dt;
    dtlast_results(k).half_max  = half_max;
end

%% filter breathing signal for hilbert
d = designfilt("lowpassfir", ...
    PassbandFrequency=0.2,StopbandFrequency=0.25, ...
    PassbandRipple=1,StopbandAttenuation=60, ...
    DesignMethod="equiripple");
breath_filt = filtfilt(d,breath);

%% hilbert transform of breathing to extract phase
y = hilbert(breath_filt);
phi = angle(y);
%%
% figure;
% subplot(3,1,1);
% plot(t_breath, breath);
%
% subplot(3,1,2);
% plot(t_breath, breath_filt);
%
% subplot(3,1,3);
% plot(t_breath, phi);

%% Figure 6: Phase histogram per spiking ROI
for k = 1:nSpkROI
    roi_id = roiSpk_id(k);
    spk_frames = ca_spk_id{k};
    spk_frames = spk_frames(spk_frames >= 1 & spk_frames <= numel(phi));
    if isempty(spk_frames), continue; end

    ca_phi = phi(spk_frames);

    figure('Color','white');
    histogram(ca_phi, 40, 'Normalization','percentage');
    xlim([-pi pi]);
    ylim([0 25]);
    axis square;
    xlabel('Phase (rad)');
    ylabel('% events');
    title(sprintf('ROI#%02d  phase preference (n=%d)', roi_id, numel(ca_phi)));

    outFig = fullfile(folderPath, sprintf('phase_hist_ROI%02d.png', roi_id));
    exportgraphics(gcf, outFig, 'ContentType','vector');
    fprintf('Saved phase histogram for ROI#%02d:\n  %s\n', roi_id, outFig);
end

%% ---- Save master file ----
outMat_master = fullfile(folderPath, 'breath_master_working.mat');
save(outMat_master, ...
    'folderPath', ...
    'fps_img', 'fps_breath', 'fs', ...
    'breath_peak', ...
    'roiSpk_id', 'nSpkROI', 'ca_spk_data', 'ca_spk_id', 't_spk', ...
    'w_sec', 'M_shuf', 'w_ms', ...
    'win', 't_evt', 'win5', ...
    'breath', 't_img', 'dFF', 'Sb', ...
    'Delta_k', 'b_frames', ...
    'dtlast_results', ...
    'dff_evt_seg_all');
fprintf('Saved master file:\n  %s\n', outMat_master);

%% ---- Local functions ----

function S = breath_scan_stat(t, w)
% BREATH_SCAN_STAT  Max number of events in any sliding window of width w.
%   t : vector of t_pre values (s)
%   w : window width (s)
%   S : scan statistic (integer count)
%
%   Uses a sorted two-pointer sweep: O(n log n).
%   Window definition: [tau, tau+w), equivalent to the centred form
%   |t_j - tau| < w/2 when tau is swept over all data-anchored positions.

t = sort(t(:));
n = numel(t);
if n == 0,  S = 0;  return;  end

S    = 1;
left = 1;
for right = 2:n
    while t(right) - t(left) >= w
        left = left + 1;
    end
    S = max(S, right - left + 1);
end
end
