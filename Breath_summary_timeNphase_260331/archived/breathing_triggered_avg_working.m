%%
close all; clear; clc;

folderPath = "D:\batch_dffQC_test_260325\260323_vgat_g8s\phys\processed\-550-933-z-50_4x_3000f_00001";

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
%F(1:30, :) = [];

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
bp_file = dir(fullfile(folderPath, '*breath_peak_data.mat'));
assert(~isempty(bp_file), 'No *breath_peak_data.mat found in %s', folderPath);
breath_peak = load(fullfile(bp_file(1).folder, bp_file(1).name));

insp_onsets = breath_peak.insp_onset_idx;
Sb = breath_peak.insp_onsets_train;

insp_onsets(insp_onsets<30)=[];
insp_onsets = insp_onsets - 30;   % adjust for tossed frames
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
figure('Color','w');
plot(fk_raw, Sk_raw,                'k',  'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([0 10]);
grid('on');
axis square;

% Save breathing spectrum
outFig = fullfile(folderPath, 'breathing_spectrum_raw.pdf');
exportgraphics(gcf, outFig, 'ContentType','vector', 'BackgroundColor','none');
fprintf('Saved breathing spectrum:\n  %s\n', outFig);
%% ------------------- FIGURE 1: BREATH + STACKED dFF (Top 10 ROIs) --------
T_img = numel(breath);
xl2 = 60;
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

figure('Position',[200 200 900 700], 'Color','w');
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
%xlim(ax1, [t_img(1) t_img(end)])
xlim(ax1, [t_img(1) xl2])
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

%xlim(ax2, [t_img(1) t_img(end)]);
xlim(ax2, [t_img(1) xl2]);
ylim(ax2, [-gap, gap*(N_roi-1)+gap]);

yticks(ax2, yt);
yticklabels(ax2, ytl);

xlabel(ax2, 'Time (s)');
ylabel(ax2, 'ROI#');
title(ax2, 'dF/F');
grid(ax2,'off');

% SCALE BAR
dFF_scale = 0.4;     % 0.4 dF/F vertical bar

% position near bottom-right of visible area
xl = xlim(ax2);
x_sb = xl(2) - 2;                         % 2 s from right edge
y_sb = -0.3 * gap;                         % below bottom trace

plot(ax2, [x_sb x_sb], [y_sb y_sb+dFF_scale], 'k', 'LineWidth', 2);
text(ax2, x_sb + 0.5, y_sb + dFF_scale/2, ...
    sprintf('%.1f dF/F', dFF_scale), ...
    'VerticalAlignment','middle','HorizontalAlignment','left', ...
    'FontSize', 9);
hold(ax2,'off');

% Save breathing + dFF
outFig = fullfile(folderPath, 'breath_N_dFF.pdf');
exportgraphics(gcf, outFig, 'ContentType','vector', 'BackgroundColor','none');
fprintf('Saved breathing and dFF traces:\n  %s\n', outFig);

% % ---- Figure 1b: middle 10 s zoom ----
% mid_fr    = round(T_img / 2);
% mid_start = mid_fr - 5 * fps_img;
% mid_end   = mid_fr + 5 * fps_img - 1;
% mid_start = max(mid_start, 1);
% mid_end   = min(mid_end,   T_img);
% 
% t_mid = t_img(mid_start:mid_end);
% 
% figure('Position',[200 200 900 700], 'Color','none');
% tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
% 
% ax1b = nexttile(1);
% plot(ax1b, t_img, -breath, 'k', 'LineWidth', 1); hold(ax1b,'on');
% if ~isempty(idx_b)
%     plot(ax1b, t_b, -y_b, 'ro', 'MarkerFaceColor','r', 'MarkerSize', 3);
% end
% xlim(ax1b, [t_mid(1) t_mid(end)]);
% xlabel(ax1b, 'Time (s)');
% ylabel(ax1b, 'Breathing (z)');
% title(ax1b, 'Breathing trace (DLC) with detected inspirations');
% grid(ax1b,'off');
% hold(ax1b,'off');
% 
% ax2b = nexttile(2, [2 1]);
% hold(ax2b,'on');
% 
% for kk = 1:N_roi
%     if ismember(kk, roiSpk_id)
%         trace_color = [1 0 0];
%     else
%         trace_color = [0 0 0];
%     end
%     plot(ax2b, t_img, Y(:,kk), 'Color', trace_color, 'LineWidth', 1);
% end
% 
% xline(ax2b, t_b, 'Color', [1 0.1 0.1], 'LineWidth', 0.3);
% 
% xlim(ax2b, [t_mid(1) t_mid(end)]);
% ylim(ax2b, [-gap, gap*(N_roi-1)+gap]);
% yticks(ax2b, yt);
% yticklabels(ax2b, ytl);
% xlabel(ax2b, 'Time (s)');
% ylabel(ax2b, 'ROI#');
% title(ax2b, 'dF/F');
% grid(ax2b,'off');
% 
% plot(ax2b, [t_mid(end) t_mid(end)], [dFF_scale 2*dFF_scale], 'k', 'LineWidth', 2);
% text(ax2b, t_mid(end)+0.1, dFF_scale + 0.05*dFF_scale, ...
%     sprintf('%.2f dF/F', dFF_scale), ...
%     'VerticalAlignment','bottom','HorizontalAlignment','left','FontSize', 10);
% 
% hold(ax2b,'off');
% 
% set(findall(gcf,'Type','axes'),'Color','none');
% outFig = fullfile(folderPath, 'breath_N_dFF_mid10s.pdf');
% exportgraphics(gcf, outFig, 'ContentType','vector', 'BackgroundColor','none');
% fprintf('Saved middle 10 s zoom:\n  %s\n', outFig);

%% Figure 2: Plot breath and detected calcium spk
fs   = fps_img;

% build stacked matrix (highest at top)
Y = zeros(T_img, nSpkROI);
for k = 1:nSpkROI
    i = roiSpk_id(k);
    Y(:,k) = dFF(:,i) + 1.5*gap*(nSpkROI - k);
end
xl2 = 60;

% y-tick positions and labels for spiking ROIs (ascending position, bottom=last ROI)
yt_spk  = fliplr(1.5*gap*(nSpkROI - (1:nSpkROI)));
ytl_spk = arrayfun(@(id) sprintf('ROI#%02d', id), roiSpk_id(end:-1:1), 'UniformOutput', false);

figure('Position',[200 200 900 700], 'Color','w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1);
plot(ax1, t_img, -breath, 'Color', 'k', 'LineWidth',1.2); xlim([0 xl2]);
xline(ax1, t_b, 'Color', [1 0.1 0.1], 'LineWidth', 0.3);

ax2 = nexttile(2, [2 1]);
plot(ax2, t_img, Y, 'Color', 'k', 'LineWidth',0.8); xlim([0 xl2]);
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

% SCALE BAR
dFF_sb = 0.4;
x_sb2 = xl2 - 1;
y_sb2 = -0.5 * 1.5*gap;
plot(ax2, [x_sb2 x_sb2], [y_sb2 y_sb2+dFF_sb], 'k', 'LineWidth', 2);
text(ax2, x_sb2 + 0.3, y_sb2 + dFF_sb/2, ...
    sprintf('%.1f dF/F', dFF_sb), ...
    'VerticalAlignment','middle','HorizontalAlignment','left', 'FontSize', 9);

yticks(ax2, yt_spk);
yticklabels(ax2, ytl_spk);
ylabel(ax2, 'ROI#');
xlabel('Time (s)')

outFig = fullfile(folderPath, sprintf('breath_ca_spk_%ds.pdf', xl2));
exportgraphics(gcf, outFig, 'ContentType','vector', 'BackgroundColor','none');
fprintf('Saved breath+ca spk (%ds):\n  %s\n', xl2, outFig);

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

figure('Name', (sprintf('Event-triggered dFF (n=%d breaths)', nEv)), 'Color','w');

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
    ylabel('dF/F');
    title(sprintf('ROI#%d', i));

    hold off;
    axis square;

end

outFig = fullfile(folderPath, 'triggered_avg.pdf');
exportgraphics(gcf, outFig, 'ContentType','vector', 'BackgroundColor','none');
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
    ylabel('dF/F');
    title(sprintf('ROI#%02d  breath-trig dFF (n=%d)', roi_id, nEv));
    axis square;
    hold off;

    % RIGHT subplot: calcium-triggered breathing
    subplot(1,2,2);
    hold on;
    if nCaValid > 0
        plot(t_ca_evt, -seg_ca', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.7);
        plot(t_ca_evt, -breath_avg, 'Color', 'k', 'LineWidth', 2);
    end
    xline(0, '--', 'Color', [0.4 0.4 0.4]);
    yline(0, '--', 'Color', [0.6 0.6 0.6]);
    xlabel('Time from Ca^{2+} spike (s)');
    ylabel('Breathing (z-score)');
    title(sprintf('ROI#%02d  Ca-trig breath (n=%d)', roi_id, nCaValid));
    axis square;
    hold off;

    % Save
    outPdf = fullfile(folderPath, sprintf('ROI%02d_triggered_avg.pdf', roi_id));
    exportgraphics(fig4, outPdf, 'ContentType','vector', 'BackgroundColor','none');
    outFigFile = fullfile(folderPath, sprintf('ROI%02d_triggered_avg.fig', roi_id));
    savefig(fig4, outFigFile);
    fprintf('Saved triggered avg for ROI#%02d:\n  %s\n  %s\n', roi_id, outPdf, outFigFile);
end

% %% ---- Figure 4b: Permutation test sweep over w_ms values ----
% w_ms_sweep = [34, 67, 100, 151, 301];           % window widths to sweep (ms)
% M_shuf_sweep = 2000;
% 
% permtest_dir = fullfile(folderPath, 'permutation_test');
% if ~exist(permtest_dir, 'dir'), mkdir(permtest_dir); end
% 
% b_frames_sw = sort(insp_onsets(:));
% Delta_k_sw  = diff(b_frames_sw) / fs;
% 
% permtest_results = struct();
% 
% for k = 1:nSpkROI
%     roi_id   = roiSpk_id(k);
%     ca_ev    = ca_spk_id{k}(:);
%     if isempty(ca_ev), continue; end
% 
%     [~, ~, bin_idx] = histcounts(ca_ev, b_frames_sw);
%     keep_ca = (bin_idx >= 1) & (bin_idx <= numel(b_frames_sw) - 1);
%     ca_in   = ca_ev(keep_ca);
%     bin_in  = bin_idx(keep_ca);
%     nEv_in  = numel(ca_in);
% 
%     if nEv_in == 0
%         fprintf('ROI#%02d : no Ca events in breath intervals, skipping sweep.\n', roi_id);
%         continue;
%     end
% 
%     dt_last_sw = (ca_in - b_frames_sw(bin_in)) / fs;
% 
%     for wi = 1:numel(w_ms_sweep)
%         wm  = w_ms_sweep(wi);
%         ws  = wm / 1000;
% 
%         % Observed scan statistic
%         S_obs_sw = breath_scan_stat(dt_last_sw, ws);
% 
%         % Shuffle
%         S_shuf_sw = nan(M_shuf_sweep, 1);
%         for m = 1:M_shuf_sweep
%             dt_sh = rand(nEv_in, 1) .* Delta_k_sw(bin_in);
%             S_shuf_sw(m) = breath_scan_stat(dt_sh, ws);
%         end
%         p_val_sw = (1 + sum(S_shuf_sw >= S_obs_sw)) / (M_shuf_sweep + 1);
% 
%         fprintf('ROI#%02d  w=%dms : S_obs=%d  p=%.4f  (n=%d)\n', ...
%             roi_id, wm, S_obs_sw, p_val_sw, nEv_in);
% 
%         % Color by significance
%         if p_val_sw < 0.01
%             obs_color_sw = [1 0.1 0.1];
%         else
%             obs_color_sw = [0.1 0.1 0.1];
%         end
% 
%         % Figure: scan statistic null distribution
%         fig_sw = figure('Color', 'w', ...
%             'Name', sprintf('Permutation ROI#%02d w=%dms', roi_id, wm), ...
%             'Position', [100 100 500 450]);
% 
%         S_shuf_pct_sw = S_shuf_sw / nEv_in * 100;
%         S_obs_pct_sw  = S_obs_sw  / nEv_in * 100;
% 
%         histogram(S_shuf_pct_sw, 'Normalization', 'percentage', ...
%             'FaceColor', [0.65 0.65 0.65], 'EdgeColor', 'none', 'FaceAlpha', 0.85);
%         hold on;
%         xline(S_obs_pct_sw, 'Color', obs_color_sw, 'LineWidth', 2);
% 
%         xlabel(sprintf('Max %%events in a %d-ms sliding window', wm));
%         ylabel('% Shuffles');
%         title(sprintf('ROI#%02d  w=%dms  (p=%.3f)', roi_id, wm, p_val_sw));
% 
%         xlim([0 100]);
%         ylim([0 30]);
%         yl_sw = ylim;  xl_sw = xlim;
%         text(mean(xl_sw), yl_sw(2) + 2, ...
%             sprintf('%d/%d  (p=%.3f)', S_obs_sw, nEv_in, p_val_sw), ...
%             'Color', obs_color_sw, 'FontSize', 11, ...
%             'HorizontalAlignment','center', 'VerticalAlignment','middle');
%         axis square;
%         box on;
%         hold off;
% 
%         % Save
%         outPdf_sw = fullfile(permtest_dir, sprintf('permtest_ROI%02d_w%dms.pdf', roi_id, wm));
%         exportgraphics(fig_sw, outPdf_sw, 'ContentType','vector', 'BackgroundColor','none');
%         outFig_sw = fullfile(permtest_dir, sprintf('permtest_ROI%02d_w%dms.fig', roi_id, wm));
%         savefig(fig_sw, outFig_sw);
%         fprintf('  Saved: %s\n', outPdf_sw);
% 
%         % Store results
%         permtest_results(k).roi_id = roi_id;
%         permtest_results(k).w_ms(wi) = wm;
%         permtest_results(k).S_obs(wi) = S_obs_sw;
%         permtest_results(k).p_val(wi) = p_val_sw;
%         permtest_results(k).S_shuf{wi} = S_shuf_sw;
%         permtest_results(k).nEv_in = nEv_in;
%     end
% end

%% ---- Figure 5: dt_last analysis (raster + histogram) ----
% For each spiking ROI, compute time from LAST inspiration (at or before)
% each calcium spike.  dt_last >= 0 always.
b_frames = sort(insp_onsets(:));
Delta_k  = diff(b_frames) / fs;       % breath interval durations (s)

win_raster = mean(Delta_k)*2; % s
win5   = round(win_raster * fs);

% Collect per-ROI results
dtlast_results = struct();

for k = 1:nSpkROI
    roi_id   = roiSpk_id(k);
    ca_ev    = ca_spk_id{k}(:);        % calcium spike frames for this ROI
    nCaTotal = numel(ca_ev);

    if nCaTotal == 0, continue; end

    % ---- Assign each Ca event to its breath interval ----
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
    H_obs        = histcounts(dt_last, bin_edges_dt, 'Normalization','probability') * 100;

    fprintf('ROI#%02d : %d Ca events in %d breath intervals\n', roi_id, nEv_in, numel(Delta_k));

    % ---- Figure: raster (2 cols) + histogram (1 col) ----
    fig5 = figure('Color', 'w', ...
        'Name', sprintf('dt_last: ROI#%02d', roi_id), ...
        'Position', [0 0 1100 400]);

    % ---- Subplot 1: Raster sorted by dt_last ascending (spans 2 of 3 cols) ----
    ax_rast = subplot(1,3,[1 2]);
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

    xline(ax_rast, 0, 'k--', 'LineWidth', 1);
    xlim(ax_rast, [-win_raster, win_raster]);
    ylim(ax_rast, [0 nEv_in+1]);
    set(ax_rast, 'YDir', 'reverse');
    xlabel(ax_rast, 'Time from Ca^{2+} event (s)');
    ylabel(ax_rast, sprintf('Ca^{2+} event # \n (sorted by dt_{last})'), 'HorizontalAlignment', 'center');
    title(ax_rast, sprintf('ROI#%02d  (%d events in %d breaths)', roi_id, nEv_in, sum(Sb)));
    hold(ax_rast, 'off');

    % ---- Subplot 2: dt_last histogram ----
    ax_hist = subplot(1,3,3);
    hold(ax_hist, 'on');

    stairs(ax_hist, bin_edges_dt, [H_obs H_obs(end)], 'Color', 'k', 'LineWidth', 2);
    histogram(ax_hist, 'BinEdges', bin_edges_dt, 'BinCounts', H_obs, ...
        'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    
    xlabel(ax_hist, 'Time from last inspiration (s)');
    ylabel(ax_hist, 'Proportion of Ca^{2+} events (%)');
    xlim(ax_hist, [0 half_max]);
    ylim(ax_hist, [0 45]);
    axis(ax_hist, 'square');
    hold(ax_hist, 'off');

    % ---- Save figure ----
    outFig_fig = fullfile(folderPath, sprintf('calcium_breath_dtlast_ROI%02d.fig', roi_id));
    savefig(fig5, outFig_fig);
    outFig_pdf = fullfile(folderPath, sprintf('calcium_breath_dtlast_ROI%02d.pdf', roi_id));
    exportgraphics(fig5, outFig_pdf, 'ContentType', 'vector', 'BackgroundColor','none');
    fprintf('Saved dt_last figure for ROI#%02d:\n  %s\n  %s\n', roi_id, outFig_fig, outFig_pdf);

    % ---- Store per-ROI results ----
    dtlast_results(k).roi_id    = roi_id;
    dtlast_results(k).dt_last   = dt_last;
    dtlast_results(k).ca_in     = ca_in;
    dtlast_results(k).bin_in    = bin_in;
    dtlast_results(k).nEv_in    = nEv_in;
    dtlast_results(k).H_obs     = H_obs;
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

%% Figure 6: Phase histogram per spiking ROI (KDE + cosine fit, side by side)
nBins_phase = 40;
bin_edges_phase = linspace(-pi, pi, nBins_phase+1);
bin_ctrs_phase  = (bin_edges_phase(1:end-1) + bin_edges_phase(2:end)) / 2;

% KDE evaluation grid
nKDE = 200;
xi_kde = linspace(-pi, pi, nKDE);

for k = 1:nSpkROI
    roi_id = roiSpk_id(k);
    spk_frames = ca_spk_id{k};
    spk_frames = spk_frames(spk_frames >= 1 & spk_frames <= numel(phi));
    if isempty(spk_frames), continue; end

    ca_phi = phi(spk_frames);
    nPh    = numel(ca_phi);

    % ---- Circular statistics ----
    mu_vec  = mean(exp(1i * ca_phi));
    theta0  = angle(mu_vec);            % preferred phase (rad)
    R       = abs(mu_vec);              % mean resultant length (0-1)

    % ---- Histogram (percentage) ----
    H_pct = histcounts(ca_phi, bin_edges_phase, 'Normalization','percentage');

    % ---- KDE (circular): wrap data, but DO NOT let wrap affect bandwidth ----
    binW   = 2*pi/nBins_phase;
    xi_kde = linspace(-pi, pi, nKDE).';   % force column
    
    % Get MATLAB's automatic bandwidth from ORIGINAL (unwrapped) data
    [~,~,bw0] = ksdensity(ca_phi(:));     % <- key line (bw from ca_phi only)
    bw0 = min(max(bw0, 0.05), 1.0);       % clamp (tune limits if you want)
    
    % Now do wrapped KDE with fixed bandwidth
    ca_phi_wrap = [ca_phi(:); ca_phi(:)+2*pi; ca_phi(:)-2*pi];
    f_kde = ksdensity(ca_phi_wrap, xi_kde, 'Bandwidth', bw0);
    
    % Convert density -> "% per bin" (match histcounts(...,'percentage'))
    f_kde_pct = f_kde / trapz(xi_kde, f_kde) * 100 * binW;

    % ---- Cosine fit on histogram ----
    % Model: y(theta) = A * cos(theta - theta0) + B
    %   B = baseline (% per bin, = 100/nBins if uniform)
    %   A = modulation amplitude (% per bin)
    %   Peak = A + B (% per bin at preferred phase)
    X_cos = [cos(bin_ctrs_phase(:) - theta0), ones(nBins_phase, 1)];
    coeff = X_cos \ H_pct(:);
    A_fit = coeff(1);                   % modulation amplitude (%)
    B_fit = coeff(2);                   % baseline (%)

    th_fine = linspace(-pi, pi, 200);
    y_fine  = A_fit * cos(th_fine - theta0) + B_fit;

    % ---- Figure: 2 panels side by side ----
    fig6 = figure('Color','w', 'Position', [100 100 900 420]);

    % ---- LEFT: KDE ----
    ax_kde = subplot(1,2,1);
    hold(ax_kde, 'on');

    % Light histogram behind KDE
    histogram(ax_kde, 'BinEdges', bin_edges_phase, 'BinCounts', H_pct, ...
        'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.6);

    % KDE curve
    plot(ax_kde, xi_kde, f_kde_pct, 'k', 'LineWidth', 2);

    % find local maxima of KDE
    [kde_pks, kde_locs] = findpeaks(f_kde_pct, xi_kde, 'MinPeakProminence', 0.3);
    xline(ax_kde, 0, '--', 'Color', [0.6 0.6 0.6]);
    xlim(ax_kde, [-pi pi]);
    ylim([0 25]);
    xticks(ax_kde, [-pi, -pi/2, 0, pi/2, pi]);
    xticklabels(ax_kde, {'-pi', '-pi/2', '0', 'pi/2', 'pi'});
    xlabel(ax_kde, 'Phase (rad)');
    ylabel(ax_kde, '% events per bin');
    nKDE_peaks = numel(kde_pks);
    if nKDE_peaks > 1
        title(ax_kde, sprintf('ROI#%02d  KDE (%d peaks)', roi_id, nKDE_peaks));
    else
        title(ax_kde, sprintf('ROI#%02d  KDE', roi_id));
    end
    axis(ax_kde, 'square');
    hold(ax_kde, 'off');

    % ---- RIGHT: Cosine fit ----
    ax_cos = subplot(1,2,2);
    hold(ax_cos, 'on');

    histogram(ax_cos, 'BinEdges', bin_edges_phase, 'BinCounts', H_pct, ...
        'FaceColor', [0.7 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

    plot(ax_cos, th_fine, y_fine, 'r', 'LineWidth', 2);

    xline(ax_cos, 0, '--', 'Color', [0.6 0.6 0.6]);
    xlim(ax_cos, [-pi pi]);
    ylim([0 25]);
    xticks(ax_cos, [-pi, -pi/2, 0, pi/2, pi]);
    xticklabels(ax_cos, {'-pi', '-pi/2', '0', 'pi/2', 'pi'});
    xlabel(ax_cos, 'Phase (rad)');
    ylabel(ax_cos, '% events per bin');
    title(ax_cos, sprintf('ROI#%02d  Cosine fit', roi_id));
    axis(ax_cos, 'square');
    hold(ax_cos, 'off');
    
    % collect fitting parameter
    % ---- Pack parameters (per ROI) ----
    % KDE peak = max of KDE curve (your "gaussian_peak")
    [gaussian_peak_pct, idxMax] = max(f_kde_pct);
    gaussian_peak_phase = xi_kde(idxMax);
    
    cos_baseline  = B_fit;
    cos_amplitude = A_fit;
    cos_peak      = A_fit + B_fit;
    
    phase_params = struct();
    phase_params.roi_id = roi_id;
    phase_params.n      = nPh;
    phase_params.theta0 = theta0;
    phase_params.R      = R;
    
    phase_params.gaussian_peak_phase = gaussian_peak_phase;   % rad
    phase_params.gaussian_peak_pct   = gaussian_peak_pct;     % % per bin
    
    phase_params.kde_pks  = kde_pks(:);   % all detected KDE peaks (% per bin)
    phase_params.kde_locs = kde_locs(:);  % their locations (rad)
    
    phase_params.cos_baseline  = cos_baseline;   % %
    phase_params.cos_amplitude = cos_amplitude;  % %
    phase_params.cos_peak      = cos_peak;       % %
    
    % ---- Save ----
    outPdf6 = fullfile(folderPath, sprintf('phase_hist_ROI%02d.pdf', roi_id));
    exportgraphics(fig6, outPdf6, 'ContentType','vector', 'BackgroundColor','none');
    outFig6 = fullfile(folderPath, sprintf('phase_hist_ROI%02d.fig', roi_id));
    savefig(fig6, outFig6);
    fprintf('ROI#%02d : theta0=%.2f rad, R=%.3f, %d KDE peak(s)\n', ...
        roi_id, theta0, R, nKDE_peaks);
    outMat6 = fullfile(folderPath, sprintf('phase_params_ROI%02d.mat', roi_id));
    save(outMat6, 'phase_params');
    fprintf('ROI#%02d : theta0=%.2f rad, R=%.3f, %d KDE peak(s)\n', ...
        roi_id, theta0, R, nKDE_peaks);
end

%% ---- Figure 7: Average projection with ROI mask overlay ----
%avg_hits = dir(fullfile(folderPath, '*_minusDark_MC_AVG_for_CP.tif'));
avg_hits = dir(fullfile(folderPath, '*_minusDark_MC_AVG_for_CP.tif'));
avgImg   = double(imread(fullfile(folderPath, avg_hits(1).name)));

% Clip at 1%/99% percentile and normalise to [0 1]
lo = prctile(avgImg(:), 1);
hi = prctile(avgImg(:), 99);
avgImg = (avgImg - lo) / (hi - lo);
avgImg(avgImg < 0) = 0;
avgImg(avgImg > 1) = 1;

% Grayscale -> RGB
avgRGB = repmat(avgImg, [1 1 3]);

% Burn yellow ROI outlines into the image
maskL    = SAMload.maskL;
outlines = bwperim(maskL > 0);
yellow   = [1 1 0];
for ch = 1:3
    plane = avgRGB(:,:,ch);
    plane(outlines) = yellow(ch);
    avgRGB(:,:,ch) = plane;
end

% --- Save 1: no labels (just outlines) ---
outNoLabel = fullfile(folderPath, 'avg_projection_ROImask.png');
imwrite(avgRGB, outNoLabel);
fprintf('Saved avg projection (no labels):\n  %s\n', outNoLabel);

% --- Save 2: with ROI labels ---
[imgH, imgW, ~] = size(avgRGB);
dpi = 150;
fig7 = figure('Color','k', ...
    'Units','inches', 'Position',[1 1 imgW/dpi imgH/dpi]);
ax7 = axes(fig7, 'Units','normalized', 'Position',[0 0 1 1]);
imshow(avgRGB, 'Parent', ax7, 'Border','tight');
hold(ax7, 'on');

% Label only spiking ROIs — placed outside ROI, no overlap
offset_px  = 35;
dirs       = [0 -1; 1 -1; 1 0; 1 1;
              0  1;-1  1;-1 0;-1 -1];
dirs       = dirs ./ vecnorm(dirs, 2, 2);
txt_hw     = [20 12];
placed_boxes = zeros(0, 4);

for k = 1:nSpkROI
    roi_id = roiSpk_id(k);
    [r, c] = find(maskL == roi_id);
    if isempty(r), continue; end
    cy = mean(r);  cx = mean(c);

    outward = [cx - imgW/2, cy - imgH/2];
    outward = outward / (norm(outward) + eps);

    best_pos   = [cx, cy];
    best_score = -Inf;

    for d = 1:size(dirs, 1)
        tx = cx + offset_px * dirs(d, 1);
        ty = cy + offset_px * dirs(d, 2);
        if tx < 1 || tx > imgW || ty < 1 || ty > imgH, continue; end

        bx = [tx - txt_hw(1)/2, ty - txt_hw(2)/2, ...
              tx + txt_hw(1)/2, ty + txt_hw(2)/2];

        ov = 0;
        for p = 1:size(placed_boxes, 1)
            pb = placed_boxes(p, :);
            if bx(1)<pb(3) && bx(3)>pb(1) && bx(2)<pb(4) && bx(4)>pb(2)
                ov = ov + 1;
            end
        end

        txi = round(max(1, min(imgW, tx)));
        tyi = round(max(1, min(imgH, ty)));
        on_roi = maskL(tyi, txi) > 0;

        score = -ov*100 - on_roi*50 + dot(dirs(d,:), outward)*10;
        if score > best_score
            best_score = score;
            best_pos   = [tx, ty];
        end
    end

    text(ax7, best_pos(1), best_pos(2), sprintf('%d', roi_id), ...
        'Color','y', 'FontSize', 8, 'FontWeight','bold', ...
        'HorizontalAlignment','center', 'VerticalAlignment','middle');

    placed_boxes(end+1, :) = [best_pos(1)-txt_hw(1)/2, best_pos(2)-txt_hw(2)/2, ...
                              best_pos(1)+txt_hw(1)/2, best_pos(2)+txt_hw(2)/2]; %#ok<AGROW>
end
hold(ax7, 'off');

outLabeled = fullfile(folderPath, 'avg_projection_ROImask_labeled.png');
exportgraphics(fig7, outLabeled, 'ContentType','image', 'Resolution', 600);
fprintf('Saved avg projection (labeled):\n  %s\n', outLabeled);

%% ---- Save master file ----
outMat_master = fullfile(folderPath, 'breath_master_working.mat');
save(outMat_master, ...
    'folderPath', ...
    'fps_img', 'fps_breath', 'fs', ...
    'breath_peak', ...
    'roiSpk_id', 'nSpkROI', 'ca_spk_data', 'ca_spk_id', 't_spk', ...
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
