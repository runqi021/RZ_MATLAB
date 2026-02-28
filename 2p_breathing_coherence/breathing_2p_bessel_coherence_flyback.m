%%
clear; clc;

% load
SAM = "D:\RUNQI\phys\processed\breathing\pFN_roi4_z0_512x512_5x_6000f_00001\pFN_roi4_z0_512x512_5x_6000f_00001_ch1_minusDark_MC_cpSAM_output.mat";
dlc_csv = "D:\RUNQI\phys\processed\breathing\pFN_roi4_z0_512x512_5x_6000f_00001\Basler_acA1300-200um__22607866__20260224_023927400DLC_Resnet101_260224_vglut2_breathingFeb25shuffle1_snapshot_best-20.csv";

[folderPath, baseName, ~] = fileparts(SAM);

%% Load F
SAMload = load(SAM);
F_raw = SAMload.F;
% toss first second (system steady state）
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

%% ------------------- DLC BREATHING INPUT --------------------
% DLC CSV
use_dot = 'mean';    % 'dot1', 'dot2', or 'mean'

% DLC CSV has 3 header lines: scorer / bodyparts / coords
dlc_raw = readmatrix(dlc_csv, 'NumHeaderLines', 3);

data_dlc = dlc_raw;
data_dlc(1:30, :) = [];
frame_idx = data_dlc(:,1);

x1 = data_dlc(:,2);  y1 = data_dlc(:,3);  p1 = data_dlc(:,4); s1 = hypot(x1, y1);
x2 = data_dlc(:,5);  y2 = data_dlc(:,6);  p2 = data_dlc(:,7); s2 = hypot(x2, y2);
x3 = data_dlc(:,8);  y3 = data_dlc(:,9);  p3 = data_dlc(:,10); s3 = hypot(x3, y3);

% Raw time for breathing video (before alignment)
t_breath = frame_idx / fps_breath;

% % Choose breathing signal
% switch lower(use_dot)
%     case 'dot1'
%         breath = s1;
%     case 'dot2'
%         breath = s2;
%     otherwise
%         breath = (s1 + s2)/2;   % average of both dots
% end

breath = (x1 + x3)/2;

% Optional: use only high-likelihood samples, interpolate others
good = (p1 > 0.6) & (p2 > 0.6);
if any(good)
    breath(~good) = interp1(t_breath(good), breath(good), ...
                                 t_breath(~good), 'linear', 'extrap');
end

% detrend + z-score breathing
breath = detrend(breath);
breath = (breath - mean(breath)) / std(breath);

fprintf('Breathing trace: %d samples after dropping first %d frames\n', ...
        numel(breath), nDrop);

%% ---- Breathing event detection: inspiration onsets from derivative ----

fs = fps_img;              % imaging / breathing sampling rateke
T_img = numel(breath);     % length of breathing trace (after resample)

f_band = [0.1 2];

% Detect strong negative-going segments:
%    inspiration onset ≈ large NEGATIVE slope, so we look for peaks in -breath_vel

peak = breath;
% Choose thresholds
minProm = 2 * std(peak);      % tweak factor 0.3–1 if too many/few
minDist = round(fs / f_band(2));         % at least one cycle at the max breathing freq

[insp_amp, insp_onsets] = findpeaks(peak, ...
    'MinPeakProminence', minProm);

%[pk_vel, loc_vel] = findpeaks(peak, ...
%    'MinPeakProminence', minProm, ...
%    'MinPeakDistance',   minDist);

% --- binary spike/event trains S(t, ROI) ---
Sb = zeros(T_img,1);   % 0/1 events

Sb(insp_onsets) = 1;

% frame shift (real in_onset is at peak rise)
preFrames = 0;
insp_onsets = insp_onsets - preFrames;
insp_onsets(insp_onsets < 1) = 1;        % clamp to valid range

% 5) Remove events too close to edges for windowing later
win = round(1 * fs);                   % e.g. ±2.5 s window later
valid = (insp_onsets > win) & (insp_onsets < (T_img - win));
insp_onsets = insp_onsets(valid);
insp_amp = insp_amp(valid);

t_insp = t_img(insp_onsets);             % event times in seconds

fprintf('Detected %d inspiration-onset events.\n', numel(insp_onsets));

% 6) QC plot
figure('Name','Breathing event detection','Position',[200 200 800 600]);

plot(t_img, breath, 'k'); hold on;
plot(t_img(insp_onsets), breath(insp_onsets), 'ro', 'MarkerFaceColor','r');
xline(t_insp, 'r:');
xlabel('Time (s)');
ylabel('d/dt (breath)');
title('Breathing derivative with detected negative-slope peaks');
grid on;

%% Chronux spectrum
TW_spec = 4;
params_spec.Fs     = fps_img;
params_spec.tapers = [TW_spec, 2*TW_spec - 1];
params_spec.pad    = 0;
params_spec.fpass  = [0.1 15];
params_spec.err    = [2 0.05];

[Sk_raw, fk_raw, Sconfk_raw] = mtspectrumc(breath, params_spec);
figure('Color','White');
plot(fk_raw, Sk_raw,                'k',  'LineWidth', 1.2);
%hold on;
%plot(fk_raw, Sconfk_raw(1,:), 'k--','LineWidth', 0.6);
%hold on;
%plot(fk_raw, Sconfk_raw(2,:), 'k--','LineWidth', 0.6);
%hold on;
xlabel('Frequency (Hz)');
ylabel('Power');
xlim([0 10]);  
grid('on');
axis square;

% Save breathing spectrum
outFig = fullfile(folderPath, 'breathing_spectrum_raw.png')
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved breathing spectrum:\n  %s\n', outFig);

% %%
% [Sk_Sb, fk_Sb, Sconfk_Sb] = mtspectrumc(Sb, params_spec);
% 
% figure;
% plot(fk_Sb, Sk_Sb,                'k',  'LineWidth', 1.2);
% %hold on;
% %plot(fk_Sb, Sconfk_Sb(1,:), 'k--','LineWidth', 0.6);
% %hold on;
% %plot(fk_Sb, Sconfk_Sb(2,:), 'k--','LineWidth', 0.6);
% %hold on;
% xlabel('Frequency (Hz)');
% ylabel('Power');
% xlim([0 10]); 
% grid('on');

%%
keepIdx = [3 4 5 7 11];
F_keep = F(:, keepIdx);
dFF_keep = dFF(:, keepIdx);

%%
g = 0.93;
[F_keep_deconv, dFF_keep_deconv, spikes_oasis, baseline_oasis] = oasis_deconv_and_dff_AR1(F_keep, g);

helper.stackDFF(dFF_keep, dFF_keep_deconv);

%% 2p event detection
% assume: dFF is T x N, fs = 30 Hz
[T,N] = size(dFF_keep);
t = (1:1:T)/30;

% --- peak detection (per ROI) ---
pks  = cell(N,1);
locs = cell(N,1);

for i = 1:N
    
    [pks{i}, locs{i}] = findpeaks(dFF_keep(:, i), ...
        'MinPeakProminence', 0.25, 'MinPeakHeight', 0.25, 'MinPeakDistance', 8);
end

% --- binary spike/event trains S(t, ROI) ---
S = zeros(T, N);   % 0/1 events

for i = 1:N
    if ~isempty(locs{i})
        S(locs{i}, i) = 1;
    end
end

figure('Name','Peak detection (good ROIs only)','Color','w');
tiledlayout('flow');  % automatically fills panels

for i = 1:N
    
    nexttile;
    plot(t, dFF_keep(:, i), 'k'); hold on;

    if ~isempty(locs{i})
        plot(t(locs{i}), dFF_keep(locs{i}, i), 'ro', 'MarkerFaceColor','r', 'MarkerSize',3);
    end

    title(sprintf('ROI %d', i));   % keep original ROI index in title
    xlabel('Time (s)');
    ylabel('\DeltaF/F');
end

% Save peak-detection figure
%outFig = fullfile(tdir, sprintf('%s_peakDetection.pdf', tname));
%exportgraphics(gcf, outFig, 'ContentType','vector');
%fprintf('Saved peak-detection figure:\n  %s\n', outFig);

%% ------------------- CHRONUX PARAMS -------------------------
% Make sure Chronux is on path:
%   addpath(genpath('C:\path\to\chronux'));

f_band = [0.3 0.4];   % breathing band

TW_coh = 4;
params_coh.Fs     = fps_img;
params_coh.tapers = [TW_coh, 2*TW_coh-1];
params_coh.pad    = 0;
params_coh.fpass  = f_band;
params_coh.err    = [2 0.05];

%% ------------------- LOOP ROIs: PLV & PHASE -----------------

X = S;

plv_vals   = nan(N,1);
plv_ci_lo  = nan(N,1);
plv_ci_hi  = nan(N,1);
phase_vals = nan(N,1);
phase_ci   = nan(N,1);

for k = 1:N
    roi_trace = X(:,k);
    roi_trace = detrend(roi_trace);
    
    valid = ~isnan(roi_trace) & ~isnan(breath);
    roi_use    = roi_trace(valid);
    breath_use = Sb(valid);
    
    if numel(roi_use) < 100
        warning('ROI %s: too few valid samples, skipping.', roi_ids{k});
        continue;
    end
    
    [C12, C, phi, ~, ~, ~, f_coh, ~, ~, Cerr] = coherencyc(roi_use, breath_use, params_coh);
    idx_band = (f_coh >= f_band(1)) & (f_coh <= f_band(2));
    if ~any(idx_band)
        warning('ROI %s: no bins in [%.1f–%.1f] Hz.', roi_ids{k}, f_band(1), f_band(2));
        continue;
    end
    
    % PLV ≈ mean coherence magnitude in band
    plv_vals(k)  = mean(C(idx_band));
    plv_ci_lo(k) = mean(Cerr(1,idx_band));
    plv_ci_hi(k) = mean(Cerr(2,idx_band));
    
    % phase (jackknife across tapers)
    [mean_phase, phase_CI] = helper.phase_jackknife(C12);
    phase_vals(k) = mean_phase;
    phase_ci(k)   = phase_CI;
end

%% ------------------- NULL THRESHOLD -------------------------

alpha  = 0.001;
nROIs  = N;
K      = params_coh.tapers(2);
alphaf = 1 - (1-alpha)^nROIs;
Cth    = sqrt(1 - alphaf^(1/(K-1)));

fprintf('Null PLV threshold (alpha=%.3f, N=%d, K=%d): %.3f\n', alpha, nROIs, K, Cth);
is_sig = plv_vals > Cth;

% Color for NON-significant
col_nonsig = [0 0 0];   % black

% Color list for significant ROIs (MATLAB default)
defaultColors = get(groot, 'defaultAxesColorOrder');   % Nx3 matrix

sigIdx = find(is_sig);
nSig   = numel(sigIdx);

% Map each significant ROI to a color
% (cycling through MATLAB default colors if more than 7)
col_sig_list = zeros(nSig, 3);
for ii = 1:nSig
    col_sig_list(ii,:) = defaultColors( mod(ii-1, size(defaultColors,1)) + 1 , : );
end

%
roi_color = zeros(N, 3);  % RGB per ROI

for ii = 1:nSig
    roi_color(sigIdx(ii), :) = col_sig_list(ii,:);
end

% non-significant ROIs remain black
for k = 1:N
    if ~is_sig(k)
        roi_color(k,:) = col_nonsig;
    end
end

%% ------------------- FIGURE 1: POLAR PLOT -------------------

figure('Position', [200 200 500 500]);
pax = polaraxes; hold(pax,'on');

th = linspace(0,2*pi,200);
polarplot(pax, th, Cth*ones(size(th)), 'r--','LineWidth',1);  % null circle
set(gca,'ThetaZeroLocation','top');%,'ThetaDir','clockwise');

for k = 1:N
    if isnan(plv_vals(k)) || isnan(phase_vals(k))
        continue;
    end

    phi0    = phase_vals(k);
    plv_val = plv_vals(k);
    rlo     = max(0, plv_ci_lo(k));
    rhi     = min(1, plv_ci_hi(k));
    dphi    = phase_ci(k);

    phi_arc = linspace(phi0 - dphi, phi0 + dphi, 40);

    if is_sig(k)
        % -------- significant: use roi_color (MATLAB default palette) --------
        col = roi_color(k,:);   % assigned earlier for sig ROIs

        % arcs
        polarplot(pax, phi_arc, plv_val*ones(size(phi_arc)), '-', ...
                  'Color', col, 'LineWidth', 1.5);

        % radial bar
        polarplot(pax, [phi0 phi0], [rlo rhi], '-', ...
                  'Color', col, 'LineWidth', 2);

        % point
        polarplot(pax, phi0, plv_val, 'o', ...
                  'MarkerFaceColor', col, ...
                  'MarkerEdgeColor', 'k', ...
                  'MarkerSize', 7);
    else
        % -------- NON-significant: original gray style --------
        col = [0.6 0.6 0.6];

        polarplot(pax, phi_arc, plv_val*ones(size(phi_arc)), '-', ...
                  'Color', col, 'LineWidth', 0.7);

        polarplot(pax, [phi0 phi0], [rlo rhi], '-', ...
                  'Color', col, 'LineWidth', 1);

        polarplot(pax, phi0, plv_val, 'o', ...
                  'MarkerFaceColor', [0.85 0.85 0.85], ...
                  'MarkerEdgeColor', [0.3 0.3 0.3], ...
                  'MarkerSize', 6);
    end
end

rlim(pax,[0 1]);
thetalim([0 360]);
title(pax, sprintf('Coherence (%.2f–%.2f Hz)', f_band(1), f_band(2)));
hold(pax,'off');

%% ------------------- FIGURE 2: BREATH + STACKED dFF (Top 10 ROIs) --------

plv_vals(isnan(plv_vals)) = 0;

% sort PLV in DESCENDING order (highest first)
[~, idx_sort] = sort(plv_vals, 'descend');

% ---- choose only top 10 ROIs ----
nTop = min(10, N);
idx_top = idx_sort(1:nTop);

% extract only these ROIs
X_top = dFF_keep(:, idx_top);     % [T_img x nTop]

% global amplitude
amp = max(X_top(:)) - min(X_top(:));
if amp == 0, amp = 1; end
gap = 0.8 * amp;

% build stacked matrix (highest at top)
Y = zeros(T_img, nTop);
for k = 1:nTop
    Y(:,k) = X_top(:,k) + gap*(nTop - k);
end

% y-ticks and labels
yt  = 0:gap:gap*(nTop-1);                    % positions (bottom → top)
ytl = roi_ids(idx_top(end:-1:1));            % labels (bottom = lowest PLV of top 10)

figure('Position',[200 200 900 700]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% ---- TOP PANEL: breathing + breathing events (Sb) ----
ax1 = nexttile(1);
plot(ax1, t_img, breath, 'k', 'LineWidth', 1); hold(ax1,'on');

% --- breathing spike overlay (Sb) ---
% assume Sb is [T_img x 1] or [1 x T_img], binary 0/1
idx_b = find(Sb(:) ~= 0);
if ~isempty(idx_b)
    t_b = t_img(idx_b);
    y_b = breath(idx_b);
    plot(ax1, t_b, y_b, 'ro', 'MarkerFaceColor','r', ...
         'MarkerSize', 3);   % breathing event markers
end

xlabel(ax1, 'Time (s)');
ylabel(ax1, 'Breathing (z)');
title(ax1, 'Breathing trace (DLC) with detected events');
grid(ax1,'on');
hold(ax1,'off');

% ---- BOTTOM PANEL: stacked top-10 ROIs + ROI spike events (S) ----
ax2 = nexttile(2, [2 1]); 
hold(ax2,'on');

for kk = 1:nTop
    roi_idx = idx_top(kk);
    col = roi_color(roi_idx,:);   % ★ consistent color

    % main trace
    plot(ax2, t_img, Y(:,kk), ...
         'Color', col, ...
         'LineWidth', 0.5);

    % --- ROI spike overlay using S ---
    idx_evt = find(S(:, roi_idx) ~= 0);   % S is [T_img x N]
    if ~isempty(idx_evt)
        t_ev = t_img(idx_evt);
        y_ev = Y(idx_evt, kk)+0.03;           % already offset version
        plot(ax2, t_ev, y_ev, '.', ...
             'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', 'r', ...
             'MarkerSize', 8);
    end
end

xlim(ax2, [t_img(1) t_img(end)]);
ylim(ax2, [-gap, gap*(nTop-1)+gap]);

yticks(ax2, yt);
yticklabels(ax2, ytl);

xlabel(ax2, 'Time (s)');
ylabel(ax2, 'ROI#');
title(ax2, sprintf('Top %d ROI traces (\\DeltaF/F, PLV-sorted) with events', nTop));
grid(ax2,'on');

% SCALE BAR
dFF_scale = 0.2;     % 0.2 ΔF/F vertical bar

% position vertical bar near right
x0 = 1*(t_img(end));
y0 = gap*(nTop-1) + 0.2*gap;

plot(ax2, [x0 x0], [dFF_scale 2*dFF_scale], 'k', 'LineWidth', 2);
text(ax2, x0+0.5, dFF_scale + 0.05*dFF_scale, ...
    sprintf('%.2f \\DeltaF/F', dFF_scale), ...
    'VerticalAlignment','bottom','HorizontalAlignment','center', ...
    'FontSize', 10);

hold(ax2,'off');

%% ---- Event-triggered segments per significant ROI (RAW dFF) ----

fs   = fps_img;
win  = round(2.5 * fs);          % ±2.5 s
t_evt = (-win:win) / fs;

% Make sure insp_onsets are valid for this window
T_img = size(dFF_keep,1);
valid = (insp_onsets > win) & (insp_onsets < (T_img - win));
insp_onsets = insp_onsets(valid);
nEv = numel(insp_onsets);

fprintf('Event-triggered : using %d inspiration onsets.\n', nEv);

sigIdx = find(is_sig);
nSig   = numel(sigIdx);

if nSig == 0
    warning('No significant ROIs; skipping event-triggered overlays.');
else
    % dff_evt_seg_all{ii} : nEv x (2*win+1) for ROI sigIdx(ii)
    dff_evt_seg_all = cell(nSig,1);

    for ii = 1:nSig
        rr = sigIdx(ii);
        roi_trace = dFF_keep(:, rr);          % RAW ΔF/F
        roi_trace = fillmissing(roi_trace,'linear');

        seg = zeros(nEv, 2*win+1);
        for e = 1:nEv
            c = insp_onsets(e);
            seg(e,:) = roi_trace(c-win : c+win);
        end

        dff_evt_seg_all{ii} = seg;
    end
end

% ---------- Per-ROI event-triggered overlays (all trials + mean, RAW dFF) ----------

if nSig > 0
    if ~exist('cmap','var') || size(cmap,1) < nSig
        cmap = lines(nSig);
    end

    for ii = 1:nSig
        rr  = sigIdx(ii);             % original ROI index
        seg = dff_evt_seg_all{ii};    % nEv x (2*win+1)
        if isempty(seg)
            continue;
        end

        % average for this ROI (RAW ΔF/F)
        dff_avg = mean(seg, 1);

        % lighter color for single trials
        col_mean  = cmap(ii,:);
        col_light = 0.7*col_mean + 0.3*[1 1 1];

        figure('Name', sprintf('Event-triggered dFF: ROI# %s', roi_ids{rr}), ...
               'Position',[450 200 700 400]);
        hold on;

        % all individual trials
        plot(t_evt, seg', 'Color', col_light, 'LineWidth', 0.7);

        % overlay mean
        plot(t_evt, dff_avg, 'Color', col_mean, 'LineWidth', 2);

        % reference lines
        xline(0, '--', 'Color', [0.4 0.4 0.4]);  % inspiration onset at 0
        yline(0, ':',  'Color', [0.6 0.6 0.6]);

        xlabel('Time from inspiration onset (s)');
        ylabel('\DeltaF/F (raw units)');
        title(sprintf('Event-triggered dFF: ROI# %s (n=%d breaths)', ...
                      roi_ids{rr}, size(seg,1)));

        grid on;
        hold off;
    end
end
