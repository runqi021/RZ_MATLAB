%%
close all; clear; clc;

folderPath = "D:\RUNQI\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi2_z5_512x512_6x_6000f_00001";

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

%breath = (x1 + x3)/2;
breath = (x1+x2+x3)/3;
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

%% ---- Breathing event detection: fine-tune parameters ----

% =========================================================
% USER-EDITABLE PEAK DETECTION PARAMETERS
% =========================================================
minProm_factor = 1.5;      % Prominence = factor × std(breath).
                            %   Increase → fewer peaks.  Decrease → more peaks.
minDist_sec    = 0.3;      % Min seconds between consecutive peaks.
minHeight      = -Inf;     % Min absolute peak height (z-score). -Inf = disabled.
minWidth_sec   = 0;        % Min peak width in seconds. 0 = disabled.
snap_win_sec   = 0.2;     % Manual-add snaps to nearest local max within this window (s).
preFrames      = 0;        % Shift all detected peaks earlier by this many frames.
% =========================================================

fs    = fps_img;
T_img = numel(breath);

peak     = breath;
minProm  = minProm_factor * std(peak);
minDist  = round(minDist_sec * fs);
minWidth = round(minWidth_sec * fs);
snap_win = round(snap_win_sec * fs);

% Build findpeaks argument list dynamically
fp_args = {'MinPeakProminence', minProm, 'MinPeakDistance', minDist};
if isfinite(minHeight)
    fp_args = [fp_args, {'MinPeakHeight', minHeight}];
end
if minWidth > 0
    fp_args = [fp_args, {'MinPeakWidth', minWidth}];
end

[insp_amp, insp_onsets] = findpeaks(peak, fp_args{:});

% Frame shift
insp_onsets = insp_onsets - preFrames;
insp_onsets(insp_onsets < 1) = 1;

fprintf('Auto-detected %d inspiration-onset events.\n', numel(insp_onsets));

% Freeze auto result (used for saving and audit trail)
insp_onsets_auto = sort(insp_onsets(:));

%% ---- Initial auto-detection QC plot ----

figure('Name', 'Breathing event detection (auto)', 'Position', [200 200 1000 400]);
plot(t_img, breath, 'k'); hold on;
plot(t_img(insp_onsets_auto), breath(insp_onsets_auto), 'ro', 'MarkerFaceColor', 'r');
xline(t_img(insp_onsets_auto), 'r:');
xlabel('Time (s)'); ylabel('Breath (z)');
title(sprintf('Auto-detected peaks: n=%d  |  minProm=%.2f  minDist=%.2f s', ...
    numel(insp_onsets_auto), minProm, minDist_sec));
grid on; hold off;

%% ---- Interactive peak editor ----
% Left-click  = add a peak (snaps to nearest local max within snap_win_sec)
% Right-click = delete nearest existing peak
% Press Enter = finish editing

insp_onsets_added   = [];            % frame indices manually added
insp_onsets_deleted = [];            % frame indices removed (from auto set)
insp_onsets_edit    = insp_onsets_auto;   % working copy

fig_edit = figure('Name', 'Peak Editor | Left-click=Add  Right-click=Delete  Enter=Done', ...
    'Position', [100 80 1100 480]);
ax_edit = axes('Parent', fig_edit);

fprintf('\n=== INTERACTIVE PEAK EDITOR ===\n');
fprintf('  Left-click  : add peak (snaps to local max within %.2f s)\n', snap_win_sec);
fprintf('  Right-click : delete nearest existing peak\n');
fprintf('  Press Enter : finish editing\n\n');

while ishandle(fig_edit)
    % Refresh plot
    cla(ax_edit);
    plot(ax_edit, t_img, breath, 'k', 'LineWidth', 1);
    hold(ax_edit, 'on');
    if ~isempty(insp_onsets_edit)
        plot(ax_edit, t_img(insp_onsets_edit), breath(insp_onsets_edit), ...
            'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    end
    xlabel(ax_edit, 'Time (s)'); ylabel(ax_edit, 'Breath (z)');
    title(ax_edit, sprintf('Peak Editor  |  n=%d  (+%d / -%d from auto)  |  L=Add  R=Delete  Enter=Done', ...
        numel(insp_onsets_edit), numel(insp_onsets_added), numel(insp_onsets_deleted)));
    grid(ax_edit, 'on'); hold(ax_edit, 'off');
    drawnow;

    % Wait for user click
    try
        [x_click, ~, button] = ginput(1);
    catch
        break;
    end
    if isempty(x_click) || ~ishandle(fig_edit)
        break;   % Enter pressed → done
    end

    % Convert clicked x (time in s) → nearest frame index
    [~, fc] = min(abs(t_img - x_click));

    if button == 1   % Left-click: add peak
        lo = max(1, fc - snap_win);
        hi = min(numel(breath), fc + snap_win);
        [~, rel_idx] = max(breath(lo:hi));
        new_frame = lo + rel_idx - 1;
        if ~ismember(new_frame, insp_onsets_edit)
            insp_onsets_edit = sort([insp_onsets_edit(:); new_frame]);
            insp_onsets_added(end+1) = new_frame;
            insp_onsets_deleted(insp_onsets_deleted == new_frame) = [];
            fprintf('Added   frame %5d  (t = %.2f s)\n', new_frame, t_img(new_frame));
        else
            fprintf('Peak at frame %d already exists — skipped.\n', new_frame);
        end

    elseif button == 3   % Right-click: delete nearest
        if ~isempty(insp_onsets_edit)
            [~, ni] = min(abs(t_img(insp_onsets_edit) - x_click));
            del_frame = insp_onsets_edit(ni);
            insp_onsets_edit(ni) = [];
            if ismember(del_frame, insp_onsets_auto)
                insp_onsets_deleted(end+1) = del_frame;
            end
            insp_onsets_added(insp_onsets_added == del_frame) = [];
            fprintf('Deleted frame %5d  (t = %.2f s)\n', del_frame, t_img(del_frame));
        end
    end
end

% Finalise
insp_onsets_final   = sort(insp_onsets_edit(:));
insp_onsets_added   = sort(insp_onsets_added(:));
insp_onsets_deleted = sort(insp_onsets_deleted(:));

fprintf('\nEditing complete: %d final peaks  (+%d added, -%d deleted from %d auto)\n', ...
    numel(insp_onsets_final), numel(insp_onsets_added), ...
    numel(insp_onsets_deleted), numel(insp_onsets_auto));

% Overwrite insp_onsets so all downstream code is unaffected
insp_onsets = insp_onsets_final;
t_insp      = t_img(insp_onsets);

% Recompute binary event train from final peaks
Sb = zeros(T_img, 1);
Sb(insp_onsets) = 1;

%% ---- Save: final peak detection figure + .mat ----

fig_qc_final = figure('Name', 'Breathing peak detection (final)', ...
    'Position', [200 200 1000 400]);
plot(t_img, breath, 'k', 'LineWidth', 1); hold on;
if ~isempty(insp_onsets_auto)
    plot(t_img(insp_onsets_auto), breath(insp_onsets_auto), ...
        'o', 'Color', [0.65 0.65 0.65], 'MarkerFaceColor', [0.65 0.65 0.65], ...
        'MarkerSize', 4, 'DisplayName', 'Auto peaks');
end
if ~isempty(insp_onsets_added)
    plot(t_img(insp_onsets_added), breath(insp_onsets_added), ...
        'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'Manually added');
end
if ~isempty(insp_onsets_final)
    plot(t_img(insp_onsets_final), breath(insp_onsets_final), ...
        'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4, 'DisplayName', 'Final accepted');
    xline(t_img(insp_onsets_final), 'r:');
end
legend('Location', 'best');
xlabel('Time (s)'); ylabel('Breath (z)');
title(sprintf('Final: %d peaks  (auto=%d, +%d added, -%d deleted)', ...
    numel(insp_onsets_final), numel(insp_onsets_auto), ...
    numel(insp_onsets_added), numel(insp_onsets_deleted)));
grid on; hold off;

outFig_peaks = fullfile(folderPath, 'breath_peak_detection_final.png');
exportgraphics(fig_qc_final, outFig_peaks, 'ContentType', 'vector');
fprintf('Saved final peak detection figure:\n  %s\n', outFig_peaks);

% findpeaks parameter struct
findpeak_params = struct( ...
    'minProm_factor',  minProm_factor, ...
    'minProm',         minProm, ...
    'minDist_sec',     minDist_sec, ...
    'minDist_frames',  minDist, ...
    'minHeight',       minHeight, ...
    'minWidth_sec',    minWidth_sec, ...
    'minWidth_frames', minWidth, ...
    'snap_win_sec',    snap_win_sec, ...
    'preFrames',       preFrames, ...
    'fps',             fs ...
);

outMat_peaks = fullfile(folderPath, 'breath_peak_data.mat');
save(outMat_peaks, ...
    'findpeak_params', ...
    'insp_onsets_auto', ...
    'insp_onsets_added', ...
    'insp_onsets_deleted', ...
    'insp_onsets_final');
fprintf('Saved peak detection data:\n  %s\n', outMat_peaks);

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

%% ------------------- FIGURE 1: BREATH + STACKED dFF (Top 10 ROIs) --------

% plv_vals(isnan(plv_vals)) = 0;
% % sort PLV in DESCENDING order (highest first)
% [~, idx_sort] = sort(plv_vals, 'descend');
% % ---- choose only top 10 ROIs ----
% nTop = min(10, N);
% idx_top = idx_sort(1:nTop);
% % extract only these ROIs
% X_top = dFF(:, idx_top);           % [T_img x nTop]

% global amplitude
amp = max(dFF(:)) - min(dFF(:));
if amp == 0, amp = 1; end; gap = 0.8 * amp;

% build stacked matrix (highest at top)
Y = zeros(T_img, N_roi);
for k = 1:N_roi
    Y(:,k) = dFF(:,k) + gap*(N_roi - k);
end

% y-ticks and labels
yt  = 0:gap:gap*(N_roi-1);                    % positions (bottom → top)
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

    % main trace
    plot(ax2, t_img, Y(:,kk), 'Color', 'k', 'LineWidth', 0.5);

    % % --- ROI spike overlay using S ---
    % idx_evt = find(S(:, roi_idx) ~= 0);   % S is [T_img x N]
    % if ~isempty(idx_evt)
    %     t_ev = t_img(idx_evt);
    %     y_ev = Y(idx_evt, kk)+0.03;           % already offset version
    %     plot(ax2, t_ev, y_ev, '.', ...
    %          'MarkerFaceColor', 'w', ...
    %          'MarkerEdgeColor', 'r', ...
    %          'MarkerSize', 8);
    % end
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
dFF_scale = 0.2;     % 0.2 ΔF/F vertical bar

% position vertical bar near right
x0 = 1*(t_img(end));
y0 = gap*(N_roi-1) + 0.2*gap;

plot(ax2, [x0 x0], [dFF_scale 2*dFF_scale], 'k', 'LineWidth', 2);
text(ax2, x0+0.5, dFF_scale + 0.05*dFF_scale, ...
    sprintf('%.2f \\DeltaF/F', dFF_scale), ...
    'VerticalAlignment','bottom','HorizontalAlignment','center', ...
    'FontSize', 10);
%xline(t_b, 'LineWidth',0.1);
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
    plot(ax2b, t_img, Y(:,kk), 'Color', 'k', 'LineWidth', 0.5);
end

xline(ax2b, t_b, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.1);

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


%% ---- Figure 3. Event-triggered segments per significant ROI (RAW dFF) ----

fs   = fps_img;
win  = round(0.5 * fs);          % ±2.5 s
t_evt = (-win:win) / fs;

T_img = size(dFF,1);
insp_ev = insp_onsets(insp_onsets > win & insp_onsets < T_img - win);  % local only
nEv = numel(insp_ev);

fprintf('Event-triggered : using %d inspiration onsets.\n', nEv);

dff_evt_seg_all = cell(N_roi,1);

for i = 1:N_roi
    roi_trace = dFF(:, i);               % RAW ΔF/F

    seg = zeros(nEv, 2*win+1);
    for e = 1:nEv
        c = insp_ev(e);
        seg(e,:) = roi_trace(c-win : c+win);
    end

    dff_evt_seg_all{i} = seg;
end

% ---------- Per-ROI event-triggered overlays (all trials + mean, RAW dFF) ----------

figure('Name', (sprintf('Event-triggered dFF (n=%d breaths)', size(seg,1))));

for i = 1:N_roi
    seg = dff_evt_seg_all{i};    % nEv x (2*win+1)
    if isempty(seg)
        continue;
    end

    % average for this ROI (RAW ΔF/F)
    dff_avg = mean(seg, 1);

    hold on;

    subplot(ceil(N_roi/10), 10, i);
    % all individual trials
    plot(t_evt, seg', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.7); hold on;
    % overlay mean
    plot(t_evt, dff_avg, 'Color', 'k', 'LineWidth', 2);

    % reference lines
    xline(0, '--', 'Color', [0.4 0.4 0.4]);  % inspiration onset at 0
    yline(0, '--',  'Color', [0.6 0.6 0.6]);

    xlabel('Time from inspiration onset (s)');
    ylabel('\DeltaF/F');
    title(sprintf('ROI#%d', i));

    %grid on;
    hold off;
    axis square;

end

outFig = fullfile(folderPath, 'triggered_avg.png')
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved triggered average:\n  %s\n', outFig);

%% ---- Figure 4. Calcium-triggered average of breathing (selected ROIs) ----

% Select ROIs to run peak detection on (only well-isolated cells)
keepIdx  = [6];          % edit to choose which ROIs
dFF_keep = dFF(:, keepIdx);
nKeep    = numel(keepIdx);

win_ca   = round(1.0 * fs);          % ±1 s window around calcium peak
t_ca_evt = (-win_ca:win_ca) / fs;

% Peak detection params — tweak per dataset
minProm_ca = 0.2;                     % prominence threshold in dF/F units
minDist_ca = round(fs * 0.2);         % min 1 s between peaks (refractory)

breath_ca_seg_all = cell(nKeep, 1);
ca_onsets_all     = cell(nKeep, 1);

for i = 1:nKeep
    roi_trace = dFF_keep(:, i);

    [~, ca_onsets] = findpeaks(roi_trace, ...
        'MinPeakProminence', minProm_ca, ...
        'MinPeakDistance',   minDist_ca);

    ca_onsets_all{i} = ca_onsets;   % save full — used by Fig 5, 6, 7

    % local filter only for windowed segment extraction
    ca_onsets_ev = ca_onsets(ca_onsets > win_ca & ca_onsets < T_img - win_ca);
    nEv_ca = numel(ca_onsets_ev);

    fprintf('ROI#%02d : %d calcium peaks detected.\n', keepIdx(i), nEv_ca);

    if nEv_ca == 0
        breath_ca_seg_all{i} = [];
        continue;
    end

    seg_ca = zeros(nEv_ca, 2*win_ca+1);
    for e = 1:nEv_ca
        c = ca_onsets_ev(e);
        seg_ca(e,:) = breath(c-win_ca : c+win_ca);
    end
    breath_ca_seg_all{i} = seg_ca;
end

% ---------- QC: raw dF/F with detected calcium peaks ----------

figure('Name', 'Calcium peak detection QC', 'Color', 'w');
tiledlayout('flow');

for i = 1:nKeep
    ca_onsets = ca_onsets_all{i};
    roi_trace = dFF_keep(:, i);

    nexttile;
    plot(t_img, roi_trace, 'k'); hold on;
    if ~isempty(ca_onsets)
        plot(t_img(ca_onsets), roi_trace(ca_onsets), 'ro', ...
            'MarkerFaceColor', 'r', 'MarkerSize', 3);
    end
    title(sprintf('ROI#%d  (n=%d)', keepIdx(i), numel(ca_onsets)));
    xlabel('Time (s)');
    ylabel('\DeltaF/F');
end

outFig = fullfile(folderPath, sprintf('ca_peak_detection_ROI%02d.png',keepIdx(i)));
exportgraphics(gcf, outFig, 'ContentType', 'vector');
fprintf('Saved calcium peak detection:\n  %s\n', outFig);

%%
% ---------- Per-ROI calcium-triggered breathing overlays ----------

figure('Name', sprintf('Calcium-triggered breathing avg (n=%d ROIs)', nKeep));

for i = 1:nKeep
    seg_ca = breath_ca_seg_all{i};
    if isempty(seg_ca)
        continue;
    end

    breath_avg = mean(seg_ca, 1);
    nEv_ca     = size(seg_ca, 1);

    subplot(ceil(nKeep/6), 6, i);
    plot(t_ca_evt, seg_ca',   'Color', [0.8 0.8 0.8], 'LineWidth', 0.7); hold on;
    plot(t_ca_evt, breath_avg, 'Color', 'k',           'LineWidth', 2);

    xline(0, '--', 'Color', [0.4 0.4 0.4]);   % calcium peak at t = 0
    yline(0, '--', 'Color', [0.6 0.6 0.6]);

    xlabel('Time from Ca^{2+} peak (s)');
    ylabel('Breathing (z)');
    title(sprintf('ROI#%d  (n=%d)', keepIdx(i), nEv_ca));

    hold off;
    axis square;
end


outFig = fullfile(folderPath, 'ca_triggered_breath_avg.png');
exportgraphics(gcf, outFig, 'ContentType', 'vector');
outFig = fullfile(folderPath, 'ca_triggered_breath_avg.fig');
savefig(gcf, outFig)
fprintf('Saved calcium-triggered breathing average:\n  %s\n', outFig);

%% Figure 5: Breathing raster around calcium events (±5 s)

win5    = round(5 * fs);            % ±5 s = 150 frames at 30 Hz
t_rast  = (-win5:win5) / fs;        % time axis for raster (-5 → +5 s)
dt_latency_all = cell(nKeep, 1);   % accumulate per-ROI latency tables

for i = 1:nKeep

    ca_ev = ca_onsets_all{i};
    nEv_r = numel(ca_ev);

    if nEv_r == 0
        fprintf('ROI#%02d : no calcium events with valid ±5 s window, skipping raster.\n', keepIdx(i));
        continue;
    end

    % Find inspiration onsets that fall inside each window (relative time)
    insp_rel = cell(nEv_r, 1);
    for e = 1:nEv_r
        c      = ca_ev(e);
        nearby = insp_onsets(insp_onsets >= (c-win5) & insp_onsets <= (c+win5));
        insp_rel{e} = (nearby - c) / fs;   % seconds relative to Ca2+ event
    end

    % ---- Raster plot: inspiration ticks around each Ca2+ event ----
    figure('Name', sprintf('Breathing raster: ROI#%d (n=%d events)', keepIdx(i), nEv_r), ...
           'Color', 'w');
    hold on;

    for e = 1:nEv_r
        t_ticks = insp_rel{e};             % inspiration times relative to Ca2+ event
        for tt = 1:numel(t_ticks)
            plot([t_ticks(tt) t_ticks(tt)], [e-0.4  e+0.4], ...
                'k', 'LineWidth', 1);
        end
    end

    xline(0, 'r--', 'LineWidth', 1.5);    % calcium event at t = 0

    xlim([-5 5]);
    ylim([0 nEv_r+1]);
    xlabel('Time from Ca^{2+} event (s)');
    ylabel('Ca^{2+} event #');
    title(sprintf('Inspiration raster: ROI#%d  (n=%d Ca^{2+} events)', keepIdx(i), nEv_r));
    set(gca, 'YDir', 'reverse');          % event 1 at top
    hold off;

    outFig = fullfile(folderPath, sprintf('ca_breath_raster_ROI%02d.png', keepIdx(i)));
    exportgraphics(gcf, outFig, 'ContentType', 'vector');
    fprintf('Saved breathing raster:\n  %s\n', outFig);

    % ---- Sorted raster: sort rows by latency to first following breath ----

    % For each event, find the first inspiration onset after t=0
    first_insp_lat = nan(nEv_r, 1);
    for e = 1:nEv_r
        post = insp_rel{e}(insp_rel{e} > 0);   % onsets after calcium event
        if ~isempty(post)
            first_insp_lat(e) = min(post);
        end
    end

    % Sort by latency (NaN — no following breath in window — go to bottom)
    [~, sort_idx] = sort(first_insp_lat, 'ascend', 'MissingPlacement', 'last');

    figure('Name', sprintf('Breathing raster (sorted): ROI#%d (n=%d events)', keepIdx(i), nEv_r), ...
           'Color', 'w');
    hold on;

    for row = 1:nEv_r
        e      = sort_idx(row);
        t_ticks = insp_rel{e};
        for tt = 1:numel(t_ticks)
            plot([t_ticks(tt) t_ticks(tt)], [row-0.4  row+0.4], ...
                'k', 'LineWidth', 1);
        end
    end

    xline(0, 'r--', 'LineWidth', 1.5);

    xlim([-5 5]);
    ylim([0 nEv_r+1]);
    xlabel('Time from Ca^{2+} event (s)');
    ylabel('Ca^{2+} event # (sorted by latency to first breath)');
    title(sprintf('Inspiration raster (sorted): ROI#%d  (n=%d Ca^{2+} events)', keepIdx(i), nEv_r));
    set(gca, 'YDir', 'reverse');
    hold off;

    outFig = fullfile(folderPath, sprintf('ca_breath_raster_sorted_ROI%02d.png', keepIdx(i)));
    exportgraphics(gcf, outFig, 'ContentType', 'vector');
    fprintf('Saved sorted breathing raster:\n  %s\n', outFig);

    % ---- Sorted raster: sort by latency of last preceding breath ----

    % For each event, find the last inspiration onset before t=0
    last_insp_lat = nan(nEv_r, 1);
    for e = 1:nEv_r
        pre = insp_rel{e}(insp_rel{e} < 0);    % onsets before calcium event
        if ~isempty(pre)
            last_insp_lat(e) = max(pre);        % closest to 0 from the left
        end
    end

    % Sort descending: least negative (closest preceding breath) at top
    [~, sort_idx_pre] = sort(last_insp_lat, 'descend', 'MissingPlacement', 'last');

    figure('Name', sprintf('Breathing raster (pre-sorted): ROI#%d (n=%d events)', keepIdx(i), nEv_r), ...
           'Color', 'w');
    hold on;

    for row = 1:nEv_r
        e       = sort_idx_pre(row);
        t_ticks = insp_rel{e};
        for tt = 1:numel(t_ticks)
            plot([t_ticks(tt) t_ticks(tt)], [row-0.4  row+0.4], ...
                'k', 'LineWidth', 1);
        end
    end

    xline(0, 'r--', 'LineWidth', 1.5);

    xlim([-5 5]);
    ylim([0 nEv_r+1]);
    xlabel('Time from Ca^{2+} event (s)');
    ylabel('Ca^{2+} event # (sorted by latency of preceding breath)');
    title(sprintf('Inspiration raster (pre-sorted): ROI#%d  (n=%d Ca^{2+} events)', keepIdx(i), nEv_r));
    set(gca, 'YDir', 'reverse');
    hold off;

    outFig = fullfile(folderPath, sprintf('ca_breath_raster_presorted_ROI%02d.png', keepIdx(i)));
    exportgraphics(gcf, outFig, 'ContentType', 'vector');
    fprintf('Saved pre-sorted breathing raster:\n  %s\n', outFig);

    % ---- Save latency table: [nEv x 3] = [dt_pre, dt_post, dt_closest] ----
    % dt_pre     : last preceding breath (s, negative)
    % dt_post    : first following breath (s, positive)
    % dt_closest : breath spike closest to calcium event (signed, pre or post)
    dt_closest = nan(nEv_r, 1);
    for e = 1:nEv_r
        all_t = insp_rel{e};
        if ~isempty(all_t)
            [~, idx_min] = min(abs(all_t));
            dt_closest(e) = all_t(idx_min);
        end
    end

    dt_latency = [last_insp_lat, first_insp_lat, dt_closest];
    dt_latency_all{i} = dt_latency;

    % ---- Sorted raster: sort by dt_closest (signed) ----
    [~, sort_idx_cl] = sort(dt_closest, 'ascend', 'MissingPlacement', 'last');

    figure('Name', sprintf('Breathing raster (closest-sorted): ROI#%d (n=%d events)', keepIdx(i), nEv_r), ...
           'Color', 'w');
    hold on;

    for row = 1:nEv_r
        e       = sort_idx_cl(row);
        t_ticks = insp_rel{e};
        for tt = 1:numel(t_ticks)
            plot([t_ticks(tt) t_ticks(tt)], [row-0.4  row+0.4], ...
                'k', 'LineWidth', 1);
        end
    end

    xline(0, 'r--', 'LineWidth', 1.5);

    xlim([-5 5]);
    ylim([0 nEv_r+1]);
    xlabel('Time from Ca^{2+} event (s)');
    ylabel('Ca^{2+} event # (sorted by closest breath)');
    title(sprintf('Inspiration raster (closest-sorted): ROI#%d  (n=%d Ca^{2+} events)', keepIdx(i), nEv_r));
    set(gca, 'YDir', 'reverse');
    hold off;

    outFig = fullfile(folderPath, sprintf('ca_breath_raster_closest_ROI%02d.png', keepIdx(i)));
    exportgraphics(gcf, outFig, 'ContentType', 'vector');
    fprintf('Saved closest-sorted breathing raster:\n  %s\n', outFig);

    outMat = fullfile(folderPath, sprintf('ca_breath_latency_ROI%02d.mat', keepIdx(i)));
    save(outMat, 'dt_latency');
    fprintf('Saved latency table (%d events x 3):\n  %s\n', nEv_r, outMat);

end


%% Figure 7: t_pre shuffle test (scan statistic)
%
%  t_pre,j = c_j - b_k  where [b_k, b_{k+1}) is the breath interval
%  containing calcium event c_j.  Units: seconds, range [0, Delta_k).
%
%  Scan statistic: S(w) = max_tau #{j : |t_pre,j - tau| < w/2}
%  — max events clustered within any ±w/2 window along the t_pre axis.
%
%  Null: for each interval [b_k, b_{k+1}) with n_k real events, draw n_k
%  surrogates ~ Uniform(0, Delta_k). Preserves count-per-interval; destroys
%  any preferred latency.
%
%  p = (1 + #{S'_m >= S_obs}) / (M + 1)   [Phipson & Smyth 2010]

w_sec  = 0.15;   % scan window half-width (s). Use 0.10 for ~100 ms jitter,
                 % 0.20 for ~200 ms. Tune to your expected coupling timescale.
M_shuf = 3000;   % shuffle count (2000-5000 for publication)
roi_stats = struct();              % accumulate per-ROI permutation results

for i = 1:nKeep

    ca_ev    = ca_onsets_all{i};
    b_frames = sort(insp_onsets(:));

    if isempty(ca_ev)
        fprintf('ROI#%02d : no Ca2+ events, skipping t_pre shuffle test.\n', keepIdx(i));
        continue;
    end

    % ---- Assign each Ca2+ event to a breath interval ----
    nBreaths = numel(b_frames);

    % histcounts: bin_idx = k  =>  ca_ev(j) in [b_frames(k), b_frames(k+1))
    % bin_idx = 0 for events before b_frames(1) or after b_frames(end)
    [~, ~, bin_idx] = histcounts(ca_ev, b_frames);

    keep_ca = (bin_idx >= 1) & (bin_idx <= nBreaths - 1);
    ca_in   = ca_ev(keep_ca);
    bin_in  = bin_idx(keep_ca);
    nEv_in  = numel(ca_in);

    if nEv_in == 0
        fprintf('ROI#%02d : no Ca2+ events inside breath intervals, skipping.\n', keepIdx(i));
        continue;
    end

    Delta_k = diff(b_frames) / fs;             % interval durations (s) [nBreaths-1 x 1]
    t_pre   = (ca_in - b_frames(bin_in)) / fs; % observed t_pre (s)    [nEv_in x 1]
    t_post  = Delta_k(bin_in) - t_pre;         % time to next breath    [nEv_in x 1]

    n_pre_w    = sum(t_pre  <= w_sec);          % events with t_pre  within fixed threshold
    n_post_w   = sum(t_post <= w_sec);          % events with t_post within fixed threshold

    % ---- Observed scan statistic ----
    S_obs = breath_scan_stat(t_pre, w_sec);

    % ---- Histogram grid for visual envelope (set once before shuffle loop) ----
    max_Delta = max(Delta_k);
    bin_edges = linspace(0, max_Delta, 25);
    bin_ctrs  = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    nBins     = numel(bin_ctrs);
    H_obs     = histcounts(t_pre, bin_edges, 'Normalization', 'probability') * 100;

    % ---- Single shuffle loop: random reassignment within each breath interval ----
    % For each Ca2+ event, keep it in the same breath interval [b_k, b_{k+1})
    % but draw its position uniformly within that interval.
    % t_shuf ~ Uniform(0, Delta_k)  for the interval containing that event.
    % This destroys any preferred latency while preserving which interval
    % each event belongs to (i.e. the breath-interval count structure).
    S_shuf = nan(M_shuf, 1);
    H_shuf = nan(M_shuf, nBins);

    for m = 1:M_shuf
        t_shuf = rand(nEv_in, 1) .* Delta_k(bin_in);   % Uniform(0, Delta_k) per event
        S_shuf(m)    = breath_scan_stat(t_shuf, w_sec);
        H_shuf(m, :) = histcounts(t_shuf, bin_edges, 'Normalization', 'probability') * 100;
    end

    % ---- p-value (Phipson & Smyth: conservative, avoids p=0) ----
    p_val = (1 + sum(S_shuf >= S_obs)) / (M_shuf + 1);

    fprintf('ROI#%02d : S_obs=%d  median(S_null)=%.1f  p=%.4f  (w=%.2fs, M=%d, n=%d)\n', ...
            keepIdx(i), S_obs, median(S_shuf), p_val, w_sec, M_shuf, nEv_in);

    % ---- Plot A: null distribution of scan statistic (% of total Ca2+ spikes) ----
    S_shuf_pct = S_shuf / nEv_in * 100;   % convert counts → % of total events
    S_obs_pct  = S_obs  / nEv_in * 100;

    % figure('Color','w', 'Name', sprintf('t_pre scan stat: ROI#%d', keepIdx(i)));
    % hold on;
    % 
    % histogram(S_shuf_pct, 'BinWidth', 100/nEv_in, 'Normalization','probability', ...
    %           'FaceColor',[0.65 0.65 0.65], 'EdgeColor','none', 'FaceAlpha',0.85);
    % xline(S_obs_pct, 'r-', 'LineWidth', 2.5);
    % yl = ylim;
    % text(S_obs_pct + 1, yl(2)*0.88, ...
    %      sprintf(' S_{obs} = %.1f%%\n (n=%d/%d)\n p = %.3f', ...
    %              S_obs_pct, S_obs, nEv_in, p_val), ...
    %      'Color','r', 'FontSize',11, 'VerticalAlignment','top');
    % 
    % xlabel('Proportion of calcium event within 100-ms lag');
    % ylabel('Proportion of shuffles');
    % title(sprintf('ROI#%d  —  t_{pre} scan statistic  (w = %.2f s,  M = %d)', ...
    %               keepIdx(i), w_sec, M_shuf));
    % xlim([0 100]);
    % axis square;  hold off;
    % 
    % outFig = fullfile(folderPath, sprintf('ca_breath_scan_stat_ROI%02d.png', keepIdx(i)));
    % exportgraphics(gcf, outFig, 'ContentType','vector');
    % fprintf('Saved scan statistic null distribution:\n  %s\n', outFig);
    % 
    % % ---- Plot B: t_pre histogram observed vs shuffled 95% envelope ----
    % H_mean = mean(H_shuf, 1);
    % H_lo   = prctile(H_shuf,  2.5, 1);
    % H_hi   = prctile(H_shuf, 97.5, 1);
    % 
    % % Use left bin edges so shuffle and observed share the same x=0 origin
    % be = bin_edges(1:end-1);
    % [xs_hi, ys_hi] = stairs(be, H_hi);
    % [xs_lo, ys_lo] = stairs(be, H_lo);
    % 
    % figure('Color','w', 'Name', sprintf('t_pre PDF vs null: ROI#%d', keepIdx(i)));
    % hold on;
    % 
    % fill([xs_hi; flipud(xs_lo)], [ys_hi; flipud(ys_lo)], ...
    %      [0.75 0.75 0.75], 'EdgeColor','none', 'FaceAlpha',0.6);
    % stairs(be, H_mean, 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);
    % stairs(be, H_obs,  'Color',[0.15 0.55 0.35], 'LineWidth',2.5);
    % 
    % xlabel('time lagging inspiration (sec)');
    % ylabel('% Ca^{2+} spike');
    % title(sprintf('ROI#%d  —  t_{pre} vs shuffle null  (p = %.3f,  n = %d events)', ...
    %               keepIdx(i), p_val, nEv_in));
    % legend('95% shuffle envelope', 'Shuffle mean', 'Observed', 'Location','best');
    % xlim([0 max_Delta]);
    % axis square;  hold off;
    % 
    % outFig = fullfile(folderPath, sprintf('ca_breath_tpre_vs_null_ROI%02d.png', keepIdx(i)));
    % exportgraphics(gcf, outFig, 'ContentType','vector');
    % fprintf('Saved t_pre PDF vs null envelope:\n  %s\n', outFig);

    % ---- Shared envelope variables for Plot C ----
    H_mean = mean(H_shuf, 1);
    H_lo   = prctile(H_shuf,  2.5, 1);
    H_hi   = prctile(H_shuf, 97.5, 1);

    be = bin_edges(1:end-1);
    [xs_hi, ys_hi] = stairs(be, H_hi);
    [xs_lo, ys_lo] = stairs(be, H_lo);

    % ---- Plot C: t_pre distribution (main) + scan statistic test (inset) ----
    sig_color = [1 0.1 0.1];   % red  — significant
    ns_color  = [0.1 0.1 0.1];   % black — not significant
    obs_color = double(p_val < 0.05) * sig_color + double(p_val >= 0.05) * ns_color;

    figure('Color','w', 'Name', sprintf('t_pre + shuffle: ROI#%d', keepIdx(i)), ...
           'Position', [100 100 1000 450]);

    % Left subplot: t_pre distribution vs shuffle envelope
    ax_main = subplot(1,2,1);
    hold(ax_main, 'on');

    fill([xs_hi; flipud(xs_lo)], [ys_hi; flipud(ys_lo)], ...
         [0.75 0.75 0.75], 'EdgeColor','none', 'FaceAlpha',0.6);
    stairs(be, H_mean, 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);
    stairs(be, H_obs,  'Color', obs_color, 'LineWidth',2.5);

    xlabel(ax_main, 'Time lagging inspiration (s)');
    ylabel(ax_main, 'Proportion of Ca^{2+} events (%)');
    title(ax_main, sprintf('ROI#%d  —  time to prior insipiration  (n=%d events)', keepIdx(i), nEv_in));
    legend(ax_main, '95% shuffle envelope','Shuffle mean','Observed','Location','northeast');
    xlim(ax_main, [0 max_Delta]);
    ylim([0 40]);
    axis(ax_main, 'square');
    hold(ax_main, 'off');

    % Right subplot: scan statistic null distribution
    ax_in = subplot(1,2,2);
    hold(ax_in, 'on');

    histogram(ax_in, S_shuf_pct, 'BinWidth', 100/nEv_in, 'Normalization','percentage', ...
              'FaceColor',[0.65 0.65 0.65], 'EdgeColor','none', 'FaceAlpha',0.85);
    xline(ax_in, S_obs_pct, 'Color', obs_color, 'LineWidth', 2);
    xlabel(ax_in, '% Ca^{2+} spike within 100-ms lag', 'FontSize', 11);
    ylabel(ax_in, '% Shuffles', 'FontSize', 11);
    ax_in.FontSize = 10;
    ax_in.XLim    = [0 100];
    ylim([0 30]);
    yl_in = ylim(ax_in);  xl_in = xlim(ax_in);
    text(ax_in, mean(xl_in), yl_in(2)+2, ...
         sprintf('%d/%d (p-val=%.3f)', n_pre_w, nEv_in, p_val), ...
         'Color', obs_color, 'FontSize',11, 'VerticalAlignment','top', 'HorizontalAlignment','center');
    axis(ax_in, 'square');
    box(ax_in, 'on');
    hold(ax_in, 'off');

    outFig = fullfile(folderPath, sprintf('ca_breath_combined_ROI%02d.png', keepIdx(i)));
    exportgraphics(gcf, outFig, 'ContentType','vector');
    fprintf('Saved combined t_pre + inset figure:\n  %s\n', outFig);

    % ---- Figure: t_post permutation test (independent of t_pre) ----
    % t_post = time from Ca event to NEXT inspiration onset.
    % Question: do Ca events cluster within w_sec BEFORE the next breath?
    % Test stat: n_post_w = sum(t_post <= w_sec)  (fixed window at t_post=0).
    % This is independent of t_pre scan stat, which asks "any cluster after current breath?"
    % Null: shuffle Ca event position uniformly within its interval,
    %       compute t_post_null = Delta_k - shuffled_position.

    max_Delta_post = max(Delta_k);
    bin_edges_post = linspace(0, max_Delta_post, 25);
    nBins_post     = numel(bin_edges_post) - 1;
    be_post        = bin_edges_post(1:end-1);

    % Observed: t_post histogram and fixed-window count (n_post_w computed above)
    H_obs_post    = histcounts(t_post, bin_edges_post, 'Normalization', 'probability') * 100;
    n_post_pct    = n_post_w / nEv_in * 100;

    % Shuffle loop: generate t_post null and compute fixed-window count + histogram
    N_shuf_post_w = nan(M_shuf, 1);
    H_shuf_post   = nan(M_shuf, nBins_post);

    for m = 1:M_shuf
        pos_shuf          = rand(nEv_in, 1) .* Delta_k(bin_in);   % random position in interval
        t_shuf_post       = Delta_k(bin_in) - pos_shuf;            % t_post from that position
        N_shuf_post_w(m)  = sum(t_shuf_post <= w_sec);
        H_shuf_post(m, :) = histcounts(t_shuf_post, bin_edges_post, 'Normalization', 'probability') * 100;
    end

    p_val_post = (1 + sum(N_shuf_post_w >= n_post_w)) / (M_shuf + 1);

    fprintf('ROI#%02d (t_post): n_obs=%d  median(n_null)=%.1f  p=%.4f  (w=%.2fs, M=%d, n=%d)\n', ...
            keepIdx(i), n_post_w, median(N_shuf_post_w), p_val_post, w_sec, M_shuf, nEv_in);

    N_shuf_post_pct = N_shuf_post_w / nEv_in * 100;

    H_mean_post = mean(H_shuf_post, 1);
    H_lo_post   = prctile(H_shuf_post,  2.5, 1);
    H_hi_post   = prctile(H_shuf_post, 97.5, 1);

    [xs_hi_post, ys_hi_post] = stairs(be_post, H_hi_post);
    [xs_lo_post, ys_lo_post] = stairs(be_post, H_lo_post);

    obs_color_post = double(p_val_post < 0.05) * [0.2 0.4 0.9] + ...
                     double(p_val_post >= 0.05) * [0.1 0.1 0.1];

    figure('Color','w', 'Name', sprintf('t_post + shuffle: ROI#%d', keepIdx(i)), ...
           'Position', [100 100 1000 450]);

    % Left subplot: t_post distribution vs shuffle envelope
    ax_main_post = subplot(1,2,1);
    hold(ax_main_post, 'on');

    fill([xs_hi_post; flipud(xs_lo_post)], [ys_hi_post; flipud(ys_lo_post)], ...
         [0.75 0.75 0.75], 'EdgeColor','none', 'FaceAlpha',0.6);
    stairs(be_post, H_mean_post, 'Color',[0.4 0.4 0.4], 'LineWidth',1.5);
    stairs(be_post, H_obs_post,  'Color', obs_color_post, 'LineWidth',2.5);

    xlabel(ax_main_post, 'Time leading inspiration (s)');
    ylabel(ax_main_post, 'Proportion of Ca^{2+} events (%)');
    title(ax_main_post, sprintf('ROI#%d  —  time to next inspiration  (n=%d events)', keepIdx(i), nEv_in));
    legend(ax_main_post, '95% shuffle envelope','Shuffle mean','Observed','Location','northeast');
    xlim(ax_main_post, [0 max_Delta_post]);
    ylim([0 40]);
    axis(ax_main_post, 'square');
    hold(ax_main_post, 'off');

    % Right subplot: fixed-window null (coupling to next inspiration onset)
    ax_in_post = subplot(1,2,2);
    hold(ax_in_post, 'on');

    histogram(ax_in_post, N_shuf_post_pct, 'BinWidth', 100/nEv_in, 'Normalization','percentage', ...
              'FaceColor',[0.65 0.65 0.65], 'EdgeColor','none', 'FaceAlpha',0.85);
    xline(ax_in_post, n_post_pct, 'Color', obs_color_post, 'LineWidth', 2);
    xlabel(ax_in_post, '% Ca^{2+} spike within 100-ms of next breath', 'FontSize', 11);
    ylabel(ax_in_post, '% Shuffles', 'FontSize', 11);
    ax_in_post.FontSize = 10;
    ax_in_post.XLim    = [0 100];
    ylim([0 30]);
    yl_in_post = ylim(ax_in_post);  xl_in_post = xlim(ax_in_post);
    text(ax_in_post, mean(xl_in_post), yl_in_post(2)+2, ...
         sprintf('%d/%d (p-val=%.3f)', n_post_w, nEv_in, p_val_post), ...
         'Color', obs_color_post, 'FontSize',11, 'VerticalAlignment','top', 'HorizontalAlignment','center');
    axis(ax_in_post, 'square');
    box(ax_in_post, 'on');
    hold(ax_in_post, 'off');

    outFig = fullfile(folderPath, sprintf('ca_breath_tpost_combined_ROI%02d.png', keepIdx(i)));
    exportgraphics(gcf, outFig, 'ContentType','vector');
    fprintf('Saved t_post combined figure:\n  %s\n', outFig);

    % ---- Accumulate per-ROI results ----
    roi_stats(i).roi_id         = keepIdx(i);
    roi_stats(i).ca_in          = ca_in;
    roi_stats(i).bin_in         = bin_in;
    roi_stats(i).nEv_in         = nEv_in;
    roi_stats(i).b_frames       = b_frames;
    roi_stats(i).Delta_k        = Delta_k;
    roi_stats(i).t_pre          = t_pre;
    roi_stats(i).t_post         = t_post;
    roi_stats(i).n_pre_w        = n_pre_w;
    roi_stats(i).n_post_w       = n_post_w;
    roi_stats(i).S_obs          = S_obs;
    roi_stats(i).S_shuf         = S_shuf;
    roi_stats(i).p_val          = p_val;
    roi_stats(i).H_obs          = H_obs;
    roi_stats(i).H_shuf         = H_shuf;
    roi_stats(i).bin_edges      = bin_edges;
    roi_stats(i).max_Delta      = max_Delta;
    roi_stats(i).N_shuf_post_w  = N_shuf_post_w;
    roi_stats(i).p_val_post     = p_val_post;
    roi_stats(i).H_obs_post     = H_obs_post;
    roi_stats(i).H_shuf_post    = H_shuf_post;
    roi_stats(i).bin_edges_post = bin_edges_post;
    roi_stats(i).max_Delta_post = max_Delta_post;

end

%% ---- Save master file ----
% Load this file + raw signals to reproduce all figures without re-running.
%
% Contents:
%   PARAMETERS      — all user-set analysis parameters
%   SIGNALS         — raw/processed traces needed for plotting
%   PEAK DETECTION  — breath peak audit trail (auto / added / deleted / final)
%   PER-ROI CACHE   — event-triggered segments and Ca peak locations
%   PERMUTATION     — per-ROI shuffle test inputs, outputs, and histograms

outMat_master = fullfile(folderPath, 'breath_master.mat');
save(outMat_master, ...
    ... % ---- folder ----
    'folderPath', ...
    ... % ---- acquisition parameters ----
    'fps_img', 'fps_breath', 'fs', ...
    ... % ---- breath peak detection parameters ----
    'findpeak_params', ...
    'insp_onsets_auto', 'insp_onsets_added', 'insp_onsets_deleted', 'insp_onsets_final', ...
    ... % ---- Ca detection parameters ----
    'keepIdx', 'nKeep', 'minProm_ca', 'minDist_ca', 'win_ca', 't_ca_evt', ...
    ... % ---- permutation test parameters ----
    'w_sec', 'M_shuf', ...
    ... % ---- event-triggered window ----
    'win', 't_evt', 'win5', 't_rast', ...
    ... % ---- signals ----
    'breath', 't_img', 'dFF', 'Sb', 't_insp', ...
    ... % ---- per-ROI caches (cell arrays, one entry per keepIdx) ----
    'ca_onsets_all', 'dff_evt_seg_all', 'breath_ca_seg_all', 'dt_latency_all', ...
    ... % ---- per-ROI permutation results (struct array, one entry per keepIdx) ----
    'roi_stats');
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
