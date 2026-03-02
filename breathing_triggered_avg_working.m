%%
close all; clear; clc;

folderPath = "D:\RUNQI\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi4_z0_512x512_5x_6000f_00001";

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

%% load breath peak
breath_peak = load(fullfile(folderPath, 'breath_peak_data.mat'));

insp_onsets = breath_peak.insp_onset_idx;
Sb = breath_peak.insp_onsets_train;

insp_onsets(insp_onsets<30)=[];
Sb(1:30)=[];

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
    fprintf('  %s → %s\n', fp.dot_selection{d}, fp.coord_types{di});
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

% %% ---- Breathing event detection: fine-tune parameters ----
%
% % =========================================================
% % USER-EDITABLE PEAK DETECTION PARAMETERS
% % =========================================================
% minProm_factor = 1.5;      % Prominence = factor × std(breath).
%                             %   Increase → fewer peaks.  Decrease → more peaks.
% minDist_sec    = 0.3;      % Min seconds between consecutive peaks.
% minHeight      = -Inf;     % Min absolute peak height (z-score). -Inf = disabled.
% minWidth_sec   = 0;        % Min peak width in seconds. 0 = disabled.
% snap_win_sec   = 0.2;     % Manual-add snaps to nearest local max within this window (s).
% preFrames      = 0;        % Shift all detected peaks earlier by this many frames.
% % =========================================================
%
% fs    = fps_img;
% T_img = numel(breath);
%
% peak     = breath;
% minProm  = minProm_factor * std(peak);
% minDist  = round(minDist_sec * fs);
% minWidth = round(minWidth_sec * fs);
% snap_win = round(snap_win_sec * fs);
%
% % Build findpeaks argument list dynamically
% fp_args = {'MinPeakProminence', minProm, 'MinPeakDistance', minDist};
% if isfinite(minHeight)
%     fp_args = [fp_args, {'MinPeakHeight', minHeight}];
% end
% if minWidth > 0
%     fp_args = [fp_args, {'MinPeakWidth', minWidth}];
% end
%
% [insp_amp, insp_onsets] = findpeaks(peak, fp_args{:});
%
% % Frame shift
% insp_onsets = insp_onsets - preFrames;
% insp_onsets(insp_onsets < 1) = 1;
%
% fprintf('Auto-detected %d inspiration-onset events.\n', numel(insp_onsets));
%
% % Freeze auto result (used for saving and audit trail)
% insp_onsets_auto = sort(insp_onsets(:));
%
% %% ---- Initial auto-detection QC plot ----
%
% figure('Name', 'Breathing event detection (auto)', 'Position', [200 200 1000 400]);
% plot(t_img, breath, 'k'); hold on;
% plot(t_img(insp_onsets_auto), breath(insp_onsets_auto), 'ro', 'MarkerFaceColor', 'r');
% xline(t_img(insp_onsets_auto), 'r:');
% xlabel('Time (s)'); ylabel('Breath (z)');
% title(sprintf('Auto-detected peaks: n=%d  |  minProm=%.2f  minDist=%.2f s', ...
%     numel(insp_onsets_auto), minProm, minDist_sec));
% grid on; hold off;
%
% %% ---- Interactive peak editor ----
% % Left-click  = add a peak (snaps to nearest local max within snap_win_sec)
% % Right-click = delete nearest existing peak
% % Press Enter = finish editing
%
% insp_onsets_added   = [];            % frame indices manually added
% insp_onsets_deleted = [];            % frame indices removed (from auto set)
% insp_onsets_edit    = insp_onsets_auto;   % working copy
%
% fig_edit = figure('Name', 'Peak Editor | Left-click=Add  Right-click=Delete  Enter=Done', ...
%     'Position', [100 80 1100 480]);
% ax_edit = axes('Parent', fig_edit);
%
% fprintf('\n=== INTERACTIVE PEAK EDITOR ===\n');
% fprintf('  Left-click  : add peak (snaps to local max within %.2f s)\n', snap_win_sec);
% fprintf('  Right-click : delete nearest existing peak\n');
% fprintf('  Press Enter : finish editing\n\n');
%
% while ishandle(fig_edit)
%     % Refresh plot
%     cla(ax_edit);
%     plot(ax_edit, t_img, breath, 'k', 'LineWidth', 1);
%     hold(ax_edit, 'on');
%     if ~isempty(insp_onsets_edit)
%         plot(ax_edit, t_img(insp_onsets_edit), breath(insp_onsets_edit), ...
%             'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
%     end
%     xlabel(ax_edit, 'Time (s)'); ylabel(ax_edit, 'Breath (z)');
%     title(ax_edit, sprintf('Peak Editor  |  n=%d  (+%d / -%d from auto)  |  L=Add  R=Delete  Enter=Done', ...
%         numel(insp_onsets_edit), numel(insp_onsets_added), numel(insp_onsets_deleted)));
%     grid(ax_edit, 'on'); hold(ax_edit, 'off');
%     drawnow;
%
%     % Wait for user click
%     try
%         [x_click, ~, button] = ginput(1);
%     catch
%         break;
%     end
%     if isempty(x_click) || ~ishandle(fig_edit)
%         break;   % Enter pressed → done
%     end
%
%     % Convert clicked x (time in s) → nearest frame index
%     [~, fc] = min(abs(t_img - x_click));
%
%     if button == 1   % Left-click: add peak
%         lo = max(1, fc - snap_win);
%         hi = min(numel(breath), fc + snap_win);
%         [~, rel_idx] = max(breath(lo:hi));
%         new_frame = lo + rel_idx - 1;
%         if ~ismember(new_frame, insp_onsets_edit)
%             insp_onsets_edit = sort([insp_onsets_edit(:); new_frame]);
%             insp_onsets_added(end+1) = new_frame;
%             insp_onsets_deleted(insp_onsets_deleted == new_frame) = [];
%             fprintf('Added   frame %5d  (t = %.2f s)\n', new_frame, t_img(new_frame));
%         else
%             fprintf('Peak at frame %d already exists — skipped.\n', new_frame);
%         end
%
%     elseif button == 3   % Right-click: delete nearest
%         if ~isempty(insp_onsets_edit)
%             [~, ni] = min(abs(t_img(insp_onsets_edit) - x_click));
%             del_frame = insp_onsets_edit(ni);
%             insp_onsets_edit(ni) = [];
%             if ismember(del_frame, insp_onsets_auto)
%                 insp_onsets_deleted(end+1) = del_frame;
%             end
%             insp_onsets_added(insp_onsets_added == del_frame) = [];
%             fprintf('Deleted frame %5d  (t = %.2f s)\n', del_frame, t_img(del_frame));
%         end
%     end
% end
%
% % Finalise
% insp_onsets_final   = sort(insp_onsets_edit(:));
% insp_onsets_added   = sort(insp_onsets_added(:));
% insp_onsets_deleted = sort(insp_onsets_deleted(:));
%
% fprintf('\nEditing complete: %d final peaks  (+%d added, -%d deleted from %d auto)\n', ...
%     numel(insp_onsets_final), numel(insp_onsets_added), ...
%     numel(insp_onsets_deleted), numel(insp_onsets_auto));
%
% % Overwrite insp_onsets so all downstream code is unaffected
% insp_onsets = insp_onsets_final;
% t_insp      = t_img(insp_onsets);
%
% % Recompute binary event train from final peaks
% Sb = zeros(T_img, 1);
% Sb(insp_onsets) = 1;
%
% %% ---- Save: final peak detection figure + .mat ----
%
% fig_qc_final = figure('Name', 'Breathing peak detection (final)', ...
%     'Position', [200 200 1000 400]);
% plot(t_img, breath, 'k', 'LineWidth', 1); hold on;
% if ~isempty(insp_onsets_auto)
%     plot(t_img(insp_onsets_auto), breath(insp_onsets_auto), ...
%         'o', 'Color', [0.65 0.65 0.65], 'MarkerFaceColor', [0.65 0.65 0.65], ...
%         'MarkerSize', 4, 'DisplayName', 'Auto peaks');
% end
% if ~isempty(insp_onsets_added)
%     plot(t_img(insp_onsets_added), breath(insp_onsets_added), ...
%         'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 7, 'DisplayName', 'Manually added');
% end
% if ~isempty(insp_onsets_final)
%     plot(t_img(insp_onsets_final), breath(insp_onsets_final), ...
%         'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4, 'DisplayName', 'Final accepted');
%     xline(t_img(insp_onsets_final), 'r:');
% end
% legend('Location', 'best');
% xlabel('Time (s)'); ylabel('Breath (z)');
% title(sprintf('Final: %d peaks  (auto=%d, +%d added, -%d deleted)', ...
%     numel(insp_onsets_final), numel(insp_onsets_auto), ...
%     numel(insp_onsets_added), numel(insp_onsets_deleted)));
% grid on; hold off;
%
% outFig_peaks = fullfile(folderPath, 'breath_peak_detection_final.png');
% exportgraphics(fig_qc_final, outFig_peaks, 'ContentType', 'vector');
% fprintf('Saved final peak detection figure:\n  %s\n', outFig_peaks);
%
% % findpeaks parameter struct
% findpeak_params = struct( ...
%     'minProm_factor',  minProm_factor, ...
%     'minProm',         minProm, ...
%     'minDist_sec',     minDist_sec, ...
%     'minDist_frames',  minDist, ...
%     'minHeight',       minHeight, ...
%     'minWidth_sec',    minWidth_sec, ...
%     'minWidth_frames', minWidth, ...
%     'snap_win_sec',    snap_win_sec, ...
%     'preFrames',       preFrames, ...
%     'fps',             fs ...
% );
%
% outMat_peaks = fullfile(folderPath, 'breath_peak_data.mat');
% save(outMat_peaks, ...
%     'findpeak_params', ...
%     'insp_onsets_auto', ...
%     'insp_onsets_added', ...
%     'insp_onsets_deleted', ...
%     'insp_onsets_final');
% fprintf('Saved peak detection data:\n  %s\n', outMat_peaks);

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
outFig = fullfile(folderPath, 'breathing_spectrum_raw.png');
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


%% ---- Figure 3. Inspiration-triggered dFF ----

fs   = fps_img;
win  = round(0.5 * fs);          % ±0.5 s
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

outFig = fullfile(folderPath, 'triggered_avg.png');
exportgraphics(gcf, outFig, 'ContentType','vector');
fprintf('Saved triggered average:\n  %s\n', outFig);

%% ---- Figure 4. Calcium-triggered average of breathing (selected ROIs) ----

% Select ROIs to run peak detection on (only well-isolated cells)
keepIdx  = 6;          % edit to choose which ROIs
dFF_keep = dFF(:, keepIdx);
nKeep    = numel(keepIdx);

win_ca   = round(1.0 * fs);          % ±1 s window around calcium peak
t_ca_evt = (-win_ca:win_ca) / fs;

% Peak detection params — tweak per dataset
minProm_ca = 0.1;                     % prominence threshold in dF/F units
minDist_ca = round(fs * 0.2);         % min 1 s between peaks (refractory)

[~, ca_onsets] = findpeaks(dFF_keep, ...
    'MinPeakProminence', minProm_ca, ...
    'MinPeakDistance',   minDist_ca);

% local filter only for windowed segment extraction
ca_onsets_ev = ca_onsets(ca_onsets > win_ca & ca_onsets < T_img - win_ca);
nEv_ca = numel(ca_onsets_ev);

fprintf('ROI#%02d : %d calcium peaks detected.\n', keepIdx, nEv_ca);

seg_ca = zeros(nEv_ca, 2*win_ca+1);
for e = 1:nEv_ca
    c = ca_onsets_ev(e);
    seg_ca(e,:) = breath(c-win_ca : c+win_ca);
end

% ---------- QC: raw dF/F with detected calcium peaks ----------
figure('Name', 'Calcium peak detection QC', 'Color', 'w');
tiledlayout('flow');

nexttile;
plot(t_img, dFF_keep, 'k'); hold on;
if ~isempty(ca_onsets)
    plot(t_img(ca_onsets), dFF_keep(ca_onsets), 'ro', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 3);
end
title(sprintf('ROI#%d  (n=%d)', keepIdx, numel(ca_onsets)));
xlabel('Time (s)');
ylabel('\DeltaF/F');

% outFig = fullfile(folderPath, sprintf('ca_peak_detection_ROI%02d.png',keepIdx));
% exportgraphics(gcf, outFig, 'ContentType', 'vector');
% fprintf('Saved calcium peak detection:\n  %s\n', outFig);

%%
% ---------- Per-ROI calcium-triggered breathing overlays ----------

figure('Name', sprintf('Calcium-triggered breathing avg (n=%d ROIs)', nKeep));

breath_avg = mean(seg_ca, 1);
nEv_ca     = size(seg_ca, 1);

plot(t_ca_evt, seg_ca',   'Color', [0.8 0.8 0.8], 'LineWidth', 0.7); hold on;
plot(t_ca_evt, breath_avg, 'Color', 'k',           'LineWidth', 2);

xline(0, '--', 'Color', [0.4 0.4 0.4]);   % calcium peak at t = 0
yline(0, '--', 'Color', [0.6 0.6 0.6]);

xlabel('Time from Ca^{2+} peak (s)');
ylabel('Breathing (z-score)');
title(sprintf('ROI#%d  (n=%d)', keepIdx, nEv_ca));

axis square;

outFig = fullfile(folderPath, 'ca_triggered_breath_avg.png');
exportgraphics(gcf, outFig, 'ContentType', 'vector');
outFig = fullfile(folderPath, 'ca_triggered_breath_avg.fig');
savefig(gcf, outFig);
fprintf('Saved calcium-triggered breathing average:\n  %s\n', outFig);

%% ---- Figure 5: dt_closest analysis (raster + histogram + permutation) ----

win5   = round(2.5 * fs);
w_ms   = 100;
w_sec  = w_ms / 1000;
M_shuf = 3000;

% ---- Assign each Ca2+ event to a breath interval ----
ca_ev    = ca_onsets;
b_frames = sort(insp_onsets(:));
[~, ~, bin_idx] = histcounts(ca_ev, b_frames);
keep_ca = (bin_idx >= 1) & (bin_idx <= numel(b_frames) - 1);
ca_in   = ca_ev(keep_ca);
bin_in  = bin_idx(keep_ca);
nEv_in  = numel(ca_in);
Delta_k = diff(b_frames) / fs;
t_pre   = (ca_in - b_frames(bin_in)) / fs;
t_post  = Delta_k(bin_in) - t_pre;

% ---- Compute signed dt_closest ----
dt_closest = t_pre;
closer_to_next = t_post < t_pre;
dt_closest(closer_to_next) = -t_post(closer_to_next);

% ---- Build insp_rel only for ca_in events ----
insp_rel_in = cell(nEv_in, 1);
for e = 1:nEv_in
    c      = ca_in(e);
    nearby = insp_onsets(insp_onsets >= (c-win5) & insp_onsets <= (c+win5));
    insp_rel_in{e} = (nearby - c) / fs;
end

% ---- Histogram bins centred on 0 ----
half_max  = max(Delta_k) / 2;
bin_edges = linspace(-half_max, half_max, 25);
bin_ctrs  = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
nBins     = numel(bin_ctrs);
H_obs     = histcounts(dt_closest, bin_edges, 'Normalization', 'probability') * 100;

% ---- Observed scan statistic ----
S_obs = breath_scan_stat(dt_closest, w_sec);

% ---- Shuffle (interval-based, recompute dt_closest from shuffled position) ----
S_shuf = nan(M_shuf, 1);
H_shuf = nan(M_shuf, nBins);
for m = 1:M_shuf
    t_pre_sh  = rand(nEv_in, 1) .* Delta_k(bin_in);
    t_post_sh = Delta_k(bin_in) - t_pre_sh;
    dt_sh     = t_pre_sh;
    dt_sh(t_post_sh < t_pre_sh) = -t_post_sh(t_post_sh < t_pre_sh);
    S_shuf(m)    = breath_scan_stat(dt_sh, w_sec);
    H_shuf(m,:)  = histcounts(dt_sh, bin_edges, 'Normalization', 'probability') * 100;
end
p_val   = (1 + sum(S_shuf >= S_obs)) / (M_shuf + 1);
n_close = sum(abs(dt_closest) <= w_sec);
fprintf('ROI#%02d : S_obs=%d  median(S_null)=%.1f  p=%.4f  (w=%.2fs, M=%d, n=%d)\n', ...
    keepIdx, S_obs, median(S_shuf), p_val, w_sec, M_shuf, nEv_in);

% ---- Color by significance ----
sig_color = [1 0.1 0.1];
ns_color  = [0.1 0.1 0.1];
obs_color = double(p_val < 0.05) * sig_color + double(p_val >= 0.05) * ns_color;

% ---- Figure: 1x3 layout ----
fig5 = figure('Color', 'w', ...
    'Name', sprintf('dt_closest: ROI#%d', keepIdx), ...
    'Position', [0 0 1100 400]);

% ---- Subplot 1: Raster sorted by dt_closest ascending ----
ax_rast = subplot(1,3,1);
hold(ax_rast, 'on');

[~, sort_idx] = sort(dt_closest, 'ascend');

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
ylabel(ax_rast, 'Ca^{2+} event #  (sorted by dt_{closest})');
title(ax_rast, sprintf('ROI%02d  (%d events in %d breaths)', keepIdx, nEv_in, sum(Sb)));
hold(ax_rast, 'off');

% ---- Subplot 2: Histogram vs shuffle envelope ----
ax_hist = subplot(1,3,2);
hold(ax_hist, 'on');

% Shuffle envelope (2.5–97.5 percentile band)
H_mean = mean(H_shuf, 1);
H_lo   = prctile(H_shuf,  2.5, 1);
H_hi   = prctile(H_shuf, 97.5, 1);

be = bin_edges(1:end-1);
[xs_hi, ys_hi] = stairs(be, H_hi);
[xs_lo, ys_lo] = stairs(be, H_lo);

fill(ax_hist, [xs_hi; flipud(xs_lo)], [ys_hi; flipud(ys_lo)], ...
     [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
stairs(ax_hist, be, H_mean, 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
stairs(ax_hist, be, H_obs,  'Color', obs_color, 'LineWidth', 2.5);

xline(ax_hist, 0, 'r--', 'LineWidth', 1);

xlabel(ax_hist, 'dt_{closest}  (s)');
ylabel(ax_hist, 'Proportion of Ca^{2+} events (%)');
legend(ax_hist, {'95% shuffle envelope', 'Shuffle mean', 'Observed'}, 'Location', 'northeast');
xlim(ax_hist, [-half_max half_max]);
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

xlabel(ax_scan, sprintf('%% events within +/-%d ms of nearest inspiration', w_ms));
ylabel(ax_scan, '% Shuffles');

yl_sc = ylim(ax_scan);  xl_sc = xlim(ax_scan);
text(ax_scan, mean(xl_sc), yl_sc(2) + 1, ...
     sprintf('%d/%d  (p=%.3f)', n_close, nEv_in, p_val), ...
     'Color', obs_color, 'FontSize', 11, ...
     'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center');

xlim(ax_scan, [0 100]);
ylim(ax_scan, [0 30]);
axis(ax_scan, 'square');
box(ax_scan, 'on');
hold(ax_scan, 'off');

% ---- Save figure ----
outFig_fig = fullfile(folderPath, sprintf('calcium_breath_dtclosest_ROI%02d.fig', keepIdx));
savefig(fig5, outFig_fig);
outFig_png = fullfile(folderPath, sprintf('calcium_breath_dtclosest_ROI%02d.png', keepIdx));
exportgraphics(fig5, outFig_png, 'ContentType', 'vector');
fprintf('Saved dt_closest figure:\n  %s\n  %s\n', outFig_fig, outFig_png);

%% ---- Save master file ----
outMat_master = fullfile(folderPath, 'breath_master_working.mat');
save(outMat_master, ...
    'folderPath', ...
    'fps_img', 'fps_breath', 'fs', ...
    'breath_peak', ...
    'keepIdx', 'nKeep', 'minProm_ca', 'minDist_ca', 'win_ca', 't_ca_evt', ...
    'w_sec', 'M_shuf', 'w_ms', ...
    'win', 't_evt', 'win5', ...
    'breath', 't_img', 'dFF', 'Sb', ...
    'ca_onsets', 'ca_in', 'bin_in', 'nEv_in', ...
    'dt_closest', 't_pre', 't_post', 'Delta_k', 'b_frames', ...
    'S_obs', 'S_shuf', 'p_val', 'n_close', ...
    'H_obs', 'H_shuf', 'bin_edges', 'half_max', ...
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
