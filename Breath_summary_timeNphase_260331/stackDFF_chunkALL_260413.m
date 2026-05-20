% breath_stacked_dFF_chunks_260413.m
%   Standalone: stacked dFF chunk plot for sig ROIs across all sessions,
%   pooled across multiple inputPaths.
%
%   Splits each sig ROI's full dFF into chunk_sec windows; drops chunks
%   with zero spikes; drops chunks listed in `exclude_chunks` (FOV idx is
%   global — increments across inputPaths in listed order).
%
%   Reads (already produced by breath_combined_summary_260329.m):
%     <each inputPath>/breath_combined_summary/fov_map.mat
%     per-session: *cpSAM_output.mat, ca_spike_data.mat
%
%   Writes:
%     <inputPaths(1)>/breath_combined_summary/stacked_dFF_chunks.pdf

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir  = fileparts(mfilename('fullpath'));
repoRoot   = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));

%% ============ PARAMETERS ============================================
inputPaths = [ ...
    "D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing", ...
    "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing" ...
];

fps_img         = 30;            % default imaging fps if detect fails
BaselineWinSec  = 20;            % dFF sliding-median baseline window (s)
chunk_sec       = 32;            % stacked-dFF chunk duration (seconds)
chunk_min_n     = 3;             % require >= chunk_min_n full chunks per FOV
align_xcorr     = true;          % iterative xcorr alignment to rolling template
max_shift_sec   = 5;             % max allowable shift per chunk (seconds)
n_align_passes  = 3;             % template->align->retemplate iterations
dFF_scale       = 0.3;           % scale-bar height in ΔF/F (label + bar length)
spacingFrac     = 0.15;          % row spacing fraction (lower => more overlap/invasion)
shear_from_vert_deg = 30;        % parallelogram lean: 30° = each row shifts right by (row_height * tan(30°)) on the rendered figure
ax_w_cm         = 36;            % axes width  (cm) — fixed so pixel/cm geometry is deterministic
ax_h_cm_per_row = 0.35;          % axes height (cm) per stacked row; min 10cm

% Each row = [global_fov_idx, roi_id, chunk_idx] to exclude.
% global_fov_idx counts across inputPaths in the order listed above.
exclude_chunks = [
    4, 29, 9;
    4, 29, 4;
];

doSave          = true;
% =====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 7);

%% Common chunk grid (target fps)
L_target       = round(chunk_sec * fps_img);
t_chunk_common = (0:L_target-1)' / fps_img;

stack_dff   = zeros(L_target, 0);
stack_evt   = zeros(L_target, 0);
stack_label = {};
sig_roi_count = 0;
fov_global    = 0;

%% Loop over inputPaths -> FOVs -> sig ROIs -> chunks
for ip = 1:numel(inputPaths)
    inputPath = inputPaths(ip);
    fprintf('\n============ inputPath %d/%d: %s ============\n', ...
        ip, numel(inputPaths), inputPath);

    fovMapPath = fullfile(inputPath, 'breath_combined_summary', 'fov_map.mat');
    if exist(fovMapPath, 'file') ~= 2
        warning('fov_map.mat not found at %s — skip.', fovMapPath); continue;
    end
    S_fm    = load(fovMapPath);
    fov_map = S_fm.fov_map;
    nFov    = size(fov_map, 1);

    allMat = dir(fullfile(inputPath, '**', 'ca_spike_data.mat'));
    sessNamesAll = cell(numel(allMat), 1);
    for i = 1:numel(allMat)
        [~, sessNamesAll{i}] = fileparts(allMat(i).folder);
    end

    for ff_local = 1:nFov
        fov_global  = fov_global + 1;
        sessName_fm = fov_map{ff_local, 1};
        sig_rois    = fov_map{ff_local, 2}(:)';
        if isempty(sig_rois)
            fprintf('F#%02d %s: no sig ROIs, skip.\n', fov_global, sessName_fm);
            continue;
        end

        match = find(strcmp(sessNamesAll, sessName_fm), 1);
        if isempty(match)
            warning('F#%02d %s: no folder match, skip.', fov_global, sessName_fm); continue;
        end
        folderPath = allMat(match).folder;

        fs = detect_session_fps(folderPath, fps_img);

        sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
        if isempty(sam_hits)
            warning('F#%02d %s: no SAM output.', fov_global, sessName_fm); continue;
        end
        SAMload = load(fullfile(folderPath, sam_hits(1).name));
        F_raw   = SAMload.F;
        dFFout  = helper.dFF_RZ(F_raw, 'FPS', fs, 'BaselineWinSec', BaselineWinSec);
        dFF     = double(dFFout.dFF);
        [T, ~]  = size(dFF);

        caLoad = load(fullfile(folderPath, 'ca_spike_data.mat'));
        caData = caLoad.roi_spikes;

        L_r  = round(chunk_sec * fs);
        N_ch = floor(T / L_r);
        fprintf('F#%02d %s: T=%d, fs=%.2f, L=%d, N_ch=%d, sig ROIs=%d\n', ...
            fov_global, sessName_fm, T, fs, L_r, N_ch, numel(sig_rois));
        if N_ch < chunk_min_n
            fprintf('  N_ch < %d, skip FOV.\n', chunk_min_n);
            continue;
        end

        for roi_id = sig_rois
            sig_roi_count = sig_roi_count + 1;
            dff_r   = dFF(:, roi_id);
            spk_idx = caData(roi_id).spike_idx(:);
            spk_idx = spk_idx(spk_idx >= 1 & spk_idx <= T);

            for c = 1:N_ch
                is_excluded = any(exclude_chunks(:,1) == fov_global & ...
                                  exclude_chunks(:,2) == roi_id & ...
                                  exclude_chunks(:,3) == c);
                if is_excluded
                    fprintf('  excluded F#%02dR#%02d c%02d\n', fov_global, roi_id, c);
                    continue;
                end

                i0 = (c-1)*L_r + 1;
                i1 = c*L_r;
                spk_in_chunk = spk_idx(spk_idx >= i0 & spk_idx <= i1);
                if isempty(spk_in_chunk), continue; end

                chunk = dff_r(i0:i1);
                if abs(fs - fps_img) < 0.1 && numel(chunk) == L_target
                    chunk_rs = chunk;
                    spk_rel  = spk_in_chunk - i0 + 1;
                else
                    tt_native = (0:numel(chunk)-1)' / fs;
                    chunk_rs  = interp1(tt_native, chunk, t_chunk_common, 'linear', 'extrap');
                    spk_rel   = round((double(spk_in_chunk) - i0) * fps_img / fs) + 1;
                end
                spk_rel = spk_rel(spk_rel >= 1 & spk_rel <= L_target);

                evt_col = zeros(L_target, 1);
                evt_col(spk_rel) = 1;

                stack_dff(:, end+1)   = chunk_rs;     %#ok<SAGROW>
                stack_evt(:, end+1)   = evt_col;      %#ok<SAGROW>
                stack_label{end+1, 1} = sprintf('F#%02dR#%02d c%02d', fov_global, roi_id, c); %#ok<SAGROW>
            end
        end
    end
end

n_kept = size(stack_dff, 2);
fprintf('\nStacked dFF chunks: %d kept (%g s windows; %d sig ROIs scanned).\n', ...
    n_kept, chunk_sec, sig_roi_count);
if n_kept == 0
    fprintf('No chunks left — nothing to plot.\n');
    return;
end

%% Iterative xcorr alignment (template-matching to common breath phase)
if align_xcorr && n_kept >= 2
    max_shift = min(round(max_shift_sec * fps_img), floor(L_target/2));

    % zero-mean each column for xcorr input
    Z = stack_dff - mean(stack_dff, 1, 'omitnan');
    Z(isnan(Z)) = 0;

    template = mean(Z, 2);

    for pass = 1:n_align_passes
        shifts = zeros(1, n_kept);
        for k = 1:n_kept
            [c, lags] = xcorr(Z(:,k), template, max_shift, 'coeff');
            [~, mi]   = max(c);
            shifts(k) = lags(mi);
        end

        % Circular shift (wraps cycle continuation from the other end —
        % appropriate for ~periodic breathing signals; no edge artifacts)
        for k = 1:n_kept
            s = shifts(k);
            Z(:, k)         = circshift(Z(:, k),         -s);
            stack_dff(:, k) = circshift(stack_dff(:, k), -s);
        end
        template  = mean(Z, 2);

        fprintf('Align pass %d/%d: |shift| mean=%.2f, max=%d frames (%.2fs)\n', ...
            pass, n_align_passes, mean(abs(shifts)), max(abs(shifts)), ...
            max(abs(shifts))/fps_img);
    end
end

%% Custom stacked plot with parallelogram shear (each row shifts right)
amp_g    = max(stack_dff(:)) - min(stack_dff(:));
if amp_g == 0, amp_g = 1; end
gap_g    = spacingFrac * amp_g;
offsets  = gap_g * (0:n_kept-1);
y_lo     = -gap_g;
y_hi     = offsets(end) + amp_g + gap_g;
y_range  = y_hi - y_lo;

% Fixed-cm axes for deterministic pixel geometry
ax_h_cm  = max(10, ax_h_cm_per_row * n_kept);
fig_w_cm = ax_w_cm + 5;     % +margins for labels/scale bar
fig_h_cm = ax_h_cm + 3;

% gap row height in cm on the rendered figure
gap_cm   = (gap_g / y_range) * ax_h_cm;
% per-row x shift in cm: tan(angle from vertical) * row height
xshift_per_row_cm = gap_cm * tand(shear_from_vert_deg);

% Convert cm -> seconds. xlim spans (chunk_sec + total_xshift_sec); plot
% pixels/sec = ax_w_cm / (chunk_sec + total_xshift_sec). Solving for it:
%   xshift_per_row_sec * ax_w_cm
%       = xshift_per_row_cm * (chunk_sec + (n_kept-1)*xshift_per_row_sec)
denom = ax_w_cm - xshift_per_row_cm * max(n_kept-1, 0);
if denom <= 0
    warning('Shear too steep for %g cm axes width; clamping.', ax_w_cm);
    xshift_per_row_cm = 0.8 * ax_w_cm / max(n_kept-1, 1);
    denom = ax_w_cm - xshift_per_row_cm * max(n_kept-1, 0);
end
xshift_per_row_sec = xshift_per_row_cm * chunk_sec / denom;
total_xshift_sec   = xshift_per_row_sec * max(n_kept-1, 0);

fprintf('Shear: %.2f cm/row -> %.3f s/row (total %.2f s, %.2f cm)\n', ...
    xshift_per_row_cm, xshift_per_row_sec, total_xshift_sec, ...
    xshift_per_row_cm * (n_kept-1));

%% Build figure
fig_stk = figure('Name','Stacked dFF chunks (parallelogram)','Color','w', ...
                 'Units','centimeters', 'Position',[2 1 fig_w_cm fig_h_cm]);
ax_stk  = axes(fig_stk, 'Units','centimeters', 'Position',[2.5 2 ax_w_cm ax_h_cm]);
hold(ax_stk, 'on');

t_base = (0:L_target-1)' / fps_img;
for k = 1:n_kept
    t_k = t_base + (k-1) * xshift_per_row_sec;
    y_k = stack_dff(:, k) + offsets(k);
    plot(ax_stk, t_k, y_k, 'k', 'LineWidth', 0.5);
end

% Scale bar (just past the rightmost shifted trace)
yTop = offsets(end) + 0.5 * amp_g;
y0_b = yTop - dFF_scale/2;
x_bar = chunk_sec + total_xshift_sec;
plot(ax_stk, [x_bar x_bar], [y0_b, y0_b + dFF_scale], 'b', 'LineWidth', 2);
text(ax_stk, x_bar, y0_b + dFF_scale + 0.05*dFF_scale, ...
     sprintf('%.2f \\DeltaF/F', dFF_scale), ...
     'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 9);

% Y-tick labels at bottom row offsets (label = ROI/chunk identity, tied to base x=0)
ytk    = 1:n_kept;
keep_t = (mod(ytk,5)==1) | (ytk==n_kept);
set(ax_stk, 'YTick', offsets(keep_t), 'YTickLabel', stack_label(keep_t));

xlim(ax_stk, [0, x_bar + 1]);   % small pad past scale bar
ylim(ax_stk, [y_lo, y_hi]);
xlabel(ax_stk, 'Time (s)');
ylabel(ax_stk, 'sig ROI \times chunk');
title(ax_stk, sprintf('Stacked dFF — %d chunks (%g s), shear %g° from vertical, %d FOVs (%d inputPaths)', ...
    n_kept, chunk_sec, shear_from_vert_deg, fov_global, numel(inputPaths)));
box(ax_stk, 'on');
hold(ax_stk, 'off');

%% Save PDF (in first inputPath's combined-summary folder)
if doSave
    outDir = fullfile(inputPaths(1), 'breath_combined_summary');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    pdfPath = fullfile(outDir, 'stacked_dFF_chunks.pdf');
    exportgraphics(fig_stk, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'none');
    fprintf('Saved: %s\n', pdfPath);
end
