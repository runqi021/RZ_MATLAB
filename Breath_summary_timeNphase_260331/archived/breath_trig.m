%% breath_trig.m  --  Single-folder dt raster, histogram & phase analysis
%  For each spiking ROI, produces THREE figures:
%
%  Figure 1 — Ca-centered (zero = Ca spike, breath ticks as marks):
%    [Raster sorted by dt_last]   [Raster sorted by dt_next]
%    [Histogram dt_last]           [Histogram dt_next]
%
%  Figure 2 — Breath-centered (zero = breath peak, Ca spikes as marks):
%    [Raster sorted by dt_last_Ca] [Raster sorted by dt_next_Ca]
%    [Histogram dt_last_Ca]         [Histogram dt_next_Ca]
%
%  Figure 3 — Phase analysis:
%    [Polar histogram]  [Linear histogram + Rayleigh test]
%    phase = dt_last / cycle_duration → 0 = insp onset, 1 = next insp
%
%  Saves .fig + .pdf per ROI into folderPath.
%  ── USER PARAMETERS ─────────────────────────────────────────────────────
clear; close all; clc;
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

folderPath = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\left_pFN_roi1_z0_3x_00001";

fps_img        = 30;
nDrop          = 30;
PixelSizeBase  = 1.7778;   % um/px at 1x zoom (ScanImage default)
crop_um        = 25;        % ROI crop box size (um)
gamma_val  = 0.6;           % gamma correction
clip_lo    = 0.5;           % brightness clip percentile (low)
clip_hi    = 99.5;          % brightness clip percentile (high)
% ──────────────────────────────────────────────────────────────────────────

%% ---- Auto-find files ----
csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
bp_file  = dir(fullfile(folderPath, '*breath_peak_data.mat'));

if isempty(csv_hits) || isempty(sam_hits) || isempty(bp_file)
    error('Missing DLC CSV, SAM output, or breath_peak_data.mat in:\n  %s', folderPath);
end

dlc_csv = fullfile(folderPath, csv_hits(1).name);
SAM     = fullfile(folderPath, sam_hits(1).name);

fprintf('DLC CSV : %s\n', csv_hits(1).name);
fprintf('SAM     : %s\n', sam_hits(1).name);

%% ---- Load SAM -> dFF ----
SAMload = load(SAM);
F = SAMload.F;
F(1:nDrop, :) = [];
dFFout = helper.dFF_RZ(F);
dFF = dFFout.dFF;
maskL = SAMload.maskL;

[T, N_roi] = size(dFF);

%% ---- Read zoom factor from source TIFF → um_per_px ----
tif_hits = dir(fullfile(folderPath, '*.tif'));
um_per_px = PixelSizeBase;   % fallback: 1x zoom
for tt = 1:numel(tif_hits)
    try
        info_t = imfinfo(fullfile(folderPath, tif_hits(tt).name));
        meta_str = "";
        if isfield(info_t(1),'Software') && ~isempty(info_t(1).Software)
            meta_str = meta_str + string(info_t(1).Software) + newline;
        end
        if isfield(info_t(1),'ImageDescription') && ~isempty(info_t(1).ImageDescription)
            meta_str = meta_str + string(info_t(1).ImageDescription) + newline;
        end
        tok = regexp(meta_str, "SI\.hRoiManager\.scanZoomFactor\s*=\s*([-+]?\d+\.?\d*)", "tokens", "once");
        if ~isempty(tok)
            zf = str2double(tok{1});
            if isfinite(zf) && zf > 0
                um_per_px = PixelSizeBase / zf;
                fprintf('Zoom = %.2fx  →  um_per_px = %.4f\n', zf, um_per_px);
                break;
            end
        end
    catch
    end
end

%% ---- Load avg projection ----
avg_hits = dir(fullfile(folderPath, '*_AVG_for_CP.tif'));
if ~isempty(avg_hits)
    avgImg = double(imread(fullfile(folderPath, avg_hits(1).name)));
else
    avgImg = [];
    warning('No avg projection TIF found, will skip ROI thumbnail.');
end
fs = fps_img;
t_img = (0:T-1)' / fps_img;

fprintf('dFF loaded: [%d x %d] (T x N_roi)\n', T, N_roi);

%% ---- Load calcium spikes ----
in = load(fullfile(folderPath, 'ca_spike_data.mat'));
ifSpk      = in.ifSpike;
roiSpk_id  = find(ifSpk == 1);
ca_spk_data = in.roi_spikes;
nSpkROI    = numel(roiSpk_id);

ca_spk_id = cell(1, nSpkROI);
for k = 1:nSpkROI
    i = roiSpk_id(k);
    ca_spk_id{k} = ca_spk_data(i).spike_idx;
end

fprintf('Spiking ROIs: %d of %d\n', nSpkROI, N_roi);
assert(nSpkROI > 0, 'No spiking ROIs found in ca_spike_data.mat');

%% ---- Load breath peaks ----
breath_peak  = load(fullfile(folderPath, 'breath_peak_data.mat'));
insp_onsets  = breath_peak.insp_onset_idx;
Sb           = breath_peak.insp_onsets_train;

insp_onsets(insp_onsets < nDrop) = [];
insp_onsets = insp_onsets - nDrop;
Sb(1:nDrop) = [];

assert(numel(insp_onsets) >= 2, 'Fewer than 2 inspiration onsets after dropping frames.');

%% ---- Reconstruct breathing signal & Hilbert phase ----
dlc_raw  = readmatrix(dlc_csv, 'NumHeaderLines', 3);
data_dlc = dlc_raw;
data_dlc(1:nDrop, :) = [];

dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
dot_idx  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);
fp       = breath_peak.findpeak_params;
nF_dlc   = size(data_dlc, 1);
nActiveDots = numel(fp.dot_selection);
traces   = NaN(nF_dlc, nActiveDots);

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
end

switch fp.combine_method
    case 'sum',  breath = sum(traces,  2, 'omitnan');
    otherwise,   breath = mean(traces, 2, 'omitnan');
end
if fp.inverted, breath = -breath; end
breath = detrend(breath);
breath = (breath - mean(breath)) / std(breath);

% Hilbert instantaneous phase: 0 at peak, +/-pi at troughs
breath_analytic = hilbert(breath);
breath_phase    = angle(breath_analytic);   % [-pi, pi]

fprintf('Breathing signal: %d samples, Hilbert phase computed.\n', numel(breath));

%% ---- Breath intervals ----
b_frames = sort(insp_onsets(:));
Delta_k  = diff(b_frames) / fs;       % breath interval durations (s)

fprintf('Breath intervals: n=%d, mean=%.3f s, range=[%.3f, %.3f] s\n', ...
    numel(Delta_k), mean(Delta_k), min(Delta_k), max(Delta_k));

%% ---- Raster / histogram parameters ----
win_raster = mean(Delta_k) * 2.5;      % seconds
win5       = round(win_raster * fs);    % frames

% Histogram bins: 1 frame per bin, same window as raster
bin_dt     = 1 / fs;                    % 1 frame
bin_edges_dt = 0 : bin_dt : win_raster; % [0, win_raster]

fprintf('Raster window: +/- %.2f s  |  Histogram bin: %.4f s (1 frame)\n', ...
    win_raster, bin_dt);

%% ---- Breathing autocorrelogram (shared across ROIs) ----
max_lag_frames = round(win_raster * fs);
[breath_acf, breath_acf_lags] = xcorr(breath, max_lag_frames, 'coeff');
breath_acf_t = breath_acf_lags / fs;

%% ---- One-cycle mapping for overlay square ----
one_cycle  = mean(Delta_k);                    % seconds
half_cycle = one_cycle / 2;
time_to_phase = @(t) t * 2*pi / one_cycle;     % seconds → radians

% Trim autocorrelation to one cycle
acf_mask       = abs(breath_acf_t) <= half_cycle;
acf_1cyc_phase = time_to_phase(breath_acf_t(acf_mask));
acf_1cyc       = breath_acf(acf_mask);

%% ---- Loop over spiking ROIs ----
for k = 1:nSpkROI
    roi_id   = roiSpk_id(k);
    ca_ev    = ca_spk_id{k}(:);
    nCaTotal = numel(ca_ev);

    if nCaTotal == 0
        fprintf('ROI#%02d : no calcium events, skipping.\n', roi_id);
        continue;
    end

    % Assign each Ca event to its breath interval
    [~, ~, bin_idx] = histcounts(ca_ev, b_frames);
    keep_ca = (bin_idx >= 1) & (bin_idx <= numel(b_frames) - 1);
    ca_in   = ca_ev(keep_ca);
    bin_in  = bin_idx(keep_ca);
    nEv_in  = numel(ca_in);

    if nEv_in == 0
        fprintf('ROI#%02d : no Ca events within breath intervals, skipping.\n', roi_id);
        continue;
    end

    % dt_last / dt_next for each Ca event
    dt_last = (ca_in - b_frames(bin_in)) / fs;
    dt_next = (b_frames(bin_in + 1) - ca_in) / fs;
    dt_nearest = min(dt_last, dt_next);

    % Build insp_rel for Ca-centered raster (breath ticks relative to each Ca event)
    insp_rel_in = cell(nEv_in, 1);
    for e = 1:nEv_in
        c      = ca_in(e);
        nearby = insp_onsets(insp_onsets >= (c - win5) & insp_onsets <= (c + win5));
        insp_rel_in{e} = (nearby - c) / fs;
    end

    % Breath-centered data: Ca spikes relative to each breath peak
    nBreath = numel(b_frames);
    ca_rel_breath = cell(nBreath, 1);
    for b = 1:nBreath
        bp = b_frames(b);
        nearby_ca = ca_ev(ca_ev >= (bp - win5) & ca_ev <= (bp + win5));
        ca_rel_breath{b} = (nearby_ca - bp) / fs;
    end

    fprintf('ROI#%02d : %d Ca events in %d breath intervals\n', ...
        roi_id, nEv_in, numel(Delta_k));

    % ---- Pre-compute dF/F snippets & PETH (needed by overlay square) ----
    win_frames = round(win_raster * fs);
    dff_roi    = dFF(:, roi_id);
    snippets   = nan(nBreath, 2*win_frames + 1);
    for b = 1:nBreath
        bp = b_frames(b);
        i1 = bp - win_frames;
        i2 = bp + win_frames;
        if i1 >= 1 && i2 <= T
            snippets(b, :) = dff_roi(i1:i2);
        end
    end
    valid_snip = ~all(isnan(snippets), 2);
    nValid     = sum(valid_snip);
    t_snip     = (-win_frames:win_frames) / fs;
    mean_dff   = mean(snippets(valid_snip, :), 1, 'omitnan');

    % Filter raster & PETH to same valid breaths as dF/F
    valid_idx       = find(valid_snip);
    b_frames_valid  = b_frames(valid_idx);
    ca_rel_valid    = ca_rel_breath(valid_idx);

    % PETH computation
    all_ca_rel   = cell2mat(ca_rel_valid);
    peth_edges   = -win_raster : bin_dt : win_raster;
    peth_counts  = histcounts(all_ca_rel, peth_edges);
    peth_rate    = peth_counts / (nValid * bin_dt);
    peth_centers = peth_edges(1:end-1) + bin_dt/2;

    % ============================================================
    % Figure layout (cm, for Illustrator):
    %   Row 1: [Avg projection (square)] [Overlay summary (square)]
    %   Row 2: Breath-triggered dF/F (2:1)
    %   Row 3: Ca spike raster sorted by dt_nearest_ca (2:1)
    %   Row 4: PETH histogram (2:1)
    %   Row 5: Breathing autocorrelogram (2:1, half-height)
    % ============================================================
    % --- Figure dimensions (cm) — all shrunk by 1.5x ---
    sc = 1.5;                           % global scale factor
    fig_w = 18 / sc;                    % cm
    margin_l = 2.2 / sc;  margin_r = 0.8 / sc;
    margin_b = 1.5 / sc;  margin_t = 0.8 / sc;
    vgap     = 1.5 / sc;

    % Row 1: two squares
    sq = 5.5 / sc;  sq_gap = 3.0 / sc;
    sq_x1 = margin_l;
    sq_x2 = margin_l + sq + sq_gap;

    % Rows 2-3: 2:1 aspect
    plot_w = fig_w - margin_l - margin_r;
    plot_h = plot_w / 2;
    % Row 4 (PETH): half the height of rows 2-3
    plot_h4 = plot_h / 2;

    fig_h = margin_b + plot_h4 + vgap + plot_h4 + vgap + plot_h + vgap + plot_h + vgap + sq + margin_t;

    fig = figure('Color', 'w', ...
        'Name', sprintf('breath_trig ROI#%02d', roi_id), ...
        'Units', 'centimeters', 'Position', [1 1 fig_w fig_h], ...
        'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
        'PaperPosition', [0 0 fig_w fig_h]);
    set(fig, 'DefaultAxesFontSize', 7, 'DefaultTextFontSize', 8);

    % Compute y positions (bottom-up)
    row5_y = margin_b;
    row4_y = row5_y + plot_h4 + vgap;
    row3_y = row4_y + plot_h4 + vgap;
    row2_y = row3_y + plot_h + vgap;
    row1_y_cm = row2_y + plot_h + vgap;

    % --- (1a) Avg projection with ROI crop ---
    ax0 = axes(fig, 'Units','centimeters', ...
        'Position', [sq_x1 row1_y_cm sq sq]);
    if ~isempty(avgImg)
        [ry, rx] = find(maskL == roi_id);
        if ~isempty(ry)
            cy = round(mean(ry));  cx = round(mean(rx));
            crop_px = round(crop_um / um_per_px / 2);   % half-width in pixels
            [imgH, imgW] = size(avgImg);
            r1 = max(1, cy - crop_px);  r2 = min(imgH, cy + crop_px);
            c1 = max(1, cx - crop_px);  c2 = min(imgW, cx + crop_px);

            roi_patch = avgImg(r1:r2, c1:c2);
            lo = prctile(roi_patch(:), clip_lo);
            hi = prctile(roi_patch(:), clip_hi);
            roi_patch = (roi_patch - lo) / (hi - lo);
            roi_patch = max(0, min(1, roi_patch));
            roi_patch = roi_patch .^ gamma_val;

            % ROI outline in crop coordinates
            mask_crop = maskL(r1:r2, c1:c2);
            roi_perim = bwperim(mask_crop == roi_id);

            patchRGB = repmat(roi_patch, [1 1 3]);
            yellow = [1 1 0];
            for ch = 1:3
                plane = patchRGB(:,:,ch);
                plane(roi_perim) = yellow(ch);
                patchRGB(:,:,ch) = plane;
            end

            imshow(patchRGB, 'Parent', ax0);

            % Scale bar: 5 um
            sb_px = round(5 / um_per_px);
            hold(ax0, 'on');
            sb_y = size(roi_patch,1) - 3;
            plot(ax0, [3 3+sb_px], [sb_y sb_y], 'w-', 'LineWidth', 2);
            text(ax0, 3+sb_px/2, sb_y-3, '5 \mum', 'Color','w', ...
                'FontSize',8, 'HorizontalAlignment','center');
            hold(ax0, 'off');
        else
            text(ax0, 0.5, 0.5, 'ROI not in mask', 'HorizontalAlignment','center');
            axis(ax0, 'off');
        end
    else
        text(ax0, 0.5, 0.5, 'No avg projection', 'HorizontalAlignment','center');
        axis(ax0, 'off');
    end
    title(ax0, sprintf('ROI#%02d', roi_id));

    % --- (1b) Overlay summary: avg dF/F, PETH KDE, autocorr ---
    ax_ov = axes(fig, 'Units','centimeters', ...
        'Position', [sq_x2 row1_y_cm sq sq]);

    % Avg dF/F (green, left y-axis) — trim to one cycle
    dff_mask = abs(t_snip) <= half_cycle;
    dff_trim = mean_dff(dff_mask);
    dff_phase = time_to_phase(t_snip(dff_mask));

    yyaxis(ax_ov, 'left');
    plot(ax_ov, dff_phase, dff_trim, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
    ylabel(ax_ov, 'dF/F');
    ax_ov.YColor = [0.2 0.7 0.2];
    hold(ax_ov, 'on');

    % Autocorrelation (black, no axis) — normalize to dF/F range
    dff_range = max(dff_trim) - min(dff_trim);
    acf_scaled = min(dff_trim) + (acf_1cyc - min(acf_1cyc)) / (max(acf_1cyc) - min(acf_1cyc) + eps) * dff_range;
    plot(ax_ov, acf_1cyc_phase, acf_scaled, 'k-', 'LineWidth', 1);

    % PETH KDE (red, right y-axis) — KDE of Ca spike times within one cycle
    all_ca_rel_trim = all_ca_rel(abs(all_ca_rel) <= half_cycle);
    yyaxis(ax_ov, 'right');
    if numel(all_ca_rel_trim) >= 2
        peth_phase_q = linspace(-pi, pi, 200);
        peth_time_q  = peth_phase_q * one_cycle / (2*pi);
        [peth_kde_f, ~] = ksdensity(all_ca_rel_trim, peth_time_q);
        peth_kde_rate = peth_kde_f * numel(all_ca_rel_trim) / nValid;
        plot(ax_ov, peth_phase_q, peth_kde_rate, 'r-', 'LineWidth', 1.5);
    end
    ylabel(ax_ov, 'Ca^{2+} spike rate');
    ax_ov.YColor = 'r';

    % Bottom x-axis: phase (radians)
    xlim(ax_ov, [-pi pi]);
    set(ax_ov, 'XTick', [-pi -pi/2 0 pi/2 pi], ...
        'XTickLabel', {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
    xlabel(ax_ov, 'Phase (rad)');

    % Top x-axis: time (overlaid transparent axes)
    ax_top = axes(fig, 'Units','centimeters', ...
        'Position', [sq_x2 row1_y_cm sq sq], ...
        'Color', 'none', 'XAxisLocation', 'top', ...
        'YTick', [], 'Box', 'off');
    xlim(ax_top, [-half_cycle half_cycle]);
    xlabel(ax_top, 'Time (s)');
    linkprop([ax_ov ax_top], 'Position');

    title(ax_ov, sprintf('#spike = %d', nEv_in));
    axis(ax_ov, 'square');
    hold(ax_ov, 'off');

    % --- (2) Breath-triggered average dF/F ---
    ax1 = axes(fig, 'Units','centimeters', ...
        'Position', [margin_l row2_y plot_w plot_h]);
    hold(ax1, 'on');

    % Individual traces (thin, light gray)
    snip_valid = snippets(valid_snip, :);
    for s = 1:size(snip_valid, 1)
        plot(ax1, t_snip, snip_valid(s, :), 'Color', [0.8 0.8 0.8 0.3], 'LineWidth', 0.5);
    end
    % Mean on top (thick black)
    plot(ax1, t_snip, mean_dff, 'k', 'LineWidth', 2);
    xline(ax1, 0, 'k--', 'LineWidth', 1);
    xlim(ax1, [-win_raster, win_raster]);
    ylabel(ax1, 'dF/F');
    set(ax1, 'XTickLabel', []);
    hold(ax1, 'off');

    % --- (3) Breath-triggered Ca spike raster sorted by dt_nearest_ca ---
    ax2 = axes(fig, 'Units','centimeters', ...
        'Position', [margin_l row3_y plot_w plot_h]);
    hold(ax2, 'on');

    % Compute dt_nearest_ca per valid breath peak
    dt_nearest_ca = nan(nValid, 1);
    for b = 1:nValid
        bp = b_frames_valid(b);
        before = ca_ev(ca_ev <= bp);
        after  = ca_ev(ca_ev >= bp);
        dists  = [];
        if ~isempty(before), dists(end+1) = (bp - before(end)) / fs; end
        if ~isempty(after),  dists(end+1) = (after(1) - bp) / fs;   end
        if ~isempty(dists),  dt_nearest_ca(b) = min(dists);          end
    end

    % Sort: NaN → bottom
    dt_near_sort = dt_nearest_ca;
    dt_near_sort(isnan(dt_near_sort)) = Inf;
    [~, sort_near_b] = sort(dt_near_sort, 'ascend');

    for row = 1:nValid
        b       = sort_near_b(row);
        t_ticks = ca_rel_valid{b};
        for tt = 1:numel(t_ticks)
            plot(ax2, [t_ticks(tt) t_ticks(tt)], [row-0.4 row+0.4], ...
                'k', 'LineWidth', 1);
        end
    end

    xline(ax2, 0, 'k--', 'LineWidth', 1);
    xlim(ax2, [-win_raster, win_raster]);
    ylim(ax2, [0 nValid+1]);
    set(ax2, 'YDir', 'reverse', 'XTickLabel', []);
    ylabel(ax2, sprintf('Breath cycle #\n(sorted by dt_{nearest Ca})'));
    hold(ax2, 'off');

    % --- (4) PETH histogram [-win_raster, +win_raster], 1-frame bins ---
    ax3 = axes(fig, 'Units','centimeters', ...
        'Position', [margin_l row4_y plot_w plot_h4]);
    hold(ax3, 'on');

    bar(ax3, peth_centers, peth_rate, 1, ...
        'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    xline(ax3, 0, 'k--', 'LineWidth', 1);
    xlim(ax3, [-win_raster, win_raster]);
    xlabel(ax3, 'Time from breath peak (s)');
    ylabel(ax3, 'Ca^{2+} spike rate (events/s)');
    hold(ax3, 'off');

    % --- (5) Breathing autocorrelogram ---
    ax5 = axes(fig, 'Units','centimeters', ...
        'Position', [margin_l row5_y plot_w plot_h4]);
    hold(ax5, 'on');
    plot(ax5, breath_acf_t, breath_acf, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    xline(ax5, 0, 'k--', 'LineWidth', 1);
    xlim(ax5, [-win_raster, win_raster]);
    ylim(ax5, [-1 1]);
    xlabel(ax5, 'Lag (s)');
    ylabel(ax5, 'Autocorrelation');
    hold(ax5, 'off');

    % ---- Save into breath_trig subfolder ----
    outDir = fullfile(folderPath, 'breath_trig');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    outFig = fullfile(outDir, sprintf('breath_trig_ROI%02d.fig', roi_id));
    savefig(fig, outFig);

    outPdf = fullfile(outDir, sprintf('breath_trig_ROI%02d.pdf', roi_id));
    exportgraphics(fig, outPdf, 'ContentType', 'vector', 'BackgroundColor', 'none');

    fprintf('Saved ROI#%02d:\n  %s\n  %s\n', roi_id, outFig, outPdf);

%{
    % ---- COMMENTED OUT: Ca-centered 2x2 (dt_last / dt_next) ----
    H_last = histcounts(dt_last, bin_edges_dt, 'Normalization', 'probability') * 100;
    H_next = histcounts(dt_next, bin_edges_dt, 'Normalization', 'probability') * 100;
    fig_ca = figure('Color','w','Name',sprintf('dt_last/dt_next: ROI#%02d',roi_id),'Position',[0 0 1000 800]);
    % ... (dt_last raster, dt_next raster, dt_last histogram, dt_next histogram)
    close(fig_ca);

    % ---- COMMENTED OUT: Breath-centered 2x2 (dt_last_ca / dt_next_ca) ----
    dt_last_ca = nan(nBreath,1); dt_next_ca = nan(nBreath,1);
    for b = 1:nBreath
        bp = b_frames(b);
        before = ca_ev(ca_ev <= bp); if ~isempty(before), dt_last_ca(b) = (bp - before(end))/fs; end
        after  = ca_ev(ca_ev >= bp); if ~isempty(after),  dt_next_ca(b) = (after(1) - bp)/fs;   end
    end
    % ... (sorted rasters, histograms)

    % ---- COMMENTED OUT: dt_nearest_ca (breath-centered) ----
    % dt_nearest_ca = min(dt_last_ca, dt_next_ca); ...

    % ---- COMMENTED OUT: Phase analysis (polar + Rayleigh) ----
    % phase_frac = dt_last ./ (dt_last + dt_next); ...
%}
end

fprintf('\nDone. Processed %d spiking ROIs.\n', nSpkROI);
