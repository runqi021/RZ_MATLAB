%% breath_trig_combined_260327.m
%  Per-ROI combined figure: 5 rows x 2 columns (10 panels).
%
%  Merges breath_trig.m, breath_trig_hmap.m, and calcium_breath_combined.m.
%
%  Layout per ROI:
%    (1,1) Avg projection crop       (1,2) Overlay summary (dF/F+KDE+autocorr)
%    (2,1) Breath-trig dF/F+autocorr (2,2) dt_last vs dt_next scatter
%    (3,1) Sorted dF/F heatmap       (3,2) PSD (ROI + breathing) + coh inset
%    (4,1) Sorted Ca spike raster    (4,2) Phase histogram + cosine fit
%    (5,1) PETH histogram            (5,2) Polar (phase / coherence)
%
%  Saves .fig + .pdf per ROI into folderPath/breath_trig_combined/.
%  ── USER PARAMETERS ─────────────────────────────────────────────────────
clear; close all; clc;
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

folderPath      = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi4_z0_512x512_5x_2000f_00001";
FPS             = 30;       % fallback imaging fps
nDrop           = 30;       % frames to drop
TW              = 3;        % multitaper time-bandwidth (K = 2*TW-1 tapers)
fmin            = 0.01;     % Hz  (PSD lower bound)
fmax            = 14;       % Hz  (PSD upper bound; [] → fs/2)
f_breath_search = [1.5, 3]; % Hz  (search band for breathing peak + coherence)
PixelSizeBase   = 1.7778;   % um/px at 1x zoom (ScanImage default)
crop_um         = 30;       % ROI crop box size (um)
gamma_val       = 0.6;      % gamma correction
clip_lo         = 0.5;      % brightness clip percentile (low)
clip_hi         = 99.5;     % brightness clip percentile (high)
prcLim          = [0.1 99.9]; % heatmap color percentile limits
signedSort      = true;    % false = sort by |dt|, true = signed dt (neg=Ca before breath)
saveFigs        = true;
% ──────────────────────────────────────────────────────────────────────────

%% ============================== LOAD ==============================

% --- *_dFF.mat ---
dffEntry = dir(fullfile(folderPath, '*_dFF.mat'));
assert(numel(dffEntry) == 1, "Expected 1 *_dFF.mat, found %d.", numel(dffEntry));
dffMatPath = fullfile(dffEntry.folder, dffEntry.name);
[~, stem]  = fileparts(dffMatPath);
D   = load(dffMatPath);
dFF = double(D.dFF);   % [T × N]
[T, N] = size(dFF);
if isfield(D,'params') && isfield(D.params,'FPS')
    fs = double(D.params.FPS);
else
    fs = FPS;
    warning("params.FPS not found — using fallback FPS=%g", FPS);
end
if isempty(fmax), fmax = fs / 2; end

% --- *cpSAM_output.mat (for maskL) ---
sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
assert(~isempty(sam_hits), "No *cpSAM_output.mat found.");
SAMload = load(fullfile(folderPath, sam_hits(1).name));
maskL   = SAMload.maskL;

% --- ca_spike_data.mat ---
spikeFile = fullfile(folderPath, 'ca_spike_data.mat');
assert(isfile(spikeFile), "ca_spike_data.mat not found.");
Ks      = load(spikeFile);
roiList = find(logical(Ks.ifSpike(:)'));
caData  = Ks.roi_spikes;
assert(~isempty(roiList), "No spiking ROIs found.");

% --- breath_peak_data.mat ---
bpFile = fullfile(folderPath, 'breath_peak_data.mat');
assert(isfile(bpFile), "breath_peak_data.mat not found.");
BP          = load(bpFile);
insp_onsets = sort(BP.insp_onset_idx(:));
fp          = BP.findpeak_params;
if isfield(fp,'camera_fps'), fs_breath = double(fp.camera_fps);
else,                        fs_breath = FPS; end

% --- DLC CSV → breath signal [T × 1] ---
csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*snapshot_200.csv')); end
if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-20.csv')); end
if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-10.csv')); end
assert(~isempty(csv_hits), "No DLC CSV found.");
dlc_raw  = readmatrix(fullfile(folderPath, csv_hits(1).name), 'NumHeaderLines', 3);
dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
dot_map  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);
nDots    = numel(fp.dot_selection);
traces   = NaN(size(dlc_raw,1), nDots);
for d = 1:nDots
    di  = dot_map.(fp.dot_selection{d});
    xc  = dlc_raw(:, dot_cols(di,1));
    yc  = dlc_raw(:, dot_cols(di,2));
    pc  = dlc_raw(:, dot_cols(di,3));
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
        sig = fillmissing(sig, 'linear', 'EndValues','nearest');
    end
    traces(:,d) = sig;
end
switch fp.combine_method
    case 'sum',  breath_dlc = sum(traces,  2, 'omitnan');
    otherwise,   breath_dlc = mean(traces, 2, 'omitnan');
end
if fp.inverted, breath_dlc = -breath_dlc; end
breath_dlc = detrend(breath_dlc);
breath_dlc = (breath_dlc - mean(breath_dlc)) / std(breath_dlc);
nB = numel(breath_dlc);
if nB >= T, breath_dlc = breath_dlc(1:T);
else,       breath_dlc(end+1:T) = 0; end

% --- Avg projection ---
avg_hits = dir(fullfile(folderPath, '*_AVG_for_CP.tif'));
if ~isempty(avg_hits)
    avgImg = double(imread(fullfile(folderPath, avg_hits(1).name)));
else
    avgImg = [];
end

% --- Zoom factor → um_per_px ---
tif_hits = dir(fullfile(folderPath, '*.tif'));
um_per_px = PixelSizeBase;
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

% --- Chronux path ---
chronuxDir = fullfile(fileparts(mfilename('fullpath')), 'chronux_2_12');
if isfolder(chronuxDir) && ~contains(path, 'chronux')
    addpath(genpath(chronuxDir));
end

fprintf("T=%d  N=%d  %d spiking ROIs  %d insp onsets\n", ...
    T, N, numel(roiList), numel(insp_onsets));

%% ============================== BREATHING PSD + PEAK + FWHM ==============================
params_breath.Fs     = fs_breath;
params_breath.tapers = [TW, 2*TW - 1];
params_breath.pad    = 0;
params_breath.fpass  = [fmin, min(fmax, fs_breath/2)];
params_breath.err    = [2, 0.05];
[Sb, fb, Serr_b] = mtspectrumc(breath_dlc, params_breath);
fbC = fb(:);  Sb = Sb(:);

mask_bp      = fbC >= f_breath_search(1) & fbC <= f_breath_search(2);
[~, rel_idx] = max(Sb(mask_bp));
idx_bp       = find(mask_bp, 1) + rel_idx - 1;
f_breath_peak = fbC(idx_bp);
half_max      = Sb(idx_bp) / 2;

idx_lo = idx_bp;
while idx_lo > 1       && Sb(idx_lo) > half_max, idx_lo = idx_lo - 1; end
idx_hi = idx_bp;
while idx_hi < numel(fbC) && Sb(idx_hi) > half_max, idx_hi = idx_hi + 1; end
f_fwhm = [fbC(idx_lo), fbC(idx_hi)];
fprintf("Breathing peak: %.3f Hz  FWHM: [%.3f, %.3f] Hz\n", ...
    f_breath_peak, f_fwhm(1), f_fwhm(2));

%% ============================== INSPIRATION TIMING ==============================
b_frames = sort(insp_onsets(:));
Delta_k  = diff(b_frames) / fs;
win_sec  = mean(Delta_k) * 2.5;
win      = round(win_sec * fs);
t_w      = (-win:win)' / fs;
nW       = 2*win + 1;

insp_valid = insp_onsets(insp_onsets > win & insp_onsets <= T - win);
nInsp      = numel(insp_valid);
fprintf("Breath-trig window: ±%.3f s  (%d valid onsets)\n", win_sec, nInsp);

% pre-compute breath-triggered dFF for all ROIs
dff_trig = cell(N, 1);
for i = 1:N
    seg = zeros(nInsp, nW);
    for e = 1:nInsp
        c = insp_valid(e);
        seg(e,:) = dFF(c-win : c+win, i);
    end
    dff_trig{i} = seg;
end

% breath-triggered autocorrelogram
breath_trig_seg = zeros(nInsp, nW);
for e = 1:nInsp
    c = insp_valid(e);
    breath_trig_seg(e,:) = breath_dlc(c-win : c+win);
end
breath_autocorr = mean(breath_trig_seg, 1);
if max(abs(breath_autocorr)) > 0
    breath_autocorr = breath_autocorr / max(abs(breath_autocorr));
end

% xcorr-based autocorrelogram (for row 1 overlay)
max_lag_frames = round(win_sec * fs);
[breath_acf, breath_acf_lags] = xcorr(breath_dlc, max_lag_frames, 'coeff');
breath_acf_t = breath_acf_lags / fs;

% one-cycle mapping
one_cycle  = mean(Delta_k);
half_cycle = one_cycle / 2;
time_to_phase = @(t) t * 2*pi / one_cycle;
acf_mask       = abs(breath_acf_t) <= half_cycle;
acf_1cyc_phase = time_to_phase(breath_acf_t(acf_mask));
acf_1cyc       = breath_acf(acf_mask);

% PETH bins
bin_dt     = 1 / fs;
peth_edges = -win_sec : bin_dt : win_sec;

%% ============================== CHRONUX PARAMS ==============================
params_mt.Fs     = fs;
params_mt.tapers = [TW, 2*TW - 1];
params_mt.pad    = 0;
params_mt.fpass  = [fmin, fmax];
params_mt.err    = [2, 0.05];

params_coh.Fs     = fs;
params_coh.tapers = [TW, 2*TW - 1];
params_coh.pad    = 0;
params_coh.fpass  = f_breath_search;
params_coh.err    = [2, 0.01];

% resample breath to calcium fs if needed
if abs(fs_breath - fs) > 0.1
    breath_rs = resample(breath_dlc, round(fs), round(fs_breath));
else
    breath_rs = breath_dlc;
end
T_coh     = min(T, numel(breath_rs));
breath_rs = breath_rs(1:T_coh);

% confC: same for every ROI
sig0 = dFF(1:T_coh, roiList(1)) - mean(dFF(1:T_coh, roiList(1)));
[~, ~, ~, ~, ~, ~, ~, confC_val] = coherencyc(breath_rs, sig0, params_coh);

cmap_roi = lines(numel(roiList));

%% ============================== FIGURE LAYOUT (cm) ==============================
% 4 rows × 2 columns, all panels square.
sq       = 3.2;       % square side (cm)
margin_l = 1.8;
margin_r = 0.6;
margin_b = 1.2;
margin_t = 1.8;
vgap     = 1.2;
hgap     = 2.2;

fig_w = margin_l + 2*sq + hgap + margin_r;
fig_h = margin_b + 4*sq + 3*vgap + margin_t;

% y positions (bottom-up)
row4_y = margin_b;
row3_y = row4_y + sq + vgap;
row2_y = row3_y + sq + vgap;
row1_y = row2_y + sq + vgap;

% x positions
col1_x = margin_l;
col2_x = col1_x + sq + hgap;

%% ============================== PER-ROI FIGURES ==============================
for k = 1:numel(roiList)
    roi_id = roiList(k);
    col    = cmap_roi(k,:);
    sig_c  = dFF(:, roi_id) - mean(dFF(:, roi_id));
    ca_ev  = caData(roi_id).spike_idx(:);

    % --- PSD ---
    [Smt, fmt, Serr_mt] = mtspectrumc(sig_c, params_mt);
    fmt = fmt(:);  Smt = Smt(:);

    % --- coherence ---
    [~, C_coh, phi_coh, ~, ~, ~, f_coh, ~, phistd_coh, Cerr_coh] = ...
        coherencyc(breath_rs, sig_c(1:T_coh), params_coh);
    f_coh = f_coh(:);  C_coh = C_coh(:);  phi_coh = phi_coh(:);

    mask_fwhm = f_coh >= f_fwhm(1) & f_coh <= f_fwhm(2);
    if ~any(mask_fwhm), mask_fwhm = true(size(f_coh)); end
    r_k  = mean(C_coh(mask_fwhm));
    th_k = angle(mean(exp(1i * (-phi_coh(mask_fwhm)))));
    dphi = mean(phistd_coh(mask_fwhm));
    rlo  = max(0,   mean(Cerr_coh(1, mask_fwhm)));
    rhi  = min(1,   mean(Cerr_coh(2, mask_fwhm)));

    % --- breath-trig dF/F snippets ---
    seg_dff  = dff_trig{roi_id};
    dff_mean = mean(seg_dff, 1);

    % --- dF/F snippets around each breath onset (for heatmap) ---
    nBreath    = numel(b_frames);
    win_frames = round(win_sec * fs);
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
    mean_dff_snip = mean(snippets(valid_snip, :), 1, 'omitnan');

    valid_idx      = find(valid_snip);
    b_frames_valid = b_frames(valid_idx);
    win5           = round(win_sec * fs);
    ca_rel_valid   = cell(nValid, 1);
    for b = 1:nValid
        bp = b_frames_valid(b);
        nearby_ca = ca_ev(ca_ev >= (bp - win5) & ca_ev <= (bp + win5));
        ca_rel_valid{b} = (nearby_ca - bp) / fs;
    end

    % --- Ca events in breath intervals (for dt scatter) ---
    [~, ~, bin_idx] = histcounts(ca_ev, b_frames);
    keep_ca = (bin_idx >= 1) & (bin_idx <= numel(b_frames) - 1);
    ca_in   = ca_ev(keep_ca);
    bin_in  = bin_idx(keep_ca);
    nEv_in  = numel(ca_in);

    dt_last_ev = (ca_in - b_frames(bin_in)) / fs;
    dt_next_ev = (b_frames(bin_in + 1) - ca_in) / fs;

    % dt_last / dt_next per spike (for phase scatter)
    dt_last_spk = nan(numel(ca_ev), 1);
    dt_next_spk = nan(numel(ca_ev), 1);
    for e = 1:numel(ca_ev)
        c   = ca_ev(e);
        bef = insp_onsets(insp_onsets <= c);
        aft = insp_onsets(insp_onsets >  c);
        if ~isempty(bef), dt_last_spk(e) = (c - bef(end)) / fs; end
        if ~isempty(aft), dt_next_spk(e) = (aft(1) - c)   / fs; end
    end
    keep3 = ~isnan(dt_last_spk) & ~isnan(dt_next_spk);
    xl3   = dt_last_spk(keep3);
    xn3   = dt_next_spk(keep3);

    % --- PETH ---
    all_ca_rel  = cell2mat(ca_rel_valid);
    peth_counts = histcounts(all_ca_rel, peth_edges);
    peth_rate   = peth_counts / (nValid * bin_dt);
    peth_centers = peth_edges(1:end-1) + bin_dt/2;

    % --- Sort by dt_nearest_ca (signed or unsigned) ---
    dt_nearest_ca = nan(nValid, 1);
    for b = 1:nValid
        bp = b_frames_valid(b);
        before = ca_ev(ca_ev <  bp);
        after  = ca_ev(ca_ev >= bp);
        dt_bef = Inf;  dt_aft = Inf;
        if ~isempty(before), dt_bef = (bp - before(end)) / fs; end
        if ~isempty(after),  dt_aft = (after(1) - bp)    / fs; end
        if signedSort
            if dt_bef <= dt_aft
                dt_nearest_ca(b) = -dt_bef;
            else
                dt_nearest_ca(b) =  dt_aft;
            end
        else
            dt_nearest_ca(b) = min(dt_bef, dt_aft);
        end
    end
    dt_near_sort = dt_nearest_ca;
    dt_near_sort(isnan(dt_near_sort)) = Inf;
    [~, sort_near_b] = sort(dt_near_sort, 'ascend');

    % sorted snippet matrix
    snip_valid  = snippets(valid_idx, :);
    snip_sorted = snip_valid(sort_near_b, :);
    cHigh = prctile(snip_sorted(:), prcLim(2));

    % ================================================================
    % FIGURE  — 2 rows × 4 columns (8 panels)
    %   Row 1: avg crop | overlay | trig dF/F+autocorr | PSD+inset
    %   Row 2: heatmap+raster | PETH | phase histogram | polar
    % ================================================================
    fig = figure('Color', 'w', ...
        'Name', sprintf('breath_trig_combined ROI#%02d', roi_id), ...
        'Units', 'centimeters', 'Position', [1 0.5 fig_w fig_h], ...
        'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
        'PaperPosition', [0 0 fig_w fig_h]);
    set(fig, 'DefaultAxesFontSize', 7, 'DefaultTextFontSize', 8);
    sgtitle(fig, sprintf('%s  —  ROI #%d', strrep(stem,'_','\_'), roi_id), 'FontSize', 9);

    % ---- (1,1) Avg projection crop ----
    ax_11 = axes(fig, 'Units','centimeters', 'Position', [col1_x row1_y sq sq]);
    if ~isempty(avgImg)
        [ry, rx] = find(maskL == roi_id);
        if ~isempty(ry)
            cy = round(mean(ry));  cx = round(mean(rx));
            crop_px = round(crop_um / um_per_px / 2);
            [imgH, imgW] = size(avgImg);
            r1 = max(1, cy - crop_px);  r2 = min(imgH, cy + crop_px);
            c1 = max(1, cx - crop_px);  c2 = min(imgW, cx + crop_px);
            roi_patch = avgImg(r1:r2, c1:c2);
            lo = prctile(roi_patch(:), clip_lo);
            hi = prctile(roi_patch(:), clip_hi);
            roi_patch = (roi_patch - lo) / (hi - lo);
            roi_patch = max(0, min(1, roi_patch));
            roi_patch = roi_patch .^ gamma_val;
            mask_crop = maskL(r1:r2, c1:c2);
            roi_perim = bwperim(mask_crop == roi_id);
            patchRGB = repmat(roi_patch, [1 1 3]);
            yellow = [1 1 0];
            for ch = 1:3
                plane = patchRGB(:,:,ch);
                plane(roi_perim) = yellow(ch);
                patchRGB(:,:,ch) = plane;
            end
            imshow(patchRGB, 'Parent', ax_11);
            sb_px = round(5 / um_per_px);
            hold(ax_11, 'on');
            sb_y = size(roi_patch,1) - 3;
            plot(ax_11, [3 3+sb_px], [sb_y sb_y], 'w-', 'LineWidth', 2);
            text(ax_11, 3+sb_px/2, sb_y-3, '5 \mum', 'Color','w', ...
                'FontSize',7, 'HorizontalAlignment','center');
            hold(ax_11, 'off');
        else
            text(ax_11, 0.5, 0.5, 'ROI not in mask', 'HorizontalAlignment','center');
            axis(ax_11, 'off');
        end
    else
        text(ax_11, 0.5, 0.5, 'No avg projection', 'HorizontalAlignment','center');
        axis(ax_11, 'off');
    end
    title(ax_11, sprintf('ROI#%02d', roi_id));

    % ---- (1,2) Overlay summary: dF/F + autocorr ----
    ax_12 = axes(fig, 'Units','centimeters', 'Position', [col2_x row1_y sq sq]);
    dff_mask  = abs(t_snip) <= half_cycle;
    dff_trim  = mean_dff_snip(dff_mask);
    t_trim    = t_snip(dff_mask);
    hold(ax_12, 'on');
    plot(ax_12, t_trim, dff_trim, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
    ylabel(ax_12, 'dF/F');
    dff_range = max(dff_trim) - min(dff_trim);
    acf_1cyc_t = breath_acf_t(acf_mask);
    acf_scaled = min(dff_trim) + (acf_1cyc - min(acf_1cyc)) / (max(acf_1cyc) - min(acf_1cyc) + eps) * dff_range;
    plot(ax_12, acf_1cyc_t, acf_scaled, 'k-', 'LineWidth', 1);
    xlim(ax_12, [-half_cycle half_cycle]);
    xlabel(ax_12, 'Time (s)');
    title(ax_12, sprintf('#spike = %d', nEv_in));
    axis(ax_12, 'square');
    hold(ax_12, 'off');

    % ---- (1,3) Breath-trig dF/F + yyaxis autocorr ----
    ax_13 = axes(fig, 'Units','centimeters', 'Position', [col1_x row2_y sq sq]);
    hold(ax_13,'on');
    plot(ax_13, t_w, seg_dff', 'Color',[0.82 0.82 0.82], 'LineWidth',0.5);
    plot(ax_13, t_w, dff_mean, 'Color', col, 'LineWidth',2);
    xline(ax_13, 0, 'k--', 'LineWidth',1);
    xlim(ax_13, [-win_sec, win_sec]);
    ylabel(ax_13, '\DeltaF/F');
    set(ax_13, 'XTickLabel', []);
    yyaxis(ax_13, 'right');
    plot(ax_13, t_w, breath_autocorr, 'k--', 'LineWidth',1.5);
    ylabel(ax_13, 'Autocorr');
    ylim(ax_13, [min(breath_autocorr), max(breath_autocorr)]);
    ax_13.YAxis(2).Color = [0 0 0];
    yyaxis(ax_13, 'left');
    ax_13.YAxis(1).Color = [0 0 0];
    hold(ax_13,'off');  box(ax_13,'on');

    % ---- (1,4) PSD + coherence inset ----
    ax_14 = axes(fig, 'Units','centimeters', 'Position', [col2_x row2_y sq sq]);
    hold(ax_14,'on');
    yLo = Serr_mt(1,:)';  yHi = Serr_mt(2,:)';
    patch(ax_14, [fmt; flipud(fmt)], [yLo; flipud(yHi)], col, ...
          'EdgeColor','none', 'FaceAlpha',0.2);
    plot(ax_14, fmt, Smt, 'Color', col, 'LineWidth',1.5);
    yLo_b = Serr_b(1,:)';  yHi_b = Serr_b(2,:)';
    patch(ax_14, [fbC; flipud(fbC)], [yLo_b; flipud(yHi_b)], [0.3 0.3 0.3], ...
          'EdgeColor','none', 'FaceAlpha',0.15);
    plot(ax_14, fbC, Sb, 'k', 'LineWidth',2);
    ax_14.XScale = 'log';  ax_14.YScale = 'log';
    psd_xmin = min(fmt(1), fbC(1));
    ax_14.XLim   = [psd_xmin, fmax];
    all_ticks = [0.01 0.1 1 10];
    all_labels = {'10^{-2}','10^{-1}','10^0','10^1'};
    keep_t = all_ticks >= psd_xmin & all_ticks <= fmax;
    ax_14.XTick  = all_ticks(keep_t);
    ax_14.XTickLabel = all_labels(keep_t);
    grid(ax_14,'on');
    xlabel(ax_14, 'Freq (Hz)');  ylabel(ax_14, 'PSD');
    hold(ax_14,'off');  box(ax_14,'on');
    pbaspect(ax_14, [1 1 1]);

    % ---- (2,1) Sorted dF/F heatmap + Ca raster overlay ----
    ax_21 = axes(fig, 'Units','centimeters', 'Position', [col1_x row3_y sq sq]);
    imagesc(ax_21, t_snip, 1:nValid, snip_sorted);
    set(ax_21, 'YDir', 'reverse');
    colormap(ax_21, flipud(gray(256)));
    caxis(ax_21, [0 cHigh]);
    hold(ax_21, 'on');
    xline(ax_21, 0, 'k--', 'LineWidth', 1);
    % overlay Ca spike raster
    for row = 1:nValid
        b       = sort_near_b(row);
        t_ticks = ca_rel_valid{b};
        for tt = 1:numel(t_ticks)
            plot(ax_21, [t_ticks(tt) t_ticks(tt)], [row-0.4 row+0.4], ...
                'Color', col, 'LineWidth', 2.3);
        end
    end
    hold(ax_21, 'off');
    xlim(ax_21, [-win_sec, win_sec]);
    ylim(ax_21, [0.5, nValid + 0.5]);
    set(ax_21, 'XTickLabel', []);
    ylabel(ax_21, 'Breath # (sorted)');
    cb = colorbar(ax_21);
    ylabel(cb, 'dF/F');
    set(ax_21, 'Units','centimeters', 'Position', [col1_x row3_y sq sq]);
    cb.Units = 'centimeters';
    cb.Position = [col1_x + sq + 0.1, row3_y, 0.15, sq];

    % ---- Row 3, Col 2: Phase histogram + cosine fit ----
    phi_spikes = 2*pi * xl3 ./ (xl3 + xn3);
    nBins    = 20;
    edges    = linspace(0, 2*pi, nBins+1);
    ctrs     = (edges(1:end-1) + edges(2:end)) / 2;
    counts   = histcounts(phi_spikes, edges);
    nCycles  = numel(b_frames) - 1;
    counts_pc = counts / nCycles;          % spikes per cycle per bin
    % duplicate to [0, 4pi]
    ctrs_ext   = [ctrs, ctrs + 2*pi];
    cnts_ext   = [counts_pc, counts_pc];
    r_vec = mean(exp(1i * phi_spikes));
    phi0  = angle(r_vec);
    Rlen  = abs(r_vec);
    B     = numel(phi_spikes) / (nBins * nCycles);
    A     = 2 * B * Rlen;
    th_fit  = linspace(0, 4*pi, 300);
    fit_ext = B + A * cos(th_fit - phi0);

    ax_22 = axes(fig, 'Units','centimeters', 'Position', [col2_x row3_y sq sq]);
    hold(ax_22,'on');
    bar(ax_22, ctrs_ext, cnts_ext, 1, 'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none');
    plot(ax_22, th_fit, fit_ext, 'Color', col, 'LineWidth',2);
    xline(ax_22, 0, 'k--', 'LineWidth',1.2);
    ax_22.XLim   = [0, 4*pi];
    ax_22.YLim   = [0 0.2];
    ax_22.XTick  = [0 pi 2*pi 3*pi 4*pi];
    ax_22.XTickLabel = {'0', '\pi', '2\pi', '3\pi', '4\pi'};
    xlabel(ax_22, 'Phase (rad)');
    ylabel(ax_22, 'Spikes/cycle');
    hold(ax_22,'off');  box(ax_22,'on');
    pbaspect(ax_22, [1 1 1]);

    % ---- Row 4, Col 1: PETH histogram ----
    ax_23 = axes(fig, 'Units','centimeters', 'Position', [col1_x row4_y sq sq]);
    hold(ax_23, 'on');
    bar(ax_23, peth_centers, peth_rate, 1, ...
        'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
    xline(ax_23, 0, 'k--', 'LineWidth', 1);
    xlim(ax_23, [-win_sec, win_sec]);
    ylim(ax_23, [0 10]);
    xlabel(ax_23, 'Time from breath peak (s)');
    ylabel(ax_23, 'Ca^{2+} rate (ev/s)');
    hold(ax_23, 'off');  box(ax_23, 'on');
    pbaspect(ax_23, [1 1 1]);

    % ---- (2,4) Polar plot ----
    fig_pos = fig.Position;
    polar_pos_norm = [col2_x/fig_pos(3), row4_y/fig_pos(4), ...
                      sq/fig_pos(3),     sq/fig_pos(4)];
    ax_24 = polaraxes(fig, 'Position', polar_pos_norm);
    hold(ax_24,'on');
    th_circ = linspace(0, 2*pi, 360);
    polarplot(ax_24, th_circ, repmat(confC_val, 1, 360), 'k--', 'LineWidth',1.2);
    polarplot(ax_24, [th_k th_k], [rlo rhi], '-', 'Color', col, 'LineWidth',2);
    phi_arc = linspace(th_k - dphi, th_k + dphi, 40);
    polarplot(ax_24, phi_arc, r_k * ones(1,40), '-', 'Color', col, 'LineWidth',1.5);
    ph = polarplot(ax_24, th_k, r_k, 'o', ...
                   'MarkerFaceColor', col, 'MarkerEdgeColor','k', 'MarkerSize',8);
    ph.UserData = struct('roi', roi_id, 'phistd', dphi, 'rlo', rlo, 'rhi', rhi);
    ax_24.RLim              = [0, 1];
    ax_24.ThetaZeroLocation = 'top';
    ax_24.ThetaDir          = 'counterclockwise';
    ax_24.FontSize          = 7;
    hold(ax_24,'off');

    dcm = datacursormode(fig);
    dcm.UpdateFcn = @(~, info) polarTip(info);

    drawnow;

    % ---- Coherence inset (bottom-left of PSD panel) ----
    rp   = ax_14.Position;
    inW  = 0.38 * rp(3);
    inH  = 0.38 * rp(4);
    inX  = rp(1) + 0.12 * rp(3);
    inY  = rp(2) + 0.12 * rp(4);
    ax_ins = axes(fig, 'Units','centimeters', 'Position',[inX inY inW inH]);
    hold(ax_ins,'on');
    plot(ax_ins, f_coh, C_coh, 'Color', col, 'LineWidth',1.2);
    patch(ax_ins, [f_fwhm(1) f_fwhm(2) f_fwhm(2) f_fwhm(1)],[0 0 1 1], ...
          [1 0.6 0.6],'EdgeColor','none','FaceAlpha',0.25);
    xline(ax_ins, f_breath_peak,'r:','LineWidth',1.2);
    yline(ax_ins, confC_val,'k--','LineWidth',1);
    hw = diff(f_breath_search) * 0.6;
    ins_fmin = f_breath_peak - hw;
    ins_fmax = f_breath_peak + hw;
    ax_ins.XLim  = [ins_fmin, ins_fmax];
    ax_ins.YLim  = [0, 1];
    ax_ins.XTick = [ins_fmin, ins_fmax];
    ax_ins.XTickLabel = {sprintf('%.1f',ins_fmin), sprintf('%.1f',ins_fmax)};
    ax_ins.YTick = [0, 1];
    ylabel(ax_ins, '|C|');
    ax_ins.XScale = 'linear';
    ax_ins.FontSize = 6;
    ax_ins.Color = [1 1 1 0.85];
    grid(ax_ins,'on');
    hold(ax_ins,'off');  box(ax_ins,'on');
    axis(ax_ins,'square');

    % ---- Save ----
    if saveFigs
        outDir = fullfile(folderPath, 'breath_trig_combined');
        if ~exist(outDir, 'dir'), mkdir(outDir); end

        outFig = fullfile(outDir, sprintf('breath_trig_combined_ROI%02d.fig', roi_id));
        savefig(fig, outFig);

        outPdf = fullfile(outDir, sprintf('breath_trig_combined_ROI%02d.pdf', roi_id));
        exportgraphics(fig, outPdf, 'ContentType','vector', 'BackgroundColor','none');

        fprintf("Saved ROI#%02d:\n  %s\n  %s\n", roi_id, outFig, outPdf);
    end
end

fprintf("\nDone. Processed %d spiking ROIs.\n", numel(roiList));

%% --- datacursor callback ---
function txt = polarTip(info)
    src = info.Target;
    ud  = [];
    if isprop(src,'UserData'), ud = src.UserData; end
    th = info.Position(1);
    r  = info.Position(2);
    if isstruct(ud) && isfield(ud,'roi')
        txt = sprintf('ROI %d\nPhase = %.3g rad (%.1f\xB0) \xB1 %.3g rad\nCoh = %.3g  [%.3g, %.3g]', ...
                      ud.roi, th, rad2deg(th), ud.phistd, r, ud.rlo, ud.rhi);
    else
        txt = sprintf('Phase = %.3g rad (%.1f\xB0)\nCoh = %.3g', th, rad2deg(th), r);
    end
end
