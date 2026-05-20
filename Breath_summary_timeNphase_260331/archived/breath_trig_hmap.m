%% breath_trig_hmap.m  --  Breath-triggered dF/F heatmap per spiking ROI
%  Replaces the tick-mark raster with a colour-coded dF/F heatmap.
%  Rows = breath onsets, sorted by dt_nearest_ca (closest Ca spike first).
%  Same window, data loading, and figure geometry as breath_trig.m.
%
%  Layout per ROI:
%    Row 1: [Avg projection (square)]  [Overlay summary (square)]
%    Row 2: Breath-triggered avg dF/F trace (2:1)
%    Row 3: dF/F heatmap sorted by dt_nearest_ca (2:1)
%    Row 4: PETH histogram (2:1, half-height)
%    Row 5: Breathing autocorrelogram (2:1, half-height)
%
%  Saves .fig + .pdf per ROI into folderPath/breath_trig_hmap/.
%  ── USER PARAMETERS ─────────────────────────────────────────────────────
%clear; close all; clc;
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

folderPath = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi2_z0_25x_512x512_3x_2000f_00001";

fps_img        = 30;
nDrop          = 30;
PixelSizeBase  = 1.7778;   % um/px at 1x zoom (ScanImage default)
crop_um        = 25;        % ROI crop box size (um)
gamma_val      = 0.6;       % gamma correction
clip_lo        = 0.5;       % brightness clip percentile (low)
clip_hi        = 99.5;      % brightness clip percentile (high)
prcLim         = [0.1 99.9]; % heatmap color percentile limits
signedSort     = true;     % false = sort by |dt|, true = sort by signed dt (neg=Ca before breath)
% ──────────────────────────────────────────────────────────────────────────

%% ---- Auto-find files ----
csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*snapshot_200.csv')); end
if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-20.csv')); end
if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-10.csv')); end
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

%% ---- Reconstruct breathing signal ----
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

fprintf('Breathing signal: %d samples.\n', numel(breath));

%% ---- Breath intervals ----
b_frames = sort(insp_onsets(:));
Delta_k  = diff(b_frames) / fs;

fprintf('Breath intervals: n=%d, mean=%.3f s, range=[%.3f, %.3f] s\n', ...
    numel(Delta_k), mean(Delta_k), min(Delta_k), max(Delta_k));

%% ---- Raster / histogram parameters ----
win_raster = mean(Delta_k) * 2.5;
win5       = round(win_raster * fs);
bin_dt     = 1 / fs;
bin_edges_dt = 0 : bin_dt : win_raster;

fprintf('Window: +/- %.2f s  |  Histogram bin: %.4f s (1 frame)\n', ...
    win_raster, bin_dt);

%% ---- Breathing autocorrelogram ----
max_lag_frames = round(win_raster * fs);
[breath_acf, breath_acf_lags] = xcorr(breath, max_lag_frames, 'coeff');
breath_acf_t = breath_acf_lags / fs;

%% ---- One-cycle mapping for overlay square ----
one_cycle  = mean(Delta_k);
half_cycle = one_cycle / 2;
time_to_phase = @(t) t * 2*pi / one_cycle;

acf_mask       = abs(breath_acf_t) <= half_cycle;
acf_1cyc_phase = time_to_phase(breath_acf_t(acf_mask));
acf_1cyc       = breath_acf(acf_mask);

%% ---- Bluewhitered colormap ----
cmap = bluewhitered;

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

    fprintf('ROI#%02d : %d Ca events in %d breath intervals\n', ...
        roi_id, nEv_in, numel(Delta_k));

    % ---- dF/F snippets around each breath onset ----
    nBreath    = numel(b_frames);
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

    valid_idx       = find(valid_snip);
    b_frames_valid  = b_frames(valid_idx);
    ca_rel_valid    = cell(nValid, 1);
    for b = 1:nValid
        bp = b_frames_valid(b);
        nearby_ca = ca_ev(ca_ev >= (bp - win5) & ca_ev <= (bp + win5));
        ca_rel_valid{b} = (nearby_ca - bp) / fs;
    end

    % ---- PETH ----
    all_ca_rel   = cell2mat(ca_rel_valid);
    peth_edges   = -win_raster : bin_dt : win_raster;
    peth_counts  = histcounts(all_ca_rel, peth_edges);
    peth_rate    = peth_counts / (nValid * bin_dt);
    peth_centers = peth_edges(1:end-1) + bin_dt/2;

    % ---- Sort by dt_nearest_ca ----
    % signedSort=false: |dt| ascending (closest Ca first)
    % signedSort=true:  signed dt ascending (Ca before breath at top, after at bottom)
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

    % Sorted snippet matrix for heatmap
    snip_valid  = snippets(valid_idx, :);       % [nValid x nTime]
    snip_sorted = snip_valid(sort_near_b, :);   % sorted by dt_nearest_ca

    % Percentile-based color limits
    cHigh = prctile(snip_sorted(:), prcLim(2));

    % ============================================================
    % Two square panels: (L) BW heatmap, (R) Ca spike raster
    % Both sorted by dt_nearest_ca, same window
    % ============================================================
    cmap_bw = flipud(gray(256));               % black = high dF/F

    %{
    % --- Blue–White–Red ---
    cmap_bwr = bluewhitered;
    cLow = prctile(snip_sorted(:), prcLim(1));
    cAbs = max(abs(cLow), abs(cHigh));

    % --- Blue–Yellow (yellow = active) ---
    n256 = 256;
    cmap_by = [linspace(0,1,n256)', linspace(0,1,n256)', linspace(0.5,0,n256)'];
    %}

    sq_sz    = 6;       % cm, square side
    hgap     = 1.8;     % cm between panels
    margin_l = 1.8;  margin_b = 1.5;  margin_t = 1.2;
    fig_w = margin_l + 2*sq_sz + hgap + 1.5;
    fig_h = margin_b + sq_sz + margin_t;

    fig = figure('Color', 'w', ...
        'Name', sprintf('breath_trig_hmap ROI#%02d', roi_id), ...
        'Units', 'centimeters', 'Position', [1 1 fig_w fig_h], ...
        'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
        'PaperPosition', [0 0 fig_w fig_h]);
    set(fig, 'DefaultAxesFontSize', 7, 'DefaultTextFontSize', 8);

    sgtitle(fig, sprintf('ROI #%d  —  sorted by dt_{nearest Ca}', roi_id), ...
        'FontSize', 9);

    % --- (L) Black–White dF/F heatmap ---
    ax_hm = axes(fig, 'Units','centimeters', ...
        'Position', [margin_l margin_b sq_sz sq_sz]);
    imagesc(ax_hm, t_snip, 1:nValid, snip_sorted);
    set(ax_hm, 'YDir', 'reverse');
    colormap(ax_hm, cmap_bw);
    caxis(ax_hm, [0 cHigh]);
    hold(ax_hm, 'on');
    xline(ax_hm, 0, 'r--', 'LineWidth', 1);
    hold(ax_hm, 'off');
    xlim(ax_hm, [-win_raster, win_raster]);
    ylim(ax_hm, [0.5, nValid + 0.5]);
    xlabel(ax_hm, 'Time from breath peak (s)');
    ylabel(ax_hm, 'Breath # (sorted)');
    title(ax_hm, 'dF/F heatmap');
    cb = colorbar(ax_hm);
    ylabel(cb, 'dF/F');
    pbaspect(ax_hm, [1 1 1]);

    % --- (R) Ca spike raster (same sorting) ---
    ax_rs = axes(fig, 'Units','centimeters', ...
        'Position', [margin_l + sq_sz + hgap, margin_b, sq_sz, sq_sz]);
    hold(ax_rs, 'on');
    for row = 1:nValid
        b       = sort_near_b(row);
        t_ticks = ca_rel_valid{b};
        for tt = 1:numel(t_ticks)
            plot(ax_rs, [t_ticks(tt) t_ticks(tt)], [row-0.4 row+0.4], ...
                'k', 'LineWidth', 1);
        end
    end
    xline(ax_rs, 0, 'r--', 'LineWidth', 1);
    xlim(ax_rs, [-win_raster, win_raster]);
    ylim(ax_rs, [0.5, nValid + 0.5]);
    set(ax_rs, 'YDir', 'reverse');
    xlabel(ax_rs, 'Time from breath peak (s)');
    set(ax_rs, 'YTickLabel', []);
    title(ax_rs, 'Ca^{2+} spike raster');
    pbaspect(ax_rs, [1 1 1]);

    % ---- Save ----
    outDir = fullfile(folderPath, 'breath_trig_hmap');
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    outFig = fullfile(outDir, sprintf('breath_trig_hmap_ROI%02d.fig', roi_id));
    savefig(fig, outFig);

    outPdf = fullfile(outDir, sprintf('breath_trig_hmap_ROI%02d.pdf', roi_id));
    exportgraphics(fig, outPdf, 'ContentType', 'vector', 'BackgroundColor', 'none');

    fprintf('Saved ROI#%02d:\n  %s\n  %s\n', roi_id, outFig, outPdf);
end

fprintf('\nDone. Processed %d spiking ROIs.\n', nSpkROI);
