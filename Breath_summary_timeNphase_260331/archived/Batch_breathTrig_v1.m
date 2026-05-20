clear; close all; clc;
%% batch_breathTrig_v1.m
%  Batch breath-triggered analysis with:
%    Row 1: [Avg projection (square)] [Overlay: dF/F + PETH KDE + breath autocorr]
%    Row 2: Breath-triggered dF/F (individual + mean)
%    Row 3: Ca spike raster sorted by dt_nearest_ca
%    Row 4: PETH histogram (Ca spike rate)
%    Row 5: Breathing autocorrelogram
%
%  Also creates:
%    1) Per-ROI .fig and .pdf
%    2) Per-FOV merged PDF
%    3) Giant merged PDF across all FOVs
%
%  Dependencies:
%    - helper.dFF_RZ
%    - ghostscript on MATLAB path (used by append_pdfs_RZ)

%% ========================= USER PARAMETERS =========================
masterFolder   = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";
skipIdx        = [];          % session indices to skip
fps_img        = 30;          % imaging frame rate (Hz)
nDrop          = 30;          % frames to drop from start
PixelSizeBase  = 1.7778;      % um/px at 1x zoom (ScanImage default)
crop_um        = 25;          % ROI crop box size (um)
gamma_val      = 0.6;         % gamma correction
clip_lo        = 0.5;         % brightness clip percentile (low)
clip_hi        = 99.5;        % brightness clip percentile (high)
mergePerFOV    = true;        % create per-FOV merged PDF
mergeAllFOV    = true;        % create giant merged PDF across all FOVs
% ==================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

%% 1 ── Discover sessions with ca_spike_data.mat
allMat = dir(fullfile(masterFolder, '**', 'ca_spike_data.mat'));

skipIdx = flip(skipIdx);
keep = true(numel(allMat),1);
skipIdx = skipIdx(skipIdx >= 1 & skipIdx <= numel(allMat));
keep(skipIdx) = false;
allMat = allMat(keep);

nSess = numel(allMat);
fprintf('Found %d sessions with ca_spike_data.mat\n', nSess);
assert(nSess > 0, 'No ca_spike_data.mat found under %s', masterFolder);

nProcessed = 0;
nSkipped   = 0;
nFailed    = 0;
failedSessions = {};

% ---- Summary PDF accumulators ----
allFovSummaryPdfs = {};
allRoiMeta = struct('sessIdx', {}, 'sessName', {}, 'roi_id', {}, 'pdf', {});

%% 2 ── Main loop
for kk_sess = 1:nSess
    folderPath = allMat(kk_sess).folder;
    [~, sessName] = fileparts(folderPath);

    fprintf('\n========== [%d/%d] %s ==========\n', kk_sess, nSess, sessName);
    close all;

    sessRoiPdfs = {};

    try
        %% ---- Auto-find required files ----
        csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
        sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
        bp_file  = dir(fullfile(folderPath, '*breath_peak_data.mat'));

        if isempty(csv_hits)
            warning('No DLC CSV found in %s, skipping.', sessName);
            nSkipped = nSkipped + 1;
            continue;
        end
        if isempty(sam_hits)
            warning('No SAM output found in %s, skipping.', sessName);
            nSkipped = nSkipped + 1;
            continue;
        end
        if isempty(bp_file)
            warning('No breath_peak_data.mat found in %s, skipping.', sessName);
            nSkipped = nSkipped + 1;
            continue;
        end

        dlc_csv = fullfile(folderPath, csv_hits(1).name);
        SAM     = fullfile(folderPath, sam_hits(1).name);

        fprintf('DLC CSV : %s\n', csv_hits(1).name);
        fprintf('SAM     : %s\n', sam_hits(1).name);

        %% ---- Auto-detect fps, load SAM → dFF ----
        fs = detect_session_fps(folderPath, fps_img);

        SAMload = load(SAM);
        F = SAMload.F;
        F(1:min(nDrop,size(F,1)), :) = [];
        dFFout = helper.dFF_RZ(F, 'FPS', fs);
        dFF = dFFout.dFF;
        maskL = SAMload.maskL;

        [T, N_roi] = size(dFF);

        fprintf('dFF: [%d x %d] (T x N_roi)\n', T, N_roi);

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

                tok = regexp(meta_str, ...
                    "SI\.hRoiManager\.scanZoomFactor\s*=\s*([-+]?\d+\.?\d*)", ...
                    "tokens", "once");

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
        end

        %% ---- Load calcium spikes ----
        in = load(fullfile(folderPath, 'ca_spike_data.mat'));
        ifSpk       = in.ifSpike;
        roiSpk_id   = find(ifSpk == 1);
        ca_spk_data = in.roi_spikes;
        nSpkROI     = numel(roiSpk_id);

        ca_spk_id = cell(1, nSpkROI);
        for kk = 1:nSpkROI
            i = roiSpk_id(kk);
            ca_spk_id{kk} = ca_spk_data(i).spike_idx;
        end

        if nSpkROI == 0
            fprintf('No spiking ROIs, skipping session.\n');
            nSkipped = nSkipped + 1;
            continue;
        end
        fprintf('Spiking ROIs: %d of %d\n', nSpkROI, N_roi);

        %% ---- Load breath peaks ----
        breath_peak = load(fullfile(bp_file(1).folder, bp_file(1).name));
        insp_onsets = breath_peak.insp_onset_idx;

        insp_onsets(insp_onsets < nDrop) = [];
        insp_onsets = insp_onsets - nDrop;
        insp_onsets = insp_onsets(:);

        % Exclude peaks within nDrop frames of either end
        insp_onsets = insp_onsets(insp_onsets > nDrop & insp_onsets <= T - nDrop);
        if numel(insp_onsets) < 1
            warning('No valid breath peaks, skipping %s.', sessName);
            nSkipped = nSkipped + 1;
            continue;
        end

        %% ---- Reconstruct breathing signal & Hilbert phase ----
        dlc_raw  = readmatrix(dlc_csv, 'NumHeaderLines', 3);
        data_dlc = dlc_raw;
        data_dlc(1:min(nDrop,size(data_dlc,1)), :) = [];

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
                case 'x'
                    sig = xc;
                case 'y'
                    sig = yc;
                case '-x'
                    sig = -xc;
                case '-y'
                    sig = -yc;
                case 'magnitude'
                    sig = sqrt(xc.^2 + yc.^2);
                otherwise
                    sig = xc;
            end

            sig(bad) = NaN;
            if sum(~isnan(sig)) >= 2
                sig = fillmissing(sig, 'linear', 'EndValues', 'nearest');
            end
            traces(:,d) = sig;
        end

        switch fp.combine_method
            case 'sum'
                breath = sum(traces, 2, 'omitnan');
            otherwise
                breath = mean(traces, 2, 'omitnan');
        end

        if fp.inverted
            breath = -breath;
        end

        breath = detrend(breath);
        breath = (breath - mean(breath)) / std(breath);

        % Hilbert instantaneous phase
        breath_analytic = hilbert(breath);
        breath_phase    = angle(breath_analytic); %#ok<NASGU>

        %% ---- Breath intervals & window ----
        b_frames = sort(insp_onsets(:));
        Delta_k  = diff(b_frames) / fs;

        win_raster = mean(Delta_k) * 2.5;   % seconds
        win5       = round(win_raster * fs);
        bin_dt     = 1 / fs;

        fprintf('Breath intervals: n=%d, mean=%.3f s | window=+/-%.2f s\n', ...
            numel(Delta_k), mean(Delta_k), win_raster);

        %% ---- Breathing autocorrelogram (shared across ROIs) ----
        max_lag_frames = round(win_raster * fs);
        [breath_acf, breath_acf_lags] = xcorr(breath, max_lag_frames, 'coeff');
        breath_acf_t = breath_acf_lags / fs;

        %% ---- One-cycle mapping for overlay square ----
        one_cycle  = mean(Delta_k);
        half_cycle = one_cycle / 2 * 1.6;
        time_to_phase = @(t) t * 2*pi / one_cycle;

        acf_mask       = abs(breath_acf_t) <= half_cycle;
        acf_1cyc_phase = time_to_phase(breath_acf_t(acf_mask));
        acf_1cyc       = breath_acf(acf_mask);

        %% ---- Loop over spiking ROIs ----
        for kk_roi = 1:nSpkROI
            roi_id = roiSpk_id(kk_roi);
            ca_ev  = ca_spk_id{kk_roi}(:);

            if isempty(ca_ev)
                continue;
            end

            % Assign each Ca event to its breath interval
            [~, ~, bin_idx] = histcounts(ca_ev, b_frames);
            keep_ca = (bin_idx >= 1) & (bin_idx <= numel(b_frames) - 1);
            ca_in   = ca_ev(keep_ca);
            bin_in  = bin_idx(keep_ca);
            nEv_in  = numel(ca_in);

            if nEv_in == 0
                continue;
            end

            % dt_last / dt_next for each Ca event
            dt_last    = (ca_in - b_frames(bin_in)) / fs; %#ok<NASGU>
            dt_next    = (b_frames(bin_in + 1) - ca_in) / fs; %#ok<NASGU>
            dt_nearest = min((ca_in - b_frames(bin_in)) / fs, ...
                             (b_frames(bin_in + 1) - ca_in) / fs); %#ok<NASGU>

            % Breath-centered data: Ca spikes relative to each breath peak
            nBreath = numel(b_frames);
            ca_rel_breath = cell(nBreath, 1);
            for b = 1:nBreath
                bp = b_frames(b);
                nearby_ca = ca_ev(ca_ev >= (bp - win5) & ca_ev <= (bp + win5));
                ca_rel_breath{b} = (nearby_ca - bp) / fs;
            end

            % Breath-triggered dF/F snippets
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
            valid_idx      = find(valid_snip);
            b_frames_valid = b_frames(valid_idx);
            ca_rel_valid   = ca_rel_breath(valid_idx);

            % PETH computation
            all_ca_rel = cell2mat(ca_rel_valid);
            peth_edges   = -win_raster : bin_dt : win_raster;
            peth_counts  = histcounts(all_ca_rel, peth_edges);
            peth_rate    = peth_counts / (max(nValid,1) * bin_dt);
            peth_centers = peth_edges(1:end-1) + bin_dt/2;

            if nValid < 1 || isempty(mean_dff) || isempty(all_ca_rel)
                fprintf('  ROI#%02d : no valid breath-aligned snippets/PETH, skipping.\n', roi_id);
                continue;
            end

            fprintf('  ROI#%02d : %d Ca events, %d valid breaths\n', ...
                roi_id, nEv_in, nValid);

            % ============================================================
            % Figure layout (cm, Illustrator-friendly)
            % ============================================================
            sc = 1.5;
            fig_w = 18 / sc;
            margin_l = 2.2 / sc;
            margin_r = 0.8 / sc;
            margin_b = 1.5 / sc;
            margin_t = 1.4 / sc;
            vgap     = 1.5 / sc;

            sq = 5.5 / sc;
            sq_gap = 3.0 / sc;
            sq_x1 = margin_l;
            sq_x2 = margin_l + sq + sq_gap;

            plot_w  = fig_w - margin_l - margin_r;
            plot_h  = plot_w / 2;
            plot_h4 = plot_h / 2;

            fig_h = margin_b + ...
                    plot_h4 + vgap + ...
                    plot_h4 + vgap + ...
                    plot_h  + vgap + ...
                    plot_h  + vgap + ...
                    sq + margin_t;

            fig = figure('Color', 'w', 'Visible', 'off', ...
                'Name', sprintf('FOV#%02d | %s | ROI#%02d', kk_sess, sessName, roi_id), ...
                'Units', 'centimeters', 'Position', [1 1 fig_w fig_h], ...
                'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
                'PaperPosition', [0 0 fig_w fig_h]);

            set(fig, 'DefaultAxesFontSize', 7, 'DefaultTextFontSize', 8);

            pageLabel = sprintf('FOV#%02d | %s | ROI#%02d', kk_sess, sessName, roi_id);
            annotation(fig, 'textbox', [0.02 0.975 0.96 0.02], ...
                'String', pageLabel, ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'top', ...
                'Interpreter', 'none', ...
                'FontName', 'Arial', ...
                'FontSize', 8, ...
                'Color', [0 0 0]);

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
                    cy = round(mean(ry));
                    cx = round(mean(rx));
                    crop_px = round(crop_um / um_per_px / 2);

                    [imgH, imgW] = size(avgImg);
                    r1 = max(1, cy - crop_px);
                    r2 = min(imgH, cy + crop_px);
                    c1 = max(1, cx - crop_px);
                    c2 = min(imgW, cx + crop_px);

                    roi_patch = avgImg(r1:r2, c1:c2);

                    lo = prctile(roi_patch(:), clip_lo);
                    hi = prctile(roi_patch(:), clip_hi);
                    if hi <= lo
                        hi = lo + eps;
                    end

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

                    imshow(patchRGB, 'Parent', ax0);

                    sb_px = round(5 / um_per_px);
                    hold(ax0, 'on');
                    sb_y = size(roi_patch,1) - 3;
                    plot(ax0, [3 3+sb_px], [sb_y sb_y], 'w-', 'LineWidth', 2);
                    text(ax0, 3+sb_px/2, sb_y-3, '5 \mum', 'Color', 'w', ...
                        'FontSize', 8, 'HorizontalAlignment', 'center');
                    hold(ax0, 'off');
                else
                    text(ax0, 0.5, 0.5, 'ROI not in mask', ...
                        'HorizontalAlignment', 'center');
                    axis(ax0, 'off');
                end
            else
                text(ax0, 0.5, 0.5, 'No avg projection', ...
                    'HorizontalAlignment', 'center');
                axis(ax0, 'off');
            end

            title(ax0, sprintf('FOV#%02d  ROI#%02d', kk_sess, roi_id), ...
                'Interpreter', 'none');

            % --- (1b) Overlay summary: avg dF/F, PETH KDE, autocorr ---
            ax_ov = axes(fig, 'Units','centimeters', ...
                'Position', [sq_x2 row1_y_cm sq sq]);
            
            % One-cycle time window
            dff_mask = abs(t_snip) <= half_cycle;
            t_1cyc   = t_snip(dff_mask);
            dff_trim = mean_dff(dff_mask);
            
            %yyaxis(ax_ov, 'left');
            plot(ax_ov, t_1cyc, dff_trim, '-', ...
                'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
            ylabel(ax_ov, 'dF/F');
            %ax_ov.YColor = [0.2 0.7 0.2];
            hold(ax_ov, 'on');
            
            % Autocorrelation (black), one cycle, scaled into dF/F range
            acf_t_1cyc = breath_acf_t(acf_mask);
            acf_1cyc_v = breath_acf(acf_mask);
            
            dff_range = max(dff_trim) - min(dff_trim);
            acf_span  = max(acf_1cyc_v) - min(acf_1cyc_v);
            
            if dff_range < eps
                dff_range = 1;
            end
            
            if acf_span < eps
                acf_scaled = ones(size(acf_1cyc_v)) * mean(dff_trim);
            else
                acf_scaled = min(dff_trim) + ...
                    (acf_1cyc_v - min(acf_1cyc_v)) / acf_span * dff_range;
            end
            
            plot(ax_ov, acf_t_1cyc, acf_scaled, 'k-', 'LineWidth', 1);
            
            % % PETH KDE (red): fit on full window, display one cycle
            % all_ca_rel_use = all_ca_rel;   % full +/- win_raster data
            % 
            % yyaxis(ax_ov, 'right');
            % if numel(all_ca_rel_use) >= 2
            %     t_kde = linspace(-win_raster, win_raster, 400);   % fit/evaluate on full window
            %     [peth_kde_f, ~] = ksdensity(all_ca_rel_use, t_kde);
            %     peth_kde_rate = peth_kde_f * numel(all_ca_rel_use) / nValid;
            %     plot(ax_ov, t_kde, peth_kde_rate, 'r-', 'LineWidth', 1.5);
            % end
            % ylabel(ax_ov, 'Ca^{2+} spike rate');
            % ax_ov.YColor = 'r';
            
            xlim(ax_ov, [-half_cycle half_cycle]);   % only show one cycle
            xlabel(ax_ov, 'Time from breath peak (s)');
            
            title(ax_ov, sprintf('#spike = %d', nEv_in));
            axis(ax_ov, 'square');
            hold(ax_ov, 'off');

            % --- (2) Breath-triggered average dF/F ---
            ax1 = axes(fig, 'Units','centimeters', ...
                'Position', [margin_l row2_y plot_w plot_h]);
            hold(ax1, 'on');

            snip_valid = snippets(valid_snip, :);
            for s = 1:size(snip_valid, 1)
                plot(ax1, t_snip, snip_valid(s, :), ...
                    'Color', [0.8 0.8 0.8 0.4], 'LineWidth', 0.5);
            end

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

            dt_nearest_ca = nan(nValid, 1);
            for b = 1:nValid
                bp = b_frames_valid(b);
                before = ca_ev(ca_ev <= bp);
                after  = ca_ev(ca_ev >= bp);

                dists = [];
                if ~isempty(before)
                    dists(end+1) = (bp - before(end)) / fs; %#ok<SAGROW>
                end
                if ~isempty(after)
                    dists(end+1) = (after(1) - bp) / fs; %#ok<SAGROW>
                end
                if ~isempty(dists)
                    dt_nearest_ca(b) = min(dists);
                end
            end

            dt_near_sort = dt_nearest_ca;
            dt_near_sort(isnan(dt_near_sort)) = Inf;
            [~, sort_near_b] = sort(dt_near_sort, 'ascend');

            for row = 1:nValid
                b = sort_near_b(row);
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

            % --- (4) PETH histogram ---
            ax3 = axes(fig, 'Units','centimeters', ...
                'Position', [margin_l row4_y plot_w plot_h4]);
            hold(ax3, 'on');

            bar(ax3, peth_centers, peth_rate, 1, ...
                'FaceColor', [0.85 0.85 0.85], ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.7);

            xline(ax3, 0, 'k--', 'LineWidth', 1);
            xlim(ax3, [-win_raster, win_raster]);
            xlabel(ax3, 'Time from breath peak (s)');
            ylabel(ax3, 'Ca^{2+} spike rate (events/s)');
            hold(ax3, 'off');

            % --- (5) Breathing autocorrelogram ---
            ax5 = axes(fig, 'Units','centimeters', ...
                'Position', [margin_l row5_y plot_w plot_h4]);
            hold(ax5, 'on');

            plot(ax5, breath_acf_t, breath_acf, ...
                'Color', [0.5 0.5 0.5], 'LineWidth', 1);
            xline(ax5, 0, 'k--', 'LineWidth', 1);
            xlim(ax5, [-win_raster, win_raster]);
            ylim(ax5, [-1 1]);
            xlabel(ax5, 'Lag (s)');
            ylabel(ax5, 'Autocorrelation');
            hold(ax5, 'off');

            % ---- Save into breath_trig subfolder ----
            outDir = fullfile(folderPath, 'breath_trig');
            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end

            outFig = fullfile(outDir, sprintf('breath_trig_ROI%02d.fig', roi_id));
            savefig(fig, outFig);

            outPdf = fullfile(outDir, sprintf('breath_trig_ROI%02d.pdf', roi_id));
            exportgraphics(fig, outPdf, 'ContentType', 'vector', 'BackgroundColor', 'none');

            fprintf('  Saved: %s\n', outPdf);

            sessRoiPdfs{end+1} = outPdf; %#ok<SAGROW>
            allRoiMeta(end+1).sessIdx  = kk_sess; %#ok<SAGROW>
            allRoiMeta(end).sessName   = sessName;
            allRoiMeta(end).roi_id     = roi_id;
            allRoiMeta(end).pdf        = outPdf;

            close(fig);
        end

        % ---- Per-FOV merged summary PDF ----
        if mergePerFOV && ~isempty(sessRoiPdfs)
            roiNums = nan(size(sessRoiPdfs));
            for ii = 1:numel(sessRoiPdfs)
                tok = regexp(sessRoiPdfs{ii}, 'ROI(\d+)\.pdf$', 'tokens', 'once');
                if ~isempty(tok)
                    roiNums(ii) = str2double(tok{1});
                end
            end
            [~, ord] = sort(roiNums);
            sessRoiPdfs = sessRoiPdfs(ord);

            outDir = fullfile(folderPath, 'breath_trig');
            sessSummaryPdf = fullfile(outDir, ...
                sprintf('FOV%02d_%s_breath_trig_ALL_ROI.pdf', kk_sess, sessName));

            if exist(sessSummaryPdf, 'file') == 2
                delete(sessSummaryPdf);
            end

            append_pdfs_RZ(sessSummaryPdf, sessRoiPdfs{:});
            fprintf('  Session summary saved: %s\n', sessSummaryPdf);

            allFovSummaryPdfs{end+1} = sessSummaryPdf; %#ok<SAGROW>
        end

        nProcessed = nProcessed + 1;

    catch ME
        nFailed = nFailed + 1;
        failedSessions{end+1} = sessName; %#ok<SAGROW>
        fprintf('ERROR in %s: %s\n', sessName, ME.message);
    end
end

%% 3 ── Giant summary across all FOVs
if mergeAllFOV && ~isempty(allFovSummaryPdfs)
    giantPdf = fullfile(masterFolder, 'breath_trig_GIANT_SUMMARY_allFOV.pdf');

    if exist(giantPdf, 'file') == 2
        delete(giantPdf);
    end

    append_pdfs_RZ(giantPdf, allFovSummaryPdfs{:});
    fprintf('Giant all-FOV summary saved:\n  %s\n', giantPdf);
end

%% 4 ── Summary
fprintf('\n====== BATCH COMPLETE ======\n');
fprintf('Processed: %d | Skipped: %d | Failed: %d\n', nProcessed, nSkipped, nFailed);

if ~isempty(allFovSummaryPdfs)
    fprintf('Per-FOV merged PDFs created: %d\n', numel(allFovSummaryPdfs));
end

if nFailed > 0
    fprintf('Failed sessions:\n');
    for f = 1:numel(failedSessions)
        fprintf('  - %s\n', failedSessions{f});
    end
end

%% ========================================================================
function append_pdfs_RZ(output, varargin)
% APPEND_PDFS_RZ  Concatenate multiple PDFs into one, silently.
% Includes a dummy last-page hack to avoid Ghostscript dropping
% the final input PDF in some installations.

    if nargin < 2
        error('append_pdfs_RZ:NotEnoughInputs', ...
              'Usage: append_pdfs_RZ(output, input1, input2, ...)');
    end

    % Remove empties / non-existing files
    keep = false(size(varargin));
    for i = 1:numel(varargin)
        keep(i) = ischar(varargin{i}) || isstring(varargin{i});
        if keep(i)
            keep(i) = exist(varargin{i}, 'file') == 2;
        end
    end
    inputFiles = varargin(keep);

    if isempty(inputFiles)
        error('append_pdfs_RZ:NoValidInputs', 'No valid input PDFs found.');
    end

    % Determine append vs create
    if exist(output, 'file') == 2
        tmpOutput  = [tempname '.pdf'];
        inputFiles = [{output} inputFiles];
    else
        tmpOutput  = output;
    end

    % ----- Create dummy page to protect the last real file -----
    dummyPDF = [tempname '_dummy.pdf'];
    hf = figure('Visible','off');
    axis off;
    set(hf, 'PaperPosition', [0 0 8.5 11]);
    print(hf, dummyPDF, '-dpdf', '-fillpage');
    close(hf);

    inputFiles = [inputFiles {dummyPDF}];
    % -----------------------------------------------------------

    % Ghostscript command file
    cmdfile = [tempname '.txt'];
    fh = fopen(cmdfile, 'w');

    fprintf(fh, ...
        '-q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="%s" -f', ...
        tmpOutput);

    for k = 1:numel(inputFiles)
        fprintf(fh, ' "%s"', inputFiles{k});
    end
    fclose(fh);

    % ---- Silent Ghostscript call ----
    try
        evalc('ghostscript(''@"'' + string(cmdfile) + ''"'' )');
    catch ME
        if exist(cmdfile,'file'), delete(cmdfile); end
        if exist(dummyPDF,'file'), delete(dummyPDF); end
        error('append_pdfs_RZ:GhostscriptFailed', ...
            'Ghostscript PDF append failed: %s', ME.message);
    end
    % ---------------------------------

    % Cleanup
    if exist(cmdfile,'file'), delete(cmdfile); end
    if exist(dummyPDF,'file'), delete(dummyPDF); end

    % Move temp into final if needed
    if ~strcmp(tmpOutput, output)
        movefile(tmpOutput, output, 'f');
    end
end