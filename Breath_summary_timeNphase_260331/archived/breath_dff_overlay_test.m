% breath_dff_overlay_test.m
%  Per-ROI: breath overlay + dF/F in TIME and PHASE domains.
%  Per-FOV: full breathing trace + all spiking dF/F traces.
%
%  Auto-detect mode:
%    - inputPath contains ca_spike_data.mat -> sub mode (single folder)
%    - Otherwise -> master mode (discover all subfolders with ca_spike_data.mat)
%
%  Per spiking ROI (1 row x 2 columns):
%    Left  : breath snippets (gray) + mean dF/F (green) vs time [-win_sec, +win_sec]
%    Right : breath phase cycles (gray) + mean dF/F (green) vs phase [0, 4pi]
%
%  Per FOV (2 rows x 1 column):
%    Top    : full breathing trace vs time
%    Bottom : all spiking dF/F traces vs time (stacked/offset)
%
%  Dependencies:
%    - helper.dFF_RZ
%    - detect_session_fps.m

clear; close all; clc;

%% ========================= USER PARAMETERS =========================
inputPath       = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\pFN_roi4_z0_512x512_5x_2000f_00001";
fps_img         = 30;             % imaging frame rate (Hz) — fallback
nDrop           = 30;             % frames dropped from DLC/breath peaks (NOT from F_raw)
win_sec         = 0.5;            % half-window for breath-trig avg (seconds)
nPhaseBins      = 200;            % bins for phase-normalized breath cycles
BaselineWinSec  = 20;             % dF/F sliding median baseline window (seconds)
doSave          = false;          % just view, don't save
% ====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 8);

%% 1 -- Auto-detect mode & discover sessions
if exist(fullfile(inputPath, 'ca_spike_data.mat'), 'file') == 2
    fprintf('Sub mode: processing single folder.\n');
    allMat = dir(fullfile(inputPath, 'ca_spike_data.mat'));
else
    allMat = dir(fullfile(inputPath, '**', 'ca_spike_data.mat'));
end

nSess = numel(allMat);
fprintf('Found %d session(s) with ca_spike_data.mat\n', nSess);
assert(nSess > 0, 'No ca_spike_data.mat found under %s', inputPath);

%% 2 -- Main session loop
for kk_sess = 1:nSess
    folderPath = allMat(kk_sess).folder;
    [~, sessName] = fileparts(folderPath);

    fprintf('\n========== [%d/%d] %s ==========\n', kk_sess, nSess, sessName);
    close all;

    try
        %% ---- Auto-find required files ----
        csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
        if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*snapshot_200.csv')); end
        if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-20.csv')); end
        if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-10.csv')); end
        sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
        bp_file  = dir(fullfile(folderPath, '*breath_peak_data.mat'));

        if isempty(csv_hits)
            warning('No DLC CSV found in %s, skipping.', sessName);
            continue;
        end
        if isempty(sam_hits)
            warning('No SAM output found in %s, skipping.', sessName);
            continue;
        end
        if isempty(bp_file)
            warning('No breath_peak_data.mat found in %s, skipping.', sessName);
            continue;
        end

        dlc_csv  = fullfile(folderPath, csv_hits(1).name);
        SAM_path = fullfile(folderPath, sam_hits(1).name);

        fprintf('DLC CSV : %s\n', csv_hits(1).name);
        fprintf('SAM     : %s\n', sam_hits(1).name);

        %% ---- Auto-detect fps, load SAM -> dFF ----
        fs = detect_session_fps(folderPath, fps_img);
        fs_breath = fs;  % breathing camera is frame-triggered 1:1

        SAMload = load(SAM_path);
        F_raw = SAMload.F;
        % Do NOT drop nDrop frames from F_raw — batch QC already saved with tossFrames=0
        dFFout = helper.dFF_RZ(F_raw, 'FPS', fs, 'BaselineWinSec', BaselineWinSec);
        dFF = double(dFFout.dFF);     % [T x N]
        [T, N_roi] = size(dFF);

        fprintf('dFF: [%d x %d] (T x N_roi)\n', T, N_roi);

        %% ---- Load calcium spikes ----
        caLoad    = load(fullfile(folderPath, 'ca_spike_data.mat'));
        ifSpk     = caLoad.ifSpike;
        roiSpk_id = find(ifSpk == 1);
        caData    = caLoad.roi_spikes;
        nSpkROI   = numel(roiSpk_id);

        if nSpkROI == 0
            fprintf('No spiking ROIs, skipping session.\n');
            continue;
        end
        fprintf('Spiking ROIs: %d of %d\n', nSpkROI, N_roi);

        %% ---- Load breath peaks ----
        breath_peak = load(fullfile(bp_file(1).folder, bp_file(1).name));
        insp_onsets = sort(breath_peak.insp_onset_idx(:));
        fp          = breath_peak.findpeak_params;

        % Adjust for dropped frames
        insp_onsets(insp_onsets <= nDrop) = [];
        insp_onsets = insp_onsets - nDrop;

        if numel(insp_onsets) < 2
            warning('Fewer than 2 breath peaks after drop, skipping %s.', sessName);
            continue;
        end

        %% ---- Reconstruct breathing signal from DLC ----
        dlc_raw  = readmatrix(dlc_csv, 'NumHeaderLines', 3);
        data_dlc = dlc_raw;
        data_dlc(1:min(nDrop, size(data_dlc,1)), :) = [];

        dot_cols = [2 3 4; 5 6 7; 8 9 10; 11 12 13];
        dot_idx  = struct('dot1',1,'dot2',2,'dot3',3,'dot4',4);
        nActiveDots = numel(fp.dot_selection);
        traces   = NaN(size(data_dlc,1), nActiveDots);

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
            case 'sum',  breath_dlc = sum(traces,  2, 'omitnan');
            otherwise,   breath_dlc = mean(traces, 2, 'omitnan');
        end
        if fp.inverted, breath_dlc = -breath_dlc; end
        breath_dlc = detrend(breath_dlc);
        breath_dlc = breath_dlc - mean(breath_dlc);

        nB = numel(breath_dlc);
        if nB >= T, breath_dlc = breath_dlc(1:T);
        else,       breath_dlc(end+1:T) = 0; end

        %% ---- Compute valid breath peaks & windows ----
        b_frames     = sort(insp_onsets(:));
        win_frames   = round(win_sec * fs);
        nW           = 2 * win_frames + 1;
        t_w          = (-win_frames:win_frames)' / fs;

        edge_margin  = max(win_frames, nDrop);
        insp_valid   = b_frames(b_frames > edge_margin & b_frames <= T - edge_margin);
        nInsp        = numel(insp_valid);

        fprintf('Breath-trig window: +/-%.3f s (%d frames @ %dHz), %d valid onsets of %d total\n', ...
            win_sec, win_frames, round(fs), nInsp, numel(b_frames));

        if nInsp < 2
            warning('Fewer than 2 valid breath peaks, skipping %s.', sessName);
            continue;
        end

        %% ---- Extract breath snippets (per-cycle [0,1] normalization) ----
        breath_snip = zeros(nInsp, nW);
        for e = 1:nInsp
            c    = insp_valid(e);
            snip = breath_dlc(c - win_frames : c + win_frames);
            mn   = min(snip);  mx = max(snip);
            if mx > mn, snip = (snip - mn) / (mx - mn); end
            breath_snip(e, :) = snip;
        end
        mean_breath_snip = mean(breath_snip, 1, 'omitnan');

        %% ---- Phase-normalize each breathing cycle [0, 2pi) ----
        phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1);
        phase_1cyc = phase_1cyc(1:end-1);   % [0, 2pi)
        nCyc       = numel(insp_valid) - 1;
        cycle_mat  = nan(nCyc, nPhaseBins);
        for c = 1:nCyc
            i1 = insp_valid(c);
            i2 = insp_valid(c + 1);
            if i1 >= 1 && i2 <= T
                chunk = breath_dlc(i1:i2);
                orig_phase = linspace(0, 2*pi, numel(chunk));
                cycle_mat(c, :) = interp1(orig_phase, chunk, phase_1cyc, 'linear');
            end
        end
        valid_cyc = ~all(isnan(cycle_mat), 2);
        cycle_mat = cycle_mat(valid_cyc, :);
        nValidCyc = size(cycle_mat, 1);

        % Normalize each cycle to [0,1] for display
        for c = 1:nValidCyc
            mn = min(cycle_mat(c,:));  mx = max(cycle_mat(c,:));
            if mx > mn
                cycle_mat(c,:) = (cycle_mat(c,:) - mn) / (mx - mn);
            end
        end
        mean_cycle = mean(cycle_mat, 1, 'omitnan');

        % Duplicate to [0, 4pi] for display
        phase_2cyc    = [phase_1cyc, phase_1cyc + 2*pi];
        cycle_mat_2x  = [cycle_mat, cycle_mat];           % [nValidCyc x 2*nPhaseBins]
        mean_cycle_2x = [mean_cycle, mean_cycle];

        fprintf('Breathing cycles: %d valid of %d total\n', nValidCyc, nCyc);

        %% ---- Output directory ----
        outDir = fullfile(folderPath, 'breath_dff_overlay');
        if ~exist(outDir, 'dir'), mkdir(outDir); end

        %% ---- Time axis for full traces ----
        t_full = (0:T-1)' / fs;   % [T x 1] seconds

        %% ==============================================================
        %%  PER-ROI FIGURES (1 row x 2 columns)
        %% ==============================================================
        cmap_roi = lines(nSpkROI);

        for kk_roi = 1:nSpkROI
            roi_id = roiSpk_id(kk_roi);
            col    = cmap_roi(kk_roi, :);
            dff_roi = dFF(:, roi_id);

            ca_ev = caData(roi_id).spike_idx(:);
            if isempty(ca_ev)
                fprintf('  ROI#%02d : no Ca events, skipping.\n', roi_id);
                continue;
            end

            %% ---- TIME: breath-triggered dF/F snippets ----
            dff_snips = zeros(nInsp, nW);
            for e = 1:nInsp
                c = insp_valid(e);
                dff_snips(e, :) = dff_roi(c - win_frames : c + win_frames);
            end
            mean_dff_time = mean(dff_snips, 1, 'omitnan');

            %% ---- PHASE: dF/F per cycle, interpolated to phase bins ----
            nCyc_roi = numel(insp_valid) - 1;
            dff_phase_roi = nan(nCyc_roi, nPhaseBins);
            for c_ = 1:nCyc_roi
                i1_ = insp_valid(c_); i2_ = insp_valid(c_ + 1);
                if i1_ >= 1 && i2_ <= T
                    chunk_ = dff_roi(i1_:i2_);
                    orig_ph_ = linspace(0, 2*pi, numel(chunk_));
                    dff_phase_roi(c_, :) = interp1(orig_ph_, chunk_, phase_1cyc, 'linear');
                end
            end
            valid_ph = ~all(isnan(dff_phase_roi), 2);
            dff_phase_valid = dff_phase_roi(valid_ph, :);
            mean_dff_phase  = mean(dff_phase_valid, 1, 'omitnan');

            % Duplicate to [0, 4pi]
            mean_dff_phase_2x = [mean_dff_phase, mean_dff_phase];

            %% ---- FIGURE: 1 row x 2 columns ----
            fig_w = 20;   % cm
            fig_h = 8;    % cm
            fig = figure('Color', 'w', 'Visible', 'on', ...
                'Name', sprintf('FOV#%02d | %s | ROI#%02d', kk_sess, sessName, roi_id), ...
                'Units', 'centimeters', 'Position', [1 1 fig_w fig_h], ...
                'PaperUnits', 'centimeters', 'PaperSize', [fig_w fig_h], ...
                'PaperPosition', [0 0 fig_w fig_h]);

            tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            title(tl, sprintf('%s  |  ROI #%02d  (%d Ca events)', ...
                strrep(sessName, '_', ' '), roi_id, numel(ca_ev)), ...
                'Interpreter', 'none', 'FontSize', 9);

            % ---- Left: TIME overlay ----
            ax1 = nexttile(tl, 1);
            hold(ax1, 'on');
            % Breath snippets (gray)
            for e = 1:nInsp
                plot(ax1, t_w, breath_snip(e, :), 'Color', [0.8 0.8 0.8 0.3], 'LineWidth', 0.3);
            end
            % Mean breath (black)
            plot(ax1, t_w, mean_breath_snip, 'k-', 'LineWidth', 1.5);

            % Mean dF/F on right y-axis (green)
            yyaxis(ax1, 'right');
            plot(ax1, t_w, mean_dff_time, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 2);
            ylabel(ax1, 'Mean dF/F');
            ax1.YAxis(2).Color = [0.2 0.7 0.2];

            yyaxis(ax1, 'left');
            ylabel(ax1, 'Breath (norm.)');
            ax1.YAxis(1).Color = [0 0 0];

            xline(ax1, 0, 'k--', 'LineWidth', 0.8);
            xlim(ax1, [-win_sec, win_sec]);
            xlabel(ax1, 'Time from breath peak (s)');
            title(ax1, 'Time domain', 'Interpreter', 'none');
            hold(ax1, 'off');
            box(ax1, 'on');

            % ---- Right: PHASE overlay ----
            ax2 = nexttile(tl, 2);
            hold(ax2, 'on');
            % Breath cycles (gray)
            for c = 1:nValidCyc
                plot(ax2, phase_2cyc, cycle_mat_2x(c, :), 'Color', [0.8 0.8 0.8 0.3], 'LineWidth', 0.3);
            end
            % Mean breath cycle (black)
            plot(ax2, phase_2cyc, mean_cycle_2x, 'k-', 'LineWidth', 1.5);

            % Mean dF/F (green) on right y-axis
            yyaxis(ax2, 'right');
            plot(ax2, phase_2cyc, mean_dff_phase_2x, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 2);
            ylabel(ax2, 'Mean dF/F');
            ax2.YAxis(2).Color = [0.2 0.7 0.2];

            yyaxis(ax2, 'left');
            ylabel(ax2, 'Breath (norm.)');
            ax2.YAxis(1).Color = [0 0 0];

            xlim(ax2, [0, 4*pi]);
            ax2.XTick      = [0 pi 2*pi 3*pi 4*pi];
            ax2.XTickLabel = {'0', '\pi', '2\pi', '3\pi', '4\pi'};
            xlabel(ax2, 'Phase (rad)');
            title(ax2, 'Phase domain', 'Interpreter', 'none');
            hold(ax2, 'off');
            box(ax2, 'on');

            drawnow;

            %% ---- Save per-ROI figure ----
            if doSave
                outPdf = fullfile(outDir, sprintf('breath_dff_overlay_ROI%02d.pdf', roi_id));
                exportgraphics(fig, outPdf, 'ContentType', 'vector', 'BackgroundColor', 'none');
                outFig = fullfile(outDir, sprintf('breath_dff_overlay_ROI%02d.fig', roi_id));
                savefig(fig, outFig);
                fprintf('  ROI#%02d : saved %s\n', roi_id, outPdf);
                close(fig);
            end

        end  % end ROI loop

        %% ==============================================================
        %%  PER-FOV FIGURE (2 rows x 1 column)
        %% ==============================================================
        fig_fov_w = 24;   % cm
        fig_fov_h = 14;   % cm
        fig_fov = figure('Color', 'w', 'Visible', 'on', ...
            'Name', sprintf('FOV#%02d | %s | overview', kk_sess, sessName), ...
            'Units', 'centimeters', 'Position', [1 1 fig_fov_w fig_fov_h], ...
            'PaperUnits', 'centimeters', 'PaperSize', [fig_fov_w fig_fov_h], ...
            'PaperPosition', [0 0 fig_fov_w fig_fov_h]);

        tl_fov = tiledlayout(fig_fov, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl_fov, sprintf('%s  (FOV overview)', strrep(sessName, '_', ' ')), ...
            'Interpreter', 'none', 'FontSize', 10);

        % ---- Row 1: full breathing trace ----
        ax_b = nexttile(tl_fov, 1);
        plot(ax_b, t_full, breath_dlc, 'k-', 'LineWidth', 0.4);
        hold(ax_b, 'on');
        % Mark breath peaks
        peak_idx_valid = insp_valid(insp_valid >= 1 & insp_valid <= T);
        plot(ax_b, t_full(peak_idx_valid), breath_dlc(peak_idx_valid), ...
            'rv', 'MarkerSize', 3, 'MarkerFaceColor', 'r');
        hold(ax_b, 'off');
        xlabel(ax_b, 'Time (s)');
        ylabel(ax_b, 'Breathing (a.u.)');
        title(ax_b, sprintf('Breathing trace  (%d peaks)', numel(peak_idx_valid)), ...
            'Interpreter', 'none');
        xlim(ax_b, [0, t_full(end)]);
        box(ax_b, 'on');

        % ---- Row 2: all spiking dF/F traces (offset-stacked) ----
        ax_d = nexttile(tl_fov, 2);
        hold(ax_d, 'on');

        % Compute offset spacing: use 95th percentile range per ROI
        offsets = zeros(nSpkROI, 1);
        spacing = 0;
        for jj = 1:nSpkROI
            roi_id = roiSpk_id(jj);
            trace_j = dFF(:, roi_id);
            p5  = prctile(trace_j, 2.5);
            p95 = prctile(trace_j, 97.5);
            range_j = p95 - p5;
            if jj == 1
                offsets(jj) = 0;
            else
                offsets(jj) = offsets(jj-1) + spacing;
            end
            spacing = max(range_j, 0.1);  % minimum gap
        end

        ytick_pos   = zeros(nSpkROI, 1);
        ytick_label = cell(nSpkROI, 1);
        for jj = 1:nSpkROI
            roi_id = roiSpk_id(jj);
            col = cmap_roi(jj, :);
            trace_j = dFF(:, roi_id) + offsets(jj);
            plot(ax_d, t_full, trace_j, '-', 'Color', col, 'LineWidth', 0.4);
            ytick_pos(jj) = median(trace_j);
            ytick_label{jj} = sprintf('ROI#%02d', roi_id);
        end

        ax_d.YTick      = ytick_pos;
        ax_d.YTickLabel = ytick_label;
        hold(ax_d, 'off');
        xlabel(ax_d, 'Time (s)');
        ylabel(ax_d, 'dF/F (stacked)');
        title(ax_d, sprintf('%d spiking ROIs', nSpkROI), 'Interpreter', 'none');
        xlim(ax_d, [0, t_full(end)]);
        box(ax_d, 'on');

        drawnow;

        %% ==============================================================
        %%  Per-FOV figure 3: breathing + dF/F overlay (one panel per ROI)
        %% ==============================================================
        fig_ov_w = 24;
        fig_ov_h = 4 * nSpkROI + 2;
        fig_ov = figure('Color', 'w', 'Visible', 'on', ...
            'Name', sprintf('FOV#%02d | %s | breath+dFF overlay', kk_sess, sessName), ...
            'Units', 'centimeters', 'Position', [1 1 fig_ov_w fig_ov_h]);

        tl_ov = tiledlayout(fig_ov, nSpkROI, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl_ov, sprintf('%s  breath + dF/F overlay', strrep(sessName, '_', ' ')), ...
            'Interpreter', 'none', 'FontSize', 10);

        % Scale breathing to 3x its range
        breath_scaled = breath_dlc * 3;

        for jj = 1:nSpkROI
            roi_id = roiSpk_id(jj);
            dff_roi = dFF(:, roi_id);

            % Shift dF/F so its baseline sits on top of the breathing mean
            breath_top = mean(breath_scaled) + std(breath_scaled);
            dff_shifted = dff_roi - mean(dff_roi) + breath_top;

            ax_ov = nexttile(tl_ov, jj);
            hold(ax_ov, 'on');
            plot(ax_ov, t_full, breath_scaled, 'k-', 'LineWidth', 0.3);
            plot(ax_ov, t_full, dff_shifted, '-', 'Color', [0.2 0.7 0.2], 'LineWidth', 0.4);
            hold(ax_ov, 'off');
            xlim(ax_ov, [0, t_full(end)]);
            ylabel(ax_ov, sprintf('ROI#%02d', roi_id));
            if jj == nSpkROI
                xlabel(ax_ov, 'Time (s)');
            else
                set(ax_ov, 'XTickLabel', []);
            end
            box(ax_ov, 'on');
        end

        drawnow;

        %% ---- Save breath+dFF overlay as starter page ----
        outDir = fullfile(folderPath, 'breath_dff_overlay');
        if ~exist(outDir, 'dir'), mkdir(outDir); end
        starterPdf = fullfile(outDir, sprintf('FOV%02d_%s_starter.pdf', kk_sess, sessName));
        set(fig_ov, 'PaperUnits', 'centimeters', ...
            'PaperSize', [fig_ov_w fig_ov_h], ...
            'PaperPosition', [0 0 fig_ov_w fig_ov_h]);
        exportgraphics(fig_ov, starterPdf, 'ContentType', 'vector', 'BackgroundColor', 'none');
        fprintf('  Starter page saved: %s\n', starterPdf);

        %% ---- Save per-FOV figure ----
        if doSave
            outPdf_fov = fullfile(outDir, 'fov_overview.pdf');
            exportgraphics(fig_fov, outPdf_fov, 'ContentType', 'vector', 'BackgroundColor', 'none');
            outFig_fov = fullfile(outDir, 'fov_overview.fig');
            savefig(fig_fov, outFig_fov);
            fprintf('  FOV overview saved: %s\n', outPdf_fov);
            close(fig_fov);
        end

        fprintf('Session %s complete: %d ROI figures + 1 FOV overview\n', sessName, nSpkROI);

    catch ME
        fprintf('ERROR in %s:\n  %s\n', sessName, ME.message);
        for me_k = 1:numel(ME.stack)
            fprintf('  in %s (line %d)\n', ME.stack(me_k).name, ME.stack(me_k).line);
        end
    end
end  % end session loop

fprintf('\n========== Done ==========\n');
