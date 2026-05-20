% spike_vs_nonspike_dff_260330.m
%  Per spiking ROI: overlay breathing-cycle dF/F snippets in TIME,
%  colored by whether the cycle contained a calcium spike.
%
%  Green  = cycles with spike
%  Black  = cycles without spike
%
%  One figure per spiking ROI, saved as PDF.
%
%  Dependencies:
%    - helper.dFF_RZ
%    - detect_session_fps.m

clear; close all; clc;

%% ========================= USER PARAMETERS =========================
inputPath       = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing";
skipIdx         = [];             % session indices to skip (master mode)
fps_img         = 30;             % imaging frame rate (Hz)
nDrop           = 30;             % frames to drop from start
BaselineWinSec  = 20;             % dF/F sliding median baseline window (seconds)
doSave          = true;           % true = save PDF
% ====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 7);

%% 1 -- Auto-detect mode & discover sessions
if exist(fullfile(inputPath, 'ca_spike_data.mat'), 'file') == 2
    fprintf('Sub mode: processing single folder.\n');
    allMat = dir(fullfile(inputPath, 'ca_spike_data.mat'));
else
    allMat = dir(fullfile(inputPath, '**', 'ca_spike_data.mat'));

    if ~isempty(skipIdx)
        skipIdx = skipIdx(:)';
        keep = true(numel(allMat), 1);
        valid_skip = skipIdx(skipIdx >= 1 & skipIdx <= numel(allMat));
        keep(valid_skip) = false;
        allMat = allMat(keep);
    end
end

nSess = numel(allMat);
fprintf('Found %d session(s) with ca_spike_data.mat\n', nSess);
assert(nSess > 0, 'No ca_spike_data.mat found under %s', inputPath);

%% 2 -- Load fov_map from combined summary
fov_map_path = fullfile(inputPath, 'breath_combined_summary', 'fov_map.mat');
assert(isfile(fov_map_path), 'fov_map.mat not found: %s\nRun breath_combined_summary first.', fov_map_path);
tmp_fm  = load(fov_map_path);
fov_map = tmp_fm.fov_map;   % {sessName, [sig_roi_ids]; ...}
fprintf('Loaded fov_map: %d FOVs\n', size(fov_map, 1));

%% 3 -- Session loop
for kk = 1:nSess
    folderPath = allMat(kk).folder;
    [~, sessName] = fileparts(folderPath);
    fprintf('\n========== [%d/%d] %s ==========\n', kk, nSess, sessName);

    try
        %% ---- Find files ----
        sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
        bp_file  = dir(fullfile(folderPath, '*breath_peak_data.mat'));
        csv_hits = dir(fullfile(folderPath, '*snapshot_best-*.csv'));
        if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*snapshot_200.csv')); end
        if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-20.csv')); end
        if isempty(csv_hits), csv_hits = dir(fullfile(folderPath, '*best-10.csv')); end

        if isempty(csv_hits),  warning('No DLC CSV, skipping.'); continue; end
        if isempty(sam_hits),  warning('No SAM output, skipping.'); continue; end
        if isempty(bp_file),   warning('No breath_peak_data, skipping.'); continue; end

        dlc_csv  = fullfile(folderPath, csv_hits(1).name);
        SAM_path = fullfile(folderPath, sam_hits(1).name);

        %% ---- FPS, SAM, dFF ----
        fs = detect_session_fps(folderPath, fps_img);

        SAMload = load(SAM_path);
        F_raw   = SAMload.F;
        dFFout  = helper.dFF_RZ(F_raw, 'FPS', fs, 'BaselineWinSec', BaselineWinSec);
        dFF     = double(dFFout.dFF);
        [T, ~]  = size(dFF);

        %% ---- Calcium spikes ----
        caLoad    = load(fullfile(folderPath, 'ca_spike_data.mat'));
        ifSpk     = caLoad.ifSpike;
        roiSpk_id = find(ifSpk == 1);
        caData    = caLoad.roi_spikes;
        nSpkROI   = numel(roiSpk_id);

        if nSpkROI == 0
            fprintf('No spiking ROIs, skipping.\n');
            continue;
        end
        fprintf('Spiking ROIs: %d\n', nSpkROI);

        %% ---- Breath peaks ----
        breath_peak = load(fullfile(bp_file(1).folder, bp_file(1).name));
        insp_onsets = sort(breath_peak.insp_onset_idx(:));
        fp          = breath_peak.findpeak_params;

        insp_onsets(insp_onsets <= nDrop) = [];
        insp_onsets = insp_onsets - nDrop;

        if numel(insp_onsets) < 2
            warning('Fewer than 2 breath peaks, skipping.');
            continue;
        end

        %% ---- Reconstruct breathing signal ----
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
            case 'sum',  breath_dlc = sum(traces, 2, 'omitnan');
            otherwise,   breath_dlc = mean(traces, 2, 'omitnan');
        end
        if fp.inverted, breath_dlc = -breath_dlc; end
        breath_dlc = detrend(breath_dlc);
        breath_dlc = breath_dlc - mean(breath_dlc);
        nB = numel(breath_dlc);
        if nB >= T, breath_dlc = breath_dlc(1:T);
        else,       breath_dlc(end+1:T) = 0; end

        %% ---- Valid breath peaks ----
        b_frames    = sort(insp_onsets(:));
        edge_margin = nDrop;
        insp_valid  = b_frames(b_frames > edge_margin & b_frames <= T - edge_margin);
        nInsp       = numel(insp_valid);
        nCyc        = nInsp - 1;

        if nCyc < 2
            warning('Fewer than 2 valid cycles, skipping.');
            continue;
        end
        fprintf('Valid breath cycles: %d\n', nCyc);

        %% ---- Look up FOV# from fov_map ----
        fov_id = NaN;
        for fm = 1:size(fov_map, 1)
            if strcmp(fov_map{fm,1}, sessName)
                fov_id = fm;
                break;
            end
        end
        if isnan(fov_id)
            warning('Session %s not found in fov_map, skipping.', sessName);
            continue;
        end

        %% ---- Per spiking ROI ----

        for jj = 1:nSpkROI
            roi_id  = roiSpk_id(jj);
            dff_roi = dFF(:, roi_id);
            ca_ev   = caData(roi_id).spike_idx(:);
            ca_ev   = ca_ev(ca_ev >= 1 & ca_ev <= T);

            % Collect each cycle's raw dF/F trace + spike flag
            cyc_traces = cell(nCyc, 1);
            cyc_times  = cell(nCyc, 1);
            has_spike  = false(nCyc, 1);

            for c = 1:nCyc
                i1 = insp_valid(c);
                i2 = insp_valid(c + 1);
                cyc_traces{c} = dff_roi(i1:i2);
                cyc_times{c}  = (0:(i2 - i1))' / fs;   % seconds from peak
                has_spike(c)  = any(ca_ev >= i1 & ca_ev < i2);
            end

            nSpkCyc = sum(has_spike);
            nNoCyc  = sum(~has_spike);

            % ---- Plot ----
            scr = get(0, 'ScreenSize');
            fig_roi = figure('Color', 'w', ...
                'Units', 'pixels', ...
                'Position', [100, 50, scr(3)*0.4, scr(4)*0.4], ...
                'Name', sprintf('FOV#%02d ROI#%02d spike vs non-spike', fov_id, roi_id));
            ax = axes(fig_roi); hold(ax, 'on');

            % Non-spike cycles (black)
            for c = find(~has_spike)'
                plot(ax, cyc_times{c}, cyc_traces{c}, ...
                    'Color', [0 0 0 0.3], 'LineWidth', 0.5);
            end
            % Spike cycles (green)
            for c = find(has_spike)'
                plot(ax, cyc_times{c}, cyc_traces{c}, ...
                    'Color', [1 0 0 0.3], 'LineWidth', 0.5);
            end

            axis square;
            hold(ax, 'off');
            xlabel(ax, 'Time from breath peak (s)');
            ylabel(ax, 'dF/F');
            title(ax, sprintf('FOV#%02d ROI#%02d — spike (green %d) | non-spike (black %d)', ...
                fov_id, roi_id, nSpkCyc, nNoCyc), 'FontSize', 8);
            box(ax, 'on');

            fprintf('  FOV#%02d ROI#%02d: %d spike, %d non-spike cycles\n', ...
                fov_id, roi_id, nSpkCyc, nNoCyc);

            if doSave
                pdfName = sprintf('spike_vs_nonspike_FOV%02d_ROI%02d.pdf', fov_id, roi_id);
                exportgraphics(fig_roi, fullfile(folderPath, pdfName), ...
                    'ContentType', 'vector', 'BackgroundColor', 'none');
                fprintf('  Saved: %s\n', pdfName);
            end
        end

    catch ME
        fprintf('ERROR in %s:\n  %s\n', sessName, ME.message);
        for me_k = 1:numel(ME.stack)
            fprintf('  in %s (line %d)\n', ME.stack(me_k).name, ME.stack(me_k).line);
        end
    end
end

fprintf('\nDone.\n');
