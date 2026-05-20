
clear all; clc; close all;
%% Run_ROI_QC_Pass2.m
%  Two-pass ROI QC — Pass 2: Per-ROI movie + sequential review.
%
%  For each FOV: loads Pass 1 exclusion, generates a survivor montage video,
%  then provides a sequential per-ROI review interface for fine exclusion.
%
%  Two input modes:
%    masterFolder: scan all subfolders with roiQC_exclusion.mat
%    folderPath:   single FOV
%
%  Outputs (per FOV folder):
%    Fmovie_perROI/survivors_montage.avi      — montage video of kept ROIs
%    Fmovie_perROI/survivors_montage_labels.png — reference image
%    roiQC_exclusion.mat                       — updated with Pass 2 results
%
%  RZ 2026-04

%% ======================================================================
%  USER-EDITABLE PARAMETERS
%  ======================================================================
masterFolder    = "D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing";      % multi-experiment scan
folderPath      = '';      % single FOV mode

reviewMode      = "both";  % "sequential" | "video_only" | "both"
skipExisting    = false;    % skip FOVs already reviewed in Pass 2

%% ======================================================================
%  PATH SETUP
%  ======================================================================
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'mod')));
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

%% ======================================================================
%  DISCOVER FOVs (find roiQC_exclusion.mat)
%  ======================================================================
if masterFolder ~= ""
    exclHits = dir(fullfile(masterFolder, '**', 'roiQC_exclusion.mat'));
elseif folderPath ~= ""
    exclHits = dir(fullfile(folderPath, 'roiQC_exclusion.mat'));
    if isempty(exclHits)
        exclHits = dir(fullfile(folderPath, '**', 'roiQC_exclusion.mat'));
    end
else
    error('Set either masterFolder or folderPath.');
end
assert(~isempty(exclHits), 'No roiQC_exclusion.mat found. Run Pass 1 first.');

nFOV = numel(exclHits);
fprintf('\n=== Found %d FOV(s) with Pass 1 exclusion ===\n', nFOV);

%% ======================================================================
%  BATCH LOOP
%  ======================================================================
nProcessed = 0;
nFailed    = 0;
failedList = {};

for ff = 1:nFOV
    exclPath  = fullfile(exclHits(ff).folder, exclHits(ff).name);
    fovFolder = exclHits(ff).folder;
    [~, fovID] = fileparts(fovFolder);

    fprintf('\n========== [%d/%d] %s ==========\n', ff, nFOV, fovID);

    try
        %% ---- Load Pass 1 exclusion ----
        excl = load(exclPath);

        % Skip if Pass 2 already completed
        if skipExisting && isfield(excl, 'pass2_datetime') && ~isnat(excl.pass2_datetime)
            fprintf('  [SKIP] Pass 2 already completed (%s).\n', char(excl.pass2_datetime));
            continue
        end
        N_roi    = excl.N_roi;
        validROI = excl.roiLabels(:)';
        survivors = find(~excl.pass1_reject);
        nSurv = numel(survivors);
        fprintf('  Total ROIs: %d  |  Pass 1 survivors: %d\n', N_roi, nSurv);

        if nSurv == 0
            fprintf('  No survivors. Skipping.\n');
            continue
        end

        %% ---- Find existing montage AVI ----
        survROI = validROI(survivors);

        montageAVI = fullfile(fovFolder, 'Fmovie_perROI', 'all_ROIs_montage.avi');
        assert(isfile(montageAVI), ...
            'No all_ROIs_montage.avi found. Run Run_ROI_montage_video.m first.\n%s', montageAVI);

        vr = VideoReader(montageAVI);
        T   = vr.NumFrames;
        fps = vr.FrameRate;
        fprintf('  Montage AVI: %d frames @ %d fps [%dx%d]\n', T, fps, vr.Width, vr.Height);

        % Derive grid layout from video dimensions + N_roi
        border = 2;
        nCols_grid = ceil(sqrt(N_roi));
        tileW = round((vr.Width - border) / nCols_grid);
        patchSz = tileW - border;
        fprintf('  Grid: %d cols, tile %d px\n', nCols_grid, patchSz);

        %% ---- Sequential review ----
        pass2_reject = false(N_roi, 1);
        pass2_reason = repmat("", N_roi, 1);

        % Load previous Pass 2 decisions if they exist
        prev_pass2_reject = [];
        if isfield(excl, 'pass2_reject') && ~isempty(excl.pass2_reject) && ...
                isfield(excl, 'pass2_datetime') && ~isnat(excl.pass2_datetime)
            prev_pass2_reject = excl.pass2_reject;
            fprintf('  Previous Pass 2 found: %d rejected. Resuming.\n', sum(prev_pass2_reject));
        end

        if reviewMode == "sequential" || reviewMode == "both"
            % Load dFF if available
            dFF_data = [];
            t_dff    = [];
            dff_hits = dir(fullfile(fovFolder, '*_dFF.mat'));
            if ~isempty(dff_hits)
                dffTmp = load(fullfile(dff_hits(1).folder, dff_hits(1).name), 'dFF');
                if isfield(dffTmp, 'dFF')
                    dFF_data = double(dffTmp.dFF);
                    t_dff = (0:size(dFF_data,1)-1)' / fps;
                end
            end

            [pass2_reject, pass2_reason] = sequential_review( ...
                montageAVI, fps, T, nCols_grid, patchSz, border, ...
                dFF_data, t_dff, survROI, survivors, N_roi, validROI, fovID, ...
                prev_pass2_reject);

        elseif reviewMode == "video_only"
            fprintf('  Review the montage: %s\n', montageAVI);
            fprintf('  Enter ROI numbers to exclude.\n');
            excludeStr = input('  Exclude ROIs (comma-separated, Enter to keep all): ', 's');
            if ~isempty(strtrim(excludeStr))
                excludeNums = str2num(excludeStr); %#ok<ST2NM>
                for ei = 1:numel(excludeNums)
                    idx = find(validROI == excludeNums(ei), 1);
                    if ~isempty(idx)
                        pass2_reject(idx) = true;
                        pass2_reason(idx) = "manual";
                    end
                end
            end
        end

        %% ---- Update exclusion data ----
        excl.pass2_reject  = pass2_reject;
        excl.pass2_reason  = pass2_reason;
        excl.final_keep    = ~excl.pass1_reject & ~pass2_reject;
        excl.pass2_datetime = datetime('now');

        save(exclPath, '-struct', 'excl', '-v7.3');

        nFinal = sum(excl.final_keep);
        nP2Rej = sum(pass2_reject);
        fprintf('  Pass 2 exclusion saved.\n');
        fprintf('  Pass 2 rejected: %d  |  Final kept: %d/%d (%.0f%%)\n', ...
            nP2Rej, nFinal, N_roi, 100*nFinal/N_roi);

        nProcessed = nProcessed + 1;

    catch ME
        fprintf(2, '\n*** FAILED [%d/%d]: %s\n%s\n', ff, nFOV, fovID, ME.message);
        for si = 1:numel(ME.stack)
            fprintf(2, '  %s (line %d)\n', ME.stack(si).name, ME.stack(si).line);
        end
        nFailed = nFailed + 1;
        failedList{end+1} = fovID; %#ok<SAGROW>
    end
end

%% ======================================================================
%  CONSOLE SUMMARY
%  ======================================================================
fprintf('\n==================== PASS 2 COMPLETE ====================\n');
fprintf('  Processed: %d\n', nProcessed);
fprintf('  Failed:    %d\n', nFailed);
fprintf('  Total:     %d\n', nFOV);
if nFailed > 0
    fprintf('  Failed FOVs:\n');
    for i = 1:numel(failedList)
        fprintf('    %s\n', failedList{i});
    end
end
fprintf('=========================================================\n');


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function [pass2_reject, pass2_reason] = sequential_review( ...
        montageAVI, fps, T, nCols_grid, patchSz, border, ...
        dFF_data, t_dff, survROI, survivors, N_roi, validROI, fovID, ...
        prev_pass2_reject)
%SEQUENTIAL_REVIEW  Per-ROI review using existing montage AVI.
%  Reads tiles from all_ROIs_montage.avi via VideoReader (no TIFF needed).
%  Left panel: per-ROI movie tile extracted from montage.
%  Right panel: dFF trace colored by decision + frame cursor.
%  Controls: slider, play/pause button, keyboard shortcuts.
    nSurv = numel(survivors);
    pass2_reject = false(N_roi, 1);
    pass2_reason = repmat("", N_roi, 1);

    % Per-survivor decision: 0=undecided, 1=keep, -1=reject
    decision = zeros(nSurv, 1);
    if ~isempty(prev_pass2_reject)
        for j = 1:nSurv
            si = survivors(j);
            if prev_pass2_reject(si)
                decision(j) = -1;
            else
                decision(j) = 1;
            end
        end
        fprintf('  Loaded %d previous Pass 2 decisions.\n', sum(decision ~= 0));
    end

    hasDFF = ~isempty(dFF_data);
    tileW = patchSz + border;

    % Map each survivor to its tile position in the montage grid.
    % The montage uses ROI label as tile index (tile i = ROI label i).
    surv_tileIdx = zeros(nSurv, 1);
    for j = 1:nSurv
        surv_tileIdx(j) = validROI(survivors(j));  % = survROI(j)
    end

    % Pre-compute tile crop coordinates in the montage frame
    surv_y0 = zeros(nSurv, 1);
    surv_x0 = zeros(nSurv, 1);
    for j = 1:nSurv
        ti = surv_tileIdx(j);
        gi = ceil(ti / nCols_grid);
        gj = mod(ti - 1, nCols_grid) + 1;
        surv_y0(j) = (gi - 1) * tileW + border + 1;
        surv_x0(j) = (gj - 1) * tileW + border + 1;
    end

    % dFF <-> video frame alignment
    T_dff = 0;
    if hasDFF, T_dff = size(dFF_data, 1); end
    frameOffset = T - T_dff;
    if frameOffset < 0, frameOffset = 0; end

    % Open VideoReader
    vr = VideoReader(montageAVI);

    % Setup figure
    fig_r = figure('Color', 'w', 'Visible', 'on', ...
        'Name', sprintf('Pass 2 Review — %s', fovID), ...
        'Position', [100 100 1200 550]);

    fig_r.UserData = '';
    fig_r.KeyPressFcn = @(~, evt) set(gcbf, 'UserData', evt.Key);

    % Transport controls
    sliderStep1 = max(1/(T-1), 1e-6);
    sliderStep2 = max(min(round(fps)/(T-1), 1), sliderStep1);
    hPlay = uicontrol(fig_r, 'Style', 'togglebutton', ...
        'String', 'Pause', 'Value', 1, ...
        'Units', 'normalized', 'Position', [0.02 0.01 0.06 0.06], ...
        'FontSize', 9, 'FontWeight', 'bold');
    hSlider = uicontrol(fig_r, 'Style', 'slider', ...
        'Min', 1, 'Max', T, 'Value', 1, ...
        'SliderStep', [sliderStep1, sliderStep2], ...
        'Units', 'normalized', 'Position', [0.09 0.01 0.72 0.05]);
    hFrameLabel = uicontrol(fig_r, 'Style', 'text', ...
        'String', sprintf('1/%d  (0.0 s)', T), ...
        'Units', 'normalized', 'Position', [0.82 0.01 0.16 0.05], ...
        'HorizontalAlignment', 'left', 'FontSize', 8, ...
        'BackgroundColor', 'w');

    cur = 1;
    frame = 1;
    done = false;
    needSetup = true;
    lastFrame = [];  % cache last read montage frame

    while ~done && isvalid(fig_r)
        if needSetup
            rid = survROI(cur);

            if decision(cur) == 1
                traceColor = [0.2 0.7 0.2]; decStr = 'KEEP'; decColor = [0 0.6 0];
            elseif decision(cur) == -1
                traceColor = [0.8 0.1 0.1]; decStr = 'REJECT'; decColor = [0.8 0 0];
            else
                traceColor = [0 0 0]; decStr = 'undecided'; decColor = [0.4 0.4 0.4];
            end

            delete(findobj(fig_r, 'Type', 'axes'));
            delete(findall(fig_r, 'Type', 'textboxshape'));

            % Left: video tile
            ax_vid = axes(fig_r, 'Units', 'normalized', 'Position', [0.02 0.16 0.30 0.74]);
            hImg = imshow(zeros(patchSz, patchSz, 3, 'uint8'), 'Parent', ax_vid);
            title(ax_vid, sprintf('ROI #%d', rid), 'FontSize', 12);

            % Right: dFF trace
            if hasDFF && rid <= size(dFF_data, 2)
                ax_tr = axes(fig_r, 'Units', 'normalized', 'Position', [0.38 0.16 0.58 0.74]);
                plot(ax_tr, t_dff, dFF_data(:, rid), '-', 'Color', traceColor, 'LineWidth', 0.8);
                hold(ax_tr, 'on');
                hCursor = xline(ax_tr, 0, '-', 'Color', [0 0 0], 'LineWidth', 2);
                hold(ax_tr, 'off');
                xlabel(ax_tr, 'Time (s)');
                ylabel(ax_tr, 'dF/F');
                title(ax_tr, sprintf('dF/F — ROI #%d', rid), 'FontSize', 10);
                box(ax_tr, 'off');
                xlim(ax_tr, [t_dff(1), t_dff(end)]);
            else
                hCursor = [];
            end

            nDecided = sum(decision ~= 0);
            annotation(fig_r, 'textbox', [0.02 0.92 0.96 0.07], ...
                'String', sprintf('ROI#%d (%d/%d) [%s] | %d/%d decided | [K]eep [R]eject [U]ndo | left/right=ROI | Space=play/pause | Enter=save', ...
                    rid, cur, nSurv, decStr, nDecided, nSurv), ...
                'EdgeColor', 'none', 'FontSize', 10, 'FontWeight', 'bold', ...
                'Color', decColor, 'HorizontalAlignment', 'center');

            frame = 1;
            hSlider.Value = 1;
            lastFrame = [];
            needSetup = false;
            drawnow;
        end

        % ---- Play / pause logic ----
        isPlaying = (hPlay.Value == 1);
        if isPlaying
            hPlay.String = 'Pause';
            frame = mod(frame, T) + 1;
            hSlider.Value = frame;
        else
            hPlay.String = 'Play';
            frame = round(hSlider.Value);
        end

        % ---- Read montage frame + extract tile ----
        montageFrame = read(vr, frame);
        y0 = surv_y0(cur);
        x0 = surv_x0(cur);
        y1 = min(y0 + patchSz - 1, size(montageFrame, 1));
        x1 = min(x0 + patchSz - 1, size(montageFrame, 2));
        tile = montageFrame(y0:y1, x0:x1, :);
        set(hImg, 'CData', tile);

        % Update frame cursor on dFF trace
        if ~isempty(hCursor) && isvalid(hCursor)
            dff_idx = frame - frameOffset;
            if dff_idx >= 1 && dff_idx <= T_dff
                set(hCursor, 'Value', t_dff(dff_idx), 'Visible', 'on');
            else
                set(hCursor, 'Visible', 'off');
            end
        end

        hFrameLabel.String = sprintf('%d/%d  (%.1f s)', frame, T, (frame-1)/fps);
        drawnow limitrate;

        % ---- Check keypress ----
        key = fig_r.UserData;
        if ~isempty(key)
            fig_r.UserData = '';
            switch key
                case 'return'
                    done = true;
                case 'escape'
                    decision(:) = 0;
                    done = true;
                case 'space'
                    hPlay.Value = 1 - hPlay.Value;
                case 'k'
                    decision(cur) = 1;
                    if cur < nSurv, cur = cur + 1; end
                    needSetup = true;
                case 'r'
                    decision(cur) = -1;
                    if cur < nSurv, cur = cur + 1; end
                    needSetup = true;
                case 'u'
                    decision(cur) = 0;
                    needSetup = true;
                case 'rightarrow'
                    if cur < nSurv, cur = cur + 1; needSetup = true; end
                case 'leftarrow'
                    if cur > 1, cur = cur - 1; needSetup = true; end
            end
        end

        if isPlaying
            pause(1 / fps);
        else
            pause(0.05);
        end
    end

    if isvalid(fig_r), close(fig_r); end

    % Apply decisions
    for j = 1:nSurv
        si = survivors(j);
        if decision(j) == -1
            pass2_reject(si) = true;
            pass2_reason(si) = "manual";
        end
    end

    nRej = sum(decision == -1);
    fprintf('  Sequential review: %d rejected, %d kept, %d undecided (treated as keep)\n', ...
        nRej, sum(decision == 1), sum(decision == 0));
end
