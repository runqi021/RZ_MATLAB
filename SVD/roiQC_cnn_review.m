clear all; clc; close all;
%% roiQC_cnn_review.m
%  Generate PDF summary of CNN predictions: avg crop + dFF trace per ROI.
%  Green = CNN keep, Red = CNN reject. Sorted by probability (borderline first).
%
%  Outputs (in outDir):
%    roiQC_cnn_review.pdf — all FOVs, all ROIs, colored by CNN prediction
%
%  RZ 2026-04

%% ======================================================================
%  USER-EDITABLE PARAMETERS
%  ======================================================================
masterFolders = [
    "D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing"
    "D:\251124_live_vglut2_soma_g8s+cy5\phys\IO"
    "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing"
];
predictionsDir = "D:\roiQC_cnn_data";   % dir with predictions_train.mat + predictions_val.mat
outDir         = "D:\roiQC_cnn_data";   % where to save the PDF

patchSize_um   = 30;       % crop half-width (microns)
gamma_val      = 0.6;
clip_lo        = 0.4;
clip_hi        = 99.5;
maxPerPage     = 7;        % ROIs per page
scaleBar_um    = 5;        % scale bar length (microns)

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
%  LOAD CNN PREDICTIONS
%  ======================================================================
% Merge train + val predictions
pred_files = ["predictions_train.mat", "predictions_val.mat"];
all_prob     = [];
all_pred     = [];
all_fov_id   = [];
all_roi_label = [];

for pf = 1:numel(pred_files)
    fpath = fullfile(predictionsDir, pred_files(pf));
    if ~isfile(fpath), continue; end
    P = load(fpath);
    all_prob      = [all_prob;      P.prob(:)];        %#ok<AGROW>
    all_pred      = [all_pred;      P.pred_keep(:)];   %#ok<AGROW>
    all_fov_id    = [all_fov_id;    P.fov_id(:)];      %#ok<AGROW>
    all_roi_label = [all_roi_label; P.roi_label(:)];   %#ok<AGROW>
end
fprintf('Loaded %d CNN predictions.\n', numel(all_prob));

%% ======================================================================
%  DISCOVER ALL FOVs
%  ======================================================================
allExcl = [];
for mi = 1:numel(masterFolders)
    hits = dir(fullfile(masterFolders(mi), '**', 'roiQC_exclusion.mat'));
    allExcl = [allExcl; hits]; %#ok<AGROW>
end
nFOV = numel(allExcl);
fprintf('Found %d FOVs.\n', nFOV);

%% ======================================================================
%  GENERATE PDF PAGES
%  ======================================================================
tmpPDFs = {};

for ff = 1:nFOV
    fovFolder = allExcl(ff).folder;
    [~, fovID] = fileparts(fovFolder);
    fprintf('[%d/%d] %s ... ', ff, nFOV, fovID);

    try
        excl = load(fullfile(fovFolder, 'roiQC_exclusion.mat'));
        N_roi = excl.N_roi;
        validROI = excl.roiLabels(:)';

        % Match CNN predictions for this FOV
        fov_mask = (all_fov_id == ff);
        fov_prob = all_prob(fov_mask);
        fov_pred = all_pred(fov_mask);
        fov_roi  = all_roi_label(fov_mask);

        if isempty(fov_prob)
            fprintf('no predictions, skip\n');
            continue
        end

        % Load avg projection
        avgImg = load_avg(fovFolder);
        if isempty(avgImg), fprintf('no avg, skip\n'); continue; end
        [H, W] = size(avgImg);

        % Load maskL
        sam_hits = dir(fullfile(fovFolder, '*_cpSAM_output.mat'));
        if isempty(sam_hits), fprintf('no SAM, skip\n'); continue; end
        samData = load(fullfile(sam_hits(1).folder, sam_hits(1).name), 'maskL');
        maskL = samData.maskL;

        % Load dFF
        dff_hits = dir(fullfile(fovFolder, '*_dFF.mat'));
        dFF = [];
        fps = 30;
        if ~isempty(dff_hits)
            dffData = load(fullfile(dff_hits(1).folder, dff_hits(1).name), 'dFF');
            if isfield(dffData, 'dFF'), dFF = double(dffData.dFF); end
        end
        [fps_det, ~] = detect_session_fps(fovFolder);
        if fps_det > 0, fps = fps_det; end

        % Pixel size
        [~, scan_meta] = detect_session_fps(fovFolder);
        if isfield(scan_meta, 'pixelSize_um') && isfinite(scan_meta.pixelSize_um)
            pixelSize_um = scan_meta.pixelSize_um;
        else
            pixelSize_um = 1.7778;
        end
        crop_px = round(patchSize_um / pixelSize_um / 2);

        props = regionprops(maskL, 'Centroid', 'Area');

        % Sort ROIs by probability (borderline near middle, clear keeps/rejects at edges)
        [~, sortIdx] = sort(fov_prob, 'descend');  % keeps first, rejects last

        % Generate pages
        nPages = ceil(numel(fov_roi) / maxPerPage);
        for iPage = 1:nPages
            i0 = (iPage-1)*maxPerPage + 1;
            i1 = min(iPage*maxPerPage, numel(fov_roi));
            nThis = i1 - i0 + 1;

            pg_w = 28; pg_h = 3.5 * nThis + 2.5;
            fig = figure('Color', 'w', 'Visible', 'off', ...
                'Units', 'centimeters', 'Position', [1 1 pg_w pg_h], ...
                'PaperUnits', 'centimeters', 'PaperSize', [pg_w pg_h], ...
                'PaperPosition', [0 0 pg_w pg_h]);

            annotation(fig, 'textbox', [0.02 0.96 0.96 0.04], ...
                'String', sprintf('FOV#%02d  %s  (page %d/%d)  |  green=keep  red=reject', ...
                    ff, strrep(fovID, '_', ' '), iPage, nPages), ...
                'EdgeColor', 'none', 'FontSize', 10, 'Interpreter', 'none');

            ml = 0.8; mr = 0.5; mb = 0.8; mt = 1.5; vgap = 0.3;
            usable_w = pg_w - ml - mr;
            row_h = (pg_h - mb - mt - (nThis-1)*vgap) / nThis;
            crop_w = row_h;  % square crop
            trace_w = usable_w - crop_w - 0.4;

            for jj = 1:nThis
                si = sortIdx(i0 + jj - 1);
                rid = fov_roi(si);
                prob = fov_prob(si);
                pred = fov_pred(si);

                row_idx = jj - 1;
                row_y = pg_h - mt - (row_idx+1)*row_h - row_idx*vgap;

                % Color
                if pred == 1
                    traceColor = [0.2 0.7 0.2];  % green
                else
                    traceColor = [0.8 0.1 0.1];  % red
                end

                % Human label for comparison
                human_idx = find(validROI == rid, 1);
                if ~isempty(human_idx)
                    human_keep = excl.final_keep(human_idx);
                else
                    human_keep = NaN;
                end

                % ---- Crop panel ----
                ax_crop = axes(fig, 'Units', 'centimeters', ...
                    'Position', [ml row_y crop_w row_h]);
                if rid <= numel(props) && props(rid).Area > 0
                    cx = props(rid).Centroid(1);
                    cy = props(rid).Centroid(2);
                    % Full crop indices (may exceed image bounds)
                    r1 = round(cy) - crop_px;  r2 = round(cy) + crop_px;
                    c1 = round(cx) - crop_px;  c2 = round(cx) + crop_px;
                    cropSz = r2 - r1 + 1;

                    % Zero-padded crop (keeps square even at edges)
                    rp = zeros(cropSz, cropSz);
                    mc = zeros(cropSz, cropSz);
                    % Source region clipped to image bounds
                    sr1 = max(1, r1); sr2 = min(H, r2);
                    sc1 = max(1, c1); sc2 = min(W, c2);
                    % Destination offsets into padded patch
                    dr1 = sr1 - r1 + 1; dr2 = sr2 - r1 + 1;
                    dc1 = sc1 - c1 + 1; dc2 = sc2 - c1 + 1;
                    rp(dr1:dr2, dc1:dc2) = avgImg(sr1:sr2, sc1:sc2);
                    mc(dr1:dr2, dc1:dc2) = maskL(sr1:sr2, sc1:sc2);

                    lo = prctile(rp(:), clip_lo); hi = prctile(rp(:), clip_hi);
                    if hi <= lo, hi = lo + eps; end
                    rp = max(0, min(1, (rp - lo) / (hi - lo)));
                    rp = rp .^ gamma_val;

                    perim = bwperim(mc == rid);
                    rpRGB = repmat(rp, [1 1 3]);
                    for ch = 1:3
                        pl = rpRGB(:,:,ch);
                        pl(perim) = [1 1 0] * (ch == [1;2;3]);
                        rpRGB(:,:,ch) = pl;
                    end
                    imshow(rpRGB, 'Parent', ax_crop);

                    % 5 um scale bar (bottom-right corner)
                    sb_px = round(scaleBar_um / pixelSize_um);
                    hold(ax_crop, 'on');
                    sb_x0 = cropSz - sb_px - 3;
                    sb_y0 = cropSz - 4;
                    plot(ax_crop, [sb_x0 sb_x0+sb_px], [sb_y0 sb_y0], '-w', 'LineWidth', 2);
                    text(ax_crop, sb_x0 + sb_px/2, sb_y0 - 3, ...
                        sprintf('%d\\mum', scaleBar_um), ...
                        'Color', 'w', 'FontSize', 5, 'HorizontalAlignment', 'center');
                    hold(ax_crop, 'off');
                else
                    axis(ax_crop, 'off');
                end

                % Label with ROI#, probability, human agreement
                if isnan(human_keep)
                    agreeStr = '';
                elseif human_keep == pred
                    agreeStr = '';
                else
                    agreeStr = '  **DISAGREE**';
                end
                title(ax_crop, sprintf('ROI#%d  p=%.2f%s', rid, prob, agreeStr), ...
                    'FontSize', 7, 'Interpreter', 'none', 'Color', traceColor);

                % ---- Trace panel ----
                if ~isempty(dFF) && rid <= size(dFF, 2)
                    ax_tr = axes(fig, 'Units', 'centimeters', ...
                        'Position', [ml + crop_w + 0.4, row_y, trace_w, row_h]);
                    t_sec = (0:size(dFF,1)-1)' / fps;
                    plot(ax_tr, t_sec, dFF(:, rid), '-', 'Color', traceColor, 'LineWidth', 0.5);
                    xlim(ax_tr, [0 t_sec(end)]);
                    set(ax_tr, 'FontSize', 5);
                    box(ax_tr, 'off');
                    if jj == nThis
                        xlabel(ax_tr, 'Time (s)', 'FontSize', 6);
                    else
                        set(ax_tr, 'XTickLabel', []);
                    end
                    ylabel(ax_tr, 'dF/F', 'FontSize', 5);
                end
            end

            tmpPdf = fullfile(tempdir, sprintf('cnn_review_%03d_%02d.pdf', ff, iPage));
            exportgraphics(fig, tmpPdf, 'ContentType', 'vector', 'BackgroundColor', 'none');
            tmpPDFs{end+1} = tmpPdf; %#ok<SAGROW>
            close(fig);
        end

        fprintf('%d ROIs, %d pages\n', numel(fov_roi), nPages);

    catch ME
        fprintf('FAILED: %s\n', ME.message);
    end
end

%% ======================================================================
%  MERGE PDF
%  ======================================================================
if ~isempty(tmpPDFs)
    outPdf = fullfile(outDir, 'roiQC_cnn_review.pdf');
    if isfile(outPdf), delete(outPdf); end
    append_pdfs_RZ(outPdf, tmpPDFs{:});
    fprintf('\nSaved: %s\n', outPdf);

    for ii = 1:numel(tmpPDFs)
        if isfile(tmpPDFs{ii}), delete(tmpPDFs{ii}); end
    end
end

fprintf('Done.\n');


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function avgImg = load_avg(fovFolder)
    avgImg = [];
    avg_hits = dir(fullfile(fovFolder, '*_AVG_for_CP.tif'));
    if ~isempty(avg_hits)
        avgImg = double(imread(fullfile(avg_hits(1).folder, avg_hits(1).name)));
        return
    end
    sam_hits = dir(fullfile(fovFolder, '*_cpSAM_output.mat'));
    if ~isempty(sam_hits)
        S = load(fullfile(sam_hits(1).folder, sam_hits(1).name), 'avgPath');
        if isfield(S, 'avgPath') && isfile(char(S.avgPath))
            avgImg = double(imread(char(S.avgPath)));
        end
    end
end
