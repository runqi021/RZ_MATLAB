clear all; clc; close all;
%% roiQC_export_cellpose_training.m
%  Export corrected masks + avg projections for CellPose fine-tuning.
%
%  All images are rescaled to a common um/px so cell diameter is consistent
%  across FOVs with different zoom levels.
%
%  After export, fine-tune CellPose (in cellpose-gpu conda env):
%    python -m cellpose --train --dir <outDir> --pretrained_model cpsam --chan 0 --n_epochs 100 --min_train_masks 1 --diameter <target_diam_px>
%
%  The original cpsam model is NOT modified. The fine-tuned model is saved
%  separately and referenced by path.
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
outDir = "D:\cellpose_training";

target_umPerPx  = 0.45;    % rescale all images to this resolution (um/px)
                            % 0.45 um/px ≈ 4x zoom at PixelSizeBase=1.7778
cell_diameter_um = 15;      % expected soma diameter (microns)

%% ======================================================================
%  PATH SETUP
%  ======================================================================
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));

%% ======================================================================
%  DISCOVER ALL EXCLUSION FILES
%  ======================================================================
allExcl = [];
for mi = 1:numel(masterFolders)
    hits = dir(fullfile(masterFolders(mi), '**', 'roiQC_exclusion.mat'));
    allExcl = [allExcl; hits]; %#ok<AGROW>
end
nFOV = numel(allExcl);
fprintf('Found %d FOVs with exclusion data.\n', nFOV);
assert(nFOV > 0, 'No roiQC_exclusion.mat found.');

if isfolder(outDir)
    % Clear old exports
    delete(fullfile(outDir, 'img_*.tif'));
end
if ~isfolder(outDir), mkdir(char(outDir)); end

target_diam_px = round(cell_diameter_um / target_umPerPx);
fprintf('Target: %.2f um/px, cell diameter = %d um = %d px\n', ...
    target_umPerPx, cell_diameter_um, target_diam_px);

%% ======================================================================
%  EXPORT LOOP
%  ======================================================================
nExported = 0;
nSkipped  = 0;

for ff = 1:nFOV
    fovFolder = allExcl(ff).folder;
    [~, fovID] = fileparts(fovFolder);
    fprintf('[%d/%d] %s ... ', ff, nFOV, fovID);

    try
        excl = load(fullfile(fovFolder, 'roiQC_exclusion.mat'));

        % Detect pixel size
        [~, scan_meta] = detect_session_fps(fovFolder);
        if isfield(scan_meta, 'pixelSize_um') && isfinite(scan_meta.pixelSize_um)
            src_umPerPx = scan_meta.pixelSize_um;
        else
            src_umPerPx = excl.pixelSize_um;
        end
        scaleFactor = src_umPerPx / target_umPerPx;
        fprintf('%.3f um/px -> scale %.2fx ... ', src_umPerPx, scaleFactor);

        % Load maskL
        sam_hits = dir(fullfile(fovFolder, '*_cpSAM_output.mat'));
        if isempty(sam_hits), fprintf('no SAM, skip\n'); nSkipped = nSkipped+1; continue; end
        samData = load(fullfile(sam_hits(1).folder, sam_hits(1).name), 'maskL');
        maskL = samData.maskL;

        % Load avg projection
        avgImg = [];
        avg_hits = dir(fullfile(fovFolder, '*_AVG_for_CP.tif'));
        if ~isempty(avg_hits)
            avgImg = imread(fullfile(avg_hits(1).folder, avg_hits(1).name));
        elseif isfield(samData, 'avgPath') && isfile(char(samData.avgPath))
            avgImg = imread(char(samData.avgPath));
        end
        if isempty(avgImg), fprintf('no avg, skip\n'); nSkipped = nSkipped+1; continue; end

        % Zero out rejected ROIs
        rejectedLabels = excl.roiLabels(~excl.final_keep);
        maskL_corrected = maskL;
        for ri = 1:numel(rejectedLabels)
            maskL_corrected(maskL_corrected == rejectedLabels(ri)) = 0;
        end

        % Relabel to contiguous 1:N_kept
        keptLabels = excl.roiLabels(excl.final_keep);
        maskL_relabel = zeros(size(maskL_corrected), 'uint16');
        for ki = 1:numel(keptLabels)
            maskL_relabel(maskL_corrected == keptLabels(ki)) = uint16(ki);
        end

        % Rescale to target um/px
        if abs(scaleFactor - 1.0) > 0.01
            avgImg_rs     = imresize(avgImg,        scaleFactor, 'bilinear');
            maskL_rs      = imresize(maskL_relabel, scaleFactor, 'nearest');
        else
            avgImg_rs = avgImg;
            maskL_rs  = maskL_relabel;
        end

        % Save with CellPose naming
        tag = sprintf('img_%03d', ff);
        imwrite(avgImg_rs, fullfile(outDir, [tag '.tif']),       'tif', 'Compression', 'lzw');
        imwrite(maskL_rs,  fullfile(outDir, [tag '_masks.tif']), 'tif', 'Compression', 'lzw');

        nKept = sum(excl.final_keep);
        [H_rs, W_rs] = size(avgImg_rs);
        fprintf('exported (%d kept, %dx%d)\n', nKept, H_rs, W_rs);
        nExported = nExported + 1;

    catch ME
        fprintf('FAILED: %s\n', ME.message);
    end
end

%% ======================================================================
%  SUMMARY
%  ======================================================================
fprintf('\n==================== EXPORT COMPLETE ====================\n');
fprintf('  Exported: %d FOVs\n', nExported);
fprintf('  Skipped:  %d FOVs\n', nSkipped);
fprintf('  Output:   %s\n', outDir);
fprintf('  Target:   %.2f um/px, cell diameter = %d px\n', target_umPerPx, target_diam_px);
fprintf('\nTo fine-tune CellPose:\n');
fprintf('  conda activate cellpose-gpu\n');
fprintf('  python -m cellpose --train --dir "%s" --pretrained_model cpsam --chan 0 --n_epochs 100 --min_train_masks 1 --diameter %d\n', ...
    outDir, target_diam_px);
fprintf('=========================================================\n');
