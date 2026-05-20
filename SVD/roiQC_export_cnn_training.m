clear all; clc; close all;
%% roiQC_export_cnn_training.m
%  Export labeled ROI data (dFF traces + avg patches + labels) to HDF5
%  for training the CNN QC classifier in Python/PyTorch.
%
%  Scans multiple experiment folders for roiQC_exclusion.mat files.
%  Splits by FOV (not by ROI) to prevent data leakage.
%
%  Outputs (in outDir):
%    train.h5  — training set
%    val.h5    — validation set
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
outDir       = "D:\roiQC_cnn_data";
patchSize_um = 30;       % crop half-width (microns)
patchResize  = 64;       % resize all patches to NxN for CNN
traceLen     = 2000;     % resample all traces to this length
valFraction  = 0.2;      % fraction of FOVs held out for validation
gamma_val    = 0.6;      % gamma correction (match Pass 1 display)
clip_lo      = 0.4;
clip_hi      = 99.5;

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

%% ======================================================================
%  TRAIN / VAL SPLIT (by FOV)
%  ======================================================================
rng(42);  % reproducible split
nVal   = max(1, round(nFOV * valFraction));
nTrain = nFOV - nVal;
perm   = randperm(nFOV);
trainIdx = sort(perm(1:nTrain));
valIdx   = sort(perm(nTrain+1:end));
fprintf('Split: %d train FOVs, %d val FOVs\n', nTrain, nVal);

%% ======================================================================
%  EXTRACT DATA
%  ======================================================================
trainData = extract_fov_data(allExcl, trainIdx, patchSize_um, patchResize, ...
    traceLen, gamma_val, clip_lo, clip_hi);
valData   = extract_fov_data(allExcl, valIdx, patchSize_um, patchResize, ...
    traceLen, gamma_val, clip_lo, clip_hi);

fprintf('\nTrain: %d ROIs (%d kept, %d rejected)\n', ...
    numel(trainData.labels), sum(trainData.labels), sum(~trainData.labels));
fprintf('Val:   %d ROIs (%d kept, %d rejected)\n', ...
    numel(valData.labels), sum(valData.labels), sum(~valData.labels));

%% ======================================================================
%  SAVE HDF5
%  ======================================================================
if ~isfolder(outDir), mkdir(char(outDir)); end

save_h5(fullfile(outDir, 'train.h5'), trainData);
save_h5(fullfile(outDir, 'val.h5'),   valData);

fprintf('\nSaved to: %s\n', outDir);
fprintf('Done.\n');


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

function data = extract_fov_data(allExcl, fovIdx, patchSize_um, patchResize, ...
        traceLen, gamma_val, clip_lo, clip_hi)
    all_traces  = {};
    all_patches = {};
    all_masks   = {};
    all_labels  = [];
    all_fov     = [];
    all_roi     = [];

    for fi = 1:numel(fovIdx)
        ii = fovIdx(fi);
        fovFolder = allExcl(ii).folder;
        [~, fovID] = fileparts(fovFolder);
        fprintf('  [%d/%d] %s ... ', fi, numel(fovIdx), fovID);

        try
            excl = load(fullfile(fovFolder, 'roiQC_exclusion.mat'));
            N_roi = excl.N_roi;
            validROI = excl.roiLabels(:)';

            % Load dFF
            dff_hits = dir(fullfile(fovFolder, '*_dFF.mat'));
            if isempty(dff_hits), fprintf('no dFF, skip\n'); continue; end
            dffData = load(fullfile(dff_hits(1).folder, dff_hits(1).name), 'dFF');
            if ~isfield(dffData, 'dFF'), fprintf('no dFF field, skip\n'); continue; end
            dFF = double(dffData.dFF);

            % Load avg projection
            avgImg = load_avg(fovFolder);
            if isempty(avgImg), fprintf('no avg, skip\n'); continue; end
            [H, W] = size(avgImg);

            % Load maskL
            sam_hits = dir(fullfile(fovFolder, '*_cpSAM_output.mat'));
            if isempty(sam_hits), fprintf('no SAM, skip\n'); continue; end
            samData = load(fullfile(sam_hits(1).folder, sam_hits(1).name), 'maskL');
            maskL = samData.maskL;

            % Pixel size
            [~, scan_meta] = detect_session_fps(fovFolder);
            if isfield(scan_meta, 'pixelSize_um') && isfinite(scan_meta.pixelSize_um)
                pixelSize_um = scan_meta.pixelSize_um;
            else
                pixelSize_um = 1.7778;
            end
            crop_px = round(patchSize_um / pixelSize_um / 2);

            props = regionprops(maskL, 'Centroid', 'Area');

            nAdded = 0;
            for i = 1:N_roi
                rid = validROI(i);
                if rid > size(dFF, 2), continue; end
                if rid > numel(props) || props(rid).Area == 0, continue; end

                % Resample trace
                trace_raw = dFF(:, rid);
                T_orig = numel(trace_raw);
                if T_orig < 10, continue; end
                t_orig = linspace(0, 1, T_orig);
                t_new  = linspace(0, 1, traceLen);
                trace_rs = interp1(t_orig, trace_raw, t_new, 'linear')';

                % Crop patch
                cx = props(rid).Centroid(1);
                cy = props(rid).Centroid(2);
                r1 = max(1, round(cy)-crop_px); r2 = min(H, round(cy)+crop_px);
                c1 = max(1, round(cx)-crop_px); c2 = min(W, round(cx)+crop_px);

                rp = avgImg(r1:r2, c1:c2);
                lo = prctile(rp(:), clip_lo); hi = prctile(rp(:), clip_hi);
                if hi <= lo, hi = lo + eps; end
                rp = max(0, min(1, (rp - lo) / (hi - lo)));
                rp = rp .^ gamma_val;
                rp = imresize(rp, [patchResize patchResize], 'bilinear');

                % Binary mask channel
                mc = maskL(r1:r2, c1:c2);
                bm = double(mc == rid);
                bm = imresize(bm, [patchResize patchResize], 'nearest');

                % Label
                label = excl.final_keep(i);

                all_traces{end+1}  = single(trace_rs);   %#ok<AGROW>
                all_patches{end+1} = single(rp);          %#ok<AGROW>
                all_masks{end+1}   = single(bm);          %#ok<AGROW>
                all_labels(end+1)  = int8(label);          %#ok<AGROW>
                all_fov(end+1)     = int32(ii);            %#ok<AGROW>
                all_roi(end+1)     = int32(rid);           %#ok<AGROW>
                nAdded = nAdded + 1;
            end
            fprintf('%d ROIs\n', nAdded);

        catch ME
            fprintf('FAILED: %s\n', ME.message);
        end
    end

    % Stack into arrays
    N = numel(all_labels);
    data.traces  = zeros(N, traceLen, 'single');
    data.patches = zeros(N, 2, patchResize, patchResize, 'single');
    data.labels  = all_labels(:);
    data.fov_id  = all_fov(:);
    data.roi_label = all_roi(:);

    for k = 1:N
        data.traces(k, :)    = all_traces{k};
        data.patches(k, 1, :, :) = all_patches{k};
        data.patches(k, 2, :, :) = all_masks{k};
    end
end


function save_h5(filepath, data)
    filepath = char(filepath);
    if isfile(filepath), delete(filepath); end

    [N, traceLen] = size(data.traces);
    if N == 0
        fprintf('  WARNING: 0 ROIs, skipping %s\n', filepath);
        return
    end
    patchSz = size(data.patches, 3);

    h5create(filepath, '/traces',    [N traceLen],    'Datatype', 'single');
    h5create(filepath, '/patches',   [N 2 patchSz patchSz], 'Datatype', 'single');
    h5create(filepath, '/labels',    [N 1],           'Datatype', 'int8');
    h5create(filepath, '/fov_id',    [N 1],           'Datatype', 'int32');
    h5create(filepath, '/roi_label', [N 1],           'Datatype', 'int32');

    h5write(filepath, '/traces',    data.traces);
    h5write(filepath, '/patches',   data.patches);
    h5write(filepath, '/labels',    data.labels(:));
    h5write(filepath, '/fov_id',    data.fov_id(:));
    h5write(filepath, '/roi_label', data.roi_label(:));

    fprintf('  Saved: %s  (%d ROIs)\n', filepath, N);
end


function avgImg = load_avg(fovFolder)
    avgImg = [];
    % Try AVG_for_CP.tif
    avg_hits = dir(fullfile(fovFolder, '*_AVG_for_CP.tif'));
    if ~isempty(avg_hits)
        avgImg = double(imread(fullfile(avg_hits(1).folder, avg_hits(1).name)));
        return
    end
    % Try cpSAM avgPath
    sam_hits = dir(fullfile(fovFolder, '*_cpSAM_output.mat'));
    if ~isempty(sam_hits)
        S = load(fullfile(sam_hits(1).folder, sam_hits(1).name), 'avgPath');
        if isfield(S, 'avgPath') && isfile(char(S.avgPath))
            avgImg = double(imread(char(S.avgPath)));
        end
    end
end
