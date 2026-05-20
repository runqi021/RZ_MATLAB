%% stackZ_dffQC_260331.m
%  Z-stack pipeline: split multi-Z TIFF -> per-slice pipeline (2-pass MC, cpSAM, dFF)
%  Uses same pipeline stages as Batch_dffQC_260325.m but for volumetric data.
%
%  Batch mode: scans subfolders for Z-stack TIFFs and processes each.
%
%  Output: <zstack_tif>_ZSTACK_PIPELINE/z0/, z5/, ... with per-slice outputs
%     per-slice: *_ch1.tif -> *_preproc.tif -> *_MC.tif -> *_cpSAM_output.mat -> *_dFF.mat
%     master:    *_ZSTACK_index.mat
%
%  Dependencies (same as Batch_dffQC_260325.m):
%    - loadtiff, saveastiff   (TIFFStack or mod/ copies)
%    - run_rigid_mc           (NoRMCorre-master/)
%    - cp_sam_extract_F_cli   (CellPose/)
%    - helper.dFF_RZ          (2p_breathing_coherence/+helper/)
%    - detect_session_fps     (optional; metadata fallback)

clear; close all; clc;

%% ========================= USER PARAMETERS =========================
% ---- Input (set ONE, leave the other "") ----
masterFolder    = "";       % scan recursively for Z-stack TIFFs
folderPath      = "E:\260124_chat_soma_g8s+cy5\7N_stack_2";  % single folder

% ---- Z-stack definition (MANUAL) ----
ZStart          = 0;        % first slice Z position (um)
ZStep           = 5;        % Z step between slices (um)
NumSlices       = 11;       % number of Z slices in the stack
FramesPerSlice  = 3000;     % imaging timepoints per slice (per channel)
SliceOrder      = "block";  % "block" (all frames z0, then z5, ...) | "interleaved"

% ---- Channel handling ----
AutoDetectChannels  = true;   % read ScanImage TIFF tags for saved channels
KeepChannelID       = 1;      % which channel to keep (1 = green / GCaMP)
NumChannelsFallback = 1;      % assumed channel count if no metadata found

% ---- Preprocessing ----
KeepFrames = Inf;   % after deinterleave (0 or Inf = keep all)
TrimL = 0;          % columns to trim from left
TrimR = 0;          % columns to trim from right

% ---- Motion Correction (2-pass rigid) ----
MaxShiftUm      = 10;        % max physical shift (um) -- auto-converted from zoom
PixelSizeBase   = 1.7778;   % um/pixel at 1x zoom
MCpasses        = 2;         % rigid MC passes (2 = template refinement)
MCinitBatch     = 500;       % frames for initial template
MCbinWidth      = 50;        % frames per batch (~1.7s at 30fps)
MCTossFrames    = 30;        % toss first N frames before MC

% ---- cpSAM ----
FPS               = 30;      % fallback; auto-detected from metadata
Diameter          = 30;      % expected cell diameter (pixels)
UseGPU            = true;
FlowThreshold     = 0.50;
CellprobThreshold = 0;

% ---- dFF ----
TossFrames      = 0;         % already tossed before MC via MCTossFrames
TraceGain       = 2.5;
MinOffset       = 1.0;
OffsetMult      = 12;
DffScaleBar     = 0.20;      % 20% dF/F scale bar
% ====================================================================

%% ========================= PATH SETUP =========================
scriptDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(scriptDir);  % RZ_MATLAB root
addpath(parentDir);

% Chronux
chronuxDir = fullfile(parentDir, 'chronux_2_12');
if isfolder(chronuxDir) && ~contains(path, 'chronux')
    addpath(genpath(chronuxDir));
end

% Helper package (+helper)
helperParent = fullfile(parentDir, '2p_breathing_coherence');
if isfolder(helperParent) && ~contains(path, '2p_breathing_coherence')
    addpath(helperParent);
end

% NoRMCorre
normcorreDir = fullfile(parentDir, 'NoRMCorre-master');
if isfolder(normcorreDir) && ~contains(path, 'NoRMCorre')
    addpath(genpath(normcorreDir));
end

% NoRMCorre modified
normcorreMod = fullfile(parentDir, 'NoRMCorre_modified_RZ_v1');
if isfolder(normcorreMod) && ~contains(path, 'NoRMCorre_modified')
    addpath(genpath(normcorreMod));
end

% minusDark
darkDir = fullfile(parentDir, 'minusDark');
if isfolder(darkDir) && ~contains(path, 'minusDark')
    addpath(genpath(darkDir));
end

% CellPose
cpDir = fullfile(parentDir, 'CellPose');
if isfolder(cpDir) && ~contains(path, 'CellPose')
    addpath(genpath(cpDir));
end

% mod/ (loadtiff, saveastiff, etc.)
modDir = fullfile(parentDir, 'mod');
if isfolder(modDir) && ~contains(path, [filesep 'mod'])
    addpath(genpath(modDir));
end

%% ========================= PYTHON =========================
cand = string([
    fullfile(getenv("USERPROFILE"), ".conda",     "envs", "cellpose-gpu", "python.exe")
    fullfile(getenv("USERPROFILE"), "anaconda3",  "envs", "cellpose-gpu", "python.exe")
    fullfile(getenv("USERPROFILE"), "miniconda3", "envs", "cellpose-gpu", "python.exe")
]);
pyExe = cand(find(isfile(cand), 1, "first"));
assert(pyExe ~= "", "Cellpose python not found. Checked:\n%s", strjoin(cand, newline));
fprintf("PythonExe: %s\n", pyExe);

%% ========================= DISCOVER Z-STACK TIFFS =========================
zstackList = discover_zstack_tifs(masterFolder, folderPath);
fprintf("\n=== Found %d Z-stack TIF files ===\n", numel(zstackList));
for i = 1:numel(zstackList)
    fprintf("  [%d] %s\n", i, zstackList(i));
end
fprintf("\n");

%% ========================= BATCH LOOP =========================
nProcessed = 0;
nFailed    = 0;
failedList = string.empty;

for ii = 1:numel(zstackList)
    zstackTif = zstackList(ii);
    fprintf("\n========== [%d/%d] %s ==========\n", ii, numel(zstackList), zstackTif);

    try
        idx = run_zstack_dffQC(zstackTif, pyExe, ...
            ZStart, ZStep, NumSlices, FramesPerSlice, SliceOrder, ...
            AutoDetectChannels, KeepChannelID, NumChannelsFallback, ...
            KeepFrames, TrimL, TrimR, ...
            MaxShiftUm, PixelSizeBase, MCpasses, MCinitBatch, MCbinWidth, MCTossFrames, ...
            FPS, Diameter, UseGPU, FlowThreshold, CellprobThreshold, ...
            TossFrames, TraceGain, MinOffset, OffsetMult, DffScaleBar);

        disp(idx);
        nProcessed = nProcessed + 1;

    catch ME
        fprintf(2, "\n*** FAILED: %s\n%s\n", zstackTif, ME.message);
        for kk = 1:numel(ME.stack)
            fprintf(2, "  in %s (line %d)\n", ME.stack(kk).name, ME.stack(kk).line);
        end
        nFailed = nFailed + 1;
        failedList(end+1) = zstackTif; %#ok<SAGROW>
    end
end

%% ========================= SUMMARY =========================
fprintf("\n==================== BATCH COMPLETE ====================\n");
fprintf("  Processed: %d\n", nProcessed);
fprintf("  Failed:    %d\n", nFailed);
fprintf("  Total:     %d\n", numel(zstackList));
if nFailed > 0
    fprintf("  Failed files:\n");
    for i = 1:numel(failedList)
        fprintf("    %s\n", failedList(i));
    end
end
fprintf("========================================================\n");

%% ======================================================================
%% ====================== Z-STACK TIFF DISCOVERY ========================
function tifList = discover_zstack_tifs(masterFolder, folderPath)
% Discover .tif files for Z-stack processing.
%   masterFolder: scan subfolders for TIFFs
%   folderPath:   scan a single folder
tifList = string.empty;

if masterFolder ~= ""
    masterFolder = string(masterFolder);
    assert(isfolder(masterFolder), "masterFolder not found: %s", masterFolder);
    subs = dir(masterFolder);
    subs = subs([subs.isdir] & ~startsWith({subs.name}, '.'));
    for i = 1:numel(subs)
        if subs(i).name == "." || subs(i).name == "..", continue; end
        subPath = fullfile(masterFolder, subs(i).name);
        tifs = dir(fullfile(subPath, '*.tif'));
        for j = 1:numel(tifs)
            if ~tifs(j).isdir
                tifList(end+1) = string(fullfile(tifs(j).folder, tifs(j).name)); %#ok<AGROW>
            end
        end
    end

elseif folderPath ~= ""
    folderPath = string(folderPath);
    assert(isfolder(folderPath), "folderPath not found: %s", folderPath);
    tifs = dir(fullfile(folderPath, '*.tif'));
    for j = 1:numel(tifs)
        if ~tifs(j).isdir
            tifList(end+1) = string(fullfile(tifs(j).folder, tifs(j).name)); %#ok<AGROW>
        end
    end
else
    error("Set either masterFolder or folderPath.");
end

tifList = unique(tifList, 'stable');
fprintf("[discover] Found %d TIF files\n", numel(tifList));
end

%% ======================================================================
%% ====================== MAIN Z-STACK PIPELINE =========================
function idx = run_zstack_dffQC(zstackTif, pyExe, ...
    ZStart, ZStep, NumSlices, FramesPerSlice, SliceOrder, ...
    AutoDetectChannels, KeepChannelID, NumChannelsFallback, ...
    KeepFrames, TrimL, TrimR, ...
    MaxShiftUm, PixelSizeBase, MCpasses, MCinitBatch, MCbinWidth, MCTossFrames, ...
    FPS, Diameter, UseGPU, FlowThreshold, CellprobThreshold, ...
    TossFrames, TraceGain, MinOffset, OffsetMult, DffScaleBar)

    zstackTif = string(zstackTif);
    assert(isfile(zstackTif), "Z-stack TIFF not found: %s", zstackTif);

    %% ---- Z list ----
    zList = ZStart + (0:NumSlices-1) * ZStep;

    %% ---- TIFF info + channel detection ----
    info = imfinfo(zstackTif);
    nPages = numel(info);
    fprintf("[TIFF] %d pages, %dx%d\n", nPages, info(1).Width, info(1).Height);

    savedCh = [];
    if AutoDetectChannels
        savedCh = scanimage_saved_channels_from_info(info(1));
    end
    if isempty(savedCh)
        savedCh = 1:NumChannelsFallback;
    end
    savedCh = savedCh(:).';
    nSaved = numel(savedCh);

    if mod(nPages, nSaved) ~= 0
        fprintf(2, "[chanDetect] WARNING: nPages=%d not divisible by nSaved=%d. Treating as single-channel.\n", nPages, nSaved);
        savedCh = KeepChannelID;
        nSaved  = 1;
    end

    kWithin = find(savedCh == KeepChannelID, 1, "first");
    keepID_used = KeepChannelID;
    if isempty(kWithin)
        kWithin = 1;
        keepID_used = savedCh(1);
        fprintf(2, "[chanDetect] KeepChannelID=%d not in savedCh=[%s]. Using %d.\n", ...
            KeepChannelID, num2str(savedCh), keepID_used);
    end
    fprintf("[chanDetect] savedCh=[%s], nSaved=%d, keep=%d (pos=%d)\n", ...
        num2str(savedCh), nSaved, keepID_used, kWithin);

    %% ---- Timepoint accounting ----
    assert(mod(nPages, nSaved) == 0, ...
        "TIFF %d pages not divisible by %d channels.", nPages, nSaved);
    nTP = nPages / nSaved;
    expectedTP = NumSlices * FramesPerSlice;
    assert(nTP == expectedTP, ...
        "Timepoint mismatch: TIFF has %d (%d pages / %d ch), expected %d*%d = %d.", ...
        nTP, nPages, nSaved, NumSlices, FramesPerSlice, expectedTP);

    %% ---- Output root ----
    [inDir, base] = fileparts(char(zstackTif));
    base = string(base);
    outRoot = fullfile(inDir, base + "_ZSTACK_PIPELINE");
    if ~isfolder(outRoot), mkdir(char(outRoot)); end

    %% ---- Index struct ----
    idx = struct();
    idx.zstackTif       = zstackTif;
    idx.outRoot         = string(outRoot);
    idx.zList_um        = zList(:);
    idx.nPages          = nPages;
    idx.nTP             = nTP;
    idx.numSlices       = NumSlices;
    idx.framesPerSlice  = FramesPerSlice;
    idx.sliceOrder      = string(SliceOrder);
    idx.savedChannels   = savedCh;
    idx.nSavedChannels  = nSaved;
    idx.keepChannelID   = keepID_used;
    idx.slices          = repmat(struct( ...
        "z_um",[], "subdir","", "slice_tif","", ...
        "preproc_tif","", "mc_tif","", "mc_mat","", ...
        "sam_mat","", "dff_mat","", "dff_png",""), NumSlices, 1);

    %% ---- Split + process each slice ----
    sliceOrder = lower(string(SliceOrder));
    tIn = Tiff(char(zstackTif), "r");

    for s = 1:NumSlices
        z = zList(s);
        subdir = fullfile(outRoot, sprintf("z%d", round(z)));
        if ~isfolder(subdir), mkdir(char(subdir)); end

        sliceTif = fullfile(subdir, sprintf("%s_z%d_ch%d.tif", base, round(z), keepID_used));

        %% ---- Compute timepoint indices ----
        if sliceOrder == "block"
            tStart = (s-1)*FramesPerSlice + 1;
            tEnd   = s*FramesPerSlice;
            tpIdx  = (tStart:tEnd)';
        elseif sliceOrder == "interleaved"
            tpIdx = (s:NumSlices:nTP)';
            assert(numel(tpIdx) == FramesPerSlice, ...
                "Interleaved indexing produced %d timepoints (expected %d).", ...
                numel(tpIdx), FramesPerSlice);
        else
            error('SliceOrder must be "block" or "interleaved".');
        end

        %% ---- Map to TIFF pages (channel deinterleave) ----
        pageIdx = (tpIdx - 1) * nSaved + kWithin;
        assert(all(pageIdx >= 1 & pageIdx <= nPages), ...
            "Computed pageIdx out of bounds for slice %d.", s);

        %% ---- Write per-slice TIFF (skip if exists) ----
        if ~isfile(sliceTif)
            fprintf("[slice %02d/%02d z=%d] Extracting %d frames...\n", ...
                s, NumSlices, round(z), numel(pageIdx));
            for k = 1:numel(pageIdx)
                tIn.setDirectory(pageIdx(k));
                fr = tIn.read();
                if ~isa(fr, "uint16"), fr = uint16(fr); end
                if k == 1
                    imwrite(fr, char(sliceTif), "tif", "Compression","lzw");
                else
                    imwrite(fr, char(sliceTif), "tif", "WriteMode","append", "Compression","lzw");
                end
            end
            fprintf("[slice %02d/%02d z=%d] Extracted -> %s\n", s, NumSlices, round(z), sliceTif);
        else
            fprintf("[slice %02d/%02d z=%d] Slice TIFF exists, skipping extraction.\n", ...
                s, NumSlices, round(z));
        end

        %% ---- Run dffQC pipeline on this slice ----
        fprintf("[slice %02d/%02d z=%d] Running pipeline...\n", s, NumSlices, round(z));

        sliceOut = run_slice_pipeline(sliceTif, pyExe, string(subdir), ...
            KeepFrames, TrimL, TrimR, ...
            MaxShiftUm, PixelSizeBase, MCpasses, MCinitBatch, MCbinWidth, MCTossFrames, ...
            FPS, Diameter, UseGPU, FlowThreshold, CellprobThreshold, ...
            TossFrames, TraceGain, MinOffset, OffsetMult, DffScaleBar);

        %% ---- Record in index ----
        idx.slices(s).z_um        = z;
        idx.slices(s).subdir      = string(subdir);
        idx.slices(s).slice_tif   = string(sliceTif);
        idx.slices(s).preproc_tif = sliceOut.preproc_tif;
        idx.slices(s).mc_tif      = sliceOut.mc_tif;
        idx.slices(s).mc_mat      = sliceOut.mc_mat;
        idx.slices(s).sam_mat     = sliceOut.sam_mat;
        idx.slices(s).dff_mat     = sliceOut.dff_mat;
        idx.slices(s).dff_png     = sliceOut.dff_png;

        fprintf("[slice %02d/%02d z=%d] DONE.\n\n", s, NumSlices, round(z));
    end

    tIn.close();

    %% ---- Save master index ----
    indexPath = fullfile(outRoot, base + "_ZSTACK_index.mat");
    save(char(indexPath), "idx", "-v7.3");
    fprintf("\nZ-stack index saved: %s\n", indexPath);
end

%% ======================================================================
%% ====================== PER-SLICE PIPELINE ============================
function out = run_slice_pipeline(sliceTif, pyExe, outFolder, ...
    KeepFrames, TrimL, TrimR, ...
    MaxShiftUm, PixelSizeBase, MCpasses, MCinitBatch, MCbinWidth, MCTossFrames, ...
    FPS, Diameter, UseGPU, FlowThreshold, CellprobThreshold, ...
    TossFrames, TraceGain, MinOffset, OffsetMult, DffScaleBar)
% Per-slice pipeline: preproc -> 2-pass MC -> cpSAM -> dFF
% Adapted from run_dffQC_pipeline in Batch_dffQC_260325.m.
% Input sliceTif is already single-channel (deinterleaved during Z-split).

    sliceTif  = string(sliceTif);
    outFolder = string(outFolder);
    [~, baseName] = fileparts(char(sliceTif));
    baseName = string(baseName);

    out = struct('preproc_tif',"", 'mc_tif',"", 'mc_mat',"", ...
                 'sam_mat',"", 'dff_mat',"", 'dff_png',"");

    %% ---- Output file names ----
    stem = baseName;  % already includes z and ch info

    out_preproc_tif = fullfile(outFolder, stem + "_preproc.tif");
    out_mc_mat      = fullfile(outFolder, stem + "_MCinfo.mat");
    out_mc_png      = fullfile(outFolder, stem + "_MCshift.png");
    out_mc_guess    = fullfile(outFolder, stem + "_preproc_MC.tif");
    out_sam_guess   = fullfile(outFolder, stem + "_preproc_MC_cpSAM_output.mat");
    out_dff_mat     = fullfile(outFolder, stem + "_dFF.mat");
    out_stack_png   = fullfile(outFolder, stem + "_stackDFF.png");
    out_meta_mat    = fullfile(outFolder, stem + "_meta.mat");

    %% ---- Metadata: zoom, FPS ----
    info = imfinfo(sliceTif);
    nPages = numel(info);
    fprintf("[META] sliceTif: %d pages, %dx%d\n", nPages, info(1).Width, info(1).Height);

    meta = "";
    if isfield(info(1), 'Software') && ~isempty(info(1).Software)
        meta = meta + string(info(1).Software) + newline;
    end
    if isfield(info(1), 'ImageDescription') && ~isempty(info(1).ImageDescription)
        meta = meta + string(info(1).ImageDescription) + newline;
    end

    % Zoom factor
    zoomFactor = parse_scanimage_scalar(meta, "SI.hRoiManager.scanZoomFactor");
    if isempty(zoomFactor) || ~isfinite(zoomFactor) || zoomFactor <= 0
        zoomFactor = [];
    end

    if ~isempty(zoomFactor)
        pixel_size_um = PixelSizeBase / zoomFactor;
        max_shift_px  = max(round(MaxShiftUm / pixel_size_um), 1);
        fprintf("[MC] zoom=%.2fx, pixel=%.4f um/px, MaxShift=%d um -> %d px\n", ...
            zoomFactor, pixel_size_um, MaxShiftUm, max_shift_px);
    else
        max_shift_px = 20;
        pixel_size_um = NaN;
        fprintf(2, "[MC] WARNING: could not read scanZoomFactor. Using default max_shift=%d px\n", max_shift_px);
    end

    % FPS
    scanFPS_raw = parse_scanimage_scalar(meta, "SI.hRoiManager.scanFrameRate");
    if ~isempty(scanFPS_raw) && isfinite(scanFPS_raw) && scanFPS_raw > 0
        detectedFPS = round(scanFPS_raw);
        fprintf("[META] scanFrameRate=%.4f -> FPS=%d\n", scanFPS_raw, detectedFPS);
    else
        detectedFPS = FPS;
        fprintf(2, "[META] WARNING: could not read scanFrameRate. Using default FPS=%d\n", detectedFPS);
    end

    % Build and save metadata struct
    scan_meta = struct();
    scan_meta.fps              = detectedFPS;
    scan_meta.scanFrameRate_raw = scanFPS_raw;
    scan_meta.zoomFactor       = zoomFactor;
    scan_meta.pixelSize_um     = pixel_size_um;
    scan_meta.max_shift_px     = max_shift_px;
    scan_meta.laserPower_pct      = parse_scanimage_vector(meta, "SI.hBeams.powers");
    scan_meta.channelSave         = parse_scanimage_vector(meta, "SI.hChannels.channelSave");
    scan_meta.channelInputRanges  = parse_scanimage_vector(meta, "SI.hChannels.channelInputRanges");
    scan_meta.motorPosition       = parse_scanimage_vector(meta, "SI.hMotors.motorPosition");
    scan_meta.source_tif       = sliceTif;
    save(char(out_meta_mat), "-struct", "scan_meta");
    fprintf("[META] Saved: %s\n", out_meta_mat);

    %% ---- STEP 0: PREPROC (trim + keep frames) ----
    % Slice TIF is already single-channel. Apply trim and frame truncation.
    if isfile(out_preproc_tif)
        fprintf("[SKIP] preproc step (found): %s\n", out_preproc_tif);
    else
        needsPreproc = (TrimL > 0) || (TrimR > 0) || ...
                        (isfinite(KeepFrames) && KeepFrames > 0 && KeepFrames < nPages);
        if needsPreproc
            fprintf("[RUN ] preproc step (trim L=%d R=%d, keepFrames=%g)\n", TrimL, TrimR, KeepFrames);
            F_raw = loadtiff(char(sliceTif));
            if TrimL > 0, F_raw(:,1:TrimL,:) = []; end
            if TrimR > 0, F_raw(:,end-TrimR+1:end,:) = []; end
            if isfinite(KeepFrames) && KeepFrames > 0 && size(F_raw,3) > KeepFrames
                F_raw = F_raw(:,:,1:KeepFrames);
                fprintf("[data] trimmed to KeepFrames=%d\n", KeepFrames);
            end
            F_raw = uint16(max(single(F_raw), 0));
            opts_sav = struct('overwrite', true, 'message', false);
            saveastiff(F_raw, char(out_preproc_tif), opts_sav);
            fprintf("[preproc] Saved: %s\n", out_preproc_tif);
        else
            % No trimming needed: symlink or copy the slice TIF as preproc
            % (use copyfile for Windows compatibility)
            copyfile(char(sliceTif), char(out_preproc_tif));
            fprintf("[preproc] No trim needed, copied slice TIF -> %s\n", out_preproc_tif);
        end
    end
    out.preproc_tif = out_preproc_tif;

    %% ---- STEP 1: MC (2-pass rigid) ----
    mc_tif = "";
    if isfile(out_mc_mat)
        SMC = load(char(out_mc_mat));
        if isfield(SMC,"mc_tif") && (isstring(SMC.mc_tif) || ischar(SMC.mc_tif))
            mc_tif = string(SMC.mc_tif);
            if ~isfile(mc_tif), mc_tif = ""; end
        end
    end
    if mc_tif == "" && isfile(out_mc_guess)
        mc_tif = out_mc_guess;
    end
    if mc_tif == ""
        mc_tif = newest_matching_file(outFolder, stem + "_preproc*MC*.tif");
        if mc_tif ~= "" && ~isfile(mc_tif), mc_tif = ""; end
    end

    if mc_tif ~= "" && isfile(mc_tif) && isfile(out_mc_mat)
        fprintf("[SKIP] MC step (found): %s\n", mc_tif);
    else
        nPasses = round(MCpasses);
        fprintf("[RUN ] MC step (%d pass, max_shift=%d px)\n", nPasses, max_shift_px);

        mc_input = out_preproc_tif;
        mcOut = [];
        for iPass = 1:nPasses
            fprintf("[MC pass %d/%d]\n", iPass, nPasses);
            if iPass == 1
                pass_max_shift = max_shift_px;
                mcOut = run_rigid_mc(mc_input, ...
                    'MaxShift', pass_max_shift, ...
                    'InitBatch', MCinitBatch, ...
                    'BinWidth', MCbinWidth, ...
                    'TossFrames', MCTossFrames);
            else
                pass_max_shift = max(round(max_shift_px / 4), 1);
                mcOut = run_rigid_mc(mc_input, ...
                    'MaxShift', pass_max_shift, ...
                    'InitBatch', MCinitBatch, ...
                    'BinWidth', MCbinWidth, ...
                    'Template', mcOut.template);
            end
            fprintf("[MC pass %d] max_shift = %d px\n", iPass, pass_max_shift);
            mc_input = string(mcOut.mc_path);
        end

        mc_tif = string(mcOut.mc_path);
        assert(mc_tif ~= "" && isfile(mc_tif), "Could not locate MC tif output.");

        % Save MC info first (protect expensive result)
        sx = []; sy = []; r = [];
        save(char(out_mc_mat), "mcOut","mc_tif","sx","sy","r","max_shift_px","-v7.3");
        fprintf("Saved MC info: %s\n", out_mc_mat);

        % Shift plot (try-catch so failure doesn't block pipeline)
        try
            if isfield(mcOut, "shifts")
                sh = mcOut.shifts;
                Tsh = numel(sh);
                sx = zeros(Tsh,1); sy = zeros(Tsh,1);
                for t = 1:Tsh
                    sft = sh(t).shifts;
                    sx(t) = sft(1); sy(t) = sft(2);
                end
                r = hypot(sx, sy);

                fig_mc = figure('Color','w','Visible','off');
                plot(r, 'LineWidth', 1.5); grid on;
                xlabel('Frame'); ylabel('|shift| (pixel)');
                title(sprintf('|Rigid shift| per frame (pass %d)', nPasses), 'Interpreter','none');
                exportgraphics(fig_mc, char(out_mc_png), "Resolution", 200);
                close(fig_mc);

                % Re-save with shift data
                save(char(out_mc_mat), "mcOut","mc_tif","sx","sy","r","max_shift_px","-v7.3");
            end
        catch ME_mc
            fprintf(2, "WARNING: shift plot failed (%s) -- MC results are saved.\n", ME_mc.message);
        end
    end
    out.mc_tif = mc_tif;
    out.mc_mat = out_mc_mat;

    %% ---- STEP 2: cpSAM (cell detection) ----
    sam_mat = "";
    if isfile(out_sam_guess)
        sam_mat = out_sam_guess;
    end
    if sam_mat == "" || ~isfile(sam_mat)
        sam_mat = newest_matching_file(outFolder, stem + "*cpSAM_output.mat");
        if sam_mat ~= "" && ~isfile(sam_mat), sam_mat = ""; end
    end

    if sam_mat ~= "" && isfile(sam_mat)
        fprintf("[SKIP] cpSAM step (found): %s\n", sam_mat);
    else
        fprintf("[RUN ] cpSAM step\n");
        samOut = cp_sam_extract_F_cli(mc_tif, ...
            'FPS', detectedFPS, ...
            'Diameter', Diameter, ...
            'PythonExe', char(pyExe), ...
            'UseGPU', UseGPU, ...
            'FlowThreshold', FlowThreshold, ...
            'CellprobThreshold', CellprobThreshold);

        sam_mat = infer_sam_mat(samOut, outFolder, stem);
        fprintf("SAM mat: %s\n", sam_mat);
    end
    out.sam_mat = sam_mat;

    %% ---- STEP 3: dFF + stack plot ----
    rerunDFF = true;
    if isfile(out_dff_mat) && isfile(out_stack_png)
        tmp = load(char(out_dff_mat), "params");
        if isfield(tmp, "params") && isfield(tmp.params, "fps") && tmp.params.fps == detectedFPS
            fprintf("[SKIP] dFF+stack step (found, fps=%d matches): %s\n", detectedFPS, out_dff_mat);
            rerunDFF = false;
        else
            fprintf("[RERUN] dFF fps mismatch (detected=%d). Recomputing.\n", detectedFPS);
        end
    end
    if rerunDFF
        fprintf("[RUN ] dFF+stack step (FPS=%d)\n", detectedFPS);

        S = load(char(sam_mat));
        assert(isfield(S,"F"), "SAM mat does not contain field F.");
        F_roi_raw = S.F;

        F_roi = F_roi_raw;
        if TossFrames > 0
            assert(size(F_roi,1) > TossFrames, "Not enough frames to toss for dFF.");
            F_roi(1:TossFrames,:) = [];
        end

        dFFout = dFF_RZ_dispatch(F_roi, detectedFPS);
        dFF = dFFout.dFF;

        params = struct();
        params.sam_mat     = sam_mat;
        params.fps         = detectedFPS;
        params.tossFrames  = TossFrames;
        params.traceGain   = TraceGain;
        params.minOffset   = MinOffset;
        params.offsetMult  = OffsetMult;
        params.dffScaleBar = DffScaleBar;

        save(char(out_dff_mat), "dFF","dFFout","F_roi_raw","F_roi","params","-v7.3");

        fig = figure("Color","w","Visible","off");
        stackDFF_with_scalebar(dFFout.dFF, detectedFPS, DffScaleBar, TraceGain, MinOffset, OffsetMult);
        title(sprintf("stackDFF | %s | fps=%d toss=%d", stem, detectedFPS, TossFrames), "Interpreter","none");
        exportgraphics(fig, char(out_stack_png), "Resolution", 200);
        close(fig);

        fprintf("Saved dFF: %s\n", out_dff_mat);
        fprintf("Saved stackDFF: %s\n", out_stack_png);
    end
    out.dff_mat = out_dff_mat;
    out.dff_png = out_stack_png;
end

%% ======================================================================
%% ===================== dFF: helper if exists, else local ===============
function dFFout = dFF_RZ_dispatch(F, fps)
if exist("helper.dFF_RZ","file") == 2
    dFFout = helper.dFF_RZ(F, 'FPS', fps);
else
    dFFout = dFF_RZ_local(F);
end
end

function out = dFF_RZ_local(F)
F = double(F);
F0 = prctile(F, 20, 1);
F0(F0 <= 0) = eps;
dFF = (F - F0) ./ F0;
out = struct();
out.F0  = F0;
out.dFF = dFF;
end

%% ======================================================================
%% ===================== stackDFF with scale bar =========================
function stackDFF_with_scalebar(dFF, fps, barAmp, gain, minOffset, offsetMult)
    [T, N] = size(dFF);
    t = (0:T-1) / fps;

    d = dFF;
    for i = 1:N
        d(:,i) = d(:,i) - median(d(:,i), "omitnan");
    end
    d = gain * d;

    s = mad(d(:), 1);
    if ~isfinite(s) || s == 0
        s = std(d(:), "omitnan");
    end
    if ~isfinite(s) || s == 0
        s = 0.05 * gain;
    end
    offset = max(minOffset, offsetMult * s);

    hold on;
    for i = 1:N
        y0 = (i-1)*offset;
        plot(t, d(:,i) + y0, "k", "LineWidth", 0.8);
        text(t(1), y0, sprintf("%d", i), ...
            "Color","k","FontSize",9, "FontName","Arial","FontWeight","normal", ...
            "HorizontalAlignment","left", "VerticalAlignment","bottom");
    end

    yticks((0:N-1)*offset);
    yticklabels(string(1:N));
    xlabel("Time (s)");
    ylabel("ROI");
    box off;

    barH = barAmp * gain;
    x0 = t(1) + 0.92*(t(end)-t(1));
    y0_bar = 0.10*offset;
    plot([x0 x0], [y0_bar y0_bar+barH], "k", "LineWidth", 2);
    text(x0, y0_bar+barH/2, sprintf("  %d%% dFF", round(barAmp*100)), ...
        "Color","k","FontSize",10, "FontName","Arial","FontWeight","normal", ...
        "HorizontalAlignment","left", "VerticalAlignment","middle");
end

%% ======================================================================
%% ====================== ScanImage channel parsing ======================
function savedCh = scanimage_saved_channels_from_info(info1)
meta = "";
if isfield(info1,"Software") && ~isempty(info1.Software)
    meta = meta + string(info1.Software) + newline;
end
if isfield(info1,"ImageDescription") && ~isempty(info1.ImageDescription)
    meta = meta + string(info1.ImageDescription) + newline;
end
meta = string(meta);

if strlength(strtrim(meta)) == 0
    savedCh = [];
    return
end

nAvail = parse_scanimage_scalar(meta, "SI.hChannels.channelsAvailable");
if isempty(nAvail) || ~isfinite(nAvail) || nAvail < 1
    nAvail = 4;
end
nAvail = min(max(round(nAvail),1), 16);

savedCh = parse_scanimage_channels(meta, "SI.hChannels.channelSave", nAvail);
if isempty(savedCh)
    savedCh = parse_scanimage_channels(meta, "SI.hChannels.channelsActive", nAvail);
end
end

function v = parse_scanimage_scalar(metaStr, key)
v = [];
pat = key + "\s*=\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)";
tok = regexp(metaStr, pat, "tokens", "once");
if isempty(tok), return; end
vv = str2double(string(tok{1}));
if isfinite(vv), v = vv; end
end

function ch = parse_scanimage_channels(metaStr, key, nAvail)
ch = [];

pat = key + "\s*=\s*([^\r\n]+)";
tok = regexp(metaStr, pat, "tokens", "once");
if isempty(tok), return; end
rhs = strtrim(string(tok{1}));

% A) scalar integer bitmask
if ~isempty(regexp(rhs, "^\d+$", "once"))
    m = str2double(rhs);
    if isfinite(m)
        ch = find(bitget(uint32(m), 1:nAvail));
        return;
    end
end

% B) explicit numeric list (e.g. [1 2])
nums = regexp(rhs, "[-+]?\d+\.?\d*", "match");
if ~isempty(nums)
    v = str2double(string(nums));
    v = v(isfinite(v));
    if ~isempty(v) && all(mod(v,1)==0) && all(v>=1) && all(v<=nAvail)
        ch = unique(v(:).', "stable");
        return;
    end
end

% C) logical list
tf = regexp(lower(rhs), "(true|false)", "match");
if ~isempty(tf)
    mask = strcmp(tf, "true");
    mask = mask(1:min(numel(mask), nAvail));
    ch = find(mask);
    return;
end
end

function vec = parse_scanimage_vector(metaStr, key)
vec = [];
pat = key + "\s*=\s*([^\r\n]+)";
tok = regexp(metaStr, pat, "tokens", "once");
if isempty(tok), return; end
rhs = strtrim(string(tok{1}));
nums = regexp(rhs, "[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", "match");
if ~isempty(nums)
    vec = str2double(string(nums));
    vec = vec(isfinite(vec));
end
end

%% ======================================================================
%% ====================== Helpers: locate outputs ========================
function sam_mat = infer_sam_mat(samOut, folderPath, stem)
sam_mat = "";
if isstring(samOut) || ischar(samOut)
    if endsWith(string(samOut), ".mat", "IgnoreCase", true)
        sam_mat = string(samOut);
    end
elseif isstruct(samOut)
    if isfield(samOut, "sam_mat"), sam_mat = string(samOut.sam_mat); end
    if sam_mat=="" && isfield(samOut, "out_mat"), sam_mat = string(samOut.out_mat); end
end
if sam_mat ~= "" && isfile(sam_mat), return; end

sam_mat = newest_matching_file(folderPath, stem + "*cpSAM_output.mat");
assert(sam_mat ~= "" && isfile(sam_mat), "Could not locate cpSAM output mat.");
end

function fpath = newest_matching_file(folderPath, pattern)
dd = dir(fullfile(folderPath, pattern));
if isempty(dd), fpath = ""; return; end
[~,ix] = max([dd.datenum]);
fpath = string(fullfile(dd(ix).folder, dd(ix).name));
end
