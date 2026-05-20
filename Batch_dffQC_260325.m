%% Batch_dffQC_260325.m
% Batch dffQC pipeline: 2-pass rigid MC -> cpSAM -> dFF
% No background subtraction. Auto zoom-based max_shift.
%
% Two input modes:
%   masterFolder: scan multiple experiments (each with phys/ or phys/processed/)
%   folderPath:   scan a single folder for TIFs (like old batch)
% Set one, leave the other "".
%
% Outputs go into run<YYMMDD>/ subfolder under each recording folder,
% so old analysis results are never overwritten.
%
% Requires on path:
%   loadtiff, saveastiff, run_rigid_mc, cp_sam_extract_F_cli
% Optional: 
%   helper.dFF_RZ
clear; clc; delete(gcp('nocreate'));

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(fullfile(repoRoot, 'NoRMCorre-master'));
addpath(fullfile(repoRoot, 'CellPose'));
addpath(fullfile(repoRoot, 'SVD'));

%% ============================== USER ==============================
masterFolder = "";   % multi-experiment master dir
%folderPath   = "C:\Users\Admin\Desktop\260330_sst_soma_g8s\phys\IO";   % single folder (ignored if masterFolder set)
folderPath = "C:\Users\Admin\Desktop\260517_shi_cal590\phys";
% Output goes directly into the data folder alongside the TIFF

% -------- re-run options --------
RerunMC       = false;   % true = delete existing MC + cpSAM + dFF and re-run from MC
RerunCellpose = false;    % true = delete existing cpSAM + dFF and re-run from Cellpose
RerunStackDFF = false;   % true = regenerate stackDFF plot only (no dFF recompute)
RerunMontageVideo = false; % true = delete existing per-ROI videos and regenerate all
MakeMontageVideo = false; % true = generate per-ROI montage video + avg proj label image

% -------- channel handling --------
AutoDetectChannels  = true;
KeepChannelID       = 2;
NumChannelsFallback = 1;

% -------- frame selection --------
KeepFrames = Inf;   % after deinterleave (0 or Inf = keep all)

% -------- trim (columns) --------
TrimL = 0;
TrimR = 0;

% -------- MC (2-pass rigid) --------
MaxShiftUm     = 10;      % max physical shift (um) — auto-converted from zoom
PixelSizeBase  = 1.7778;  % um/px at 1x zoom
MCpasses       = 1;       % rigid MC passes (2 = template refinement)
MCinitBatch    = 500;     % frames for initial template
MCbinWidth     = 50;      % frames per batch (~1.7s at 30fps)
MCTossFrames   = 30;      % toss first N frames before MC

% -------- cpSAM --------
FPS               = 30;   % fallback if TIFF metadata missing; auto-detected per file
Diameter          = [];   % [] = Cellpose auto-estimates; or set numeric pixels (e.g. 30)
UseGPU            = true;
FlowThreshold     = 0.4;
CellprobThreshold = -0.5;

% -------- dFF --------
BaselineWinSec = 15;
TossFrames  = 0;       % already tossed before MC via MCTossFrames
TraceGain   = 2.5;
MinOffset   = 1.0;
OffsetMult  = 12;
DffScaleBar = 0.20;    % 20%

% -------- python exe --------
cand = string([
    fullfile(getenv("USERPROFILE"), ".conda",     "envs", "cellpose-gpu", "python.exe")
    fullfile(getenv("USERPROFILE"), "anaconda3",  "envs", "cellpose-gpu", "python.exe")
    fullfile(getenv("USERPROFILE"), "miniconda3", "envs", "cellpose-gpu", "python.exe")
]);
pyExe = cand(find(isfile(cand), 1, "first"));
assert(pyExe ~= "", "Cellpose python not found. Checked:\n%s", strjoin(cand, newline));
fprintf("Using PythonExe:\n%s\n", pyExe);

%% ============================== DISCOVER TIFS ==============================
tifList = discover_tifs(masterFolder, folderPath);
fprintf("\n=== Found %d TIF files to process ===\n", numel(tifList));
for i = 1:numel(tifList)
    fprintf("  [%d] %s\n", i, tifList(i));
end
fprintf("\n");

%% ============================== BATCH LOOP ==============================
nProcessed = 0;
nSkipped   = 0;
nFailed    = 0;
failedList = string.empty;

for ii = 1:numel(tifList)
    tifPath = tifList(ii);
    [recFolder, ~, ~] = fileparts(char(tifPath));
    outFolder = string(recFolder);

    fprintf("\n========== [%d/%d] %s ==========\n", ii, numel(tifList), tifPath);
    fprintf("[output] %s\n", outFolder);

    try
        out = run_dffQC_pipeline(tifPath, pyExe, outFolder, ...
            "AutoDetectChannels",AutoDetectChannels, ...
            "KeepChannelID",KeepChannelID, ...
            "NumChannelsFallback",NumChannelsFallback, ...
            "KeepFrames",KeepFrames, ...
            "TrimL",TrimL, "TrimR",TrimR, ...
            "MaxShiftUm",MaxShiftUm, "PixelSizeBase",PixelSizeBase, "MCpasses",MCpasses, ...
            "MCinitBatch",MCinitBatch, "MCbinWidth",MCbinWidth, "MCTossFrames",MCTossFrames, ...
            "FPS",FPS, "Diameter",Diameter, ...
            "UseGPU",UseGPU, "FlowThreshold",FlowThreshold, "CellprobThreshold",CellprobThreshold, ...
            "RerunMC",RerunMC, "RerunCellpose",RerunCellpose, "RerunStackDFF",RerunStackDFF, ...
            "RerunMontageVideo",RerunMontageVideo, "MakeMontageVideo",MakeMontageVideo, ...
            "TossFrames",TossFrames, "BaselineWinSec",BaselineWinSec, ...
            "TraceGain",TraceGain, "MinOffset",MinOffset, "OffsetMult",OffsetMult, ...
            "DffScaleBar",DffScaleBar);

        disp(out);
        nProcessed = nProcessed + 1;

    catch ME
        fprintf(2, "\n*** FAILED [%d/%d]: %s\n%s\n\n", ii, numel(tifList), tifPath, ME.message);
        nFailed = nFailed + 1;
        failedList(end+1) = tifPath; %#ok<SAGROW>
    end
end

%% ============================== SUMMARY ==============================
fprintf("\n==================== BATCH COMPLETE ====================\n");
fprintf("  Processed: %d\n", nProcessed);
fprintf("  Failed:    %d\n", nFailed);
fprintf("  Total:     %d\n", numel(tifList));
if nFailed > 0
    fprintf("  Failed files:\n");
    for i = 1:numel(failedList)
        fprintf("    %s\n", failedList(i));
    end
end
fprintf("========================================================\n");

%% ======================================================================
%% ====================== TIF DISCOVERY ==================================
function tifList = discover_tifs(masterFolder, folderPath)
% Discover .tif files to process.
%   masterFolder: scan multiple experiments (each with phys/ somewhere inside)
%   folderPath:   scan a single folder — smart: auto-dives into phys/processed/ if found
%
% Works at any level:
%   "D:\BATCH_DFFQC_TEST_260325"              -> finds all experiments' phys/ dirs
%   "D:\260322_sst_soma_g8s"                   -> finds phys/ inside, scans there
%   "D:\260322_sst_soma_g8s\phys"              -> finds processed/ inside, scans there
%   "D:\...\phys\processed\maybe breathing"    -> scans directly for recordings

tifList = string.empty;

if masterFolder ~= ""
    % --- masterFolder mode: collect all root dirs, then scan each ---
    masterFolder = string(masterFolder);
    assert(isfolder(masterFolder), "masterFolder not found: %s", masterFolder);

    % find all 'phys' directories recursively under master
    allDirs = dir(fullfile(masterFolder, '**', 'phys'));
    physDirs = string.empty;
    for i = 1:numel(allDirs)
        if allDirs(i).isdir
            physDirs(end+1) = string(fullfile(allDirs(i).folder, allDirs(i).name)); %#ok<AGROW>
        end
    end

    % resolve each phys dir to its scan root
    scanDirs = string.empty;
    for i = 1:numel(physDirs)
        scanDirs = [scanDirs, resolve_scan_dirs(physDirs(i))]; %#ok<AGROW>
    end

    for i = 1:numel(scanDirs)
        tifList = [tifList, find_tifs_recursive(scanDirs(i))]; %#ok<AGROW>
    end

elseif folderPath ~= ""
    % --- folderPath mode: smart — auto-dive into phys/processed/ ---
    folderPath = string(folderPath);
    assert(isfolder(folderPath), "folderPath not found: %s", folderPath);

    scanDirs = resolve_scan_dirs(folderPath);
    for i = 1:numel(scanDirs)
        tifList = [tifList, find_tifs_recursive(scanDirs(i))]; %#ok<AGROW>
    end
else
    error("Set either masterFolder or folderPath.");
end

tifList = unique(tifList, 'stable');
fprintf("[discover] Found %d TIF files\n", numel(tifList));
end

function scanDirs = resolve_scan_dirs(rootDir)
% Given a directory, figure out where recordings actually live.
%   - If phys/ exists underneath -> dive in (then check for processed/)
%   - If processed/ exists underneath -> dive in
%   - Otherwise -> scan this folder directly
scanDirs = string.empty;

physDir = fullfile(rootDir, "phys");
if isfolder(physDir)
    % found phys/, now check for processed/
    procDir = fullfile(physDir, "processed");
    if isfolder(procDir)
        scanDirs(end+1) = procDir;
    else
        scanDirs(end+1) = physDir;
    end
    return
end

procDir = fullfile(rootDir, "processed");
if isfolder(procDir)
    scanDirs(end+1) = procDir;
    return
end

% no phys/ or processed/ found — scan this folder directly
scanDirs(end+1) = rootDir;
end

function tifList = find_tifs_recursive(scanDir)
% Recursively find .tif files under scanDir.
% Filter: keep only recording_name/recording_name.tif (filename matches parent folder).
% Also handle loose .tif files (move into own subfolder).
tifList = string.empty;

dd = dir(fullfile(scanDir, '**', '*.tif'));
for j = 1:numel(dd)
    if dd(j).isdir, continue; end

    % use split instead of fileparts to avoid dots in folder names
    % (e.g., "4.2x_..." would be split into name="4" ext=".2x_..." by fileparts)
    parts = split(string(dd(j).folder), filesep);
    parName = parts(end);
    tifBase = erase(string(dd(j).name), ".tif");

    if string(tifBase) == string(parName)
        % matches convention: recording_name/recording_name.tif
        tifList(end+1) = string(fullfile(dd(j).folder, dd(j).name)); %#ok<AGROW>
    elseif string(dd(j).folder) == string(scanDir)
        % loose .tif directly in scan dir — move into own subfolder
        subDir = fullfile(scanDir, tifBase);
        if ~isfolder(subDir), mkdir(char(subDir)); end
        src = fullfile(dd(j).folder, dd(j).name);
        dst = fullfile(subDir, dd(j).name);
        if ~isfile(string(dst))
            movefile(char(src), char(dst));
        end
        tifList(end+1) = string(dst); %#ok<AGROW>
    end
    % else: .tif in a subfolder but name doesn't match parent -> skip (not a recording)
end
end

%% ======================================================================
%% ====================== PER-FILE PIPELINE ==============================
function out = run_dffQC_pipeline(tifPath, pyExe, outFolder, varargin)

p = inputParser;
p.addRequired("tifPath", @(s)isstring(s)||ischar(s));
p.addRequired("pyExe", @(s)isstring(s)||ischar(s));
p.addRequired("outFolder", @(s)isstring(s)||ischar(s));

% channel
p.addParameter("AutoDetectChannels", true, @(x)islogical(x)&&isscalar(x));
p.addParameter("KeepChannelID", 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter("NumChannelsFallback", 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% frame selection
p.addParameter("KeepFrames", 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% trim
p.addParameter("TrimL", 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("TrimR", 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% MC
p.addParameter("MaxShiftUm", 10, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("PixelSizeBase", 1.7778, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("MCpasses", 2, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter("MCinitBatch", 500, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("MCbinWidth", 50, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("MCTossFrames", 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% cpSAM
p.addParameter("FPS", 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("Diameter", [], @(x)isempty(x) || (isnumeric(x)&&isscalar(x)&&x>0));
p.addParameter("UseGPU", true, @(x)islogical(x)&&isscalar(x));
p.addParameter("FlowThreshold", 0.50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter("CellprobThreshold", 0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter("RerunMC", false, @(x)islogical(x)&&isscalar(x));
p.addParameter("RerunCellpose", false, @(x)islogical(x)&&isscalar(x));
p.addParameter("RerunStackDFF", false, @(x)islogical(x)&&isscalar(x));
p.addParameter("RerunMontageVideo", false, @(x)islogical(x)&&isscalar(x));
p.addParameter("MakeMontageVideo", true, @(x)islogical(x)&&isscalar(x));

% dFF + plot
p.addParameter("TossFrames", 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("BaselineWinSec", 15, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("TraceGain", 2.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("MinOffset", 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("OffsetMult", 12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("DffScaleBar", 0.20, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.parse(tifPath, pyExe, outFolder, varargin{:});
pr = p.Results;

tifPath   = string(pr.tifPath);
pyExe     = string(pr.pyExe);
outFolder = string(pr.outFolder);
assert(isfile(tifPath), "TIFF not found: %s", tifPath);
assert(isfile(pyExe), "PythonExe not found: %s", pyExe);

% create output subfolder
if ~isfolder(outFolder), mkdir(char(outFolder)); end
fprintf("[setup] Output folder: %s\n", outFolder);

[~, baseName, ~] = fileparts(char(tifPath));
baseName = string(baseName);

% --- detect channels from TIFF tags ---
info = imfinfo(tifPath);
nPages = numel(info);

savedCh = [];
if pr.AutoDetectChannels
    savedCh = scanimage_saved_channels_from_info(info(1));
end
if isempty(savedCh)
    savedCh = 1:pr.NumChannelsFallback;
end
savedCh = savedCh(:).';
nSaved = numel(savedCh);

if mod(nPages, nSaved) ~= 0
    fprintf(2, "[chanDetect] WARNING: nPages=%d not divisible by nSaved=%d. Treating as single-channel.\n", nPages, nSaved);
    savedCh = pr.KeepChannelID;
    nSaved  = 1;
end

kWithin = find(savedCh == pr.KeepChannelID, 1, "first");
keepID_used = pr.KeepChannelID;
if isempty(kWithin)
    kWithin = 1;
    keepID_used = savedCh(1);
    fprintf(2, "[chanDetect] KeepChannelID=%d not in savedCh=[%s]. Using %d.\n", pr.KeepChannelID, num2str(savedCh), keepID_used);
end

fprintf("[chanDetect] savedCh=[%s] (nSaved=%d), keepID=%d (kWithin=%d)\n", ...
    num2str(savedCh), nSaved, keepID_used, kWithin);

% --- auto max_shift from zoom factor ---
meta = "";
if isfield(info(1),"Software") && ~isempty(info(1).Software)
    meta = meta + string(info(1).Software) + newline;
end
if isfield(info(1),"ImageDescription") && ~isempty(info(1).ImageDescription)
    meta = meta + string(info(1).ImageDescription) + newline;
end

zoomFactor = parse_scanimage_scalar(meta, "SI.hRoiManager.scanZoomFactor");
if ~isempty(zoomFactor) && isfinite(zoomFactor) && zoomFactor > 0
    pixel_size_um = pr.PixelSizeBase / zoomFactor;
    max_shift_px  = round(pr.MaxShiftUm / pixel_size_um);
    max_shift_px  = max(max_shift_px, 1);
    fprintf("[MC] zoom=%.2fx, pixel=%.4f um/px, MaxShift=%d um -> %d px\n", ...
        zoomFactor, pixel_size_um, pr.MaxShiftUm, max_shift_px);
else
    max_shift_px = 20;
    fprintf(2, "[MC] WARNING: could not read scanZoomFactor. Using default max_shift=%d px\n", max_shift_px);
end

% --- auto-detect frame rate from TIFF metadata ---
scanFPS_raw = parse_scanimage_scalar(meta, "SI.hRoiManager.scanFrameRate");
if ~isempty(scanFPS_raw) && isfinite(scanFPS_raw) && scanFPS_raw > 0
    detectedFPS = round(scanFPS_raw);
    fprintf("[META] scanFrameRate=%.4f -> FPS=%d\n", scanFPS_raw, detectedFPS);
else
    detectedFPS = pr.FPS;
    fprintf(2, "[META] WARNING: could not read scanFrameRate. Using default FPS=%d\n", detectedFPS);
end

% --- build imaging metadata struct ---
scan_meta = struct();
scan_meta.fps = detectedFPS;
scan_meta.scanFrameRate_raw = scanFPS_raw;
scan_meta.zoomFactor = zoomFactor;
if ~isempty(zoomFactor) && isfinite(zoomFactor) && zoomFactor > 0
    scan_meta.pixelSize_um = pr.PixelSizeBase / zoomFactor;
else
    scan_meta.pixelSize_um = NaN;
end
scan_meta.laserPower_pct      = parse_scanimage_vector(meta, "SI.hBeams.powers");
scan_meta.channelSave         = parse_scanimage_vector(meta, "SI.hChannels.channelSave");
scan_meta.channelInputRanges  = parse_scanimage_vector(meta, "SI.hChannels.channelInputRanges");
scan_meta.motorPosition       = parse_scanimage_vector(meta, "SI.hMotors.motorPosition");
scan_meta.framesPerSlice      = parse_scanimage_scalar(meta, "SI.hStackManager.framesPerSlice");
scan_meta.numSlices           = parse_scanimage_scalar(meta, "SI.hStackManager.numSlices");
scan_meta.source_tif          = tifPath;

% --- define output stem ---
tag = "ch" + string(keepID_used);
if contains(lower(baseName), lower(tag))
    stem = baseName;
else
    stem = baseName + "_" + tag;
end

% output paths — all in outFolder (the run<YYMMDD> subfolder)
out_preproc_tif = outFolder + filesep + stem + "_preproc.tif";

out_mc_mat   = outFolder + filesep + stem + "_MCinfo.mat";
out_mc_png   = outFolder + filesep + stem + "_MCshift.png";
out_mc_guess = outFolder + filesep + stem + "_preproc_MC.tif";

out_sam_guess = outFolder + filesep + stem + "_preproc_MC_cpSAM_output.mat";
out_dff_mat   = outFolder + filesep + stem + "_dFF.mat";
out_stack_png = outFolder + filesep + stem + "_stackDFF.png";
out_meta_mat  = outFolder + filesep + stem + "_meta.mat";

% Save _meta.mat (always overwrite — small file)
save(char(out_meta_mat), "-struct", "scan_meta");
fprintf("[META] Saved: %s\n", out_meta_mat);

% --- STEP 1: LOAD + DEINTERLEAVE + KEEP FRAMES ---
if isfile(out_preproc_tif)
    fprintf("[SKIP] preproc step (found): %s\n", out_preproc_tif);
else
    fprintf("[RUN ] preproc step (load + deinterleave + keep frames)\n");

    F_raw = loadtiff(tifPath);
    F = F_raw;
    if pr.TrimL > 0, F(:,1:pr.TrimL,:) = []; end
    if pr.TrimR > 0, F(:,end-pr.TrimR+1:end,:) = []; end

    if nSaved > 1
        assert(mod(size(F,3), nSaved) == 0, "Pages=%d not divisible by nSaved=%d.", size(F,3), nSaved);
        F_keep = F(:,:,kWithin:nSaved:end);
    else
        F_keep = F;
    end
    fprintf("[data] raw pages=%d, after deinterleave=%d\n", size(F,3), size(F_keep,3));

    if pr.KeepFrames > 0 && isfinite(pr.KeepFrames) && size(F_keep,3) > pr.KeepFrames
        F_keep = F_keep(:,:,1:pr.KeepFrames);
        fprintf("[data] trimmed to KeepFrames=%d\n", pr.KeepFrames);
    else
        fprintf("[data] keeping all %d frames\n", size(F_keep,3));
    end

    F_keep = uint16(max(single(F_keep), 0));

    out_preproc_tif_c = char(out_preproc_tif);
    if exist(out_preproc_tif_c, "file"), delete(out_preproc_tif_c); end
    opts = struct('overwrite', true, 'message', false);
    saveastiff(F_keep, out_preproc_tif_c, opts);

    fprintf("Saved preproc TIFF:\n%s\n", out_preproc_tif);
end

% --- STEP 2: MC (2-pass rigid) ---
if pr.RerunMC
    % Delete existing MC + cpSAM + dFF so they get regenerated
    delPatterns = ["*_MC*.tif", "*_MCinfo.mat", "*_MCshift.png", ...
                   "*cpSAM_output.mat", "*_dFF.mat", "*_stackDFF.png"];
    for ip = 1:numel(delPatterns)
        hits = dir(fullfile(outFolder, delPatterns(ip)));
        for ih = 1:numel(hits)
            delPath = fullfile(hits(ih).folder, hits(ih).name);
            fprintf("[RERUN-MC] Deleting: %s\n", delPath);
            delete(delPath);
        end
    end
end
mc_tif = "";
if isfile(out_mc_mat)
    SMC = load(out_mc_mat);
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
    nPasses = round(pr.MCpasses);
    fprintf("[RUN ] MC step (%d pass, max_shift=%d px)\n", nPasses, max_shift_px);

    mc_input = out_preproc_tif;
    mcOut = [];
    for iPass = 1:nPasses
        fprintf("[MC pass %d/%d]\n", iPass, nPasses);
        if iPass == 1
            pass_max_shift = max_shift_px;
            mcOut = run_rigid_mc(mc_input, ...
                'MaxShift', pass_max_shift, ...
                'InitBatch', pr.MCinitBatch, ...
                'BinWidth', pr.MCbinWidth, ...
                'TossFrames', pr.MCTossFrames);
        else
            pass_max_shift = max(round(max_shift_px / 4), 1);
            mcOut = run_rigid_mc(mc_input, ...
                'MaxShift', pass_max_shift, ...
                'InitBatch', pr.MCinitBatch, ...
                'BinWidth', pr.MCbinWidth, ...
                'Template', mcOut.template);
        end
        fprintf("[MC pass %d] max_shift = %d px\n", iPass, pass_max_shift);
        mc_input = string(mcOut.mc_path);
    end

    mc_tif = string(mcOut.mc_path);
    assert(mc_tif ~= "" && isfile(mc_tif), "Could not locate MC tif output.");

    % save MC info first (protect expensive result)
    sx = []; sy = []; r = [];
    save(char(out_mc_mat), "mcOut","mc_tif","sx","sy","r","max_shift_px","-v7.3");
    fprintf("Saved MC info:\n%s\n", out_mc_mat);

    % shift plot (try-catch)
    try
        if isfield(mcOut,"shifts")
            sh = mcOut.shifts;
            Tsh = numel(sh);
            sx = zeros(Tsh,1); sy = zeros(Tsh,1);
            for t = 1:Tsh
                sft = sh(t).shifts;
                sx(t) = sft(1); sy(t) = sft(2);
            end
            r = hypot(sx, sy);

            fig3 = figure('Color','w');
            plot(r,'LineWidth',1.5); grid on;
            xlabel('Frame'); ylabel('|shift| (pixel)');
            title(sprintf('|Rigid shift| per frame (pass %d)', nPasses));
            exportgraphics(fig3, char(out_mc_png), "Resolution", 200);
            close(fig3);

            save(char(out_mc_mat), "mcOut","mc_tif","sx","sy","r","max_shift_px","-v7.3");
        end
    catch ME
        fprintf("WARNING: shift plot failed (%s) — MC results are saved.\n", ME.message);
    end
end

% --- STEP 3: cpSAM ---
if pr.RerunCellpose
    % Delete existing cpSAM + dFF so they get regenerated
    delPatterns = ["*cpSAM_output.mat", "*_dFF.mat", "*_stackDFF.png"];
    for ip = 1:numel(delPatterns)
        hits = dir(fullfile(outFolder, delPatterns(ip)));
        for ih = 1:numel(hits)
            delPath = fullfile(hits(ih).folder, hits(ih).name);
            fprintf("[RERUN] Deleting: %s\n", delPath);
            delete(delPath);
        end
    end
end

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
        'FPS', pr.FPS, ...
        'Diameter', pr.Diameter, ...
        'PythonExe', char(pyExe), ...
        'UseGPU', pr.UseGPU, ...
        'FlowThreshold', pr.FlowThreshold, ...
        'CellprobThreshold', pr.CellprobThreshold);

    sam_mat = infer_sam_mat(samOut, outFolder, stem);
    fprintf("SAM mat:\n%s\n", sam_mat);
end

% --- STEP 4: dFF + stack plot ---
rerunDFF = true;
if pr.RerunStackDFF && isfile(out_dff_mat)
    % Replot only — load existing dFF, delete old png, replot
    if isfile(out_stack_png), delete(char(out_stack_png)); end
    fprintf("[REPLOT] stackDFF only (FPS=%d)\n", detectedFPS);
    tmp = load(char(out_dff_mat), "dFF");
    dFF = tmp.dFF;

    fig = figure("Color","w");
    stackDFF_with_scalebar(dFF, detectedFPS, pr.DffScaleBar, pr.TraceGain, pr.MinOffset, pr.OffsetMult);
    title(sprintf("stackDFF | %s | fps=%d toss=%d", stem, detectedFPS, pr.TossFrames), "Interpreter","none");
    exportgraphics(fig, char(out_stack_png), "Resolution", 200);
    close(fig);
    fprintf("Saved stackDFF:\n%s\n", out_stack_png);
    rerunDFF = false;

elseif isfile(out_dff_mat) && isfile(out_stack_png)
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

    S = load(sam_mat);
    assert(isfield(S,"F"), "SAM mat does not contain field F.");
    F_roi_raw = S.F;

    F_roi = F_roi_raw;
    if pr.TossFrames > 0
        assert(size(F_roi,1) > pr.TossFrames, "Not enough frames to toss for dFF.");
        F_roi(1:pr.TossFrames,:) = [];
    end

    dFFout = helper.dFF_RZ(F_roi, 'FPS', detectedFPS, 'BaselineWinSec', pr.BaselineWinSec);
    dFF = dFFout.dFF;

    params = struct();
    params.sam_mat       = sam_mat;
    params.fps           = detectedFPS;
    params.tossFrames    = pr.TossFrames;
    params.baselineWinSec = pr.BaselineWinSec;
    params.traceGain     = pr.TraceGain;
    params.minOffset  = pr.MinOffset;
    params.offsetMult = pr.OffsetMult;
    params.dffScaleBar = pr.DffScaleBar;

    save(char(out_dff_mat), "dFF","dFFout","F_roi_raw","F_roi","params","-v7.3");

    fig = figure("Color","w");
    stackDFF_with_scalebar(dFF, detectedFPS, pr.DffScaleBar, pr.TraceGain, pr.MinOffset, pr.OffsetMult);
    title(sprintf("stackDFF | %s | fps=%d toss=%d", stem, detectedFPS, pr.TossFrames), "Interpreter","none");
    exportgraphics(fig, char(out_stack_png), "Resolution", 200);
    close(fig);

    fprintf("Saved dFF:\n%s\n", out_dff_mat);
    fprintf("Saved stackDFF:\n%s\n", out_stack_png);
end

% --- STEP 5: Per-ROI videos ---
if pr.MakeMontageVideo && sam_mat ~= "" && isfile(sam_mat)
    movieDir = fullfile(outFolder, "Fmovie_perROI");
    if pr.RerunMontageVideo && isfolder(movieDir)
        fprintf("[RERUN] Deleting existing per-ROI videos: %s\n", movieDir);
        rmdir(char(movieDir), 's');
    end
    try
        generate_perROI_videos(char(outFolder), char(sam_mat), 15);
        fprintf("[DONE] Per-ROI videos saved.\n");
    catch ME
        fprintf(2, "[WARN] Per-ROI videos failed: %s\n", ME.message);
    end
end

% --- output struct ---
out = struct();
out.tifPath = tifPath;
out.outFolder = outFolder;
out.savedChannels = savedCh;
out.nSaved = nSaved;
out.keepChannelID_used = keepID_used;

out.preproc_tif = out_preproc_tif;
out.mc_tif      = mc_tif;
out.mcInfo_mat  = out_mc_mat;
out.mcShift_png = out_mc_png;
out.sam_mat     = sam_mat;
out.dff_mat     = out_dff_mat;
out.stack_png   = out_stack_png;
out.meta_mat    = out_meta_mat;

end

%% ======================================================================
function stackDFF_with_scalebar(dFF, fps, barAmp, gain, minOffset, offsetMult)
    [T, N] = size(dFF);
    t = (0:T-1) / fps;

    % --- force white-on-black regardless of MATLAB theme ---
    ax = gca;
    set(ax, "Color","w", "XColor","k", "YColor","k");
    set(ax.Parent, "Color","w");

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
    y0 = 0.10*offset;
    plot([x0 x0], [y0 y0+barH], "k", "LineWidth", 2);
    text(x0, y0+barH/2, sprintf("  %d%% dFF", round(barAmp*100)), ...
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

if ~isempty(regexp(rhs, "^\d+$", "once"))
    m = str2double(rhs);
    if isfinite(m)
        ch = find(bitget(uint32(m), 1:nAvail));
        return;
    end
end

nums = regexp(rhs, "[-+]?\d+\.?\d*", "match");
if ~isempty(nums)
    v = str2double(string(nums));
    v = v(isfinite(v));
    if ~isempty(v) && all(mod(v,1)==0) && all(v>=1) && all(v<=nAvail)
        ch = unique(v(:).', "stable");
        return;
    end
end

tf = regexp(lower(rhs), "(true|false)", "match");
if ~isempty(tf)
    mask = strcmp(tf, "true");
    mask = mask(1:min(numel(mask), nAvail));
    ch = find(mask);
    return;
end
end

%% ====================== helpers: locate outputs ========================
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
