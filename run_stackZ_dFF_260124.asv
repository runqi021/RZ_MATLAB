%% run_zstack_pipeline_full.m
% Manual Z-stack split (NO reliance on SI stack metadata) + per-slice pipeline:
%   slice tif -> dark -> rigid MC -> cpSAM -> dFF (toss first N frames)
%
% NOW: auto-detect ScanImage saved channels from TIFF tags and deinterleave if needed.
%
% Output:
%   <same folder as zstack tif>/<base>_ZSTACK_PIPELINE/z0, z20, ..., z400/
%   and a master index: <base>_ZSTACK_index.mat

clear; clc;

%% ========================= USER INPUTS =========================
zstackTif = "E:\260124_chat_soma_g8s+cy5\7N_stack_2\7N_R_-1000_2050_z0-30_00003.tif";

% Z definition (MANUAL)
ZStart    = 0;
ZStep     = 5;
NumSlices = 11;

% Frames per slice (MANUAL)  (this is FRAMES PER SLICE PER CHANNEL TIMEPOINTS)
FramesPerSlice = 3000;

% If acquisition is slice blocks (all frames of z0, then all frames of zStep, ...)
SliceOrder = "block";    % "block" | "interleaved"

% dFF settings
TossFrames = 30;         % toss first N frames BEFORE dFF

% ---- channel handling (AUTO) ----
AutoDetectChannels  = true;  % reads ScanImage tags from TIFF to infer saved channels
KeepChannelID       = 1;     % keep Channel 1 (cells)
NumChannelsFallback = 1;     % used only if no SI metadata found

% cpSAM settings
FPS               = 30;
Diameter          = 30;
UseGPU            = true;
FlowThreshold     = 0.50;
CellprobThreshold = 0;

% dark / trim settings
TrimL          = 10;
TrimR          = 10;
DarkMaxFrames  = 500;
DarkSubsample  = 4;
ApplyDarkMask  = true;
DarkSigma      = 3;
SaveDiagnostics = true;

%% ========================= Python auto-detect =========================
cand = string([
    fullfile(getenv("USERPROFILE"), ".conda",     "envs", "cellpose-gpu", "python.exe")
    fullfile(getenv("USERPROFILE"), "anaconda3",  "envs", "cellpose-gpu", "python.exe")
    fullfile(getenv("USERPROFILE"), "miniconda3", "envs", "cellpose-gpu", "python.exe")
]);

pyExe = cand(find(isfile(cand), 1, "first"));
assert(pyExe ~= "", "Cellpose python not found. Checked:\n%s", strjoin(cand, newline));
fprintf("Using PythonExe:\n%s\n", pyExe);

%% ========================= RUN =========================
idx = run_zstack_pipeline_manual( ...
    zstackTif, pyExe, ...
    "ZStart",ZStart, "ZStep",ZStep, "NumSlices",NumSlices, ...
    "FramesPerSlice",FramesPerSlice, ...
    "SliceOrder",SliceOrder, ...
    "AutoDetectChannels",AutoDetectChannels, ...
    "KeepChannelID",KeepChannelID, ...
    "NumChannelsFallback",NumChannelsFallback, ...
    "TrimL",TrimL, "TrimR",TrimR, ...
    "DarkMaxFrames",DarkMaxFrames, "DarkSubsample",DarkSubsample, ...
    "ApplyDarkMask",ApplyDarkMask, "DarkSigma",DarkSigma, ...
    "SaveDiagnostics",SaveDiagnostics, ...
    "FPS",FPS, "Diameter",Diameter, ...
    "UseGPU",UseGPU, "FlowThreshold",FlowThreshold, "CellprobThreshold",CellprobThreshold, ...
    "TossFrames",TossFrames, ...
    "TraceGain",2.5, "MinOffset",1.0, "OffsetMult",12);

disp(idx);

%% ======================================================================
function idx = run_zstack_pipeline_manual(zstackTif, pyExe, varargin)
% Manual Z-stack splitter + pipeline runner.
% Channel handling:
%   - Detect saved channels from ScanImage metadata in TIFF tags.
%   - Treat TIFF pages as: (frameIgnoringChannelIndex) * nSavedChannels + channelWithinSaved
%
% Assumes TIFF is either:
%   - "block": slice1 timepoints 1..K, slice2 timepoints K+1..2K, ...
%   - "interleaved": timepoints cycle through slices: t1=z1, t2=z2, ..., repeat
%
% Where "timepoint" here means one acquired frame at some z, ignoring channel pages.

p = inputParser;
p.addRequired("zstackTif", @(s)isstring(s)||ischar(s));
p.addRequired("pyExe",     @(s)isstring(s)||ischar(s));

% manual Z spec
p.addParameter("ZStart", 0, @(x)isnumeric(x)&&isscalar(x));
p.addParameter("ZStep",  20, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("NumSlices", 21, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% manual time spec (per slice)
p.addParameter("FramesPerSlice", 1000, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% order
p.addParameter("SliceOrder", "block", @(s)isstring(s)||ischar(s)); % "block" | "interleaved"

% channel handling
p.addParameter("AutoDetectChannels", true, @(x)islogical(x)&&isscalar(x));
p.addParameter("KeepChannelID", 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter("NumChannelsFallback", 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% pass-through to getSAM
p.addParameter("TrimL",10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("TrimR",10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("DarkMaxFrames",500, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("DarkSubsample",4, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter("ApplyDarkMask",true, @(x)islogical(x)&&isscalar(x));
p.addParameter("DarkSigma",3, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("SaveDiagnostics",true, @(x)islogical(x)&&isscalar(x));

% cpSAM
p.addParameter("FPS",30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("Diameter",30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("UseGPU",true, @(x)islogical(x)&&isscalar(x));
p.addParameter("FlowThreshold",0.50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter("CellprobThreshold",0, @(x)isnumeric(x)&&isscalar(x));

% dFF
p.addParameter("TossFrames",30, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("TraceGain",2.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("MinOffset",1.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("OffsetMult",12, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.parse(zstackTif, pyExe, varargin{:});
pr = p.Results;

zstackTif = string(pr.zstackTif);
pyExe     = string(pr.pyExe);
assert(isfile(zstackTif), "Z-stack TIFF not found: %s", zstackTif);
assert(isfile(pyExe),     "PythonExe not found: %s", pyExe);

% Z list (manual)
zList = pr.ZStart + (0:pr.NumSlices-1)*pr.ZStep;

% TIFF info
info = imfinfo(zstackTif);
nPages = numel(info);

% ---- detect saved channels ----
savedCh = [];
if pr.AutoDetectChannels
    savedCh = scanimage_saved_channels_from_info(info(1));
end
if isempty(savedCh)
    savedCh = 1:pr.NumChannelsFallback;
end
savedCh = savedCh(:).';
nSaved = numel(savedCh);

% map KeepChannelID -> position within saved channels
kWithin = find(savedCh == pr.KeepChannelID, 1, "first");
keepID_used = pr.KeepChannelID;
if isempty(kWithin)
    kWithin = 1;
    keepID_used = savedCh(1);
    fprintf(2, "[chanDetect] KeepChannelID=%d not found in savedCh=[%s]. Using Channel %d.\n", ...
        pr.KeepChannelID, num2str(savedCh), keepID_used);
end

fprintf("[chanDetect] savedCh=[%s] (nSaved=%d), keepChannelID=%d (kWithin=%d)\n", ...
    num2str(savedCh), nSaved, keepID_used, kWithin);

% ---- timepoint accounting (timepoints = pages / nSaved) ----
assert(mod(nPages, nSaved) == 0, ...
    "TIFF has %d pages but nSavedChannels=%d -> not divisible. Check metadata or fallback.", nPages, nSaved);

nTimepointsTotal = nPages / nSaved;  % frames ignoring channel pages
expectedTimepoints = pr.NumSlices * pr.FramesPerSlice;

assert(nTimepointsTotal == expectedTimepoints, ...
    "Timepoint count mismatch.\nTIFF has %d timepoints (= %d pages / %d savedCh), but NumSlices*FramesPerSlice = %d*%d = %d.\nFix inputs or verify TIFF.", ...
    nTimepointsTotal, nPages, nSaved, pr.NumSlices, pr.FramesPerSlice, expectedTimepoints);

% output root
[inDir, base] = fileparts(zstackTif);
outRoot = fullfile(inDir, base + "_ZSTACK_PIPELINE");
if ~exist(outRoot, "dir"), mkdir(outRoot); end

idx = struct();
idx.zstackTif = zstackTif;
idx.outRoot   = string(outRoot);
idx.zList_um  = zList(:);
idx.nPages    = nPages;
idx.nTimepointsTotal = nTimepointsTotal;

idx.numSlices = pr.NumSlices;
idx.framesPerSlice = pr.FramesPerSlice;
idx.sliceOrder = string(pr.SliceOrder);

idx.savedChannels   = savedCh;
idx.nSavedChannels  = nSaved;
idx.keepChannelID_requested = pr.KeepChannelID;
idx.keepChannelID_used      = keepID_used;
idx.keepChannelPosWithinSaved = kWithin;

idx.slices = repmat(struct( ...
    "z_um",[], "subdir","", ...
    "timepoints",[], "pages",[], ...
    "slice_tif","", ...
    "sam_mat","", "mc_info_mat","", "dff_mat","", "dff_png",""), pr.NumSlices, 1);

sliceOrder = lower(string(pr.SliceOrder));

% open once for reading pages
tIn = Tiff(zstackTif, "r");

for s = 1:pr.NumSlices
    z = zList(s);
    subdir = fullfile(outRoot, sprintf("z%d", round(z)));
    if ~exist(subdir, "dir"), mkdir(subdir); end

    sliceTif = fullfile(subdir, sprintf("%s_z%d_ch%d.tif", base, round(z), keepID_used));

    % ---- timepoint indices (ignoring channel pages) ----
    if sliceOrder == "block"
        tStart = (s-1)*pr.FramesPerSlice + 1;
        tEnd   = s*pr.FramesPerSlice;
        tpIdx  = tStart:tEnd;
    elseif sliceOrder == "interleaved"
        tpIdx = s:pr.NumSlices:nTimepointsTotal;
        assert(numel(tpIdx) == pr.FramesPerSlice, ...
            "Interleaved indexing produced %d timepoints (expected %d).", numel(tpIdx), pr.FramesPerSlice);
    else
        error('SliceOrder must be "block" or "interleaved".');
    end
    tpIdx = tpIdx(:);

    % ---- map timepoints -> TIFF page indices for kept channel ----
    pageIdx = (tpIdx - 1) * nSaved + kWithin;
    assert(all(pageIdx >= 1 & pageIdx <= nPages), "Computed pageIdx out of bounds.");

    % write slice tif (uint16, LZW)
    if exist(sliceTif, "file"), delete(sliceTif); end
    for k = 1:numel(pageIdx)
        tIn.setDirectory(pageIdx(k));
        fr = tIn.read();
        if ~isa(fr,"uint16"), fr = uint16(fr); end
        if k == 1
            imwrite(fr, sliceTif, "tif", "Compression","lzw");
        else
            imwrite(fr, sliceTif, "tif", "WriteMode","append", "Compression","lzw");
        end
    end

    % run getSAM inside the slice folder (outputs land next to sliceTif)
    infoSAM = getSAM(sliceTif, ...
        "TrimL",pr.TrimL,"TrimR",pr.TrimR, ...
        "DarkMaxFrames",pr.DarkMaxFrames,"DarkSubsample",pr.DarkSubsample, ...
        "ApplyDarkMask",pr.ApplyDarkMask,"DarkSigma",pr.DarkSigma, ...
        "SaveDiagnostics",pr.SaveDiagnostics, ...
        "FPS",pr.FPS,"Diameter",pr.Diameter, ...
        "PythonExe",pyExe, ...
        "UseGPU",pr.UseGPU,"FlowThreshold",pr.FlowThreshold,"CellprobThreshold",pr.CellprobThreshold);

    % dFF
    outDFF = run_and_plot_dFF(infoSAM.sam_mat, ...
        "FPS",pr.FPS,"TossFrames",pr.TossFrames, ...
        "TraceGain",pr.TraceGain,"MinOffset",pr.MinOffset,"OffsetMult",pr.OffsetMult);

    % record
    idx.slices(s).z_um        = z;
    idx.slices(s).subdir      = string(subdir);

    idx.slices(s).timepoints  = tpIdx;
    idx.slices(s).pages       = pageIdx;

    idx.slices(s).slice_tif   = string(sliceTif);
    idx.slices(s).sam_mat     = infoSAM.sam_mat;
    idx.slices(s).mc_info_mat = infoSAM.mc_info_mat;
    idx.slices(s).dff_mat     = outDFF.dff_mat;
    idx.slices(s).dff_png     = outDFF.dff_png;

    fprintf("[slice %02d/%02d | z=%d | ch=%d] wrote %d frames -> %s\n", ...
        s, pr.NumSlices, round(z), keepID_used, numel(tpIdx), subdir);
end

tIn.close();

indexPath = fullfile(outRoot, base + "_ZSTACK_index.mat");
save(indexPath, "idx", "-v7.3");
fprintf("\nDONE. Index saved:\n%s\n", indexPath);
end

%% ======================================================================
function info = getSAM(tiffPath, varargin)
% getSAM
% Pipeline per slice:
%   load -> trim -> estimate dark -> subtract (safe) -> optional dark mask
%   -> save *_minusDark.tif
%   -> run_rigid_mc -> save *_MCinfo.mat + shift plot
%   -> cp_sam_extract_F_cli -> return SAM mat path

p = inputParser;
p.addRequired("tiffPath", @(s)isstring(s)||ischar(s));

p.addParameter("TrimL", 10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("TrimR", 10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

p.addParameter("DarkMaxFrames", 500, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("DarkSubsample", 4, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

p.addParameter("ApplyDarkMask", true, @(x)islogical(x)&&isscalar(x));
p.addParameter("DarkSigma", 3, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.addParameter("SaveDiagnostics", true, @(x)islogical(x)&&isscalar(x));

% cpSAM params
p.addParameter("FPS", 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("Diameter", 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("PythonExe", "", @(s)isstring(s)||ischar(s));
p.addParameter("UseGPU", true, @(x)islogical(x)&&isscalar(x));
p.addParameter("FlowThreshold", 0.50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter("CellprobThreshold", 0, @(x)isnumeric(x)&&isscalar(x));

p.parse(tiffPath, varargin{:});
pr = p.Results;

tiffPath = string(pr.tiffPath);
assert(isfile(tiffPath), "TIFF not found: %s", tiffPath);
assert(pr.PythonExe ~= "", "PythonExe is empty. Pass PythonExe=pyExe.");

[folderPath, baseName] = fileparts(tiffPath);

out_dark_tif   = fullfile(folderPath, baseName + "_minusDark.tif");
out_dark_mat   = fullfile(folderPath, baseName + "_darkInfo.mat");
out_dark_png   = fullfile(folderPath, baseName + "_darkMask_sumProj.png");
out_hist_png   = fullfile(folderPath, baseName + "_darkHist.png");

out_mc_mat     = fullfile(folderPath, baseName + "_MCinfo.mat");
out_mc_png     = fullfile(folderPath, baseName + "_MCshift.png");

%% load + trim (NOTE: loads full movie into RAM)
F = loadtiff(tiffPath);

if pr.TrimL > 0
    F(:, 1:pr.TrimL, :) = [];
end
if pr.TrimR > 0
    F(:, end-pr.TrimR+1:end, :) = [];
end

%% estimate dark
darkOut = estimate_dark_from_movie(F, ...
    'MaxFrames', pr.DarkMaxFrames, ...
    'Subsample', pr.DarkSubsample);

dark_mean = double(darkOut.dark_mean);
dark_std  = double(darkOut.dark_std);

%% subtract dark safely
Fs = single(F) - single(dark_mean);
Fs(Fs < 0) = 0;
F_dark = uint16(Fs);

%% optional per-pixel mask (NO giant 3D logical allocation)
mask = true(size(F_dark,1), size(F_dark,2));
thr  = pr.DarkSigma * dark_std;

if pr.ApplyDarkMask
    meanF = mean(single(F_dark), 3);
    mask  = meanF > thr;

    fracKeep = nnz(mask)/numel(mask)*100;
    fprintf('[darkMask] Keeping %.2f %% of pixels (mean > %.3f)\n', fracKeep, thr);

    for k = 1:size(F_dark,3)
        fr = F_dark(:,:,k);
        fr(~mask) = 0;
        F_dark(:,:,k) = fr;
    end
end

%% diagnostics
if pr.SaveDiagnostics
    sumProj = sum(single(F_dark), 3);
    fig = figure('Color','w');
    imagesc(sumProj); axis image off; colormap(gray); hold on;
    badMask = ~mask;
    [r,c] = find(badMask);
    if ~isempty(r)
        scatter(c, r, 2, 'r', 'filled');
    end
    title(sprintf('Sum projection with masked pixels (%.2f%%)', 100*nnz(badMask)/numel(badMask)));
    exportgraphics(fig, out_dark_png, "Resolution", 200);
    close(fig);

    if exist("Fhist", "file") == 2
        F_sat = 65535 - dark_mean;

        try
            Fhist(F_dark, sat=F_sat);
        catch
            Fhist(F_dark, F_sat);
        end

        figH = gcf;
        set(figH, "Color", "w");
        drawnow;

        exportgraphics(figH, char(out_hist_png), "Resolution", 200);
        close(figH);
    end
end

%% save dark tif + info
out_dark_tif_c = char(out_dark_tif);
if exist(out_dark_tif_c, "file"), delete(out_dark_tif_c); end

opts = struct('overwrite', true, 'message', false);   % avoid interactive retry prompt
saveastiff(F_dark, out_dark_tif_c, opts);

save(out_dark_mat, ...
    "tiffPath", "dark_mean", "dark_std", "thr", "mask", "darkOut", ...
    "-v7.3");

fprintf("Saved dark-corrected TIFF:\n%s\n", out_dark_tif);

%% rigid motion correction
mcOut = run_rigid_mc(out_dark_tif);

% infer MC tif path robustly
mc_tif = "";
candFields = ["outPath","outFile","mcPath","tifPath","mcTif","mcTifPath"];
for f = candFields
    if isfield(mcOut, f)
        val = mcOut.(f);
        if isstring(val) || ischar(val)
            mc_tif = string(val);
            if isfile(mc_tif), break; end
        end
    end
end

if mc_tif == "" || ~isfile(mc_tif)
    mc_guess = fullfile(folderPath, baseName + "_minusDark_MC.tif");
    if isfile(mc_guess)
        mc_tif = mc_guess;
    else
        mc_tif = newest_matching_file(folderPath, baseName + "_minusDark*MC*.tif");
    end
end
assert(mc_tif ~= "" && isfile(mc_tif), "Could not locate MC tif output for: %s", out_dark_tif);

%% save MC shifts + plot
sx = []; sy = []; r = [];
if isfield(mcOut, "shifts")
    sh = mcOut.shifts;
    T = numel(sh);
    sx = zeros(T,1); sy = zeros(T,1);
    for t = 1:T
        sft = sh(t).shifts;
        sx(t) = sft(1);
        sy(t) = sft(2);
    end
    r = hypot(sx, sy);

    fig3 = figure('Color','w');
    plot(r, 'LineWidth', 1.5); grid on;
    xlabel('Frame'); ylabel('|shift| (pixel)');
    title('|Rigid shift| per frame');
    exportgraphics(fig3, out_mc_png, "Resolution", 200);
    close(fig3);
end

save(out_mc_mat, "mcOut", "mc_tif", "sx", "sy", "r", "-v7.3");
fprintf("Saved MC info:\n%s\n", out_mc_mat);

%% cpSAM extraction
samOut = cp_sam_extract_F_cli(mc_tif, ...
    'FPS', pr.FPS, ...
    'Diameter', pr.Diameter, ...
    'PythonExe', char(pr.PythonExe), ...
    'UseGPU', pr.UseGPU, ...
    'FlowThreshold', pr.FlowThreshold, ...
    'CellprobThreshold', pr.CellprobThreshold);

% infer SAM mat path
sam_mat = "";
if isstring(samOut) || ischar(samOut)
    if endsWith(string(samOut), ".mat", "IgnoreCase", true)
        sam_mat = string(samOut);
    end
elseif isstruct(samOut)
    if isfield(samOut, "sam_mat"), sam_mat = string(samOut.sam_mat); end
    if sam_mat=="" && isfield(samOut, "out_mat"), sam_mat = string(samOut.out_mat); end
end
if sam_mat == "" || ~isfile(sam_mat)
    sam_mat = newest_matching_file(folderPath, baseName + "*cpSAM_output.mat");
end
assert(sam_mat ~= "" && isfile(sam_mat), "Could not locate cpSAM output mat for: %s", mc_tif);

info = struct();
info.input_tif     = tiffPath;
info.folder        = string(folderPath);
info.baseName      = string(baseName);

info.dark_tif      = string(out_dark_tif);
info.dark_info_mat = string(out_dark_mat);
info.dark_mask_png = string(out_dark_png);
info.dark_hist_png = string(out_hist_png);

info.mc_tif        = string(mc_tif);
info.mc_info_mat   = string(out_mc_mat);
info.mc_shift_png  = string(out_mc_png);

info.sam_mat       = string(sam_mat);
info.cp_sam_out    = samOut;

fprintf("SAM output:\n%s\n", info.sam_mat);
end

%% ======================================================================
function out = run_and_plot_dFF(samMatPath, varargin)
% Load cpSAM output (.mat with field F), toss first N frames, compute dF/F,
% save dFF .mat + stacked trace .png

p = inputParser;
p.addRequired("samMatPath", @(s)isstring(s)||ischar(s));
p.addParameter("FPS", 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("TossFrames", 30, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("TraceGain", 2.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("MinOffset", 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("OffsetMult", 12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("SavePrefix", "", @(s)isstring(s)||ischar(s));
p.parse(samMatPath, varargin{:});
pr = p.Results;

samMatPath = string(pr.samMatPath);
assert(isfile(samMatPath), "SAM mat not found: %s", samMatPath);

S = load(samMatPath);
assert(isfield(S,"F"), "SAM mat does not contain field F.");
F_raw = S.F;
F = F_raw;

if pr.TossFrames > 0
    assert(size(F,1) > pr.TossFrames, "Not enough frames to toss.");
    F(1:pr.TossFrames,:) = [];
end

dFFout = helper.dFF_RZ(F);
dFF = dFFout.dFF;

[folderPath, baseName] = fileparts(samMatPath);
stem = string(baseName);
if pr.SavePrefix ~= ""
    stem = string(pr.SavePrefix);
end

out_dff_mat = fullfile(folderPath, stem + "_dFF.mat");
out_png     = fullfile(folderPath, stem + "_stackDFF.png");

params = struct();
params.samMatPath  = samMatPath;
params.fps         = pr.FPS;
params.tossFrames  = pr.TossFrames;
params.traceGain   = pr.TraceGain;
params.minOffset   = pr.MinOffset;
params.offsetMult  = pr.OffsetMult;

save(out_dff_mat, "dFF", "dFFout", "F_raw", "F", "params", "-v7.3");

fig = figure("Color","w");
stackDFF_with_scalebar(dFF, pr.FPS, 0.20, pr.TraceGain, pr.MinOffset, pr.OffsetMult);
title(sprintf("stackDFF | %s | toss=%d frames", stem, pr.TossFrames), "Interpreter","none");
exportgraphics(fig, out_png, "Resolution", 200);
close(fig);

out = struct();
out.dff_mat  = string(out_dff_mat);
out.dff_png  = string(out_png);
out.params   = params;
out.nFrames  = size(dFF,1);
out.nROI     = size(dFF,2);
end

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
    y0 = 0.10*offset;
    plot([x0 x0], [y0 y0+barH], "k", "LineWidth", 2);
    text(x0, y0+barH/2, sprintf("  %d%% dFF", round(barAmp*100)), ...
        "Color","k","FontSize",10, "FontName","Arial","FontWeight","normal", ...
        "HorizontalAlignment","left", "VerticalAlignment","middle");
end

function fpath = newest_matching_file(folderPath, pattern)
    dd = dir(fullfile(folderPath, pattern));
    if isempty(dd)
        fpath = "";
        return
    end
    [~,ix] = max([dd.datenum]);
    fpath = fullfile(dd(ix).folder, dd(ix).name);
end

%% ======================================================================
% ====================== ScanImage channel parsing =======================
function savedCh = scanimage_saved_channels_from_info(info1)
% Return vector of saved channel IDs inferred from ScanImage metadata.
% Tries:
%   SI.hChannels.channelSave
%   SI.hChannels.channelsActive
% And interprets scalar as bitmask.

meta = string(info1.Software) + newline;

meta = string(meta);
if strlength(strtrim(meta)) == 0
    savedCh = [];
    return
end

% how many channels exist? (helps bitmask decode)
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
% parse "KEY = number"
v = [];
pat = key + "\s*=\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)";
tok = regexp(metaStr, pat, "tokens", "once");
if isempty(tok), return; end
vv = str2double(string(tok{1}));
if isfinite(vv), v = vv; end
end

function ch = parse_scanimage_channels(metaStr, key, nAvail)
% Parse ScanImage metadata line for channels.
% Supports:
%   KEY = 3            (bitmask -> [1 2] if nAvail>=2)
%   KEY = [1 2]
%   KEY = [true false true false]
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
    if ~isempty(v)
        if all(mod(v,1)==0) && all(v>=1) && all(v<=nAvail)
            ch = unique(v(:).', "stable");
            return;
        end
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
