%% single_MC_SAM_260324.m
% Single-file pipeline (no batch, no dark subtraction):
%   - AUTO ScanImage channel detection from TIFF tags
%   - Deinterleave only if needed (keep KeepChannelID)
%   - Keep only KeepFrames frames (after deinterleave) — set 0 or Inf to keep all
%   - Two-pass rigid MC -> cpSAM
%   - dFF (toss first frames) -> save *_dFF.mat + *_stackDFF.png
%   - SKIP steps if outputs exist
%
% No dark subtraction — raw data is cast to uint16 directly.
%
% Auto max_shift from zoom:
%   Pixel size = 1.7778 / scanZoomFactor (um/px).
%   Calibrated: 6x->0.2963, 4.2x->0.4233, 1.2x->1.4815 um/px.
%   MaxShiftUm (default 10 um) is converted to pixels automatically.
%
% Multi-pass MC:
%   Pass 1 registers with auto-generated template (noisy).
%   Pass 2 re-registers using the cleaner template from pass 1.
%   Both passes are rigid, integer-pixel shifts.
%
% Requires on path:
%   loadtiff, saveastiff, run_rigid_mc, cp_sam_extract_F_cli
% Optional:
%   helper.dFF_RZ
clear; clc; delete(gcp('nocreate'));

%% ============================== USER ==============================
filePath = "D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\random\7x_-850-1032-z30_3000f_00001\7x_-850-1032-z30_3000f_00001.tif";   % <<< exact TIFF path

% -------- channel handling --------
AutoDetectChannels  = true;   % parse SI tags if present
KeepChannelID       = 1;      % keep ch1
NumChannelsFallback = 1;      % if metadata missing or nonsense

% -------- frame selection --------
KeepFrames = 2000;   % after deinterleave, keep this many frames (0 or Inf = keep all)

% -------- trim (columns) --------
TrimL = 0;
TrimR = 0;

% -------- MC --------
MaxShiftUm     = 10;      % max physical shift (um) — auto-converted to pixels from zoom
PixelSizeBase  = 1.7778;  % um/px at 1x zoom (calibrated from 6x, 4.2x, 1.2x data)
MCpasses       = 2;       % number of rigid MC passes (2 = template refinement)
MCinitBatch    = 500;     % frames for initial template estimation
MCbinWidth     = 50;      % frames per batch (~1.7s at 30fps, ~2 breath cycles)
MCTossFrames   = 30;      % toss first N frames before MC (shutter/PMT artifacts)

% -------- cpSAM --------
FPS               = 30;
Diameter          = 30;
UseGPU            = true;
FlowThreshold     = 0.50;
CellprobThreshold = 0;

% -------- dFF --------
TossFrames = 0;       % already tossed before MC via MCTossFrames
TraceGain  = 2.5;
MinOffset  = 1.0;
OffsetMult = 12;
DffScaleBar = 0.20;   % 20%

% -------- python exe --------
cand = string([
    fullfile(getenv("USERPROFILE"), ".conda",     "envs", "cellpose-gpu", "python.exe")
    fullfile(getenv("USERPROFILE"), "anaconda3",  "envs", "cellpose-gpu", "python.exe")
    fullfile(getenv("USERPROFILE"), "miniconda3", "envs", "cellpose-gpu", "python.exe")
]);
pyExe = cand(find(isfile(cand), 1, "first"));
assert(pyExe ~= "", "Cellpose python not found. Checked:\n%s", strjoin(cand, newline));
fprintf("Using PythonExe:\n%s\n", pyExe);

%% ============================== RUN ==============================
out = run_single_MC_SAM(filePath, pyExe, ...
    "AutoDetectChannels",AutoDetectChannels, ...
    "KeepChannelID",KeepChannelID, ...
    "NumChannelsFallback",NumChannelsFallback, ...
    "KeepFrames",KeepFrames, ...
    "TrimL",TrimL, "TrimR",TrimR, ...
    "MaxShiftUm",MaxShiftUm, "PixelSizeBase",PixelSizeBase, "MCpasses",MCpasses, ...
    "MCinitBatch",MCinitBatch, "MCbinWidth",MCbinWidth, "MCTossFrames",MCTossFrames, ...
    "FPS",FPS, "Diameter",Diameter, ...
    "UseGPU",UseGPU, "FlowThreshold",FlowThreshold, "CellprobThreshold",CellprobThreshold, ...
    "TossFrames",TossFrames, ...
    "TraceGain",TraceGain, "MinOffset",MinOffset, "OffsetMult",OffsetMult, ...
    "DffScaleBar",DffScaleBar);

disp(out);

%% ======================================================================
function out = run_single_MC_SAM(filePath, pyExe, varargin)

p = inputParser;
p.addRequired("filePath", @(s)isstring(s)||ischar(s));
p.addRequired("pyExe", @(s)isstring(s)||ischar(s));

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
p.addParameter("Diameter", 30, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("UseGPU", true, @(x)islogical(x)&&isscalar(x));
p.addParameter("FlowThreshold", 0.50, @(x)isnumeric(x)&&isscalar(x));
p.addParameter("CellprobThreshold", 0, @(x)isnumeric(x)&&isscalar(x));

% dFF + plot
p.addParameter("TossFrames", 30, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("TraceGain", 2.5, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("MinOffset", 1.0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("OffsetMult", 12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("DffScaleBar", 0.20, @(x)isnumeric(x)&&isscalar(x)&&x>0);

p.parse(filePath, pyExe, varargin{:});
pr = p.Results;

filePath = string(pr.filePath);
pyExe    = string(pr.pyExe);
assert(isfile(filePath), "TIFF not found: %s", filePath);
assert(isfile(pyExe), "PythonExe not found: %s", pyExe);

[parentDir, baseName, ~] = fileparts(char(filePath));
parentDir = string(parentDir);
baseName  = string(baseName);

% --- create output folder for this file (in same parent directory) ---
outFolder = fullfile(parentDir, baseName);
if ~isfolder(outFolder), mkdir(char(outFolder)); end
fprintf("[setup] Output folder: %s\n", outFolder);

% --- detect channels from TIFF tags ---
info = imfinfo(filePath);
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
    fprintf(2, "[chanDetect] WARNING: nPages=%d not divisible by nSaved=%d from metadata. Treating as single-channel.\n", nPages, nSaved);
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
    max_shift_px  = max(max_shift_px, 1);  % at least 1 pixel
    fprintf("[MC] zoom=%.2fx, pixel=%.4f um/px, MaxShift=%d um -> %d px\n", ...
        zoomFactor, pixel_size_um, pr.MaxShiftUm, max_shift_px);
else
    max_shift_px = 20;
    fprintf(2, "[MC] WARNING: could not read scanZoomFactor from metadata. Using default max_shift=%d px\n", max_shift_px);
end

% --- define output stem ---
tag = "ch" + string(keepID_used);
if contains(lower(baseName), lower(tag))
    stem = baseName;
else
    stem = baseName + "_" + tag;
end

% output paths (no dark step — input TIFF goes straight to MC)
out_preproc_tif = outFolder + filesep + stem + "_preproc.tif";

out_mc_mat   = outFolder + filesep + stem + "_MCinfo.mat";
out_mc_png   = outFolder + filesep + stem + "_MCshift.png";
out_mc_guess = outFolder + filesep + stem + "_preproc_MC.tif";

out_sam_guess = outFolder + filesep + stem + "_preproc_MC_cpSAM_output.mat";
out_dff_mat   = outFolder + filesep + stem + "_dFF.mat";
out_stack_png = outFolder + filesep + stem + "_stackDFF.png";

% --- STEP 1: LOAD + DEINTERLEAVE + KEEP FRAMES (skip if preproc tif exists) ---
if isfile(out_preproc_tif)
    fprintf("[SKIP] preproc step (found): %s\n", out_preproc_tif);
else
    fprintf("[RUN ] preproc step (load + deinterleave + keep frames)\n");

    F_raw = loadtiff(filePath); % uint16 [Y X Tpages]
    F = F_raw;
    if pr.TrimL > 0, F(:,1:pr.TrimL,:) = []; end
    if pr.TrimR > 0, F(:,end-pr.TrimR+1:end,:) = []; end

    % deinterleave only if needed
    if nSaved > 1
        assert(mod(size(F,3), nSaved) == 0, "Pages=%d not divisible by nSaved=%d.", size(F,3), nSaved);
        F_keep = F(:,:,kWithin:nSaved:end);
    else
        F_keep = F;
    end
    fprintf("[data] raw pages=%d, after deinterleave=%d\n", size(F,3), size(F_keep,3));

    % keep only KeepFrames (0 or Inf = keep all)
    if pr.KeepFrames > 0 && isfinite(pr.KeepFrames) && size(F_keep,3) > pr.KeepFrames
        F_keep = F_keep(:,:,1:pr.KeepFrames);
        fprintf("[data] trimmed to KeepFrames=%d\n", pr.KeepFrames);
    else
        fprintf("[data] keeping all %d frames\n", size(F_keep,3));
    end

    % cast to uint16 (ScanImage may store as int16, which imwrite can't handle)
    F_keep = uint16(max(single(F_keep), 0));

    % save preprocessed TIFF (no dark subtraction)
    out_preproc_tif_c = char(out_preproc_tif);
    if exist(out_preproc_tif_c, "file"), delete(out_preproc_tif_c); end
    opts = struct('overwrite', true, 'message', false);
    saveastiff(F_keep, out_preproc_tif_c, opts);

    fprintf("Saved preproc TIFF:\n%s\n", out_preproc_tif);
end

% --- STEP 2: MC (skip if MC tif exists OR MCinfo exists + can infer tif) ---
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
            % refinement pass: 1/4 max_shift, with cleaner template from previous pass
            % no TossFrames on pass 2 (already tossed in pass 1 output)
            pass_max_shift = max(round(max_shift_px / 4), 1);
            mcOut = run_rigid_mc(mc_input, ...
                'MaxShift', pass_max_shift, ...
                'InitBatch', pr.MCinitBatch, ...
                'BinWidth', pr.MCbinWidth, ...
                'Template', mcOut.template);
        end
        fprintf("[MC pass %d] max_shift = %d px\n", iPass, pass_max_shift);
        mc_input = string(mcOut.mc_path);  % next pass uses this pass's output
    end

    mc_tif = string(mcOut.mc_path);
    assert(mc_tif ~= "" && isfile(mc_tif), "Could not locate MC tif output.");

    % shift plot (from final pass)
    sx = []; sy = []; r = [];
    if isfield(mcOut,"shifts")
        sh = mcOut.shifts;
        T = numel(sh);
        sx = zeros(T,1); sy = zeros(T,1);
        for t = 1:T
            sft = sh(t).shifts;
            sx(t) = sft(1); sy(t) = sft(2);
        end
        r = hypot(sx, sy);

        fig3 = figure('Color','w');
        plot(r,'LineWidth',1.5); grid on;
        xlabel('Frame'); ylabel('|shift| (pixel)');
        title(sprintf('|Rigid shift| per frame (pass %d)', nPasses));
        exportgraphics(fig3, out_mc_png, "Resolution", 200);
        close(fig3);
    end

    save(char(out_mc_mat), "mcOut","mc_tif","sx","sy","r","max_shift_px","-v7.3");
    fprintf("Saved MC info:\n%s\n", out_mc_mat);
end

% --- STEP 3: cpSAM (skip if cpSAM output exists) ---
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

% --- STEP 4: dFF + stack plot (skip if both exist) ---
if isfile(out_dff_mat) && isfile(out_stack_png)
    fprintf("[SKIP] dFF+stack step (found): %s\n", out_dff_mat);
else
    fprintf("[RUN ] dFF+stack step\n");

    S = load(sam_mat);
    assert(isfield(S,"F"), "SAM mat does not contain field F.");
    F_roi_raw = S.F;  % [T x Nroi]

    F_roi = F_roi_raw;
    if pr.TossFrames > 0
        assert(size(F_roi,1) > pr.TossFrames, "Not enough frames to toss for dFF.");
        F_roi(1:pr.TossFrames,:) = [];
    end

    % dFF
    dFFout = dFF_RZ_dispatch(F_roi);
    dFF = dFFout.dFF;

    params = struct();
    params.sam_mat    = sam_mat;
    params.fps        = pr.FPS;
    params.tossFrames = pr.TossFrames;
    params.traceGain  = pr.TraceGain;
    params.minOffset  = pr.MinOffset;
    params.offsetMult = pr.OffsetMult;
    params.dffScaleBar = pr.DffScaleBar;

    save(char(out_dff_mat), "dFF","dFFout","F_roi_raw","F_roi","params","-v7.3");

    fig = figure("Color","w");
    stackDFF_with_scalebar(dFF, pr.FPS, pr.DffScaleBar, pr.TraceGain, pr.MinOffset, pr.OffsetMult);
    title(sprintf("stackDFF | %s | toss=%d", stem, pr.TossFrames), "Interpreter","none");
    exportgraphics(fig, char(out_stack_png), "Resolution", 200);
    close(fig);

    fprintf("Saved dFF:\n%s\n", out_dff_mat);
    fprintf("Saved stackDFF:\n%s\n", out_stack_png);
end

% --- output struct ---
out = struct();
out.filePath = filePath;
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

end

%% ======================================================================
% ===================== dFF: helper if exists, else local =================
function dFFout = dFF_RZ_dispatch(F)
if exist("helper.dFF_RZ","file")==2
    dFFout = helper.dFF_RZ(F);
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

%% ======================================================================
% ====================== ScanImage channel parsing =======================
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

%% ====================== helpers: locate outputs =========================
function mc_tif = infer_mc_tif(mcOut, folderPath, stem)
mc_tif = "";
candFields = ["outPath","outFile","mcPath","tifPath","mcTif","mcTifPath"];
for f = candFields
    if isstruct(mcOut) && isfield(mcOut, f)
        val = mcOut.(f);
        if isstring(val) || ischar(val)
            mc_tif = string(val);
            if isfile(mc_tif), return; end
        end
    end
end

guess = fullfile(folderPath, stem + "_preproc_MC.tif");
if isfile(guess), mc_tif = string(guess); return; end

mc_tif = newest_matching_file(folderPath, stem + "_preproc*MC*.tif");
end

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
