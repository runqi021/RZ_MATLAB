%% run_singleZ_pipeline_full_autochan_skip_dff.m
% Single-Z pipeline with:
%   - AUTO ScanImage channel detection from TIFF tags
%   - Deinterleave only if needed (keep KeepChannelID)
%   - trim -> dark -> mask -> save *_minusDark.tif
%   - rigid MC -> cpSAM
%   - dFF (toss first frames) -> save *_dFF.mat + *_stackDFF.png
%   - SKIP steps if outputs exist
%
% Requires on path:
%   loadtiff, estimate_dark_from_movie, saveastiff, run_rigid_mc, cp_sam_extract_F_cli
% Optional:
%   Fhist, helper.dFF_RZ
clear; clc; delete(gcp('nocreate'));

%% ============================== USER ==============================
folderPath = 'D:\RUNQI\phys';
list = dir(folderPath);

% -------- channel handling --------
AutoDetectChannels  = true;   % parse SI tags if present
KeepChannelID       = 1;      % keep ch1
NumChannelsFallback = 1;      % if metadata missing or nonsense

% -------- trim --------
TrimL = 0;
TrimR = 0;

% -------- dark --------
DarkMaxFrames   = 500;
DarkSubsample   = 4;
DarkSigma       = 3;
ApplyDarkMask   = true;
SaveDiagnostics = true;
TossFramesForHist = 30; % for Fhist preview only

% -------- MC + cpSAM --------
FPS               = 30;
Diameter          = 30;
UseGPU            = true;
FlowThreshold     = 0.50;
CellprobThreshold = 0;

% -------- dFF --------
TossFrames = 30;      % toss BEFORE dFF (like your pipeline)
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

%% ============================== RUN BATCH==============================
finalPath = fullfile(folderPath, 'processed');
mkdir(finalPath);

ct = 0;
for k = 1:numel(list)
    if strcmp(list(k).name, 'processed'); continue;
    % subfolder
    elseif list(k).isdir
        if strcmp(list(k).name,'.') || strcmp(list(k).name,'..'); continue;
        else
            sublist = dir(fullfile(folderPath, list(k).name));
            for i = 1:numel(sublist)
                if strcmp(sublist(i).name,'.') || strcmp(sublist(i).name,'..'); continue;
                elseif endsWith(sublist(i).name, '.tif')   
                    filePath = fullfile(folderPath, list(k).name); 
                    filePath = fullfile(filePath, sublist(i).name); 
                    sprintf(filePath); 
                    subPath = fullfile(finalPath, list(k).name);
                    subPath = fullfile(subPath, erase(sublist(i).name, '.tif'));
                    % individual file subfolder
                    mkdir(subPath);
                    movefile(filePath, subPath);
                    target = fullfile(subPath, sublist(i).name);
    
                    %Run pipeline here
                    out = run_singleZ_pipeline_full( ...
                    target, pyExe, ...
                    "AutoDetectChannels",AutoDetectChannels, ...
                    "KeepChannelID",KeepChannelID, ...
                    "NumChannelsFallback",NumChannelsFallback, ...
                    "TrimL",TrimL, "TrimR",TrimR, ...
                    "DarkMaxFrames",DarkMaxFrames, "DarkSubsample",DarkSubsample, ...
                    "DarkSigma",DarkSigma, "ApplyDarkMask",ApplyDarkMask, ...
                    "SaveDiagnostics",SaveDiagnostics, "TossFramesForHist",TossFramesForHist, ...
                    "FPS",FPS, "Diameter",Diameter, ...
                    "UseGPU",UseGPU, "FlowThreshold",FlowThreshold, "CellprobThreshold",CellprobThreshold, ...
                    "TossFrames",TossFrames, ...
                    "TraceGain",TraceGain, "MinOffset",MinOffset, "OffsetMult",OffsetMult, ...
                    "DffScaleBar",DffScaleBar);

                    disp(out);
                    ct = ct+1;
                end
            end
            sprintf('%d videos processed under path %s', ct, sublist(i).name);
        end
    % no subfolder
    elseif ~list(k).isdir
        if strcmp(list(k).name,'.') || strcmp(list(k).name,'..'); continue;
        elseif endsWith(list(k).name, '.tif')
            filePath = fullfile(folderPath, list(k).name); 
            sprintf(filePath); 
            
            % individual file subfolder
            subPath = fullfile(finalPath, erase(list(k).name, '.tif'));
            mkdir(subPath);
            movefile(filePath, subPath);
            target = fullfile(subPath, list(k).name);

            %Run pipeline here
            out = run_singleZ_pipeline_full( ...
            target, pyExe, ...
            "AutoDetectChannels",AutoDetectChannels, ...
            "KeepChannelID",KeepChannelID, ...
            "NumChannelsFallback",NumChannelsFallback, ...
            "TrimL",TrimL, "TrimR",TrimR, ...
            "DarkMaxFrames",DarkMaxFrames, "DarkSubsample",DarkSubsample, ...
            "DarkSigma",DarkSigma, "ApplyDarkMask",ApplyDarkMask, ...
            "SaveDiagnostics",SaveDiagnostics, "TossFramesForHist",TossFramesForHist, ...
            "FPS",FPS, "Diameter",Diameter, ...
            "UseGPU",UseGPU, "FlowThreshold",FlowThreshold, "CellprobThreshold",CellprobThreshold, ...
            "TossFrames",TossFrames, ...
            "TraceGain",TraceGain, "MinOffset",MinOffset, "OffsetMult",OffsetMult, ...
            "DffScaleBar",DffScaleBar);

            disp(out);
            ct = ct+1;
        end
        sprintf('%d videos processed under path %s', ct, folderPath);
    end

end

%% ======================================================================
function out = run_singleZ_pipeline_full(filePath, pyExe, varargin)

p = inputParser;
p.addRequired("filePath", @(s)isstring(s)||ischar(s));
p.addRequired("pyExe", @(s)isstring(s)||ischar(s));

% channel
p.addParameter("AutoDetectChannels", true, @(x)islogical(x)&&isscalar(x));
p.addParameter("KeepChannelID", 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter("NumChannelsFallback", 1, @(x)isnumeric(x)&&isscalar(x)&&x>=1);

% trim
p.addParameter("TrimL", 10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter("TrimR", 10, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

% dark
p.addParameter("DarkMaxFrames", 500, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("DarkSubsample", 4, @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter("DarkSigma", 3, @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter("ApplyDarkMask", true, @(x)islogical(x)&&isscalar(x));
p.addParameter("SaveDiagnostics", true, @(x)islogical(x)&&isscalar(x));
p.addParameter("TossFramesForHist", 30, @(x)isnumeric(x)&&isscalar(x)&&x>=0);

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

[folderPath, baseName, ~] = fileparts(char(filePath));
folderPath = string(folderPath);
baseName   = string(baseName);

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

% If this TIFF is already extracted (e.g., *_ch1.tif), metadata may still say 2ch.
% Use a conservative sanity check: if pages not divisible -> treat as single-channel.
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

% --- define output stem ---
tag = "ch" + string(keepID_used);
if contains(lower(baseName), lower(tag))
    stem = baseName;
else
    stem = baseName + "_" + tag;
end

out_dark_tif = folderPath + filesep + stem + "_minusDark.tif";
out_dark_mat = folderPath + filesep + stem + "_darkInfo.mat";
out_dark_png = folderPath + filesep + stem + "_darkMask_sumProj.png";
out_hist_png = folderPath + filesep + stem + "_darkHist.png";

out_mc_mat   = folderPath + filesep + stem + "_MCinfo.mat";
out_mc_png   = folderPath + filesep + stem + "_MCshift.png";

% MC tif naming is not guaranteed; we infer it after running or from folder
out_mc_guess = folderPath + filesep + stem + "_minusDark_MC.tif"; % common guess

% cpSAM mat
out_sam_guess = folderPath + filesep + stem + "_minusDark_MC_cpSAM_output.mat"; % maybe
out_dff_mat   = folderPath + filesep + stem + "_dFF.mat";
out_stack_png = folderPath + filesep + stem + "_stackDFF.png";

% --- STEP 1: DARK (skip if exists) ---
if isfile(out_dark_tif) && isfile(out_dark_mat)
    fprintf("[SKIP] dark step (found): %s\n", out_dark_tif);
else
    fprintf("[RUN ] dark step\n");

    % load + trim
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
    fprintf("[data] raw pages=%d, kept pages=%d\n", size(F,3), size(F_keep,3));

    % estimate dark
    darkOut = estimate_dark_from_movie(F_keep, ...
        'MaxFrames', pr.DarkMaxFrames, ...
        'Subsample', pr.DarkSubsample);

    dark_mean = double(darkOut.dark_mean);
    dark_std  = double(darkOut.dark_std);

    % SAFE subtract
    Fs = single(F_keep) - single(dark_mean);
    Fs(Fs < 0) = 0;
    F_dark = uint16(Fs);

    % optional mask
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

    % diagnostics
    if pr.SaveDiagnostics
        sumProj = sum(single(F_dark), 3);
        fig = figure('Color','w');
        imagesc(sumProj); axis image off; colormap(gray); hold on;
        badMask = ~mask;
        [r,c] = find(badMask);
        if ~isempty(r), scatter(c, r, 2, 'r', 'filled'); end
        title(sprintf('Sum projection with masked pixels (%.2f%%)', 100*nnz(badMask)/numel(badMask)));
        exportgraphics(fig, out_dark_png, "Resolution", 200);
        close(fig);

        if exist("Fhist","file")==2
            if isscalar(dark_mean)
                F_sat = 65535 - dark_mean;
            else
                F_sat = 65535 - median(dark_mean(:));
            end

            F_plot = F_dark;
            if pr.TossFramesForHist > 0 && size(F_plot,3) > pr.TossFramesForHist
                F_plot(:,:,1:pr.TossFramesForHist) = [];
            end

            try
                Fhist(F_plot, sat=F_sat);
            catch
                Fhist(F_plot, F_sat);
            end
            figH = gcf; set(figH,"Color","w"); drawnow;
            exportgraphics(figH, char(out_hist_png), "Resolution", 200);
            close(figH);
        end
    end

    % save
    out_dark_tif_c = char(out_dark_tif);
    if exist(out_dark_tif_c, "file"), delete(out_dark_tif_c); end
    opts = struct('overwrite', true, 'message', false);
    saveastiff(F_dark, out_dark_tif_c, opts);

    save(char(out_dark_mat), "filePath","savedCh","nSaved","keepID_used","kWithin", ...
        "dark_mean","dark_std","thr","mask","darkOut","-v7.3");

    fprintf("Saved minusDark TIFF:\n%s\n", out_dark_tif);
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
    % last resort: look for newest *_minusDark*MC*.tif
    mc_tif = newest_matching_file(folderPath, stem + "_minusDark*MC*.tif");
    if mc_tif ~= "" && ~isfile(mc_tif), mc_tif = ""; end
end

if mc_tif ~= "" && isfile(mc_tif) && isfile(out_mc_mat)
    fprintf("[SKIP] MC step (found): %s\n", mc_tif);
else
    fprintf("[RUN ] MC step\n");
    mcOut = run_rigid_mc(out_dark_tif);

    mc_tif = infer_mc_tif(mcOut, folderPath, stem);
    assert(mc_tif ~= "" && isfile(mc_tif), "Could not locate MC tif output.");

    % shift plot
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
        title('|Rigid shift| per frame');
        exportgraphics(fig3, out_mc_png, "Resolution", 200);
        close(fig3);
    end

    save(char(out_mc_mat), "mcOut","mc_tif","sx","sy","r","-v7.3");
    fprintf("Saved MC info:\n%s\n", out_mc_mat);
end

% --- STEP 3: cpSAM (skip if cpSAM output exists) ---
sam_mat = "";
if isfile(out_sam_guess)
    sam_mat = out_sam_guess;
end
if sam_mat == "" || ~isfile(sam_mat)
    sam_mat = newest_matching_file(folderPath, stem + "*cpSAM_output.mat");
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

    sam_mat = infer_sam_mat(samOut, folderPath, stem);
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
    dFFout = dFF_RZ_dispatch(F_roi); % calls helper.dFF_RZ if available, else local
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

out.minusDark_tif = out_dark_tif;
out.darkInfo_mat  = out_dark_mat;
out.darkMask_png  = out_dark_png;
out.darkHist_png  = out_hist_png;

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
% F: [T x N] single/double ok
% If helper.dFF_RZ exists, use it. Otherwise use a robust local fallback.
if exist("helper.dFF_RZ","file")==2
    dFFout = helper.dFF_RZ(F);
else
    dFFout = dFF_RZ_local(F);
end
end

function out = dFF_RZ_local(F)
% Simple robust dF/F:
%   F0 per ROI = 20th percentile over time (robust baseline)
%   dFF = (F - F0) ./ max(F0, eps)
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

    % baseline removal per trace (visual)
    d = dFF;
    for i = 1:N
        d(:,i) = d(:,i) - median(d(:,i), "omitnan");
    end

    % gain (visual)
    d = gain * d;

    % spacing
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

    % scale bar (barAmp dFF; plotted height scales with gain)
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
% Return vector of saved channel IDs inferred from ScanImage metadata.
% Tries:
%   SI.hChannels.channelSave
%   SI.hChannels.channelsActive
% Interprets scalar as bitmask.

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
% Supports:
%   KEY = 3            (bitmask)
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

% B) explicit numeric list
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

guess = fullfile(folderPath, stem + "_minusDark_MC.tif");
if isfile(guess), mc_tif = string(guess); return; end

mc_tif = newest_matching_file(folderPath, stem + "_minusDark*MC*.tif");
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
