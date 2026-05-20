function out = cp_sam_extract_F_cli_finetuned(tiffPath, varargin)
% cp_sam_extract_F_cli_finetuned  Run fine-tuned CellPose, save one overlay PNG.

% -------------------- options --------------------
cand = string([
fullfile(getenv("USERPROFILE"), ".conda",     "envs", "cellpose-gpu", "python.exe")
fullfile(getenv("USERPROFILE"), "anaconda3",  "envs", "cellpose-gpu", "python.exe")
fullfile(getenv("USERPROFILE"), "miniconda3", "envs", "cellpose-gpu", "python.exe")
]);
pyExe = cand(find(isfile(cand), 1, "first"));
assert(pyExe ~= "", "Cellpose python not found. Checked:\n%s", strjoin(cand, newline));

p = inputParser;
addParameter(p,'ModelPath','');
addParameter(p,'DiameterUm',16);
addParameter(p,'FlowThreshold',0.5);
addParameter(p,'CellprobThreshold',-0.5);
addParameter(p,'BatchSize',8);
addParameter(p,'UseGPU',true);
addParameter(p,'PythonExe', pyExe);
addParameter(p,'Savedir',[]);
parse(p,varargin{:});
opt = p.Results;

assert(~isempty(opt.ModelPath) && isfile(opt.ModelPath), ...
    'ModelPath not found: %s', opt.ModelPath);

% -------------------- paths & info --------------------
assert(isfile(tiffPath),'File not found: %s', tiffPath);
[dirn,base,~] = fileparts(tiffPath);
if isempty(opt.Savedir), opt.Savedir = dirn; end
if ~isfolder(opt.Savedir), mkdir(opt.Savedir); end

% Temp dir for cellpose intermediate files
tmpDir = fullfile(tempdir, 'cp_ft_tmp');
if ~isfolder(tmpDir), mkdir(tmpDir); end
avgPath = fullfile(tmpDir, base + "_AVG_for_ftCP.tif");

info = imfinfo(tiffPath);
H = info(1).Height; W = info(1).Width; T = numel(info);

% -------------------- detect pixel size --------------------
[~, scan_meta] = detect_session_fps(dirn);
if isfield(scan_meta, 'pixelSize_um') && isfinite(scan_meta.pixelSize_um)
    pixelSize_um = scan_meta.pixelSize_um;
else
    pixelSize_um = 1.7778;
    fprintf(2, 'WARNING: pixelSize_um not found. Using fallback %.4f um/px.\n', pixelSize_um);
end
diameterPx = round(opt.DiameterUm / pixelSize_um);
fprintf('Diameter: %d um / %.4f um/px = %d px\n', opt.DiameterUm, pixelSize_um, diameterPx);

% -------------------- 1) build time-average --------------------
avgIm = zeros(H,W,'double');
for k = 1:T
    fr = imread(tiffPath, k, 'Info', info);
    avgIm = avgIm + double(fr);
end
avgIm = avgIm / T;
imwrite(uint16(avgIm), avgPath);

% -------------------- 2) run fine-tuned CellPose --------------------
args = sprintf(' -m cellpose --image_path "%s" --pretrained_model "%s" --save_tif --chan 0 --chan2 0 --batch_size %d --flow_threshold %g --cellprob_threshold %g --diameter %d', ...
                avgPath, opt.ModelPath, opt.BatchSize, opt.FlowThreshold, opt.CellprobThreshold, diameterPx);
if opt.UseGPU, args = sprintf('%s --use_gpu', args); end

cmd = sprintf('"%s"%s', opt.PythonExe, args);
fprintf('Running CellPose CLI:\n%s\n', cmd);
[status, cmdout] = system(cmd);
if status ~= 0
    error('CellPose CLI failed.\n%s', cmdout);
end

% -------------------- 3) load mask --------------------
d = dir(fullfile(tmpDir, "*cp_masks*.tif"));
assert(~isempty(d), 'No mask file found after CLI run.');
[~,idx] = max([d.datenum]);
maskPath = fullfile(d(idx).folder, d(idx).name);
maskL = uint16(imread(maskPath));

% Clean up temp files
delete(fullfile(tmpDir, '*'));

labels = unique(maskL); labels(labels==0) = [];
N_roi = numel(labels);

% -------------------- 4) write overlay PNG directly --------------------
% Clip 0.5-99.5%, gamma 0.6
lo = prctile(avgIm(:), 0.5);
hi = prctile(avgIm(:), 99.5);
if hi <= lo, hi = lo + eps; end
fovImg = (avgIm - lo) / (hi - lo);
fovImg = max(0, min(1, fovImg));
fovImg = fovImg .^ 0.6;

% Build RGB with yellow ROI outlines
fovRGB = repmat(fovImg, [1 1 3]);
for ri = 1:N_roi
    perim_ri = bwperim(maskL == labels(ri));
    for ch = 1:3
        plane = fovRGB(:,:,ch);
        plane(perim_ri) = [1 1 0] * (ch == [1;2;3]);
        fovRGB(:,:,ch) = plane;
    end
end

% Burn in ROI number labels
for ri = 1:N_roi
    [ry_, rx_] = find(maskL == labels(ri));
    if isempty(ry_), continue; end
    cx = mean(rx_); cy = mean(ry_);
    fovRGB = insertText(fovRGB, [cx cy], sprintf('%d',ri), ...
        'FontSize', 10, 'BoxOpacity', 0, ...
        'TextColor', 'yellow', 'AnchorPoint', 'Center');
end

% Scale bar 50 um
sb_px = round(50 / pixelSize_um);
sb_y = H - 10;
fovRGB(sb_y-1:sb_y, 10:10+sb_px-1, :) = 1;

[~, base2, ~] = fileparts(tiffPath);
outPng = fullfile(opt.Savedir, base2 + "_ftCP_overlay.png");
imwrite(fovRGB, outPng);

% -------------------- return --------------------
out.maskL      = maskL;
out.nROI       = N_roi;
out.overlayPng = outPng;
out.pixelSize_um = pixelSize_um;
out.diameterPx   = diameterPx;

fprintf('Done. %d ROIs  ->  %s\n', N_roi, outPng);
end
