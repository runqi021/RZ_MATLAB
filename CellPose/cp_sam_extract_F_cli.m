function out = cp_sam_extract_F_cli(tiffPath, varargin)
% cp_sam_extract_F_cli  Use Cellpose-SAM via CLI on avg image, then extract F(t)
% Usage:
%   out = cp_sam_extract_F_cli('movie.tif', 'FPS',30, 'Diameter',26, ...
%                              'PythonExe','C:\Users\zhang\anaconda3\envs\cellpose-gpu\python.exe', ...
%                              'UseGPU',true, 'Savedir',[]);
%
% Outputs:
%   out.maskL, out.F (T x N), out.t, out.avgPath, out.maskPath, out.overlayPath

% -------------------- options --------------------
cand = string([
fullfile(getenv("USERPROFILE"), ".conda",     "envs", "cellpose-gpu", "python.exe")
fullfile(getenv("USERPROFILE"), "anaconda3",  "envs", "cellpose-gpu", "python.exe")
fullfile(getenv("USERPROFILE"), "miniconda3", "envs", "cellpose-gpu", "python.exe")
]);
pyExe = cand(find(isfile(cand), 1, "first"));
assert(pyExe ~= "", "Cellpose python not found. Checked:\n%s", strjoin(cand, newline));
fprintf("Using PythonExe:\n%s\n", pyExe);

p = inputParser;
addParameter(p,'FPS',30);
addParameter(p,'Diameter',[]);             % pixels across soma
addParameter(p,'FlowThreshold',0.4);
addParameter(p,'CellprobThreshold',0.0);
addParameter(p,'BatchSize',8);
addParameter(p,'UseGPU',true);
addParameter(p,'PythonExe', pyExe);  %'C:\Users\zhang\anaconda3\envs\cellpose-gpu\python.exe'
addParameter(p,'Savedir',[]);
parse(p,varargin{:});
opt = p.Results;

% -------------------- paths & info --------------------
assert(isfile(tiffPath),'File not found.');
[dirn,base,~] = fileparts(tiffPath);
if isempty(opt.Savedir), opt.Savedir = dirn; end
if ~isfolder(opt.Savedir), mkdir(opt.Savedir); end
avgPath  = fullfile(opt.Savedir, base + "_AVG_for_CP.tif");

info = imfinfo(tiffPath);
H = info(1).Height; W = info(1).Width; T = numel(info);

% -------------------- 1) build time-average image (streaming) --------------------
avgIm = zeros(H,W,'double');
for k = 1:T
    fr = imread(tiffPath, k, 'Info', info);
    avgIm = avgIm + double(fr);
end
avgIm = avgIm / T;

% save average as 16-bit
if isa(fr,'uint16')
    avg16 = uint16(avgIm);
    imwrite(avg16, avgPath);
else
    s = prctile(avgIm(:), [0.1 99.9]);
    avg16 = uint16( max(0, min(65535, (avgIm - s(1)) * 65535/max(1, s(2)-s(1)))) );
    imwrite(avg16, avgPath);
end

% -------------------- 2) run Cellpose-SAM via CLI --------------------
assert(isfile(opt.PythonExe),'PythonExe not found: %s',opt.PythonExe);

% Build command
% Notes:
%   --image_path <avgPath>
%   --pretrained_model cpsam
%   --save_tif         (writes *_cp_masks.tif)
%   --chan 0 --chan2 0 (grayscale)
%   --diameter <pixels> (omit if [])
%   --use_gpu          (if desired)
args = sprintf(' -m cellpose --image_path "%s" --pretrained_model cpsam --save_tif --chan 0 --chan2 0 --batch_size %d --flow_threshold %g --cellprob_threshold %g', ...
                avgPath, opt.BatchSize, opt.FlowThreshold, opt.CellprobThreshold);
if ~isempty(opt.Diameter)
    args = sprintf('%s --diameter %g', args, opt.Diameter);
end
if opt.UseGPU
    args = sprintf('%s --use_gpu', args);
end

cmd = sprintf('"%s"%s', opt.PythonExe, args);
fprintf('Running Cellpose-SAM CLI:\n%s\n', cmd);
[status, cmdout] = system(cmd);
if status ~= 0
    error('Cellpose CLI failed.\nCommand:\n%s\n\nOutput:\n%s', cmd, cmdout);
end

% -------------------- 3) find mask file --------------------
d = dir(fullfile(opt.Savedir, base + "*cp_masks*.tif"));
if isempty(d)
    % sometimes CP writes into a subfolder named after avg or into ~/.cellpose; try broader search
    d = dir(fullfile(fileparts(avgPath), "*cp_masks*.tif"));
end
assert(~isempty(d), 'Could not find any *_cp_masks*.tif after CLI run.');
[~,idx] = max([d.datenum]);
maskPath = fullfile(d(idx).folder, d(idx).name);

% -------------------- 4) load mask & extract F(t) --------------------
maskL = uint16(imread(maskPath));
assert(all(size(maskL)==[H W]), 'Mask size (%dx%d) != movie frame size (%dx%d).', size(maskL,1), size(maskL,2), H, W);

labels = unique(maskL); labels(labels==0)=[];
N = numel(labels);

mvec   = double(maskL(:));
valid  = mvec>0;
lab    = mvec(valid);
maxLab = double(max(mvec));
npix   = accumarray(lab, 1, [maxLab 1], @sum, 0);

F = zeros(T, N, 'double');
for k = 1:T
    fr = double(imread(tiffPath, k, 'Info', info));
    F(k,:) = accumarray(lab, fr(valid), [maxLab 1], @sum, 0);
end
F = F(:,labels) ./ npix(labels)';

t = (0:T-1)'/opt.FPS;

% -------------------- 5) save overlay for QC --------------------
% Make periphery-only ROI mask movie (white outline + ROI #)
% Assumes in workspace:
%   tiffPath : path to the MC movie (e.g. roi5_..._MC.tif)
%   maskL    : [H x W] labeled mask from SAM output
% -------- A) Per-ROI periphery (avoid union-style merging) --------
labels = unique(maskL);
labels(labels == 0) = [];
N_roi = numel(labels);

edge_all = false(H, W);
for ii = 1:N_roi
    mk = (maskL == labels(ii));
    if any(mk(:))
        edge_all = edge_all | bwperim(mk);
    end
end
% -------- B) ROI centroids + static label image (once) --------
props = regionprops(maskL, 'Centroid');
centers = reshape([props.Centroid], 2, []).';
label_str = arrayfun(@(k) sprintf('%d', k), 1:N_roi, 'UniformOutput', false);

labelMask = false(H, W);

tmpRGB = zeros(H, W, 3, 'uint8');
tmpRGB = insertText(tmpRGB, centers, label_str, ...
                    'FontSize', 14, ...
                    'BoxOpacity', 0, ...
                    'TextColor', 'white');
labelGray = rgb2gray(tmpRGB) > 0;   % boolean mask of text positions
labelMask = labelGray;
% -------- C) Output path --------
[folder, base, ~] = fileparts(tiffPath);
ovPath1 = fullfile(folder, base + "_AVG_ROImask.tif");
ovPath2 = fullfile(folder, base + "_AVG_ROIlabel.tif");
%fprintf('Writing: %s\n', ovPath1);
% -------- D) SUM Proj and mask --------
WHITE = uint16(65535);
avg16(edge_all) = WHITE;         % uint16 white

imwrite(avg16, ovPath1, 'Compression', 'none');

if any(labelMask(:))
    avg16(labelMask) = WHITE;
end

imwrite(avg16, ovPath2, 'Compression', 'none');

% -------------------- save & return --------------------
csvPath = fullfile(opt.Savedir, base + "_F.csv");
hdr = "ROI_" + string(labels(:))';
fid = fopen(csvPath,'w'); fprintf(fid, '%s\n', strjoin(hdr, ',')); fclose(fid);
dlmwrite(csvPath, F, '-append');

save(fullfile(opt.Savedir, base+"_cpSAM_output.mat"), 'F','t','maskL','avgPath','maskPath','ovPath1', 'ovPath2','opt','-v7.3');

out.maskL      = maskL;
out.F          = F;
out.t          = t;
out.avgPath    = avgPath;
out.maskPath   = maskPath;
out.overlayMaskPath  = ovPath1;
out.overlayLabelPath = ovPath2;
out.params     = opt;

fprintf('Done.\n  Avg:   %s\n  Mask:  %s\n  F.csv: %s\n  QC:    %s\n', avgPath, maskPath, csvPath, ovPath1);
end
