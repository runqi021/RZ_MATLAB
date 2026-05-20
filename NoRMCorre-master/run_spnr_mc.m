%% MC NoRMCorre — Sub-Pixel Non-Rigid (piecewise rigid)
function mcOut = run_spnr_mc(tiffPath, varargin)
% run_spnr_mc  Piecewise-rigid (non-rigid) NoRMCorre on a TIFF movie.
%
% Splits the FOV into overlapping patches. Each patch gets its own rigid
% shift, then shifts are smoothly interpolated across boundaries.
% Overlap prevents neurons at tile edges from being cut.
%
% Usage:
%   mcOut = run_spnr_mc('movie.tif');
%   mcOut = run_spnr_mc('movie.tif', 'MaxShift', 7, 'PatchSize', [128 128]);
%   mcOut = run_spnr_mc('movie.tif', 'Template', T);
%
% Inputs:
%   tiffPath : path to input TIFF movie
%   'MaxShift'   (optional) : max shift in pixels per patch (default 20)
%   'PatchSize'  (optional) : [rows cols] patch size in pixels (default [128 128])
%   'Overlap'    (optional) : [rows cols] overlap in pixels (default = PatchSize/2)
%   'MaxDev'     (optional) : [dy dx] max deviation between adjacent patches (default [5 5])
%   'UsFac'      (optional) : upsampling factor for subpixel (1=integer, default 1)
%   'Template'   (optional) : pre-computed template (e.g. from previous pass)
%   'InitBatch'  (optional) : frames for initial template (default 500)
%   'BinWidth'   (optional) : frames per batch (default 50)
%   'TossFrames' (optional) : toss first N frames before MC (default 0)
%
% Outputs (struct mcOut, also saved as *_SPNR_output.mat next to movie):
%   .mc_path     : path to saved motion-corrected TIFF (*_SPNR.tif)
%   .shifts      : NoRMCorre shifts struct (per-patch)
%   .template    : final template from NoRMCorre
%   .options     : NoRMCorre options used
%   .metrics     : struct with motion metrics (cY, cM, mY, mM, vY, vM)
%

% --- parse optional args ---
ip = inputParser;
ip.addRequired('tiffPath', @(s) isstring(s)||ischar(s));
ip.addParameter('MaxShift', 20, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('PatchSize', [128 128], @(x) isnumeric(x)&&numel(x)<=2);
ip.addParameter('Overlap', [], @(x) isnumeric(x)||isempty(x));
ip.addParameter('MaxDev', [5 5], @(x) isnumeric(x)&&numel(x)<=2);
ip.addParameter('UsFac', 1, @(x) isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('Template', [], @(x) isnumeric(x)||isempty(x));
ip.addParameter('InitBatch', 500, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('BinWidth', 50, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('TossFrames', 0, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
ip.parse(tiffPath, varargin{:});

maxShift   = ip.Results.MaxShift;
patchSize  = ip.Results.PatchSize;
if isscalar(patchSize), patchSize = [patchSize patchSize]; end
overlapIn  = ip.Results.Overlap;
if isempty(overlapIn)
    overlap = round(patchSize / 2);  % default 50% overlap
else
    overlap = overlapIn;
    if isscalar(overlap), overlap = [overlap overlap]; end
end
maxDev     = ip.Results.MaxDev;
if isscalar(maxDev), maxDev = [maxDev maxDev]; end
usFac      = round(ip.Results.UsFac);
templateIn = ip.Results.Template;
initBatch  = round(ip.Results.InitBatch);
binWidth   = round(ip.Results.BinWidth);
tossFrames = round(ip.Results.TossFrames);

% --- basic checks ---
tiffPath = char(string(tiffPath));
assert(exist(tiffPath,'file') == 2, 'File not found: %s', tiffPath);

% --- start parallel pool if available ---
try
    if isempty(gcp('nocreate'))
        gcp;
    end
catch
end

% --- read movie ---
fprintf('Reading movie: %s\n', tiffPath);
tic;
Y = read_file(tiffPath);
toc;

inClass = class(Y);

% --- toss first frames (shutter/PMT artifacts) ---
if tossFrames > 0 && size(Y,3) > tossFrames
    Y(:,:,1:tossFrames) = [];
    fprintf('Tossed first %d frames before MC (%d remaining)\n', tossFrames, size(Y,3));
end

Y = single(Y);
Y = Y - min(Y(:));
T = size(Y, ndims(Y));

% --- set non-rigid MC parameters ---
fprintf('SPNR MC: max_shift=%d px, patch=[%d %d], overlap=[%d %d], max_dev=[%d %d], us_fac=%d\n', ...
    maxShift, patchSize(1), patchSize(2), overlap(1), overlap(2), maxDev(1), maxDev(2), usFac);
fprintf('bin_width=%d, init_batch=%d\n', binWidth, initBatch);

options_nr = NoRMCorreSetParms( ...
    'd1',           size(Y,1), ...
    'd2',           size(Y,2), ...
    'grid_size',    patchSize, ...
    'overlap_pre',  overlap, ...
    'overlap_post', overlap, ...
    'max_shift',    maxShift, ...
    'max_dev',      maxDev, ...
    'us_fac',       usFac, ...
    'bin_width',    binWidth, ...
    'init_batch',   initBatch, ...
    'mot_uf',       4, ...          % upsample shifts 4x for smooth interpolation
    'correct_bidir', false);        % already handled by ScanImage

% --- run NoRMCorre (non-rigid) ---
if ~isempty(templateIn)
    fprintf('Running SPNR NoRMCorre with provided template...\n');
    tic;
    [M1, shifts1, template1, options_nr] = normcorre(Y, options_nr, templateIn);
    toc;
else
    fprintf('Running SPNR NoRMCorre...\n');
    tic;
    [M1, shifts1, template1, options_nr] = normcorre(Y, options_nr);
    toc;
end

% --- save corrected movie as TIFF first (protect expensive MC result) ---
[folder, base, ~] = fileparts(tiffPath);
folder = char(folder); base = char(base);
outname = char(fullfile(folder, [base '_SPNR.tif']));

mc = M1;
switch inClass
    case {'uint8','uint16','uint32'}
        mc(mc < 0) = 0;
        maxType = double(intmax(inClass));
        mc(mc > maxType) = maxType;
        mc = cast(mc, inClass);

    case {'int8','int16','int32'}
        minType = double(intmin(inClass));
        maxType = double(intmax(inClass));
        mc(mc < minType) = minType;
        mc(mc > maxType) = maxType;
        mc = cast(mc, inClass);

    otherwise
        % leave as single
end

fprintf('Saving SPNR corrected movie to:\n%s\n', outname);

Tsave = size(mc,3);
for t = 1:Tsave
    frOut = mc(:,:,t);
    if t == 1
        imwrite(frOut, outname, 'Compression', 'none');
    else
        imwrite(frOut, outname, 'WriteMode', 'append', 'Compression', 'none');
    end
end

% --- compute motion metrics ---
fprintf('Computing motion metrics...\n');
nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);

[cY, mY, vY]   = motion_metrics(Y,  10);
[cM1, mM1, vM1]= motion_metrics(M1, 10);

% --- pack outputs and save .mat ---
metrics = struct();
metrics.cY   = cY;
metrics.cM1  = cM1;
metrics.mY   = mY;
metrics.mM1  = mM1;
metrics.vY   = vY;
metrics.vM1  = vM1;

mcOut = struct();
mcOut.mc_path  = outname;
mcOut.shifts   = shifts1;
mcOut.template = template1;
mcOut.options  = options_nr;
mcOut.metrics  = metrics;
mcOut.qc_png   = '';

outMat = char(fullfile(folder, [base '_SPNR_output.mat']));
save(outMat, 'mcOut', '-v7.3');
fprintf('Saved SPNR metadata to:\n%s\n', outMat);

% --- QC figures (try-catch so a figure error never loses MC results) ---
try
    fig1 = figure;
        ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off;
              title('mean raw data','fontsize',14,'fontweight','bold');
        ax2 = subplot(2,2,2); imagesc(mM1,[nnY,mmY]); axis equal; axis tight; axis off;
              title('mean SPNR corrected','fontsize',14,'fontweight','bold');
        subplot(2,2,3); plot(1:T, cY, 1:T, cM1);
              legend('raw','SPNR'); title('correlation coefficients','fontsize',14,'fontweight','bold');
        subplot(2,2,4);
              scatter(cY, cM1); hold on;
              plot([0.9*min(cY),1.05*max(cM1)], [0.9*min(cY),1.05*max(cM1)], '--r');
              axis square;
              xlabel('raw','fontsize',14,'fontweight','bold');
              ylabel('SPNR','fontsize',14,'fontweight','bold');
        linkaxes([ax1,ax2],'xy');

    qcPng = char(fullfile(folder, [base '_SPNR_QC.png']));
    exportgraphics(fig1, qcPng, "Resolution", 200);
    close(fig1);

    % update mcOut with QC path and re-save
    mcOut.qc_png = qcPng;
    save(outMat, 'mcOut', '-v7.3');
catch ME
    fprintf('WARNING: QC figure failed (%s) — MC results are saved.\n', ME.message);
end

fprintf('SPNR MC complete.\n');

end
