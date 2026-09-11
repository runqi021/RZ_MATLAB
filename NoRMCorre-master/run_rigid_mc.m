%% MC NoRMCorre
function mcOut = run_rigid_mc(tiffPath, varargin)
% run_rigid_mc  Rigid, integer-pixel NoRMCorre on a TIFF movie.
%
% Usage:
%   mcOut = run_rigid_mc('movie.tif');
%   mcOut = run_rigid_mc('movie.tif', 'MaxShift', 7);
%   mcOut = run_rigid_mc('movie.tif', 'MaxShift', 7, 'Template', T);
%
% Inputs:
%   tiffPath : path to input TIFF movie
%   'MaxShift' (optional) : max shift in pixels (default 20)
%   'Template' (optional) : pre-computed template for registration
%                           (e.g. from a previous pass)
%
% Outputs (struct mcOut, also saved as *_MC_output.mat next to movie):
%   .mc_path     : path to saved motion-corrected TIFF (*_MC.tif)
%   .shifts      : NoRMCorre shifts struct
%   .template    : final template from NoRMCorre
%   .options     : NoRMCorre options used
%   .metrics     : struct with motion metrics (cY, cM, mY, mM, vY, vM)
%

% --- parse optional args ---
ip = inputParser;
ip.addRequired('tiffPath', @(s) isstring(s)||ischar(s));
ip.addParameter('MaxShift', 20, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('Template', [], @(x) isnumeric(x)||isempty(x));
ip.addParameter('InitBatch', 500, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('BinWidth', 50, @(x) isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('TossFrames', 0, @(x) isnumeric(x)&&isscalar(x)&&x>=0);
ip.parse(tiffPath, varargin{:});
maxShift   = ip.Results.MaxShift;
templateIn = ip.Results.Template;
initBatch  = round(ip.Results.InitBatch);
binWidth   = round(ip.Results.BinWidth);
tossFrames = round(ip.Results.TossFrames);

% --- basic checks ---
tiffPath = char(string(tiffPath));
assert(exist(tiffPath,'file') == 2, 'File not found: %s', tiffPath);

% --- parallel pool ---
% Rigid MC runs with use_parallel=false, so NoRMCorre does NOT use the pool
% here. A 12-worker local pool needlessly holds several GB during registration
% -- the exact moment long movies (e.g. 18000 frames) hit the RAM ceiling and
% throw "Out of memory". So do NOT auto-start a pool; an existing one is left
% untouched (still unused by rigid MC).

% --- read movie ---
fprintf('Reading movie: %s\n', tiffPath);
tic;
Y = read_file(tiffPath);   % NoRMCorre helper
toc;

inClass = class(Y);

% --- toss first frames (shutter/PMT artifacts) ---
if tossFrames > 0 && size(Y,3) > tossFrames
    Y(:,:,1:tossFrames) = [];
    fprintf('Tossed first %d frames before MC (%d remaining)\n', tossFrames, size(Y,3));
end

Y   = single(Y);
Y   = Y - min(Y(:));       % shift to start at 0
T   = size(Y, ndims(Y));

% --- set rigid, integer-shift MC parameters ---
fprintf('max_shift=%d px, bin_width=%d, init_batch=%d\n', maxShift, binWidth, initBatch);
options_rigid = NoRMCorreSetParms( ...
    'd1',        size(Y,1), ...
    'd2',        size(Y,2), ...
    'bin_width', binWidth, ...
    'max_shift', maxShift, ...
    'us_fac',    1, ... %20, ...   <--- integer-pixel shifts (no subpixel)
    'init_batch',initBatch);

% --- run NoRMCorre (rigid) ---
if ~isempty(templateIn)
    fprintf('Running rigid NoRMCorre with provided template...\n');
    tic;
    [M1, shifts1, template1, options_rigid] = normcorre(Y, options_rigid, templateIn);
    toc;
else
    fprintf('Running rigid NoRMCorre...\n');
    tic;
    [M1, shifts1, template1, options_rigid] = normcorre(Y, options_rigid);
    toc;
end

% --- output path ---
[folder, base, ext] = fileparts(tiffPath);
folder = char(folder); base = char(base); ext = char(ext);
outname = char(fullfile(folder, [base '_MC' ext]));

% --- motion metrics (QC-only) -------------------------------------------
% motion_metrics runs corr() on a full [pixels x frames] matrix, which builds
% several whole-movie copies and exhausts RAM on long movies (this is where
% 12000/18000-frame runs OOMed / hung in corr). Metrics only feed the QC plot,
% so: (a) subsample frames when T is large, (b) free Y before corr(). The saved
% MC movie is unaffected (full resolution, all frames).
fprintf('Computing motion metrics...\n');
mstride = max(1, ceil(T/3000));            % feed <=~3000 frames into corr()
tcorr   = 1:mstride:T;                      % frame indices used for the QC plot

% display range from a cheap linear-index subsample (avoids quantile(Y(:)) copy)
ssv = Y(1:max(1,round(numel(Y)/5e6)):end);
nnY = quantile(ssv, 0.005);  mmY = quantile(ssv, 0.995);  clear ssv;

if mstride > 1
    Ysub = Y(:,:,tcorr);                   % small copy (~3 GB)
    clear Y;                               % free raw movie (~19 GB) before corr()
    [cY, mY, vY] = motion_metrics(Ysub, 10);
    clear Ysub;
    [cM1, mM1, vM1] = motion_metrics(M1(:,:,tcorr), 10);
else
    [cY,  mY,  vY ] = motion_metrics(Y,  10);
    clear Y;
    [cM1, mM1, vM1] = motion_metrics(M1, 10);
end

% --- save corrected movie (chunked cast + BigTIFF append) ---------------
% Cast per chunk instead of materializing a whole-movie 'mc' duplicate, and
% write via saveastiff (append+big) so long movies stream out as BigTIFF with
% bounded RAM and no classic-TIFF 4GB cap. Small movies still produce a normal
% (auto classic-vs-big) TIFF.
fprintf('Saving motion-corrected movie to:\n%s\n', outname);
if exist('saveastiff','file') == 2
    Twrite   = size(M1,3);
    wchunk   = 1000;                 % frames per write block (~1 GB at 512x512 single)
    firstBlk = true;
    for c0 = 1:wchunk:Twrite
        c1  = min(c0 + wchunk - 1, Twrite);
        blk = clip_cast_block(M1(:,:,c0:c1), inClass);
        saveastiff(blk, outname, struct('append', ~firstBlk, 'big', true, ...
            'message', false, 'overwrite', firstBlk));
        firstBlk = false;
    end
    clear blk;
else
    % Fallback: whole-array cast + classic imwrite append (4 GB-limited).
    mc = clip_cast_block(M1, inClass);
    for t = 1:size(mc,3)
        if t == 1
            imwrite(mc(:,:,t), outname, 'Compression', 'none');
        else
            imwrite(mc(:,:,t), outname, 'WriteMode', 'append', 'Compression', 'none');
        end
    end
    clear mc;
end

% --- pack outputs and save .mat ---
metrics = struct();
metrics.cY   = cY;
metrics.cM1  = cM1;
metrics.mY   = mY;
metrics.mM1  = mM1;
metrics.vY   = vY;
metrics.vM1  = vM1;

mcOut = struct();
mcOut.mc_path   = outname;
mcOut.shifts    = shifts1;
mcOut.template  = template1;
mcOut.options   = options_rigid;
mcOut.metrics   = metrics;
mcOut.qc_png    = '';
mcOut.shift_png = '';

outMat = char(fullfile(folder, [base '_MC_output.mat']));
save(outMat, 'mcOut', '-v7.3');
fprintf('Saved MC metadata to:\n%s\n', outMat);

% --- QC figures (try-catch so a figure error never loses MC results) ---
try
    fig1 = figure;
        ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off;
              title('mean raw data','fontsize',14,'fontweight','bold');
        ax2 = subplot(2,2,2); imagesc(mM1,[nnY,mmY]); axis equal; axis tight; axis off;
              title('mean rigid corrected','fontsize',14,'fontweight','bold');
        subplot(2,2,3); plot(tcorr, cY, tcorr, cM1);
              legend('raw','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold');
        subplot(2,2,4);
              scatter(cY, cM1); hold on;
              plot([0.9*min(cY),1.05*max(cM1)], [0.9*min(cY),1.05*max(cM1)], '--r');
              axis square;
              xlabel('raw','fontsize',14,'fontweight','bold');
              ylabel('rigid','fontsize',14,'fontweight','bold');
        linkaxes([ax1,ax2],'xy');

    qcPng = char(fullfile(folder, [base '_MC_QC.png']));
    exportgraphics(fig1, qcPng, "Resolution", 200);
    close(fig1);

    % shifts plot
    shifts_r = squeeze(cat(3,shifts1(:).shifts));
    fig2 = figure;
        ax1 = subplot(2,1,1); hold on;
              plot(shifts_r(:,1),'--k','linewidth',2);
              title('displacements along x','fontsize',14,'fontweight','bold');
              set(gca,'Xtick',[]);
        ax2 = subplot(2,1,2); hold on;
              plot(shifts_r(:,2),'--k','linewidth',2);
              title('displacements along y','fontsize',14,'fontweight','bold');
              xlabel('timestep','fontsize',14,'fontweight','bold');
        linkaxes([ax1,ax2],'x');

    shiftPng = char(fullfile(folder, [base '_MC_shifts.png']));
    exportgraphics(fig2, shiftPng, "Resolution", 200);
    close(fig2);

    % update mcOut with QC paths and re-save
    mcOut.qc_png    = qcPng;
    mcOut.shift_png = shiftPng;
    save(outMat, 'mcOut', '-v7.3');
catch ME
    fprintf('WARNING: QC figure failed (%s) — MC results are saved.\n', ME.message);
end

fprintf('Rigid MC complete.\n');

end

% ------------------------------------------------------------------------
function blk = clip_cast_block(blk, inClass)
% Clip out-of-range values then cast a single-precision block back to the
% original integer class (float classes are left as single).
switch inClass
    case {'uint8','uint16','uint32'}
        blk(blk < 0) = 0;
        mx = double(intmax(inClass)); blk(blk > mx) = mx;
        blk = cast(blk, inClass);
    case {'int8','int16','int32'}
        mn = double(intmin(inClass)); blk(blk < mn) = mn;
        mx = double(intmax(inClass)); blk(blk > mx) = mx;
        blk = cast(blk, inClass);
    otherwise
        % leave as single
end
end