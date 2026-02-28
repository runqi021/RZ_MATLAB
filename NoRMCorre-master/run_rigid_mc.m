%% MC NoRMCorre
function mcOut = run_rigid_mc(tiffPath)
% run_rigid_mc  Rigid, integer-pixel NoRMCorre on a TIFF movie.
%
% Usage:
%   mcOut = run_rigid_mc('movie.tif');
%
% Inputs:
%   tiffPath : path to input TIFF movie
%
% Outputs (struct mcOut, also saved as *_MC_output.mat next to movie):
%   .mc_path     : path to saved motion-corrected TIFF (*_MC.tif)
%   .shifts      : NoRMCorre shifts struct
%   .template    : final template from NoRMCorre
%   .options     : NoRMCorre options used
%   .metrics     : struct with motion metrics (cY, cM, mY, mM, vY, vM)
%

% --- basic checks ---
if nargin < 1 || isempty(tiffPath)
    error('Please provide a TIFF path, e.g. mcOut = run_rigid_mc(''movie.tif'');');
end

if isstring(tiffPath)
    tiffPath = char(tiffPath);
end

assert(exist(tiffPath,'file') == 2, 'File not found: %s', tiffPath);

% --- start parallel pool if available ---
try
    if isempty(gcp('nocreate'))
        gcp;    % start default pool
    end
catch
    % no parallel toolbox, ignore
end

% --- read movie ---
fprintf('Reading movie: %s\n', tiffPath);
tic;
Y = read_file(tiffPath);   % NoRMCorre helper
toc;

inClass = class(Y);

Y   = single(Y);
Y   = Y - min(Y(:));       % shift to start at 0
T   = size(Y, ndims(Y));

% --- set rigid, integer-shift MC parameters ---
options_rigid = NoRMCorreSetParms( ...
    'd1',        size(Y,1), ...
    'd2',        size(Y,2), ...
    'bin_width', 200, ...
    'max_shift', 20, ...
    'us_fac',    1, ... %20, ...   <--- integer-pixel shifts (no subpixel)
    'init_batch',200);

% --- run NoRMCorre (rigid) ---
fprintf('Running rigid NoRMCorre...\n');
tic;
[M1, shifts1, template1, options_rigid] = normcorre(Y, options_rigid);
toc;

% --- compute motion metrics (optional but useful) ---
fprintf('Computing motion metrics...\n');
nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);

[cY, mY, vY]   = motion_metrics(Y,  10);
[cM1, mM1, vM1]= motion_metrics(M1, 10);

% --- (optional) quick QC figures (comment out if you don't want plots) ---
figure;
    ax1 = subplot(2,2,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off;
          title('mean raw data','fontsize',14,'fontweight','bold');
    ax2 = subplot(2,2,2); imagesc(mM1,[nnY,mmY]); axis equal; axis tight; axis off;
          title('mean rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,2,3); plot(1:T, cY, 1:T, cM1);
          legend('raw','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold');
    subplot(2,2,4);
          scatter(cY, cM1); hold on;
          plot([0.9*min(cY),1.05*max(cM1)], [0.9*min(cY),1.05*max(cM1)], '--r');
          axis square;
          xlabel('raw','fontsize',14,'fontweight','bold');
          ylabel('rigid','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2],'xy');

% shifts plot
shifts_r = squeeze(cat(3,shifts1(:).shifts));
figure;
    ax1 = subplot(2,1,1); hold on;
          plot(shifts_r(:,1),'--k','linewidth',2);
          title('displacements along x','fontsize',14,'fontweight','bold');
          set(gca,'Xtick',[]);
    ax2 = subplot(2,1,2); hold on;
          plot(shifts_r(:,2),'--k','linewidth',2);
          title('displacements along y','fontsize',14,'fontweight','bold');
          xlabel('timestep','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2],'x');

% --- save rigid-corrected movie as 16-bit TIFF next to input ---
[folder, base, ext] = fileparts(tiffPath);
outname = fullfile(folder, [base '_MC' ext]);

mc = M1;
% clamp and cast to uint16
%mc(mc < 0)      = 0;
%mc(mc > 65535)  = 65535;
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
        % leave as single (or whatever) if it was float
        % If you want strictly 'single', uncomment:
        % V_masked = single(V_masked);
end

fprintf('Saving motion-corrected movie to:\n%s\n', outname);

T = size(mc,3);

for t = 1:T
    frOut = mc(:,:,t);

    if t == 1
        imwrite(frOut, outname, 'Compression', 'none');
    else
        imwrite(frOut, outname, 'WriteMode', 'append', 'Compression', 'none');
    end
end

% --- pack outputs and save MC_output .mat ---
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
mcOut.options  = options_rigid;
mcOut.metrics  = metrics;

outMat = fullfile(folder, [base '_MC_output.mat']);
save(outMat, 'mcOut', '-v7.3');

fprintf('Saved MC metadata to:\n%s\n', outMat);
fprintf('Rigid MC complete.\n');

end