%% Run_SVD_recon_260811.m
% Pixel-SVD decomposition + low-rank reconstruction VIDEO for one FOV,
% run on several versions of the same movie (raw / mcmc) so the spectra
% and the residuals are directly comparable.
%
% Updates over Run_SVD_working_251230.m + SVD_denoise_recon.m:
%   - processes a LIST of movie variants in one go, into separate output
%     folders (the old script overwrote all_modes\SVD_result.mat each run)
%   - explicit filename suffixes instead of fuzzy contains('_MC') discovery,
%     which in this folder matches the single-frame *_AVG_for_CP.tif
%   - front-drop frame alignment so every variant covers the SAME frames
%   - SVD via the temporal Gram matrix (exact, ~50x faster + far less RAM
%     than svd(Y,'econ') on a [190k x 3000] matrix)
%   - writes a real time-series MP4: original | rank-K recon | residual
%     (LINEAR contrast, shared limits) instead of 20-bin per-mode TIFFs
%
% Outputs go into <folderPath>\<outDirName>\<variant>\ .

clear; clc;

%% PATH SETUP (script lives in SVD/ subfolder)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'mod')));

%% USER PARAMS
folderPath = "D:\Ventral_surface_summary\Vgat\0730\cell1\roi1_3.2x_x1100y1050_z6_6000f_12lp_00001";

% Movie variants to run. suffix is appended to the FOV stem ('' = the raw
% ScanImage TIFF). tossFront = frames the pipeline had NOT yet removed.
variants = struct( ...
    'label',    {'raw',                    'mcmc'}, ...
    'suffix',   {'',                       '_ch1_preproc_MC_MC'}, ...
    'tossFront',{30,                       0});

outDirName     = 'svd_check_260811';

pad            = 36;      % crop border (pixels), 0 = no crop
nModes         = 200;     % SVD modes to keep / save
demeanPerPixel = true;    % subtract per-pixel mean before SVD

KreconList     = [1 10];        % rank(s): 1 = mode 1 alone, 10 = modes 1-10
videoWinSec    = [0 20];        % video time window (s); [] = whole movie
videoSpeed     = 1;             % 1 = real time
nGallery       = 12;            % modes shown in the eigenimage gallery
climPrc        = [1 99.5];      % LINEAR display percentiles (never gamma)

%% AUTO-DETECT FPS
[fps, ~] = detect_session_fps(folderPath);
fprintf('[SVD] Detected FPS = %g\n', fps);

folderPath = char(folderPath);

% NOT fileparts(): FOV folder names contain dots ('roi1_2.4x_...'), which
% fileparts would read as an extension and truncate.
pp      = strsplit(regexprep(folderPath, '[\\/]+$', ''), {'\', '/'});
fovStem = pp{end};
outRoot = fullfile(folderPath, outDirName);
if ~exist(outRoot, 'dir'); mkdir(outRoot); end

%% ---------------- LOAD ALL VARIANTS, THEN ALIGN FRAMES ----------------
nV = numel(variants);
Ystack = cell(nV,1);

for iv = 1:nV
    tifPath = fullfile(folderPath, [fovStem variants(iv).suffix '.tif']);
    assert(isfile(tifPath), 'Missing TIFF for variant "%s": %s', ...
           variants(iv).label, tifPath);

    fprintf('\n[SVD] === variant "%s" ===\n%s\n', variants(iv).label, tifPath);
    tLoad = tic;
    [Y, Hc, Wc] = load_tif_matrix(tifPath, pad);
    fprintf('[SVD] Loaded [%d px x %d frames] in %.1f s\n', ...
            size(Y,1), size(Y,2), toc(tLoad));

    assert(size(Y,2) > 100, ...
           'Variant "%s" has only %d frames — that is a projection, not a movie.', ...
           variants(iv).label, size(Y,2));

    nDrop = variants(iv).tossFront;
    if nDrop > 0
        Y = Y(:, nDrop+1:end);
        fprintf('[SVD] Dropped first %d frames -> %d\n', nDrop, size(Y,2));
    end

    variants(iv).tifPath = tifPath;
    variants(iv).Hc = Hc;
    variants(iv).Wc = Wc;
    Ystack{iv} = Y;
end

% front-drop alignment: keep the last T_common frames of every variant
Tall = cellfun(@(y) size(y,2), Ystack);
Tc   = min(Tall);
for iv = 1:nV
    if Tall(iv) > Tc
        Ystack{iv} = Ystack{iv}(:, end-Tc+1:end);
        fprintf('[SVD] "%s": truncated %d -> %d frames (front-drop align)\n', ...
                variants(iv).label, Tall(iv), Tc);
    end
end
tvec = (0:Tc-1)' / fps;
fprintf('\n[SVD] Common frame count across variants: %d (%.1f s)\n', Tc, Tc/fps);

% pixel dims must match or the comparison is meaningless
assert(numel(unique([variants.Hc])) == 1 && numel(unique([variants.Wc])) == 1, ...
       'Variants have different frame sizes after cropping.');

%% ---------------- PER-VARIANT SVD + FIGURES + VIDEOS ----------------
res = struct('label',{},'varExp',{},'singvals',{},'cumVar',{});

for iv = 1:nV
    lab = variants(iv).label;
    Hc  = variants(iv).Hc;
    Wc  = variants(iv).Wc;
    outDir = fullfile(outRoot, lab);
    if ~exist(outDir, 'dir'); mkdir(outDir); end

    Y = Ystack{iv};

    % ---- demean (in place; keep mu to add back for reconstruction) ----
    if demeanPerPixel
        mu = mean(Y, 2);
        Y  = Y - mu;
    else
        mu = zeros(size(Y,1), 1, 'single');
    end

    % ---- SVD through the temporal Gram matrix -------------------------
    % C = Y'Y is [T x T] (small). eig(C) gives V and S^2 exactly;
    % U = Y*V/S. Same answer as svd(Y,'econ') for the retained modes.
    fprintf('\n[SVD] "%s": Gram matrix [%d x %d]...\n', lab, Tc, Tc);
    tSVD = tic;
    C = Y' * Y;
    C = double((C + C')/2);
    [Vfull, Dfull] = eig(C);
    [lam, ord] = sort(max(diag(Dfull), 0), 'descend');
    K = min(nModes, sum(lam > 0));

    singvals = sqrt(lam);
    V = single(Vfull(:, ord(1:K)));
    S = single(singvals(1:K));
    U = (Y * V) ./ S(:)';                  % [P x K], unit-norm columns
    varExp = single(lam / sum(lam));
    fprintf('[SVD] "%s": %d modes in %.1f s | mode1 = %.1f%% var\n', ...
            lab, K, toc(tSVD), 100*varExp(1));

    % ---- save ---------------------------------------------------------
    SVD_result = struct();
    SVD_result.label          = lab;
    SVD_result.sourceTif      = variants(iv).tifPath;
    SVD_result.Hc             = Hc;
    SVD_result.Wc             = Wc;
    SVD_result.T              = Tc;
    SVD_result.fps            = fps;
    SVD_result.pad            = pad;
    SVD_result.tossFront      = variants(iv).tossFront;
    SVD_result.demeanPerPixel = demeanPerPixel;
    SVD_result.t              = tvec;
    SVD_result.mu             = mu;
    SVD_result.U              = U;
    SVD_result.S              = diag(S);           % kept square for back-compat
    SVD_result.V              = V;
    SVD_result.singvals       = S;
    SVD_result.varExp         = varExp(1:K);
    SVD_result.varExpFull     = single(varExp);    % full spectrum, for cumulatives

    outMat = fullfile(outDir, 'SVD_result.mat');
    save(outMat, 'SVD_result', '-v7.3');
    d = dir(outMat);
    fprintf('[SVD] Saved %s (%.0f MB)\n', outMat, d.bytes/1e6);

    % ---- singular spectrum -------------------------------------------
    plot_spectrum(varExp, lab, outDir);

    % ---- eigenimage / temporal-trace gallery --------------------------
    plot_mode_gallery(U, S, V, varExp, tvec, Hc, Wc, min(nGallery,K), lab, outDir);

    % ---- reconstruction videos ---------------------------------------
    for K0 = KreconList
        K0e = min(K0, K);
        write_recon_video(Y, mu, U, S, V, K0e, Hc, Wc, fps, ...
                          videoWinSec, videoSpeed, climPrc, lab, outDir);
    end

    res(end+1) = struct('label', lab, 'varExp', varExp, ...
                        'singvals', singvals, ...
                        'cumVar', cumsum(varExp)); %#ok<SAGROW>

    Ystack{iv} = [];   % free
    clear Y U V C Vfull Dfull;
end

%% ---------------- CROSS-VARIANT COMPARISON ----------------
plot_variant_comparison(res, outRoot, KreconList);

fprintf('\n[SVD] ---- summary ----\n');
for iv = 1:numel(res)
    fprintf('  %-6s  mode1 = %5.1f%%   cum@20 = %5.1f%%   cum@100 = %5.1f%%\n', ...
            res(iv).label, 100*res(iv).varExp(1), ...
            100*res(iv).cumVar(min(20,end)), 100*res(iv).cumVar(min(100,end)));
end
fprintf('\n[SVD] Done -> %s\n', outRoot);


%% ========================================================================
function [Y, Hc, Wc] = load_tif_matrix(tifPath, pad)
% Load a TIFF stack into [pixels x frames] single, cropped by pad.

    tw = Tiff(tifPath, 'r');
    cleaner = onCleanup(@() close(tw)); %#ok<NASGU>

    H0 = getTag(tw, 'ImageLength');
    W0 = getTag(tw, 'ImageWidth');

    % count directories
    T0 = 1;
    while ~lastDirectory(tw)
        nextDirectory(tw);
        T0 = T0 + 1;
    end

    rows = (1+pad):(H0-pad);
    cols = (1+pad):(W0-pad);
    Hc = numel(rows);
    Wc = numel(cols);

    Y = zeros(Hc*Wc, T0, 'single');
    for kk = 1:T0
        setDirectory(tw, kk);
        fr = read(tw);
        fr = fr(rows, cols);
        Y(:,kk) = single(fr(:));
    end
end


%% ========================================================================
function plot_spectrum(varExp, lab, outDir)
    nShow = min(200, numel(varExp));
    f = figure('Color','w','Visible','off','Position',[100 100 900 380]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing','compact', 'Padding','compact');

    ax1 = nexttile(tl);
    semilogy(ax1, 1:nShow, 100*varExp(1:nShow), 'o-', ...
             'LineWidth',1.1, 'MarkerSize',3, 'MarkerFaceColor','b', 'Color','b');
    xlabel(ax1,'mode index'); ylabel(ax1,'variance explained (%)');
    title(ax1,'spectrum'); grid(ax1,'on');

    ax2 = nexttile(tl);
    plot(ax2, 1:nShow, 100*cumsum(varExp(1:nShow)), 'k-', 'LineWidth',1.4);
    xlabel(ax2,'modes retained'); ylabel(ax2,'cumulative variance (%)');
    title(ax2,'cumulative'); grid(ax2,'on'); ylim(ax2,[0 100]);

    title(tl,sprintf('%s — pixel-SVD spectrum', lab), 'Interpreter','none');
    exportgraphics(f, fullfile(outDir,'singular_spectrum.png'), 'Resolution',200);
    exportgraphics(f, fullfile(outDir,'singular_spectrum.pdf'), 'ContentType','vector');
    close(f);
end


%% ========================================================================
function plot_mode_gallery(U, S, V, varExp, tvec, Hc, Wc, nG, lab, outDir)
% Rows of (eigenimage, temporal trace) pairs, 2 modes per row.

    nRows = ceil(nG/2);
    f = figure('Color','w','Visible','off','Position',[50 50 1500 260*nRows]);
    tl = tiledlayout(f, nRows, 6, 'TileSpacing','compact', 'Padding','compact');

    for k = 1:nG
        us = reshape(single(U(:,k)) * S(k), [Hc, Wc]);
        mx = prctile(abs(us(:)), 99);
        if mx == 0; mx = 1; end

        ax = nexttile(tl);
        imagesc(ax, imgaussfilt(us, 1.0));
        axis(ax, 'image'); axis(ax, 'off');
        caxis(ax, [-mx mx]); %#ok<CAXIS>  (clim() is R2022a+)
        colormap(ax, bluewhitered);
        title(ax, sprintf('mode %d  (%.2f%%)', k, 100*varExp(k)), 'FontSize',9);

        ax = nexttile(tl, [1 2]);
        plot(ax, tvec, V(:,k), 'k', 'LineWidth', 0.3);
        xlim(ax, [tvec(1) tvec(end)]);
        xlabel(ax, 'time (s)', 'FontSize',8);
        ylabel(ax, sprintf('v_{%d}', k), 'FontSize',8);
        grid(ax, 'on');
        set(ax, 'FontSize', 8);
    end

    title(tl,sprintf('%s — top %d modes: U_k S_k and v_k(t)', lab, nG), ...
            'Interpreter','none');
    exportgraphics(f, fullfile(outDir,'mode_gallery.png'), 'Resolution',180);
    close(f);
end


%% ========================================================================
function write_recon_video(Y, mu, U, S, V, K, Hc, Wc, fps, winSec, speed, ...
                           climPrc, lab, outDir)
% MP4: original | rank-K reconstruction | residual.  LINEAR contrast only.

    T = size(Y,2);
    if isempty(winSec)
        idx = 1:T;
    else
        i0  = max(1, round(winSec(1)*fps) + 1);
        i1  = min(T, round(winSec(2)*fps));
        idx = i0:i1;
    end

    Orig  = Y(:,idx) + mu;                                   % [P x n]
    Recon = U(:,1:K) * (S(1:K) .* V(idx,1:K)') + mu;
    Resid = Orig - Recon;

    % shared LINEAR limits from the original; residual symmetric about 0.
    % subsample for the percentiles — sorting 10^8 elements is pointless here.
    sub   = 1:37:numel(Orig);
    lims  = prctile(Orig(sub), climPrc);
    if lims(2) <= lims(1); lims(2) = lims(1) + 1; end
    rmax  = prctile(abs(Resid(sub)), 99);
    if rmax == 0; rmax = 1; end

    % H.264 needs even frame dimensions
    stripH = 26;
    padR   = mod(3*Wc, 2);
    padB   = mod(stripH + Hc, 2);
    strip  = make_label_strip({'original', sprintf('rank-%d recon', K), ...
                              'residual'}, Wc, stripH);
    if padR; strip = cat(2, strip, zeros(stripH,1,3,'uint8')); end

    vPath = fullfile(outDir, sprintf('recon_k%03d.mp4', K));
    vw = VideoWriter(vPath, 'MPEG-4');
    vw.FrameRate = max(1, fps*speed);
    vw.Quality   = 95;
    open(vw);
    cleaner = onCleanup(@() close(vw)); %#ok<NASGU>

    for n = 1:numel(idx)
        a = to_u8(reshape(Orig(:,n),  [Hc Wc]), lims(1), lims(2));
        b = to_u8(reshape(Recon(:,n), [Hc Wc]), lims(1), lims(2));
        c = to_u8(reshape(Resid(:,n), [Hc Wc]), -rmax,   rmax);
        frame = repmat(cat(2, a, b, c), [1 1 3]);
        if padR; frame = cat(2, frame, zeros(Hc,1,3,'uint8')); end
        frame = cat(1, strip, frame);
        if padB; frame = cat(1, frame, zeros(1,size(frame,2),3,'uint8')); end
        writeVideo(vw, frame);
    end

    fprintf('[SVD] "%s": wrote %s  (%d frames, %.1f s of data)\n', ...
            lab, vPath, numel(idx), numel(idx)/fps);
end


%% ========================================================================
function img = to_u8(x, lo, hi)
    img = uint8(255 * min(max((single(x) - lo) / (hi - lo), 0), 1));
end


%% ========================================================================
function strip = make_label_strip(labels, tileW, stripH)
% One-off rendered caption bar, concatenated on top of every video frame.

    nP = numel(labels);
    strip = zeros(stripH, nP*tileW, 3, 'uint8');   % fallback: plain black bar
    try
        f = figure('Visible','off','Color','k','Units','pixels', ...
                   'Position',[100 100 nP*tileW stripH]);
        cleaner = onCleanup(@() close(f)); %#ok<NASGU>
        ax = axes('Parent', f, 'Position',[0 0 1 1], 'Color','k', ...
                  'XLim',[0 nP], 'YLim',[0 1]);
        axis(ax, 'off');
        for i = 1:nP
            text(ax, i-0.5, 0.5, labels{i}, 'Color','w', 'FontSize',11, ...
                 'FontWeight','bold', 'HorizontalAlignment','center', ...
                 'VerticalAlignment','middle');
        end
        strip = imresize(getframe(f).cdata, [stripH, nP*tileW]);
    catch ME
        fprintf('[SVD] label strip render failed (%s) — using blank bar\n', ...
                ME.identifier);
    end
end


%% ========================================================================
function plot_variant_comparison(res, outRoot, KreconList)
    if numel(res) < 2; return; end

    f = figure('Color','w','Visible','off','Position',[100 100 1000 400]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing','compact', 'Padding','compact');
    cols = lines(numel(res));

    ax1 = nexttile(tl); hold(ax1,'on');
    ax2 = nexttile(tl); hold(ax2,'on');
    for iv = 1:numel(res)
        nShow = min(200, numel(res(iv).varExp));
        semilogy(ax1, 1:nShow, 100*res(iv).varExp(1:nShow), '-', ...
                 'LineWidth',1.4, 'Color',cols(iv,:), 'DisplayName',res(iv).label);
        plot(ax2, 1:nShow, 100*res(iv).cumVar(1:nShow), '-', ...
             'LineWidth',1.4, 'Color',cols(iv,:), 'DisplayName',res(iv).label);
    end
    set(ax1,'YScale','log');
    xlabel(ax1,'mode index'); ylabel(ax1,'variance explained (%)');
    title(ax1,'spectrum'); grid(ax1,'on'); legend(ax1,'Location','northeast');

    for K0 = KreconList
        xline(ax2, K0, 'k:', sprintf('k=%d',K0), 'HandleVisibility','off');
    end
    xlabel(ax2,'modes retained'); ylabel(ax2,'cumulative variance (%)');
    title(ax2,'cumulative'); grid(ax2,'on'); ylim(ax2,[0 100]);
    legend(ax2,'Location','southeast');

    title(tl,'pixel-SVD: variant comparison');
    exportgraphics(f, fullfile(outRoot,'variant_comparison.png'), 'Resolution',200);
    exportgraphics(f, fullfile(outRoot,'variant_comparison.pdf'), 'ContentType','vector');
    close(f);
    fprintf('[SVD] Saved variant comparison -> %s\n', ...
            fullfile(outRoot,'variant_comparison.png'));
end
