%% Run_SVD_working_251230.m
% Run SVD on a calcium imaging TIFF to visualize motion artifacts.
%
% Give a data subfolder path — the script auto-discovers the TIFF,
% reads FPS from metadata, and saves SVD results + per-mode PDFs
% into an all_modes/ subfolder alongside the data.

clear; clc;

%% PATH SETUP (script lives in SVD/ subfolder)
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));
addpath(genpath(fullfile(repoRoot, 'mod')));

%% USER PARAMS
folderPath      = "D:\251124_live_vglut2_soma_g8s+cy5\phys\IO\roi3_-1000-300-80_48lp_930_x1.4_256x512_3000f_00001";   % <-- paste data subfolder path here

demeanPerPixel  = true;
pad             = 36;        % crop border (pixels), 0 = no crop
nModes          = 200;      % number of SVD modes to save (max for scree/recon)
alpha_polarity  = 0.05;     % significance level for polarity test
savePDF         = false;     % set true to also dump per-mode PDFs + binned movies

% TIFF preference: 'mc' = motion-corrected, 'preproc' = preprocessed, 'raw' = original
tifPref         = 'mc';

%% AUTO-DETECT FPS
[fps, scan_meta] = detect_session_fps(folderPath);
fprintf('[SVD] Detected FPS = %d\n', fps);

%% DISCOVER TIFF
tifPath = discover_tif(folderPath, tifPref);
fprintf('[SVD] Using TIFF: %s\n', tifPath);

% Determine TossFrames: MC TIFFs already had frames tossed, raw/preproc may not
[~, tifStem] = fileparts(tifPath);
if contains(tifStem, '_MC', 'IgnoreCase', true)
    dropFirstN = 0;   % MC pipeline already tossed initial frames
    fprintf('[SVD] MC TIFF detected — no additional frame toss\n');
else
    dropFirstN = 30;  % raw/preproc: toss stabilization frames
    fprintf('[SVD] Non-MC TIFF — tossing first %d frames\n', dropFirstN);
end

%% OUTPUT DIRECTORY
modesDir = fullfile(folderPath, 'all_modes');
if ~exist(modesDir, 'dir')
    mkdir(modesDir);
    fprintf('[SVD] Created output folder: %s\n', modesDir);
end

%% RUN SVD
run_svd_for_file(tifPath, savePDF, modesDir, fps, ...
                 pad, dropFirstN, demeanPerPixel, nModes, alpha_polarity);

fprintf('\n[SVD] Done.\n');


%% ========================================================================
function tifPath = discover_tif(folderPath, pref)
% Find the best TIFF in folderPath based on preference.
%   pref: 'mc' | 'preproc' | 'raw'

    assert(isfolder(folderPath), 'Folder does not exist: %s', folderPath);
    allTifs = dir(fullfile(folderPath, '*.tif'));
    assert(~isempty(allTifs), 'No .tif files found in: %s', folderPath);

    names = {allTifs.name};

    % Classify TIFFs
    is_mc      = cellfun(@(n) contains(n, '_MC', 'IgnoreCase', true), names);
    is_preproc = cellfun(@(n) contains(n, '_preproc', 'IgnoreCase', true), names);
    is_raw     = ~is_mc & ~is_preproc;

    switch lower(pref)
        case 'mc'
            order = [find(is_mc), find(is_preproc), find(is_raw)];
        case 'preproc'
            order = [find(is_preproc), find(is_mc), find(is_raw)];
        case 'raw'
            order = [find(is_raw), find(is_preproc), find(is_mc)];
        otherwise
            order = 1:numel(allTifs);
    end

    assert(~isempty(order), 'No matching TIFFs for pref="%s" in: %s', pref, folderPath);
    tifPath = fullfile(allTifs(order(1)).folder, allTifs(order(1)).name);
end


%% ========================================================================
function run_svd_for_file(tiffPath, savePDF, modesDir, fps, ...
                          pad, dropFirstN, demeanPerPixel, nModes, alpha_polarity)
% Core SVD analysis on a single TIFF stack.

    info = imfinfo(tiffPath);
    H0   = info(1).Height;
    W0   = info(1).Width;
    T0   = numel(info);

    fprintf('[SVD] Loading %d frames of %dx%d from %s\n', T0, H0, W0, tiffPath);

    % crop
    rows = (1+pad):(H0-pad);
    cols = (1+pad):(W0-pad);
    Hc   = numel(rows);
    Wc   = numel(cols);

    Y_F = zeros(Hc*Wc, T0, 'single');
    for kk = 1:T0
        fr = imread(tiffPath, kk, 'Info', info);
        fr = fr(rows, cols);
        Y_F(:,kk) = single(fr(:));
    end

    % toss first N frames
    if dropFirstN > 0 && T0 > dropFirstN
        Y_F = Y_F(:, dropFirstN+1:end);
    end

    T = size(Y_F, 2);
    t = (0:T-1)' / fps;
    fprintf('[SVD] Data matrix: [%d px x %d frames]\n', Hc*Wc, T);

    % demean
    if demeanPerPixel
        Y_F_demean = Y_F - mean(Y_F, 2);
    else
        Y_F_demean = Y_F;
    end

    % SVD
    fprintf('[SVD] Computing SVD (econ)...\n');
    tic;
    [U_F,S_F,V_F] = svd(Y_F_demean, 'econ');
    fprintf('[SVD] SVD done in %.2f s\n', toc);

    singvals_F = diag(S_F);
    varExp_F   = singvals_F.^2 / sum(singvals_F.^2);

    %% Save SVD results (top K modes only)
    outSVD = fullfile(modesDir, 'SVD_result.mat');
    Ksave = min(nModes, size(V_F,2));

    SVD_result = struct();
    SVD_result.Hc = Hc;
    SVD_result.Wc = Wc;
    SVD_result.T  = T;
    SVD_result.fps = fps;
    SVD_result.pad = pad;
    SVD_result.dropFirstN = dropFirstN;
    SVD_result.demeanPerPixel = demeanPerPixel;
    SVD_result.sourceTif = tiffPath;

    SVD_result.t  = t(:);
    SVD_result.mu = single(mean(Y_F,2));
    SVD_result.U  = single(U_F(:,1:Ksave));
    SVD_result.S  = single(S_F(1:Ksave,1:Ksave));
    SVD_result.V  = single(V_F(:,1:Ksave));
    SVD_result.singvals = single(diag(S_F(1:Ksave,1:Ksave)));
    SVD_result.varExp   = single(varExp_F(1:Ksave));

    save(outSVD, 'SVD_result', '-v7.3');
    d = dir(outSVD);
    fprintf('[SVD] Saved SVD_result (top %d) -> %s (%.1f MB)\n', ...
            Ksave, outSVD, d.bytes/1e6);

    % global raw limits for binned raw video
    lims_raw = prctile(Y_F(:), [0 99]);

    %% Normalized singular value spectrum
    figSV = figure('Name','Singular value spectrum','Color','w','Visible','off');
    nShow = min(200, numel(varExp_F));
    semilogy(1:nShow, 100*varExp_F(1:nShow), 'bo-', ...
             'LineWidth',1.2, 'MarkerFaceColor','b');
    xlabel('mode index');
    ylabel('variance explained (%)');
    title('(Raw F) Top singular values — normalized');
    grid on;

    outSV = fullfile(modesDir, 'normSingularSpectrum.pdf');
    exportgraphics(figSV, outSV, 'ContentType','vector');
    close(figSV);
    fprintf('[SVD] Saved singular value spectrum -> %s\n', outSV);

    %% Raw singular value spectrum
    figSV2 = figure('Name','Eigenvalue spectrum','Color','w','Visible','off');
    nShow = min(200, numel(singvals_F));
    semilogy(1:nShow, singvals_F(1:nShow).^2, 'bo-', ...
             'LineWidth',1.2, 'MarkerFaceColor','b');
    xlabel('index of eigenvalue');
    ylabel('Eigenvalue');
    title('(Raw F) Eigenvalue spectrum');
    grid on;

    outSV2 = fullfile(modesDir, 'SingularSpectrum.pdf');
    exportgraphics(figSV2, outSV2, 'ContentType','vector');
    close(figSV2);
    fprintf('[SVD] Saved eigenvalue spectrum -> %s\n', outSV2);

    nModesEff = min(nModes, size(V_F,2));

    %% Per-mode PDFs + combined movies
    if savePDF
        for k = 1:nModesEff
            fprintf('[SVD] Mode %d/%d\n', k, nModesEff);

            u = U_F(:, k);
            v = V_F(:, k);
            s = S_F(k,k);

            Tloc = numel(v);
            tloc = t(:);

            %% Polarity test on v(t)
            v_dt = detrend(v);
            pct  = 0.015;             % ~1.5% tails
            nSel = max(1, round(pct * Tloc));

            [~, ord] = sort(v_dt, 'ascend');
            idx_low  = ord(1:nSel);
            idx_high = ord(end-nSel+1:end);

            v_lo = abs(v_dt(idx_low));
            v_hi = abs(v_dt(idx_high));

            mu_lo = mean(v_lo);
            mu_hi = mean(v_hi);
            sd_lo = std(v_lo);
            sd_hi = std(v_hi);
            n_lo = numel(v_lo);
            n_hi = numel(v_hi);

            se = sqrt(sd_lo^2/n_lo + sd_hi^2/n_hi);
            tstat = (mu_hi - mu_lo) / (se + eps);
            df = (sd_lo^2/n_lo + sd_hi^2/n_hi)^2 / ...
                 ((sd_lo^2/n_lo)^2/(n_lo-1) + (sd_hi^2/n_hi)^2/(n_hi-1));
            p = 1 - tcdf(tstat, df);

            if p >= alpha_polarity
                u = -u;
                v = -v;
                idx_l = idx_high;
                idx_h = idx_low;
            else
                idx_l = idx_low;
                idx_h = idx_high;
            end

            us = u * s;

            avg_low  = mean(Y_F(:, idx_l),  2);
            avg_high = mean(Y_F(:, idx_h), 2);

            img_low  = reshape(avg_low,  [Hc, Wc]);
            img_high = reshape(avg_high, [Hc, Wc]);
            diffMap  = img_high - img_low;

            us = reshape(us, [Hc, Wc]);

            %% Figure: trace + eigenface + raw diff
            fig = figure('Name', sprintf('Mode %d', k), ...
                         'Color', 'w', 'Visible','off');

            tl = tiledlayout(fig, 4, 6, ...
                'TileSpacing','compact', 'Padding','compact');

            % 1) temporal trace
            ax1 = nexttile(tl, [1 6]);
            plot(tloc, v, 'k', 'LineWidth', 0.2); hold on;
            plot(tloc(idx_l), v(idx_l), 'blueo', 'MarkerFaceColor','b', 'MarkerSize',4);
            plot(tloc(idx_h), v(idx_h), 'ro',   'MarkerFaceColor','r', 'MarkerSize',4);
            hold off;
            xlim([tloc(1) tloc(end)]);
            xlabel('time (s)');
            ylabel('v_k(t)');
            grid on;
            title(sprintf('Mode %d V(t)', k));

            sigma = 1.2;
            us_sm = imgaussfilt(us, sigma);
            diff_sm = imgaussfilt(diffMap, sigma);

            % 2) eigenface
            ax2 = nexttile(tl, [3 3]);
            imagesc(ax2, us_sm);
            axis(ax2, 'off');
            absVals_us = abs(us(:));
            mx_us = prctile(absVals_us, 99);
            if mx_us == 0, mx_us = 1; end
            caxis(ax2, [-mx_us mx_us]);
            colormap(ax2, bluewhitered);
            colorbar(ax2);
            title(ax2, 'U*S (eigen image)');

            % 3) raw high - low
            ax3 = nexttile(tl, [3 3]);
            imagesc(ax3, diff_sm);
            axis(ax3, 'off');
            absVals_diff = abs(diffMap(:));
            mx_diff = prctile(absVals_diff, 99);
            if mx_diff == 0, mx_diff = 1; end
            caxis(ax3, [-mx_diff mx_diff]);
            colormap(ax3, bluewhitered);
            colorbar(ax3);
            title(ax3, 'Raw (mean top v) - (mean bottom v)');

            sgtitle(tl, sprintf('Mode %d: eigenface vs raw top/bottom(v)', k));

            outPDF = fullfile(modesDir, sprintf('Mode_%02d_top_bottom.pdf', k));
            exportgraphics(fig, outPDF, 'ContentType','vector');
            close(fig);
            fprintf('  Saved %s\n', outPDF);

            %% Percentile-binned raw + recon combined movie
            Y_mode = s * (u(:) * v(:)');          % [Hc*Wc x T]
            lims_recon = prctile(Y_mode(:), [1 99]);

            stepP  = 5;
            edgesP = 0:stepP:100;
            nBins  = numel(edgesP) - 1;
            vEdges = prctile(v, edgesP);

            % raw bins
            Y_bin = zeros(Hc, Wc, nBins, 'single');
            for b = 1:nBins
                v_lo_b = vEdges(b);
                v_hi_b = vEdges(b+1);

                if b < nBins
                    idx_bin = find(v >= v_lo_b & v <  v_hi_b);
                else
                    idx_bin = find(v >= v_lo_b & v <= v_hi_b);
                end
                if isempty(idx_bin); continue; end

                img_bin = mean(Y_F(:, idx_bin), 2);
                Y_bin(:,:,b) = reshape(img_bin, [Hc, Wc]);
            end

            % scale raw
            lo = lims_raw(1);
            hi = lims_raw(2);
            Y_norm = (Y_bin - lo) / (hi - lo);
            Y_norm = min(max(Y_norm, 0), 1);
            Y_scaled_uint8 = uint8(255 * Y_norm);

            % recon bins
            Recon_bin = zeros(Hc, Wc, nBins, 'single');
            for b = 1:nBins
                v_lo_b = vEdges(b);
                v_hi_b = vEdges(b+1);

                if b < nBins
                    idx_bin = find(v >= v_lo_b & v <  v_hi_b);
                else
                    idx_bin = find(v >= v_lo_b & v <= v_hi_b);
                end
                if isempty(idx_bin); continue; end

                v_bin   = mean(v(idx_bin));
                img_bin = reshape(s * u(:) * v_bin, [Hc, Wc]);
                Recon_bin(:,:,b) = img_bin;
            end

            % scale recon
            lo_r = lims_recon(1);
            hi_r = lims_recon(2);
            Recon_norm = (Recon_bin - lo_r) / (hi_r - lo_r);
            Recon_norm = min(max(Recon_norm, 0), 1);
            Recon_scaled_uint8 = uint8(255 * Recon_norm);

            % combined movie: raw | recon
            Combined     = cat(2, Y_scaled_uint8, Recon_scaled_uint8);
            combinedPath = fullfile(modesDir, sprintf('Mode_%02d_combined.tif', k));

            nFrames_comb = size(Combined, 3);
            for b = 1:nFrames_comb
                frame = Combined(:,:,b);
                if b == 1
                    imwrite(frame, combinedPath, 'tif','Compression','none');
                else
                    imwrite(frame, combinedPath, 'tif', ...
                            'WriteMode','append','Compression','none');
                end
            end
            fprintf('  Saved combined movie: %s\n', combinedPath);
        end

        %% Merge per-mode PDFs into AllModes.pdf
        pdfList = dir(fullfile(modesDir, 'Mode_*.pdf'));
        if ~isempty(pdfList)
            modeNums = cellfun(@(s) sscanf(s, 'Mode_%d'), {pdfList.name});
            [~, idxSort] = sort(modeNums);
            pdfList = pdfList(idxSort);

            fileNames = fullfile(modesDir, {pdfList.name});
            mergedPDF = fullfile(modesDir, 'AllModes.pdf');
            append_pdfs_RZ(mergedPDF, fileNames{:});
            fprintf('[SVD] Merged %d PDFs into: %s\n', numel(fileNames), mergedPDF);
        else
            fprintf('[SVD] No per-mode PDFs found to merge.\n');
        end
    else
        fprintf('[SVD] Skipping PDF summary.\n');
    end

end
