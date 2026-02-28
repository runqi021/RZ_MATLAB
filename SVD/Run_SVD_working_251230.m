%% step2_svd_singleFolder.m
% Run SVD + PDFs on all .tif / .tiff files directly inside ONE folder
% (MC or not, doesn't matter).

clear; clc;
%gcp;  % optional parallel pool

%% USER PARAMS
filePath        = "E:\RZ\data\shi\260115_homo_cal590\test\fov2_-200_L_30lp\raw\fov2_-200_L_30lp_1100nm_5x_00003.tif";
fps             = 30;       % <-- change if needed
demeanPerPixel  = true;
pad             = 0;        % <-- NO CROP
dropFirstN      = 30;        % <-- toss FIRST 6 FRAMES
nModes          = 30;       % number of modes to visualize
alpha_polarity  = 0.05;     % for polarity test
%%
[folderPath, baseName, ~] = fileparts(char(filePath));

modesDir = fullfile(folderPath, 'all_modes');
if ~exist(modesDir, 'dir')
    mkdir(modesDir);
    fprintf('  [FS] Created modes subfolder: %s\n', modesDir);
end

savePDF=true;
run_svd_for_file(filePath, savePDF, modesDir, fps, ...
                     pad, dropFirstN, demeanPerPixel, nModes, alpha_polarity);

fprintf('\n SVD Done.\n');


%% ========================================================================
function run_svd_for_file(tiffPath, savePDF, modesDir, fps, ...
                          pad, dropFirstN, demeanPerPixel, nModes, alpha_polarity)
% run_svd_for_file
%   tiffPath   : TIFF stack
%   modesDir   : subfolder for per-mode PDFs and movies
%   fps        : sampling rate
%   pad        : crop border (pixels)
%   dropFirstN : # frames to drop from start
%   demeanPerPixel : logical
%   nModes     : # modes to visualize
%   alpha_polarity : sig level for polarity test

    info = imfinfo(tiffPath);
    H0   = info(1).Height;
    W0   = info(1).Width;

    T0 = numel(info);          % use ALL frames
    t  = (0:T0-1)'/fps;

    fprintf('  [SVD] Loading %d frames of %dx%d from %s\n', ...
            T0, H0, W0, tiffPath);

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
        t   = t(dropFirstN+1:end);
    end

    T = numel(t);
    fprintf('  [SVD] Data matrix: [%d px × %d frames]\n', Hc*Wc, T);

    % demean
    if demeanPerPixel
        Y_F_demean = Y_F - mean(Y_F, 2);
    else
        Y_F_demean = Y_F;
    end

    % SVD
    fprintf('  [SVD] Computing SVD (econ)...\n');
    tic;
    [U_F,S_F,V_F] = svd(Y_F_demean, 'econ');
    fprintf('  [SVD] SVD done in %.2f s\n', toc);

    singvals_F = diag(S_F);
    varExp_F   = singvals_F.^2 / sum(singvals_F.^2);

    %% -------- SAVE SVD EARLY (robust) --------
    % save SVD results (TOP K only, v7.3 safe)
    outSVD = fullfile(modesDir, 'SVD_result.mat');
    
    Ksave = min(nModes, size(V_F,2));   % save only top nModes
    
    SVD_result = struct();
    SVD_result.Hc = Hc;
    SVD_result.Wc = Wc;
    SVD_result.T  = T;                 % after dropFirstN
    SVD_result.fps = fps;
    SVD_result.pad = pad;
    SVD_result.dropFirstN = dropFirstN;
    SVD_result.demeanPerPixel = demeanPerPixel;
    
    SVD_result.t  = t(:);                              % [T x 1]
    SVD_result.mu = single(mean(Y_F,2));               % [px x 1] raw mean (post-drop)
    SVD_result.U  = single(U_F(:,1:Ksave));            % [px x K]
    SVD_result.S  = single(S_F(1:Ksave,1:Ksave));      % [K x K]
    SVD_result.V  = single(V_F(:,1:Ksave));            % [T x K]
    SVD_result.singvals = single(diag(S_F(1:Ksave,1:Ksave)));
    SVD_result.varExp   = single(varExp_F(1:Ksave));
    
    save(outSVD, 'SVD_result', '-v7.3');
    
    d = dir(outSVD);
    fprintf('  [SVD] Saved SVD_result (top %d) -> %s (%.1f MB)\n', ...
            Ksave, outSVD, d.bytes/1e6);


    % global raw limits for visualizing binned raw video
    lims_raw = prctile(Y_F(:), [0 99]);

    %% normalized singular value spectrum (per movie)
    figSV = figure('Name','Singular value spectrum','Color','w', ...
                   'Visible','off');

    nShow = min(200, numel(varExp_F));
    semilogy(1:nShow, 100*varExp_F(1:nShow), 'bo-', ...
             'LineWidth',1.2, 'MarkerFaceColor','b');
    xlabel('mode index');
    ylabel('variance explained (%)');
    title('(Raw F) Top 100 variance fractions');
    grid on;

    outSV = fullfile(modesDir, 'normSingularSpectrum.pdf');
    exportgraphics(figSV, outSV, 'ContentType','vector');
    close(figSV);
    fprintf('  [SVD] Saved singular value spectrum → %s\n', outSV);

    % number of modes to use
    nModesEff = min(nModes, size(V_F,2));

    %% singular value spectrum (per movie)
    figSV = figure('Name','Singular value spectrum','Color','w', ...
                   'Visible','off');

    nShow = min(200, numel(singvals_F));
    semilogy(1:nShow, singvals_F(1:nShow).^2, 'bo-', ...
             'LineWidth',1.2, 'MarkerFaceColor','b');
    xlabel('index of eigenvalue');
    ylabel('Eigenvalue');
    title('(Raw F) Top 100 variance fractions');
    grid on;

    outSV = fullfile(modesDir, 'SingularSpectrum.pdf');
    exportgraphics(figSV, outSV, 'ContentType','vector');
    close(figSV);
    fprintf('  [SVD] Saved singular value spectrum → %s\n', outSV);

    % number of modes to use
    nModesEff = min(nModes, size(V_F,2));

    %%
    if savePDF
        for k = 1:nModesEff
            fprintf('  [SVD] Mode %d/%d\n', k, nModesEff);
    
            u = U_F(:, k);
            v = V_F(:, k);
            s = S_F(k,k);
    
            Tloc = numel(v);
            tloc = t(:);
    
            %% polarity test on v(t)
            v_dt = detrend(v);        % or v - median(v)
            pct  = 0.015;             % ~1.5% tails
            nSel = max(1, round(pct * Tloc));
    
            [~, ord] = sort(v_dt, 'ascend');
            idx_low  = ord(1:nSel);              % most negative
            idx_high = ord(end-nSel+1:end);      % most positive
    
            v_lo = abs(v_dt(idx_low));
            v_hi = abs(v_dt(idx_high));
    
            mu_lo = mean(v_lo);
            mu_hi = mean(v_hi);
    
            sd_lo = std(v_lo);
            sd_hi = std(v_hi);
    
            n_lo = numel(v_lo);
            n_hi = numel(v_hi);
    
            se = sqrt(sd_lo^2/n_lo + sd_hi^2/n_hi);
    
            % H1: |high| > |low|
            tstat = (mu_hi - mu_lo) / (se + eps);
    
            df = (sd_lo^2/n_lo + sd_hi^2/n_hi)^2 / ...
                 ((sd_lo^2/n_lo)^2/(n_lo-1) + (sd_hi^2/n_hi)^2/(n_hi-1));
    
            p = 1 - tcdf(tstat, df);
    
            if p >= alpha_polarity
                % positive tail not significantly stronger -> flip
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
    
            % smooth eigenface for vis
            %sigma = 1.2;
            %us_smooth = imgaussfilt(reshape(us, [Hc, Wc]), sigma);
            us = reshape(us, [Hc, Wc]);
            %% Figure: trace + eigenface + raw diff
            fig = figure('Name', sprintf('Mode %d', k), ...
                         'Color', 'w', 'Visible','off');
    
            tl = tiledlayout(fig, 4, 6, ...
                'TileSpacing','compact', ...
                'Padding','compact');
    
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
            colormap(ax2, bluewhitered);   % needs bluewhitered on path
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
            %colormap(ax3, gray);
            colormap(ax3, bluewhitered);
            colorbar(ax3);
            title(ax3, 'Raw (mean top v) - (mean bottom v)');
    
            sgtitle(tl, sprintf('Mode %d: eigenface vs raw top/bottom(v)', k));
    
            outPDF = fullfile(modesDir, sprintf('Mode_%02d_top_bottom.pdf', k));
            exportgraphics(fig, outPDF, 'ContentType','vector');
            close(fig);
    
            fprintf('    [SVD] Saved %s\n', outPDF);
    
            %% percentile-binned raw + recon + combined movie
            % full single-mode recon for limits
            Y_mode = s * (u * v');                 % [Hc*Wc x T]
            lims_recon = prctile(Y_mode(:), [1 99]);
    
            % percentile edges on v
            stepP  = 5;                     % 1%-wide bins
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
                img_bin = reshape(s * u * v_bin, [Hc, Wc]);
                Recon_bin(:,:,b) = img_bin;
            end
    
            % scale recon
            lo_r = lims_recon(1);
            hi_r = lims_recon(2);
            Recon_norm = (Recon_bin - lo_r) / (hi_r - lo_r);
            Recon_norm = min(max(Recon_norm, 0), 1);
            Recon_scaled_uint8 = uint8(255 * Recon_norm);
    
            % combined movie: raw | recon
            Combined     = cat(2, Y_scaled_uint8, Recon_scaled_uint8);   % H x (2W) x nBins
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
    
            fprintf('    [SVD] Saved combined movie: %s\n', combinedPath);
        end
    
        %% merge per-mode PDFs into AllModes.pdf (per movie)
        pdfList = dir(fullfile(modesDir, 'Mode_*.pdf'));
        if ~isempty(pdfList)
            % numeric sort on Mode_XX
            modeNums = cellfun(@(s) sscanf(s, 'Mode_%d'), {pdfList.name});
            [~, idxSort] = sort(modeNums);
            pdfList = pdfList(idxSort);
    
            fileNames = fullfile(modesDir, {pdfList.name});
            mergedPDF = fullfile(modesDir, 'AllModes.pdf');  % <-- per movie
    
            append_pdfs_RZ(mergedPDF, fileNames{:});
    
            fprintf('  [SVD] Merged %d PDFs into:\n         %s\n', ...
                    numel(fileNames), mergedPDF);
        else
            fprintf('  [SVD] No per-mode PDFs found to merge.\n');
        end
    else fprintf('Not saving PDF summary'); end

end
