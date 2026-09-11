% breathing_population_manifold.m
%   PCA-based visualization of the population CALCIUM manifold.
%   (Breathing / inspiration-phase analysis removed -- this is now a pure
%    calcium population-manifold + neuron-space clustering tool.)
%
%   Per session: PCA on dFF [T x N] (z-scored per ROI) -> time-domain manifold,
%   plus neuron-space k-means clustering on the PC loadings. Figures:
%     Fig 5     PC pairs plot (first 10 PCs), per session
%     Fig 7     dFF stacked by cluster
%     Fig 8     per-cluster mean dFF (+/- SEM)
%     Fig 9     cluster anatomy: FOV avg projection with ROI outlines by cluster
%     COMBINED  time-space 3D | scree | PC traces | neuron-space 3D | neuron-space 2D
%
%   Auto-detects input mode:
%     - inputPath itself contains *cpSAM_output.mat -> single-session
%     - Otherwise -> recurse to find all session folders
%
%   Dependencies:
%     - helper.dFF_RZ
%     - detect_session_fps.m

clear; close all; clc;

%% ========================= USER PARAMETERS =============================
%inputPath = "D:\Ventral_surface_summary\Vglut2\1124\cell1\roi5_1400-1230-0_x4.4_15lp_6000f_00001";
%inputPath = "D:\Ventral_surface_summary\Vglut2\1124\IO\roi3_-1000-300-40_32lp_930_x1.4_256x512_3000f_00001";

%inputPath = "D:\Ventral_surface_summary\Vglut2\1124\IO\roi3_-1000-300-60_38lp_930_x1.4_256x512_3000f_00001";

inputPath = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\airpuff_redo\Right_POandDAO_z80_airpuff_00001";
fps_img        = 30;        % fallback fps if detect_session_fps fails
BaselineWinSec = 20;        % dF/F sliding-median window (s)

zscore_per_roi = true;      % z-score each ROI (column) before PCA

nPCsForCluster = 20;        % neuron-space clustering: # PCs to use
k_clusters     = 5;         % neuron-space clustering: # k-means clusters
% ========================================================================

%% PATH SETUP
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(fullfile(repoRoot, '2p_breathing_coherence'));

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 8);

%% DISCOVER SESSIONS
% Discover by *cpSAM_output.mat (the actual data we load).
inputPath = char(inputPath);
local_hits = dir(fullfile(inputPath, '*cpSAM_output.mat'));
if ~isempty(local_hits)
    allMat = local_hits;
else
    allMat = dir(fullfile(inputPath, '**', '*cpSAM_output.mat'));
end
nSess = numel(allMat);
assert(nSess > 0, 'No *cpSAM_output.mat found under %s', inputPath);
fprintf('Found %d session(s) under %s\n', nSess, inputPath);

%% PER-SESSION COLLECTION
sess = struct('name',{},'folder',{},'dFF',{},'fs',{});
for kk = 1:nSess
    folderPath = allMat(kk).folder;
    [~, sessName] = fileparts(folderPath);
    fprintf('\n[%d/%d] %s\n', kk, nSess, sessName);

    try
        sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
        if isempty(sam_hits)
            warning('  No *cpSAM_output.mat, skipping'); continue;
        end

        % ---- fps (auto-detect, fallback to fps_img) ----
        fs = fps_img;
        if exist('detect_session_fps', 'file') == 2
            try
                fs = detect_session_fps(folderPath, fps_img);
            catch
                fs = fps_img;
            end
        end

        % ---- F -> dFF ----
        SAM    = load(fullfile(sam_hits(1).folder, sam_hits(1).name));
        F_raw  = SAM.F;
        dFFout = helper.dFF_RZ(F_raw, 'FPS', fs, 'BaselineWinSec', BaselineWinSec);
        dFF    = double(dFFout.dFF);   % [T x N], cpSAM F is already nDrop-trimmed
        [T, N] = size(dFF);
        fprintf('  dFF: [%d x %d]  fs=%.2f Hz\n', T, N, fs);

        s = struct('name', sessName, 'folder', folderPath, 'dFF', dFF, 'fs', fs);
        sess(end+1) = s; %#ok<SAGROW>
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        continue;
    end
end

nSessOK = numel(sess);
assert(nSessOK > 0, 'No usable sessions.');
fprintf('\n=== %d usable session(s) ===\n', nSessOK);

[~, inputName] = fileparts(inputPath);

%% PER-SESSION TIME-DOMAIN PCA
for kk = 1:nSessOK
    D = sess(kk).dFF;
    roi_idx_used = (1:size(D, 2))';
    if size(D, 2) < 3
        sess(kk).scores_t   = [];
        sess(kk).varexp_t   = [];
        sess(kk).V_t        = [];
        sess(kk).roi_idx_t  = [];
        continue;
    end
    if zscore_per_roi, Dz = zscore(D, 0, 1); else, Dz = D - mean(D, 1); end
    [Ut, St, Vt] = svd(Dz, 'econ');
    sess(kk).scores_t  = Ut * St;             % [T x r]: time-space coords
    sess(kk).V_t       = Vt;                  % [N x r]: neuron-space loadings
    sess(kk).varexp_t  = (diag(St)).^2 / sum((diag(St)).^2);
    sess(kk).roi_idx_t = roi_idx_used(:);
end

%% Neuron-space k-means clustering (used by fig7/8/9 and COMBINED)
for kk = 1:nSessOK
    if isempty(sess(kk).V_t)
        sess(kk).cluster_id_v = [];
        continue;
    end
    Vt = sess(kk).V_t;
    [Nroi, rk] = size(Vt);
    nC = min(nPCsForCluster, rk);
    k  = min(k_clusters, max(1, Nroi - 1));
    if k > 1 && Nroi >= k
        rng(1, 'twister');
        sess(kk).cluster_id_v = kmeans(Vt(:, 1:nC), k, ...
            'Replicates', 10, 'MaxIter', 200, 'EmptyAction', 'singleton');
    else
        sess(kk).cluster_id_v = ones(Nroi, 1);
    end
end

%% ============================ FIGURE 5: pairs plot, first 10 PCs ========
%  10×10 grid of PCi-vs-PCj 2D trajectories. Diagonal = PC vs time.
%  One figure per session. Useful for spotting structure that lives in
%  axes beyond PC1-3 (e.g. orthogonal rotations in PC4-5).
nPC_pairs = 10;
fig5 = gobjects(nSessOK, 1);
for kk = 1:nSessOK
    if isempty(sess(kk).scores_t), continue; end
    sct = sess(kk).scores_t;
    ve  = sess(kk).varexp_t;
    nP  = min(nPC_pairs, size(sct, 2));
    if nP < 2, continue; end

    fs_kk = sess(kk).fs;
    t_kk  = (0:size(sct,1)-1)' / fs_kk;

    fig5(kk) = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 24 24], 'PaperPositionMode', 'auto', ...
        'Renderer', 'opengl');

    for r = 1:nP
        for c = 1:nP
            ax = subplot(nP, nP, (r-1)*nP + c);
            if r == c
                % Diagonal: PC score vs time
                plot(ax, t_kk, sct(:, r), 'k-', 'LineWidth', 0.3);
                axis(ax, 'tight');
            else
                % Off-diagonal: PCc vs PCr trajectory line
                plot(ax, sct(:, c), sct(:, r), 'k-', 'LineWidth', 0.2);
                axis(ax, 'tight'); axis(ax, 'square');
            end
            % Tick labels only on outer edges
            if r < nP, set(ax, 'XTickLabel', []); end
            if c > 1,  set(ax, 'YTickLabel', []); end
            % Axis labels only on outer edges
            if r == nP
                xlabel(ax, sprintf('PC%d (%.1f%%)', c, 100*ve(c)), 'FontSize', 6);
            end
            if c == 1
                ylabel(ax, sprintf('PC%d (%.1f%%)', r, 100*ve(r)), 'FontSize', 6);
            end
            set(ax, 'FontSize', 5);
            box(ax, 'on');
        end
    end

    sgtitle(fig5(kk), sprintf('FIG 5: PC pairs plot (first %d PCs) | %s', ...
        nP, strrep(sess(kk).name, '_', '\_')), ...
        'Interpreter', 'tex', 'FontSize', 10);
end

%% ============================ FIGURE 7: dFF stacked by cluster ===========
%  All ROIs stacked vertically, grouped by cluster (1 -> k top to bottom).
%  Each trace z-scored. Cluster colors match the neuron-space scatter.
fig7 = gobjects(nSessOK, 1);
for kk = 1:nSessOK
    if isempty(sess(kk).V_t) || ~isfield(sess(kk), 'cluster_id_v') || ...
            isempty(sess(kk).cluster_id_v)
        continue;
    end
    cluster_id = sess(kk).cluster_id_v;
    roi_idx    = sess(kk).roi_idx_t;
    fs_kk      = sess(kk).fs;
    k          = max(cluster_id);

    % Within each cluster, sort by V1 magnitude descending so the most
    % strongly participating neurons show first.
    Vt = sess(kk).V_t;
    sort_key = nan(numel(cluster_id), 2);
    sort_key(:, 1) = cluster_id;             % primary: cluster
    sort_key(:, 2) = -abs(Vt(:, 1));         % secondary: |V1| descending
    [~, ix] = sortrows(sort_key);
    roi_sorted     = roi_idx(ix);
    cluster_sorted = cluster_id(ix);

    % z-score each ROI's dFF
    D  = sess(kk).dFF(:, roi_sorted);
    Dz = (D - mean(D, 1, 'omitnan')) ./ std(D, 0, 1, 'omitnan');
    [T, N] = size(Dz);
    t_kk     = (0:T-1)' / fs_kk;
    trace_amp_scale = 2.0;    % per-trace amplitude scale (overlap is OK)
    fig_scale       = 0.64;   % overall figure size scale
    spacing = 1.2 * median(max(Dz, [], 1) - min(Dz, [], 1), 'omitnan');

    fig7(kk) = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 30*fig_scale max(15, 0.18 * N + 4)*fig_scale], ...
        'PaperPositionMode', 'auto', 'Renderer', 'painters');
    movegui(fig7(kk), 'center');   % avoid top-right clipping off-screen
    ax = axes(fig7(kk)); %#ok<LAXES>
    hold(ax, 'on');

    cmap_k = lines(k);
    cluster_centers = nan(k, 1);
    for jj = 1:N
        cc    = cluster_sorted(jj);
        y_off = (N - jj) * spacing;          % jj=1 (cluster 1, top) -> highest y
        plot(ax, t_kk, trace_amp_scale * Dz(:, jj) + y_off, '-', ...
            'Color', cmap_k(cc, :), 'LineWidth', 0.8);
    end

    % Cluster boundaries: separator lines + center labels
    for cc = 1:k
        ix_cc = find(cluster_sorted == cc);
        if isempty(ix_cc), continue; end
        cluster_centers(cc) = (N - mean(ix_cc)) * spacing;
        if cc < k
            y_sep = (N - max(ix_cc) - 0.5) * spacing;
            yline(ax, y_sep, 'k:', 'LineWidth', 0.5);
        end
    end

    [centers_sorted, ord] = sort(cluster_centers, 'ascend');
    cluster_labels = arrayfun(@(c) sprintf('cluster %d (n=%d)', c, ...
        sum(cluster_id == c)), 1:k, 'UniformOutput', false);
    set(ax, 'YTick', centers_sorted, 'YTickLabel', cluster_labels(ord), ...
        'TickLabelInterpreter', 'tex', 'FontSize', 8);
    xlabel(ax, 'time (s)');
    xlim(ax, [0 t_kk(end)]);
    % Tight y-limits
    y_top = (N-1)*spacing + trace_amp_scale * max(Dz, [], 'all', 'omitnan');
    y_bot = trace_amp_scale * min(Dz, [], 'all', 'omitnan');
    ylim(ax, [y_bot y_top]);
    box(ax, 'on'); grid(ax, 'on');
    title(ax, sprintf('FIG 7: dFF stacked by cluster (N=%d ROIs) | %s', ...
        N, strrep(sess(kk).name, '_', '\_')), ...
        'Interpreter', 'tex', 'FontSize', 10);

    fprintf('Fig 7 [%s]: %d ROIs stacked, %d clusters\n', sess(kk).name, N, k);
end

%% ============================ FIGURE 8: per-cluster mean dFF ============
%  One trace per cluster, mean dFF across that cluster's ROIs (z-scored
%  per ROI first to prevent loud cells from dominating).
fig8 = gobjects(nSessOK, 1);
for kk = 1:nSessOK
    if isempty(sess(kk).cluster_id_v), continue; end
    cluster_id = sess(kk).cluster_id_v;
    roi_idx    = sess(kk).roi_idx_t;
    fs_kk      = sess(kk).fs;
    k          = max(cluster_id);

    D  = sess(kk).dFF(:, roi_idx);
    Dz = (D - mean(D, 1, 'omitnan')) ./ std(D, 0, 1, 'omitnan');
    [T, ~] = size(Dz);
    t_kk   = (0:T-1)' / fs_kk;

    % Population mean & SEM (across all ROIs)
    pop_n    = size(Dz, 2);
    pop_mean = mean(Dz, 2, 'omitnan');
    pop_sem  = std(Dz, 0, 2, 'omitnan') ./ sqrt(pop_n);

    % Mean & SEM per cluster
    mean_per_cluster = nan(T, k);
    sem_per_cluster  = nan(T, k);
    cluster_n        = zeros(k, 1);
    for cc = 1:k
        m = cluster_id == cc;
        cluster_n(cc) = sum(m);
        if cluster_n(cc) == 0, continue; end
        mean_per_cluster(:, cc) = mean(Dz(:, m), 2, 'omitnan');
        sem_per_cluster(:, cc)  = std(Dz(:, m), 0, 2, 'omitnan') ./ ...
            sqrt(cluster_n(cc));
    end

    fig8(kk) = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 28 10 + 1.2*(k+1)], 'PaperPositionMode', 'auto', ...
        'Renderer', 'painters');
    movegui(fig8(kk), 'center');

    cmap_k = lines(k);
    ax = axes(fig8(kk)); %#ok<LAXES>
    hold(ax, 'on');

    % Stacked: population on top, then cluster 1, ..., cluster k.
    ranges_clu = max(mean_per_cluster, [], 1, 'omitnan') - ...
                 min(mean_per_cluster, [], 1, 'omitnan');
    range_pop  = max(pop_mean) - min(pop_mean);
    spacing    = 1.4 * max([ranges_clu, range_pop], [], 'omitnan');

    yticks_all = nan(k+1, 1);
    yticklbl   = cell(k+1, 1);

    % --- Population trace on top ---
    y_off_pop = k * spacing;
    valid = ~isnan(pop_mean);
    fill(ax, [t_kk(valid); flipud(t_kk(valid))], ...
        [pop_mean(valid) + pop_sem(valid); flipud(pop_mean(valid) - pop_sem(valid))] + y_off_pop, ...
        [0 0 0], 'FaceAlpha', 0.20, 'EdgeColor', 'none');
    plot(ax, t_kk, pop_mean + y_off_pop, 'k-', 'LineWidth', 1.4);
    yticks_all(1) = y_off_pop;
    yticklbl{1}   = sprintf('population (n=%d)', pop_n);

    % --- Per-cluster traces below ---
    for cc = 1:k
        y_off = (k - cc) * spacing;
        mu    = mean_per_cluster(:, cc);
        sem   = sem_per_cluster(:, cc);
        valid = ~isnan(mu);
        fill(ax, [t_kk(valid); flipud(t_kk(valid))], ...
            [mu(valid) + sem(valid); flipud(mu(valid) - sem(valid))] + y_off, ...
            cmap_k(cc, :), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        plot(ax, t_kk, mu + y_off, '-', 'Color', cmap_k(cc, :), 'LineWidth', 1.0);
        yticks_all(cc + 1) = y_off;
        yticklbl{cc + 1}   = sprintf('cluster %d (n=%d)', cc, cluster_n(cc));
    end

    [yticks_sorted, ord] = sort(yticks_all, 'ascend');
    set(ax, 'YTick', yticks_sorted, 'YTickLabel', yticklbl(ord), ...
        'TickLabelInterpreter', 'tex', 'FontSize', 8);
    xlabel(ax, 'time (s)');
    ylabel(ax, 'mean dF/F (z-scored, +/-SEM)');
    xlim(ax, [0 t_kk(end)]);
    title(ax, sprintf('FIG 8: per-cluster mean dF/F | %s', ...
        strrep(sess(kk).name, '_', '\_')), 'Interpreter', 'tex', 'FontSize', 10);
    grid(ax, 'on'); box(ax, 'on');

    fprintf('Fig 8 [%s]: %d cluster averages, sizes = [%s]\n', ...
        sess(kk).name, k, num2str(cluster_n(:)'));
end

%% ============================ FIGURE 9: cluster anatomy (FOV) ===========
%  Average projection (gamma, percentile clip) with each ROI's cellpose
%  boundary outlined in its cluster color.
fig9 = gobjects(nSessOK, 1);
gamma_disp   = 0.5;
clip_pct     = [10 99.8];
boundary_lw  = 1.2;
for kk = 1:nSessOK
    if isempty(sess(kk).cluster_id_v), continue; end
    folderPath = sess(kk).folder;
    cluster_id = sess(kk).cluster_id_v;
    roi_idx    = sess(kk).roi_idx_t;
    k          = max(cluster_id);
    cmap_k     = lines(k);

    % Locate AVG projection (cellpose input) and cellpose label masks
    avg_hits = dir(fullfile(folderPath, '*_AVG_for_CP.tif'));
    if isempty(avg_hits)
        avg_hits = dir(fullfile(folderPath, '*AVG*.tif'));
    end
    mask_hits = dir(fullfile(folderPath, '*_cp_masks.tif'));
    if isempty(mask_hits)
        mask_hits = dir(fullfile(folderPath, '*ROIlabel*.tif'));
    end
    if isempty(avg_hits) || isempty(mask_hits)
        fprintf('Fig 9 [%s]: missing AVG or masks tif, skipping\n', sess(kk).name);
        continue;
    end

    avg_img  = double(imread(fullfile(avg_hits(1).folder, avg_hits(1).name)));
    mask_img = double(imread(fullfile(mask_hits(1).folder, mask_hits(1).name)));

    % Display: percentile clip then gamma
    clim     = prctile(avg_img(:), clip_pct);
    avg_norm = (avg_img - clim(1)) / max(clim(2) - clim(1), eps);
    avg_norm = min(max(avg_norm, 0), 1);
    avg_disp = avg_norm .^ gamma_disp;

    fig9(kk) = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 18 14], 'PaperPositionMode', 'auto', ...
        'Renderer', 'painters');
    movegui(fig9(kk), 'center');
    ax = axes(fig9(kk)); %#ok<LAXES>

    imshow(avg_disp, [], 'Parent', ax, 'InitialMagnification', 'fit');
    colormap(ax, gray(256));
    hold(ax, 'on');

    % Outline each ROI's boundary in its cluster color
    h_leg = gobjects(k, 1);
    n_drawn = zeros(k, 1);
    for ii = 1:numel(roi_idx)
        roi_label_id = roi_idx(ii);
        cc = cluster_id(ii);
        roi_mask = (mask_img == roi_label_id);
        if ~any(roi_mask, 'all'), continue; end
        B = bwboundaries(roi_mask, 'noholes');
        for bb = 1:numel(B)
            h = plot(ax, B{bb}(:, 2), B{bb}(:, 1), '-', ...
                'Color', cmap_k(cc, :), 'LineWidth', boundary_lw);
            if ~isgraphics(h_leg(cc)), h_leg(cc) = h; end
        end
        n_drawn(cc) = n_drawn(cc) + 1;
    end

    legend(ax, h_leg(isgraphics(h_leg)), ...
        arrayfun(@(c) sprintf('cluster %d (n=%d)', c, n_drawn(c)), ...
            find(isgraphics(h_leg))', 'UniformOutput', false), ...
        'Location', 'northeastoutside', 'FontSize', 8, 'Box', 'off', ...
        'TextColor', 'k');
    title(ax, sprintf('FIG 9: cluster anatomy | %s | gamma=%.1f, clip=[%.1f %.1f]%%', ...
        strrep(sess(kk).name, '_', '\_'), gamma_disp, clip_pct(1), clip_pct(2)), ...
        'Interpreter', 'tex', 'FontSize', 9);
    axis(ax, 'image'); axis(ax, 'off');

    fprintf('Fig 9 [%s]: %d ROIs outlined\n', sess(kk).name, sum(n_drawn));
end

%% ============================ COMBINED: 5 panels in one row =============
%  Panels: time-space 3D | scree | PC traces | neuron-space 3D | neuron-space 2D
fig_combined = gobjects(nSessOK, 1);
nPC_traces_combined = 20;
for kk = 1:nSessOK
    if isempty(sess(kk).scores_t) || size(sess(kk).scores_t, 2) < 3, continue; end
    sct = sess(kk).scores_t;
    ve  = sess(kk).varexp_t;
    Vt         = sess(kk).V_t;
    cluster_id = sess(kk).cluster_id_v;
    Nroi       = size(Vt, 1);
    k          = max(cluster_id);
    cmap_k     = lines(k);

    fig_combined(kk) = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 56 12], 'PaperPositionMode', 'auto', ...
        'Renderer', 'opengl');
    movegui(fig_combined(kk), 'center');

    tlo = tiledlayout(fig_combined(kk), 1, 5, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    % --- Tile 1: time-space 3D manifold (colored by time) ---
    ax1 = nexttile(tlo, 1);
    x = sct(:,1); y = sct(:,2); z = sct(:,3);
    c   = (1:numel(x))' / sess(kk).fs;
    cmp = parula(256); cax = [0 max(c)];
    surface(ax1, [x x], [y y], [z z], [c c], ...
        'EdgeColor','interp','FaceColor','none','LineWidth',1.0);
    colormap(ax1, cmp); caxis(ax1, cax);
    xlabel(ax1, sprintf('PC1 (%.1f%%)', 100*ve(1)));
    ylabel(ax1, sprintf('PC2 (%.1f%%)', 100*ve(2)));
    zlabel(ax1, sprintf('PC3 (%.1f%%)', 100*ve(3)));
    title(ax1, 'time-space 3D');
    view(ax1, 30, 25); grid(ax1,'on'); box(ax1,'on'); axis(ax1,'tight');
    cb1 = colorbar(ax1); cb1.Label.String = 'time (s)';

    % --- Tile 2: scree ---
    ax2 = nexttile(tlo, 2);
    nsv = min(20, numel(ve));
    yyaxis(ax2, 'left');
    bar(ax2, 1:nsv, 100*ve(1:nsv), 'FaceColor', [0.30 0.50 0.80], 'EdgeColor', 'none');
    ylabel(ax2, '% var explained');
    yyaxis(ax2, 'right');
    plot(ax2, 1:nsv, 100*cumsum(ve(1:nsv)), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12);
    ylabel(ax2, 'cumulative %');
    ylim(ax2, [0 100]);
    xlabel(ax2, 'PC #');
    title(ax2, sprintf('Scree (top 3: %.1f%%)', 100*sum(ve(1:min(3,end)))));
    grid(ax2, 'on'); box(ax2, 'on');

    % --- Tile 3: first N PC traces stacked ---
    ax3   = nexttile(tlo, 3);
    nP    = min(nPC_traces_combined, size(sct, 2));
    fs_kk = sess(kk).fs;
    t_kk  = (0:size(sct,1)-1)' / fs_kk;
    pc_block  = sct(:, 1:nP);
    ranges_pc = max(pc_block, [], 1) - min(pc_block, [], 1);
    spacing_pc = 1.2 * median(ranges_pc);
    hold(ax3, 'on');
    yticks_pc = nan(nP, 1);
    yticklbl  = cell(nP, 1);
    for p = 1:nP
        y_off = (nP - p) * spacing_pc;
        plot(ax3, t_kk, pc_block(:, p) - mean(pc_block(:, p)) + y_off, ...
            'k-', 'LineWidth', 0.4);
        yticks_pc(p) = y_off;
        yticklbl{p}  = sprintf('PC%d (%.1f%%)', p, 100*ve(p));
    end
    set(ax3, 'YTick', flipud(yticks_pc), 'YTickLabel', flipud(yticklbl), ...
        'TickLabelInterpreter', 'tex', 'FontSize', 6);
    xlabel(ax3, 'time (s)');
    xlim(ax3, [0 t_kk(end)]);
    title(ax3, sprintf('First %d PC traces', nP));
    grid(ax3, 'on'); box(ax3, 'on');

    % --- Tile 4: neuron-space 3D scatter (V1 x V2 x V3, cluster-colored) ---
    ax4 = nexttile(tlo, 4);
    hold(ax4, 'on');
    h_leg = gobjects(k, 1);
    for cc = 1:k
        m = cluster_id == cc;
        h_leg(cc) = scatter3(ax4, Vt(m,1), Vt(m,2), Vt(m,3), 30, ...
            cmap_k(cc, :), 'filled', 'MarkerEdgeColor', 'k');
    end
    xlabel(ax4, sprintf('V1 (PC1, %.1f%%)', 100*ve(1)));
    ylabel(ax4, sprintf('V2 (PC2, %.1f%%)', 100*ve(2)));
    zlabel(ax4, sprintf('V3 (PC3, %.1f%%)', 100*ve(3)));
    title(ax4, sprintf('neuron-space 3D (N=%d, k=%d)', Nroi, k));
    view(ax4, 30, 25); grid(ax4,'on'); box(ax4,'on'); axis(ax4,'tight');
    legend(ax4, h_leg, arrayfun(@(c) sprintf('c%d', c), 1:k, ...
        'UniformOutput', false), 'Location', 'best', 'FontSize', 6, 'Box', 'off');

    % --- Tile 5: neuron-space 2D scatter (V1 vs V2) ---
    ax5 = nexttile(tlo, 5);
    hold(ax5, 'on');
    for cc = 1:k
        m = cluster_id == cc;
        scatter(ax5, Vt(m,1), Vt(m,2), 50, cmap_k(cc, :), 'filled', ...
            'MarkerEdgeColor', 'k');
    end
    xline(ax5, 0, 'k:'); yline(ax5, 0, 'k:');
    xlabel(ax5, sprintf('V1 (PC1, %.1f%%)', 100*ve(1)));
    ylabel(ax5, sprintf('V2 (PC2, %.1f%%)', 100*ve(2)));
    title(ax5, 'neuron-space 2D');
    grid(ax5,'on'); box(ax5,'on'); axis(ax5,'equal');

    title(tlo, sprintf('COMBINED | %s', strrep(sess(kk).name, '_', '\_')), ...
        'Interpreter', 'tex', 'FontSize', 10);
    fprintf('Combined fig [%s]: 5 panels (time-3D, scree, %d PCs, neuron-3D, neuron-2D)\n', ...
        sess(kk).name, nP);
end
