% breathing_population_manifold.m
%   PCA-based visualization of the population calcium manifold over breathing.
%
%   Two views are produced from the same input folder:
%     1. PHASE-BINNED MANIFOLD (Figure 1)
%        Pool ROIs across all sessions in the input folder. Build an
%        [nPhaseBins x N_total] matrix of mean dFF per inspiration-phase bin,
%        z-score per ROI, run PCA -> project the nPhaseBins phase points into
%        PC1/PC2/PC3. If the population is breath-locked, this draws a closed
%        loop; the loop's ring-closure ratio quantifies how cleanly it closes.
%
%     2. TIME-DOMAIN MANIFOLD (Figure 2)
%        Per session: PCA on dFF [T x N] (z-scored per ROI), project every
%        frame into PC1/PC2 and color by instantaneous breathing phase.
%        Reveals whether the manifold is organized by phase as a single axis
%        or whether phase is one feature among several.
%
%   Auto-detects input mode:
%     - inputPath itself contains ca_spike_data.mat -> single-session
%     - Otherwise -> recurse to find all session folders
%
%   Outputs (saved into <inputPath>/population_manifold/):
%     - pcaManifold_phase_<tag>.pdf    Figure 1 (phase-binned)
%     - pcaManifold_time_<tag>.pdf     Figure 2 (per-session time-domain)
%     - pcaManifold_<tag>.mat          PCs, scores, variance explained
%
%   Dependencies:
%     - helper.dFF_RZ
%     - detect_session_fps.m

clear; close all; clc;

%% ========================= USER PARAMETERS =============================
inputPath      = "D:\251124_live_vglut2_soma_g8s+cy5\phys\IO\roi3_-1000-300-80_48lp_930_x1.4_256x512_3000f_00001";

nPhaseBins     = 24;        % phase bins for phase-locked manifold
nDrop          = 30;        % frames dropped at start (must match upstream)
fps_img        = 30;        % fallback fps if detect_session_fps fails
BaselineWinSec = 20;        % dF/F sliding-median window (s)

zscore_per_roi = true;      % z-score each ROI (column) before PCA
sigOnly        = false;     % restrict to coherence-sig ROIs (via fov_map.mat)
nPC_show       = 3;         % # PCs to keep for ring projection

T_trunc_sec     = 60;       % Figure 4: truncate each session to this many seconds
useSig_for_fig4 = true;     % Figure 4: project only sig ROIs

nPCsForCluster   = 10;       % neuron-space clustering: # PCs to use
k_clusters  = 5;       % neuron-space clustering: # k-means clusters

doSave         = true;
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
% Discover by *cpSAM_output.mat (the actual data we load). ca_spike_data.mat
% was the old discovery anchor but it's a downstream artifact and not
% required for manifold analysis.
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

% Optional sig filter via fov_map.mat (master-mode breath_combined_summary)
fovMat   = fullfile(inputPath, 'breath_combined_summary', 'fov_map.mat');
has_fov  = isfile(fovMat);
if has_fov
    Sf = load(fovMat); fov_map = Sf.fov_map;
    fov_sess_names = string(fov_map(:,1));
else
    fov_map = {}; fov_sess_names = strings(0,1);
end
if sigOnly && ~has_fov
    warning('sigOnly requested but no fov_map.mat found; falling back to all ROIs');
    sigOnly = false;
end

%% PER-SESSION COLLECTION
sess = struct('name',{},'folder',{},'PB',{},'dFF',{},'frame_phase',{},'fs',{},'sig_mask',{});
for kk = 1:nSess
    folderPath = allMat(kk).folder;
    [~, sessName] = fileparts(folderPath);
    fprintf('\n[%d/%d] %s\n', kk, nSess, sessName);

    try
        sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
        bp_file  = dir(fullfile(folderPath, '*breath_peak_data.mat'));
        if isempty(sam_hits)
            warning('  No *cpSAM_output.mat, skipping'); continue;
        end
        has_breath = ~isempty(bp_file);
        if ~has_breath
            fprintf('  No *breath_peak_data.mat -> phase analysis disabled, time-domain only\n');
        end

        % ---- fps (auto-detect, fallback to fps_img) ----
        if exist('detect_session_fps', 'file') == 2
            try
                fs = detect_session_fps(folderPath, fps_img);
            catch
                fs = fps_img;
            end
        else
            fs = fps_img;
        end

        % ---- F -> dFF ----
        SAM    = load(fullfile(sam_hits(1).folder, sam_hits(1).name));
        F_raw  = SAM.F;
        dFFout = helper.dFF_RZ(F_raw, 'FPS', fs, 'BaselineWinSec', BaselineWinSec);
        dFF    = double(dFFout.dFF);   % [T x N], cpSAM F is already nDrop-trimmed
        [T, N] = size(dFF);
        fprintf('  dFF: [%d x %d]  fs=%.2f Hz\n', T, N, fs);

        % ---- breath peaks (raw frame coords -> post-drop coords) ----
        frame_phase = nan(T, 1);
        PB          = nan(nPhaseBins, N);
        if has_breath
            bp = load(fullfile(bp_file(1).folder, bp_file(1).name));
            insp_onsets = sort(bp.insp_onset_idx(:));
            insp_onsets(insp_onsets <= nDrop) = [];
            insp_onsets = insp_onsets - nDrop;
            insp_onsets(insp_onsets > T) = [];
            if numel(insp_onsets) < 3
                warning('  Fewer than 3 breath peaks; phase analysis disabled');
                has_breath = false;
            else
                for c = 1:numel(insp_onsets)-1
                    i1 = insp_onsets(c); i2 = insp_onsets(c+1);
                    if i1 >= 1 && i2 <= T && i2 > i1
                        L   = i2 - i1;
                        idx = i1:(i2-1);
                        frame_phase(idx) = 2*pi * ((0:L-1)' + 0.5) / L;
                    end
                end

                % phase-binned mean dFF: [nPhaseBins x N]
                edges  = linspace(0, 2*pi, nPhaseBins + 1);
                bin_id = discretize(frame_phase, edges);
                for k = 1:nPhaseBins
                    m = (bin_id == k);
                    if any(m), PB(k, :) = mean(dFF(m, :), 1, 'omitnan'); end
                end
            end
        end

        % ---- sig mask ----
        % Default: all-true (no fov_map info) OR all-false (fov_map present
        % but session not listed). This ensures sessions absent from fov_map
        % don't leak ALL their ROIs into the sig pool.
        if has_fov
            sig_mask = false(N, 1);
            mm = find(fov_sess_names == string(sessName), 1);
            if ~isempty(mm)
                sig_ids = fov_map{mm, 2};
                sig_mask(sig_ids(sig_ids >= 1 & sig_ids <= N)) = true;
            else
                fprintf('  WARNING: session not in fov_map, 0 sig ROIs\n');
            end
        else
            sig_mask = true(N, 1);
        end

        s = struct('name', sessName, 'folder', folderPath, 'PB', PB, 'dFF', dFF, ...
            'frame_phase', frame_phase, 'fs', fs, 'sig_mask', sig_mask);
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

%% POOL phase-binned matrix: ALL ROIs (Figure 1)
[PB_all, sess_of_roi_all] = pool_phase_matrix(sess, false);
[~, N_all] = size(PB_all);
fprintf('Phase-binned pool (ALL): [%d phase bins x %d ROIs]\n', nPhaseBins, N_all);
fig1 = []; all_pca = struct();
if N_all >= 3
    [fig1, all_pca] = plot_phase_manifold(PB_all, nPhaseBins, nPC_show, zscore_per_roi, ...
        sprintf('PHASE-BINNED MANIFOLD | ALL ROIs (N=%d) | %s', N_all, inputName));
    fprintf('Ring closure (ALL): %.3f\n', all_pca.ring_closure);
else
    fprintf('Skipping Figure 1: no usable phase-binned data (no breath info?)\n');
end

%% PER-SESSION TIME-DOMAIN PCA (Figure 2 — uses sigOnly toggle for D selection)
for kk = 1:nSessOK
    if sigOnly && any(sess(kk).sig_mask)
        D = sess(kk).dFF(:, sess(kk).sig_mask);
        roi_idx_used = find(sess(kk).sig_mask);
    else
        D = sess(kk).dFF;
        roi_idx_used = (1:size(sess(kk).dFF, 2))';
    end
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

%% Neuron-space k-means clustering (used by fig7 and COMBINED)
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

%% ===== Per-session time-domain PCA visualization (3 SEPARATE figures) =====
%  Figure 2 -- 3D trajectory PC1 × PC2 × PC3 (colored by phase or time)
%  Figure 3 -- scree (% variance per PC + cumulative)
%  Figure 4 -- first 20 PC score traces, stacked vertically vs time
nPC_traces = 20;
nColG = ceil(sqrt(nSessOK));
nRowG = ceil(nSessOK / nColG);
used_phase_color = false;

%% Figures 2, 3, 4 are now in the COMBINED row figure below; standalones disabled.
%{
%% ============================ FIGURE 2: 3D trajectories =================
fig2 = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 max(8, 7*nColG) max(7, 6*nRowG)], ...
    'PaperPositionMode', 'auto', 'Renderer', 'opengl');
ax_last2 = [];
for kk = 1:nSessOK
    ax = subplot(nRowG, nColG, kk);
    if isempty(sess(kk).scores_t) || size(sess(kk).scores_t, 2) < 3
        text(ax, 0.5, 0.5, '<3 PCs', 'HorizontalAlignment', 'center');
        axis(ax, 'off');
        title(ax, strrep(sess(kk).name, '_', '\_'), 'FontSize', 7);
        continue;
    end
    sct = sess(kk).scores_t;
    ph  = sess(kk).frame_phase;
    x = sct(:,1); y = sct(:,2); z = sct(:,3);
    if all(isnan(ph))
        c   = (1:numel(x))' / sess(kk).fs;
        cmp = parula(256);
        cax = [0 c(end)];
    else
        keep = ~isnan(ph);
        x = x(keep); y = y(keep); z = z(keep);
        c = ph(keep);
        cmp = hsv(256);
        cax = [0 2*pi];
        used_phase_color = true;
    end
    surface(ax, [x x], [y y], [z z], [c c], ...
        'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 1.0);
    colormap(ax, cmp); caxis(ax, cax);
    ve = sess(kk).varexp_t;
    xlabel(ax, sprintf('PC1 (%.1f%%)', 100*ve(1)));
    ylabel(ax, sprintf('PC2 (%.1f%%)', 100*ve(2)));
    zlabel(ax, sprintf('PC3 (%.1f%%)', 100*ve(3)));
    title(ax, strrep(sess(kk).name, '_', '\_'), 'FontSize', 7);
    view(ax, 30, 25); axis(ax, 'tight'); grid(ax, 'on'); box(ax, 'on');
    ax_last2 = ax;
end
if ~isempty(ax_last2)
    cb2 = colorbar(ax_last2);
    if used_phase_color
        cb2.Label.String = 'breath phase (rad)';
        cb2.Ticks = [0 pi/2 pi 3*pi/2 2*pi];
        cb2.TickLabels = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    else
        cb2.Label.String = 'time (s)';
    end
end
sgtitle(fig2, sprintf('FIG 2: 3D MANIFOLD | %s', inputName), ...
    'Interpreter', 'none', 'FontSize', 10);

%% ============================ FIGURE 3: scree (% variance) ==============
fig3 = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 max(8, 7*nColG) max(7, 5*nRowG)], ...
    'PaperPositionMode', 'auto', 'Renderer', 'painters');
for kk = 1:nSessOK
    ax = subplot(nRowG, nColG, kk);
    if isempty(sess(kk).scores_t)
        axis(ax, 'off'); continue;
    end
    ve  = sess(kk).varexp_t;
    nsv = min(20, numel(ve));
    yyaxis(ax, 'left');
    bar(ax, 1:nsv, 100*ve(1:nsv), 'FaceColor', [0.30 0.50 0.80], 'EdgeColor', 'none');
    ylabel(ax, '% var explained');
    yyaxis(ax, 'right');
    plot(ax, 1:nsv, 100*cumsum(ve(1:nsv)), 'r.-', 'LineWidth', 1.2, 'MarkerSize', 12);
    ylabel(ax, 'cumulative %');
    ylim(ax, [0 100]);
    xlabel(ax, 'PC #');
    title(ax, sprintf('%s | top 3: %.1f%%', strrep(sess(kk).name, '_', '\_'), ...
        100*sum(ve(1:min(3,end)))), 'FontSize', 7);
    grid(ax, 'on'); box(ax, 'on');
end
sgtitle(fig3, sprintf('FIG 3: PCA SCREE | %s', inputName), ...
    'Interpreter', 'none', 'FontSize', 10);

%% ============================ FIGURE 4: first N PC traces ===============
fig4 = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 max(10, 9*nColG) max(10, 8*nRowG)], ...
    'PaperPositionMode', 'auto', 'Renderer', 'painters');
for kk = 1:nSessOK
    ax = subplot(nRowG, nColG, kk);
    if isempty(sess(kk).scores_t)
        axis(ax, 'off'); continue;
    end
    sct   = sess(kk).scores_t;
    ve    = sess(kk).varexp_t;
    nP    = min(nPC_traces, size(sct, 2));
    fs_kk = sess(kk).fs;
    t_kk  = (0:size(sct,1)-1)' / fs_kk;
    pc_block  = sct(:, 1:nP);
    ranges_pc = max(pc_block, [], 1) - min(pc_block, [], 1);
    spacing_pc = 1.2 * median(ranges_pc);
    hold(ax, 'on');
    yticks_pc = nan(nP, 1);
    yticklbl  = cell(nP, 1);
    for p = 1:nP
        y_off = (nP - p) * spacing_pc;     % PC1 on top
        plot(ax, t_kk, pc_block(:, p) - mean(pc_block(:, p)) + y_off, ...
            'k-', 'LineWidth', 0.4);
        yticks_pc(p) = y_off;
        yticklbl{p}  = sprintf('PC%d (%.1f%%)', p, 100*ve(p));
    end
    set(ax, 'YTick', flipud(yticks_pc), 'YTickLabel', flipud(yticklbl), ...
        'TickLabelInterpreter', 'tex', 'FontSize', 6);
    xlabel(ax, 'time (s)');
    xlim(ax, [0 t_kk(end)]);
    title(ax, sprintf('%s | first %d PCs', strrep(sess(kk).name, '_', '\_'), nP), ...
        'FontSize', 7);
    grid(ax, 'on'); box(ax, 'on');
end
sgtitle(fig4, sprintf('FIG 4: first %d PC traces | %s', nPC_traces, inputName), ...
    'Interpreter', 'none', 'FontSize', 10);
%}

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

%% ============================ FIGURE 6: NEURON-SPACE manifold === DISABLED
%  Now folded into the COMBINED figure below (panels 4 + 5).
%{
fig6 = gobjects(nSessOK, 1);
for kk = 1:nSessOK
    if isempty(sess(kk).V_t) || size(sess(kk).V_t, 2) < 3, continue; end

    Vt = sess(kk).V_t;
    [Nroi, rk] = size(Vt);
    nC  = min(nPCsForCluster, rk);

    % k-means in the first nC loading dimensions; bound k by Nroi-1.
    k = min(k_clusters, max(1, Nroi - 1));
    if k > 1 && Nroi >= k
        rng(1, 'twister');                                  % reproducible
        cluster_id = kmeans(Vt(:, 1:nC), k, ...
            'Replicates', 10, 'MaxIter', 200, 'EmptyAction', 'singleton');
    else
        cluster_id = ones(Nroi, 1);
    end

    fig6(kk) = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 36 11], 'PaperPositionMode', 'auto', ...
        'Renderer', 'opengl');

    % --- Left: 3D scatter V1 x V2 x V3 ---
    axA = subplot(1, 3, 1);
    cmap_k = lines(k);
    hold(axA, 'on');
    h_leg = gobjects(k, 1);
    for cc = 1:k
        m = cluster_id == cc;
        h_leg(cc) = scatter3(axA, Vt(m, 1), Vt(m, 2), Vt(m, 3), 35, ...
            cmap_k(cc, :), 'filled', 'MarkerEdgeColor', 'k');
    end
    ve = sess(kk).varexp_t;
    xlabel(axA, sprintf('V1 (PC1, %.1f%%)', 100*ve(1)));
    ylabel(axA, sprintf('V2 (PC2, %.1f%%)', 100*ve(2)));
    zlabel(axA, sprintf('V3 (PC3, %.1f%%)', 100*ve(3)));
    title(axA, sprintf('3D neuron-space (N=%d ROIs, k=%d on V1..V%d)', ...
        Nroi, k, nC), 'FontSize', 8);
    view(axA, 30, 25); grid(axA, 'on'); box(axA, 'on'); axis(axA, 'tight');
    legend(axA, h_leg, arrayfun(@(c) sprintf('cluster %d', c), 1:k, ...
        'UniformOutput', false), 'Location', 'best', 'FontSize', 7, 'Box', 'off');

    % --- Middle: 2D scatter V1 vs V2 ---
    axB = subplot(1, 3, 2);
    hold(axB, 'on');
    for cc = 1:k
        m = cluster_id == cc;
        scatter(axB, Vt(m, 1), Vt(m, 2), 50, cmap_k(cc, :), 'filled', ...
            'MarkerEdgeColor', 'k');
    end
    xline(axB, 0, 'k:'); yline(axB, 0, 'k:');
    xlabel(axB, sprintf('V1 (PC1, %.1f%%)', 100*ve(1)));
    ylabel(axB, sprintf('V2 (PC2, %.1f%%)', 100*ve(2)));
    title(axB, '2D neuron-space (V1 vs V2)', 'FontSize', 8);
    grid(axB, 'on'); box(axB, 'on'); axis(axB, 'equal');

    % --- Right: scree (% variance per PC + cumulative) ---
    axC = subplot(1, 3, 3);
    nsv = min(20, numel(ve));
    yyaxis(axC, 'left');
    bar(axC, 1:nsv, 100*ve(1:nsv), 'FaceColor', [0.30 0.50 0.80], 'EdgeColor', 'none');
    ylabel(axC, '% var explained');
    yyaxis(axC, 'right');
    plot(axC, 1:nsv, 100*cumsum(ve(1:nsv)), 'r.-', ...
        'LineWidth', 1.2, 'MarkerSize', 12);
    ylabel(axC, 'cumulative %');
    ylim(axC, [0 100]);
    xlabel(axC, 'PC #');
    title(axC, sprintf('Scree (top 3: %.1f%%; clustered on V1..V%d)', ...
        100*sum(ve(1:min(3,end))), nC), 'FontSize', 8);
    grid(axC, 'on'); box(axC, 'on');
    % Mark the clustering cutoff
    xline(axC, nC + 0.5, 'k--', 'LineWidth', 0.8);

    % Stash cluster IDs back on sess for downstream reuse
    sess(kk).cluster_id_v = cluster_id;

    sgtitle(fig6(kk), sprintf('FIG 6: NEURON-SPACE MANIFOLD | %s', ...
        strrep(sess(kk).name, '_', '\_')), 'Interpreter', 'tex', 'FontSize', 10);

    % Cluster size summary
    cluster_sizes = arrayfun(@(c) sum(cluster_id == c), 1:k);
    fprintf('Fig 6 [%s]: cluster sizes = [%s]\n', sess(kk).name, ...
        num2str(cluster_sizes));
end
%}

%% ============================ FIGURE 7: dFF stacked by cluster ===========
%  All ROIs stacked vertically, grouped by Fig 6 cluster (1 -> 4 top to bottom).
%  Each trace z-scored. Cluster colors match Fig 6.
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
    fig_scale       = 0.64;   % overall figure size scale (was 0.8, now * 0.8)
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
    % Tight y-limits: just below the bottom trace's lowest excursion to just
    % above the top trace's highest excursion.
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
%  per ROI first to prevent loud cells from dominating). Same colors as
%  fig7 / combined neuron-space scatter.
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

    % Stacked: population on top, then cluster 1, ..., cluster k. Spacing
    % from the largest mean trace range (clusters or population).
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
%  Average projection (gamma 0.5, percentile clip 0.5-99.5) with each ROI's
%  cellpose boundary outlined in its cluster color.
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
%  Replaces stand-alone fig2 + fig3 + fig4 + fig6. Pairs (fig5) and
%  stacked dFF (fig7) remain as their own figures.
fig_combined = gobjects(nSessOK, 1);
nPC_traces_combined = 20;
for kk = 1:nSessOK
    if isempty(sess(kk).scores_t) || size(sess(kk).scores_t, 2) < 3, continue; end
    sct = sess(kk).scores_t;
    ve  = sess(kk).varexp_t;
    ph  = sess(kk).frame_phase;
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

    % --- Tile 1: time-space 3D manifold ---
    ax1 = nexttile(tlo, 1);
    x = sct(:,1); y = sct(:,2); z = sct(:,3);
    if all(isnan(ph))
        c   = (1:numel(x))' / sess(kk).fs;
        cmp = parula(256); cax = [0 c(end)];
        cb_label = 'time (s)'; cb_ticks = []; cb_lbls = {};
    else
        keep = ~isnan(ph);
        x = x(keep); y = y(keep); z = z(keep); c = ph(keep);
        cmp = hsv(256); cax = [0 2*pi];
        cb_label = 'breath phase (rad)';
        cb_ticks = [0 pi/2 pi 3*pi/2 2*pi];
        cb_lbls  = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    end
    surface(ax1, [x x], [y y], [z z], [c c], ...
        'EdgeColor','interp','FaceColor','none','LineWidth',1.0);
    colormap(ax1, cmp); caxis(ax1, cax);
    xlabel(ax1, sprintf('PC1 (%.1f%%)', 100*ve(1)));
    ylabel(ax1, sprintf('PC2 (%.1f%%)', 100*ve(2)));
    zlabel(ax1, sprintf('PC3 (%.1f%%)', 100*ve(3)));
    title(ax1, 'time-space 3D');
    view(ax1, 30, 25); grid(ax1,'on'); box(ax1,'on'); axis(ax1,'tight');
    cb1 = colorbar(ax1);
    cb1.Label.String = cb_label;
    if ~isempty(cb_ticks), cb1.Ticks = cb_ticks; cb1.TickLabels = cb_lbls; end

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

%% ============================ FIGURE 3 =============================== — DISABLED
%{
%  Phase-binned manifold restricted to coherence-significantly-tuned ROIs
fig3 = []; sig_pca = struct(); PB_sig = []; sess_of_roi_sig = [];
if has_fov
    [PB_sig, sess_of_roi_sig] = pool_phase_matrix(sess, true);
    [~, N_sig] = size(PB_sig);
    fprintf('Phase-binned pool (SIG): [%d phase bins x %d ROIs]\n', nPhaseBins, N_sig);
    if N_sig >= 3
        [fig3, sig_pca] = plot_phase_manifold(PB_sig, nPhaseBins, nPC_show, zscore_per_roi, ...
            sprintf('PHASE-BINNED MANIFOLD | SIG ONLY (N=%d) | %s', N_sig, inputName));
        fprintf('Ring closure (SIG): %.3f\n', sig_pca.ring_closure);
    else
        fprintf('Skipping Figure 3: only %d sig ROIs (need >=3)\n', N_sig);
    end
else
    fprintf('Skipping Figure 3: no fov_map.mat available\n');
end
%}

%% ============================ FIGURE 4 =============================== — DISABLED
%{
%  Per-session LONGITUDINAL trajectory in time-domain PC space.
%  One continuous curve per session (PC1 x PC2 x PC3), colored by time
%  (sequential parula: start=blue, end=yellow). Truncates each session to
%  the first T_trunc_sec seconds for visual comparability across sessions.

fprintf('\n--- Figure 4: longitudinal trajectories (first %.1fs per session) ---\n', T_trunc_sec);

nCol4 = ceil(sqrt(nSessOK));
nRow4 = ceil(nSessOK / nCol4);
fig4 = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 max(8, 7*nCol4) max(7, 6*nRow4)], ...
    'PaperPositionMode', 'auto', 'Renderer', 'opengl');

n_used    = 0;
n_skipped = 0;
ax_last   = [];
for kk = 1:nSessOK
    ax = subplot(nRow4, nCol4, kk);
    fs_kk     = sess(kk).fs;
    T_full    = size(sess(kk).dFF, 1);
    T_target  = round(T_trunc_sec * fs_kk);
    if T_full < T_target
        text(ax, 0.5, 0.5, sprintf('too short\n(%.1fs < %.1fs)', T_full/fs_kk, T_trunc_sec), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
        axis(ax, 'off');
        title(ax, strrep(sess(kk).name, '_', '\_'), 'FontSize', 7);
        n_skipped = n_skipped + 1;
        continue;
    end

    if useSig_for_fig4
        D = sess(kk).dFF(1:T_target, sess(kk).sig_mask);
        roi_label = 'sig';
    else
        D = sess(kk).dFF(1:T_target, :);
        roi_label = 'all';
    end
    if size(D, 2) < 3
        text(ax, 0.5, 0.5, sprintf('<3 %s ROIs\n(have %d)', roi_label, size(D,2)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
        axis(ax, 'off');
        title(ax, strrep(sess(kk).name, '_', '\_'), 'FontSize', 7);
        n_skipped = n_skipped + 1;
        continue;
    end

    if zscore_per_roi
        Dz = zscore(D, 0, 1);
    else
        Dz = D - mean(D, 1);
    end
    [Ut, St, ~] = svd(Dz, 'econ');
    sct = Ut * St;
    ve  = (diag(St)).^2 / sum((diag(St)).^2);

    x = sct(:, 1); y = sct(:, 2); z = sct(:, 3);
    t_sec = (0:numel(x)-1)' / fs_kk;

    % MATLAB has no native multi-color line, so use the surface-as-line trick.
    surface(ax, [x x], [y y], [z z], [t_sec t_sec], ...
        'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 1.2);
    colormap(ax, parula(256));
    caxis(ax, [0 T_trunc_sec]);

    xlabel(ax, sprintf('PC1 (%.1f%%)', 100*ve(1)));
    ylabel(ax, sprintf('PC2 (%.1f%%)', 100*ve(2)));
    zlabel(ax, sprintf('PC3 (%.1f%%)', 100*ve(3)));
    title(ax, sprintf('%s', strrep(sess(kk).name, '_', '\_')), 'FontSize', 7);
    view(ax, 30, 25); grid(ax, 'on'); box(ax, 'on'); axis(ax, 'tight');
    n_used  = n_used + 1;
    ax_last = ax;
end

if ~isempty(ax_last)
    cb4 = colorbar(ax_last);
    cb4.Label.String = 'time (s)';
end
sgtitle(fig4, sprintf('LONGITUDINAL MANIFOLD (per-session, first %.1fs) | %s | useSig=%d', ...
    T_trunc_sec, inputName, useSig_for_fig4), 'Interpreter', 'none', 'FontSize', 10);
fprintf('Figure 4: %d/%d sessions plotted, %d skipped\n', n_used, nSessOK, n_skipped);

%% ============================ FIGURE 5 =============================== — DISABLED
%  Sig-tuned ROIs: full-length dFF traces, one panel per sig ROI.
%  Plotted at native session length (no truncation, no PCA).
sig_traces = struct('dFF',{}, 't',{}, 'sessName',{}, 'roi_id',{}, 'fs',{});
for kk = 1:nSessOK
    if ~any(sess(kk).sig_mask), continue; end
    sig_ids = find(sess(kk).sig_mask);
    fs_kk   = sess(kk).fs;
    T_kk    = size(sess(kk).dFF, 1);
    t_vec   = (0:T_kk-1)' / fs_kk;
    for jj = 1:numel(sig_ids)
        s = struct( ...
            'dFF',      sess(kk).dFF(:, sig_ids(jj)), ...
            't',        t_vec, ...
            'sessName', sess(kk).name, ...
            'roi_id',   sig_ids(jj), ...
            'fs',       fs_kk);
        sig_traces(end+1) = s; %#ok<SAGROW>
    end
end

n_sig_traces = numel(sig_traces);
fig5 = [];
if n_sig_traces == 0
    fprintf('Figure 5 skipped: no sig ROIs\n');
else
    % Stack all sig traces on a single shared x-axis. Traces end naturally
    % at their own duration; longer ones extend further to the right.
    max_dur = max(arrayfun(@(s) s.t(end), sig_traces));
    ranges  = arrayfun(@(s) max(s.dFF) - min(s.dFF), sig_traces);
    spacing = 1.2 * median(ranges);

    fig5 = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 30 max(10, 1.6*n_sig_traces)], ...
        'PaperPositionMode', 'auto', 'Renderer', 'painters');
    ax = axes(fig5); %#ok<LAXES>
    hold(ax, 'on');

    yticks_pos = nan(n_sig_traces, 1);
    ylabels    = cell(n_sig_traces, 1);
    for jj = 1:n_sig_traces
        % Bottom-up: jj=1 at bottom, last at top
        y_off  = (jj - 1) * spacing;
        trace  = sig_traces(jj).dFF - mean(sig_traces(jj).dFF, 'omitnan');
        plot(ax, sig_traces(jj).t, trace + y_off, 'k-', 'LineWidth', 0.5);
        yticks_pos(jj) = y_off;
        ylabels{jj} = sprintf('%s | ROI %d', ...
            strrep(sig_traces(jj).sessName, '_', '\_'), sig_traces(jj).roi_id);
    end
    xlabel(ax, 'time (s)');
    ylabel(ax, 'dF/F (stacked, mean-subtracted)');
    xlim(ax, [0, max_dur]);
    set(ax, 'YTick', yticks_pos, 'YTickLabel', ylabels, 'TickLabelInterpreter', 'tex');
    grid(ax, 'on'); box(ax, 'on');
    title(ax, sprintf('SIG ROIs | full-length dFF stack (n=%d) | %s', ...
        n_sig_traces, inputName), 'Interpreter', 'none', 'FontSize', 10);
    fprintf('Figure 5: %d sig ROI traces stacked (max dur = %.1fs)\n', n_sig_traces, max_dur);
end
%}

%% SAVE — DISABLED
%{
if doSave
    outDir = fullfile(inputPath, 'population_manifold');
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    tag = sprintf('%dbins_z%d', nPhaseBins, zscore_per_roi);

    pdf1 = fullfile(outDir, sprintf('pcaManifold_phase_all_%s.pdf', tag));
    pdf2 = fullfile(outDir, sprintf('pcaManifold_time_%s.pdf',      tag));
    pdf3 = fullfile(outDir, sprintf('pcaManifold_phase_sig_%s.pdf', tag));
    matp = fullfile(outDir, sprintf('pcaManifold_%s.mat',           tag));

    exportgraphics(fig1, pdf1, 'ContentType', 'vector', 'BackgroundColor', 'none');
    exportgraphics(fig2, pdf2, 'ContentType', 'vector', 'BackgroundColor', 'none');
    if ~isempty(fig3)
        exportgraphics(fig3, pdf3, 'ContentType', 'vector', 'BackgroundColor', 'none');
    end

    sess_save = rmfield(sess, 'dFF');
    save(matp, 'PB_all', 'all_pca', 'sess_of_roi_all', ...
        'PB_sig', 'sig_pca', 'sess_of_roi_sig', ...
        'sess_save', 'nPhaseBins', 'zscore_per_roi', 'inputPath', '-v7.3');

    fprintf('\nSaved:\n  %s\n  %s\n', pdf1, pdf2);
    if ~isempty(fig3), fprintf('  %s\n', pdf3); end
    fprintf('  %s\n', matp);
end
%}

%% ============================ LOCAL FUNCTIONS ============================
function [PB_pool, sess_of_roi] = pool_phase_matrix(sess, sigOnly)
% Concatenate per-session phase-binned matrices [nPhaseBins x N_sess] into a
% single pool [nPhaseBins x N_total]. Drops ROIs with any empty (NaN) phase
% bin (typical for short sessions that don't visit all bins).
    nSessOK     = numel(sess);
    PB_blocks   = cell(nSessOK, 1);
    sess_of_roi = [];
    for kk = 1:nSessOK
        if sigOnly
            mask = sess(kk).sig_mask;
        else
            mask = true(size(sess(kk).PB, 2), 1);
        end
        if ~any(mask), continue; end
        PB_blocks{kk} = sess(kk).PB(:, mask);
        sess_of_roi   = [sess_of_roi; repmat(kk, sum(mask), 1)]; %#ok<AGROW>
    end
    PB_pool = cat(2, PB_blocks{:});
    if isempty(PB_pool), return; end

    nanCol = any(isnan(PB_pool), 1);
    if any(nanCol)
        fprintf('  Dropping %d ROIs with empty phase bins\n', sum(nanCol));
        PB_pool     = PB_pool(:, ~nanCol);
        sess_of_roi = sess_of_roi(~nanCol);
    end
end

function [fig, out] = plot_phase_manifold(PB_pool, nPhaseBins, nPC_show, zscore_per_roi, ttl)
% PCA on a phase-binned pool [nPhaseBins x N_ROIs] and 4-panel ring figure.
% Returns fig handle and a struct with PCs, scores, var explained, ring closure.
    if zscore_per_roi
        PB_in = zscore(PB_pool, 0, 1);
    else
        PB_in = PB_pool - mean(PB_pool, 1);
    end
    [U, S, V] = svd(PB_in, 'econ');
    sv     = diag(S);
    varexp = sv.^2 / sum(sv.^2);
    scores = U * S;
    nPC    = min(nPC_show, size(scores, 2));

    phase_ctrs = linspace(0, 2*pi, nPhaseBins + 1);
    phase_ctrs = phase_ctrs(1:end-1) + pi/nPhaseBins;

    ring_pts     = scores(:, 1:nPC);
    d_close      = norm(ring_pts(end, :) - ring_pts(1, :));
    d_consec     = mean(vecnorm(diff(ring_pts), 2, 2));
    ring_closure = d_close / max(d_consec, eps);

    out = struct('PB_in', PB_in, 'scores', scores, 'V', V, 'sv', sv, ...
        'varexp', varexp, 'ring_closure', ring_closure, 'phase_ctrs', phase_ctrs);

    % Periodic cubic spline through the 24 phase points -> smooth closed
    % curve. Trick: replicate points 3x along phase so the spline at the
    % wrap-around (phase=0 ~= 2pi) matches its neighbors on both sides.
    phase_ext  = [phase_ctrs(:) - 2*pi; phase_ctrs(:); phase_ctrs(:) + 2*pi];
    phase_dense = linspace(0, 2*pi, 400)';
    smooth_curve = nan(numel(phase_dense), nPC);
    for p = 1:nPC
        s_ext = repmat(scores(:, p), 3, 1);
        smooth_curve(:, p) = interp1(phase_ext, s_ext, phase_dense, 'spline');
    end

    fig = figure('Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 28 22], 'PaperPositionMode', 'auto', 'Renderer', 'painters');

    % --- A: 3D ring ---
    axA = subplot(2, 2, 1);
    plot3(axA, smooth_curve(:,1), smooth_curve(:,2), smooth_curve(:,3), ...
        'k-', 'LineWidth', 1.5);
    hold(axA, 'on');
    scatter3(axA, scores(:,1), scores(:,2), scores(:,3), 60, phase_ctrs, ...
        'filled', 'MarkerEdgeColor', 'k');
    colormap(axA, hsv(256)); caxis(axA, [0 2*pi]);
    xlabel(axA, sprintf('PC1 (%.1f%%)', 100*varexp(1)));
    ylabel(axA, sprintf('PC2 (%.1f%%)', 100*varexp(2)));
    zlabel(axA, sprintf('PC3 (%.1f%%)', 100*varexp(3)));
    title(axA, sprintf('3D ring (N=%d ROIs)', size(PB_pool, 2)));
    view(axA, 30, 25); grid(axA, 'on'); box(axA, 'on'); axis(axA, 'tight');

    % --- B: 2D ring with phase colorbar ---
    axB = subplot(2, 2, 2);
    plot(axB, smooth_curve(:,1), smooth_curve(:,2), 'k-', 'LineWidth', 1.5);
    hold(axB, 'on');
    scatter(axB, scores(:,1), scores(:,2), 90, phase_ctrs, ...
        'filled', 'MarkerEdgeColor', 'k');
    colormap(axB, hsv(256)); caxis(axB, [0 2*pi]);
    cb = colorbar(axB);
    cb.Label.String = 'breath phase (rad)';
    cb.Ticks = [0 pi/2 pi 3*pi/2 2*pi];
    cb.TickLabels = {'0','\pi/2','\pi','3\pi/2','2\pi'};
    xlabel(axB, sprintf('PC1 (%.1f%%)', 100*varexp(1)));
    ylabel(axB, sprintf('PC2 (%.1f%%)', 100*varexp(2)));
    title(axB, sprintf('PC1-PC2 ring | closure ratio = %.3f', ring_closure));
    axis(axB, 'equal'); grid(axB, 'on'); box(axB, 'on');

    % --- C: scree ---
    axC = subplot(2, 2, 3);
    nsv = min(20, numel(varexp));
    yyaxis(axC, 'left');
    bar(axC, 1:nsv, 100*varexp(1:nsv), 'FaceColor', [0.30 0.50 0.80], ...
        'EdgeColor', 'none');
    ylabel(axC, '% var explained');
    yyaxis(axC, 'right');
    plot(axC, 1:nsv, 100*cumsum(varexp(1:nsv)), 'r.-', ...
        'LineWidth', 1.2, 'MarkerSize', 12);
    ylabel(axC, 'cumulative %');
    ylim(axC, [0 100]);
    xlabel(axC, 'PC #');
    title(axC, sprintf('Scree (top 3: %.1f%%)', 100*sum(varexp(1:min(3,end)))));
    grid(axC, 'on'); box(axC, 'on');

    % --- D: PCs vs phase, two cycles ---
    axD = subplot(2, 2, 4);
    phase_2cyc = [phase_ctrs(:); phase_ctrs(:) + 2*pi];
    hold(axD, 'on');
    clrs = lines(nPC);
    hL = gobjects(nPC, 1);
    for p = 1:nPC
        y = [scores(:, p); scores(:, p)];
        hL(p) = plot(axD, phase_2cyc, y, '-', 'Color', clrs(p,:), 'LineWidth', 1.5);
    end
    xline(axD, 2*pi, 'k--');
    xlim(axD, [0 4*pi]);
    set(axD, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
        'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
    xlabel(axD, 'breath phase (rad)');
    ylabel(axD, 'PC score');
    title(axD, 'PCs vs phase (2 cycles shown)');
    legend(axD, hL, arrayfun(@(p) sprintf('PC%d', p), 1:nPC, ...
        'UniformOutput', false), 'Location', 'best', 'Box', 'off');
    grid(axD, 'on'); box(axD, 'on');

    sgtitle(fig, ttl, 'Interpreter', 'none', 'FontSize', 10);
end
