%% Phase-aligned breathing overlay
%  Hilbert-transform the breathing signal, extract instantaneous phase,
%  and plot every breath cycle (inspiration-to-inspiration) as normalised
%  position vs phase φ ∈ [0, 1].
%
%  INPUT MODES (set ONE, leave the other empty):
%    singlePath — one folder containing breath_peak_data.mat
%    masterPath — parent folder; all subfolders with breath_peak_data.mat
%                 are pooled onto a single overlay figure
%
%  OUTPUTS — figure saved as breath_phase_overlay.png / .pdf

clear; clc; close all;

%% ===== USER PARAMETERS ================================================
singlePath = '';    % <-- single recording folder
masterPath = "D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\maybe breathing";    % <-- parent folder (searches subfolders)
% ========================================================================

%% Resolve folder list
if ~isempty(masterPath) && ~isempty(singlePath)
    error('Set either singlePath or masterPath, not both.');
elseif ~isempty(singlePath)
    d = dir(fullfile(singlePath, '*breath_peak_data.mat'));
    if isempty(d)
        error('No breath_peak_data.mat found in %s', singlePath);
    end
    matFiles = fullfile({d.folder}, {d.name})';
    savePath = singlePath;
elseif ~isempty(masterPath)
    d = dir(fullfile(masterPath, '**', '*breath_peak_data.mat'));
    if isempty(d)
        error('No breath_peak_data.mat found under %s', masterPath);
    end
    matFiles = fullfile({d.folder}, {d.name})';
    savePath = masterPath;
    fprintf('Found %d recordings with breath_peak_data.mat\n', numel(matFiles));
else
    error('Set either singlePath or masterPath.');
end

%% Collect all cycles across folders
nPts = 200;
phi_grid = linspace(0, 1, nPts)';
all_cycles_global = [];
totalCycles = 0;

for iF = 1:numel(matFiles)
    bp = load(matFiles{iF});
    breath = bp.breath(:);
    onsets = bp.insp_onset_idx(:);

    %% Hilbert transform
    analytic = hilbert(breath);
    inst_phase_raw = angle(analytic);        %#ok  — kept for reference

    nCycles = numel(onsets) - 1;
    cycles_this = NaN(nPts, nCycles);

    for k = 1:nCycles
        idx = onsets(k):onsets(k+1);
        if numel(idx) < 3; continue; end

        seg = breath(idx);
        seg_min = min(seg);
        seg_max = max(seg);
        if seg_max == seg_min; continue; end
        seg_norm = (seg - seg_min) / (seg_max - seg_min);

        phi_orig = linspace(0, 1, numel(idx))';
        cycles_this(:, k) = interp1(phi_orig, seg_norm, phi_grid, 'pchip');
    end

    all_cycles_global = [all_cycles_global, cycles_this]; %#ok<AGROW>
    totalCycles = totalCycles + nCycles;

    fprintf('  [%d/%d] %s  —  %d cycles\n', iF, numel(matFiles), fileparts(matFiles{iF}), nCycles);
end

%% Plot
figure('Color','w','Position',[100 100 700 450]); hold on;

nCol = size(all_cycles_global, 2);
cmap = parula(nCol);

for k = 1:nCol
    if all(isnan(all_cycles_global(:,k))); continue; end
    plot(phi_grid, all_cycles_global(:,k), '-', ...
        'Color', [cmap(k,:) 0.15], 'LineWidth', 0.6);
end

% cycle-averaged trace
avg_trace = mean(all_cycles_global, 2, 'omitnan');
plot(phi_grid, avg_trace, 'k-', 'LineWidth', 2.5);

xlabel('Breath-cycle phase \phi (\tau)');
ylabel('Normalised breath position');
if numel(matFiles) == 1
    title(sprintf('Breath position vs phase  |  %d cycles', totalCycles));
else
    title(sprintf('Breath position vs phase  |  %d recordings, %d cycles', ...
        numel(matFiles), totalCycles));
end
xlim([0 1]);
ylim([-0.05 1.05]);
set(gca, 'FontSize', 12, 'Box', 'off', 'TickDir', 'out');

%% Save
saveas(gcf, fullfile(savePath, 'breath_phase_overlay.png'));
exportgraphics(gcf, fullfile(savePath, 'breath_phase_overlay.pdf'), ...
    'ContentType','vector','BackgroundColor','none');
fprintf('Saved to %s\n', savePath);
