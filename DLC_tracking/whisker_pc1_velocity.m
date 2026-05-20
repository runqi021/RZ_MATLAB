% whisker_pc1_velocity.m
%  Per-whisker PCA on (x,y) → PC1 position & dPC1/dt velocity.
%  Compares traces, PSD, and cross-whisker correlation for both.
%
%  Dependencies: Chronux (mtspectrumc), Statistics Toolbox (pca)

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir  = fileparts(mfilename('fullpath'));
repoRoot   = fileparts(scriptDir);
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot, 'chronux_2_12')));

%% ========================= USER PARAMETERS ===========================
dlc_csv = "C:\260407_KA_electro_dorsal_whisking\basler_20260407_164056_run003DLC_Resnet50_whisker_dorsalApr8shuffle1_snapshot_best-10.csv";
fps       = 90;
point_ids = [1 2 3 4];
pmin      = 0.5;
fill_method = 'linear';

% Multitaper PSD
TW    = 4;
fpass = [0.1 44];

% Plotting
trace_xlim = [10 20];
doSave     = true;

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 7);

%% 1 -- Load DLC CSV
data = readmatrix(dlc_csv, 'NumHeaderLines', 3);
nTriplets = floor((size(data,2)-1)/3);
frame = data(:,1);
t = frame / fps;
T = numel(t);
nPts = numel(point_ids);
fs = fps;

X_all = nan(T, nPts);
Y_all = nan(T, nPts);
for k = 1:nPts
    pid = point_ids(k);
    cx = 2 + 3*(pid-1);
    cy = 3 + 3*(pid-1);
    cp = 4 + 3*(pid-1);
    x = data(:,cx); y = data(:,cy); p = data(:,cp);
    bad = (p < pmin) | isnan(p);
    x(bad) = NaN; y(bad) = NaN;
    x = fillmissing(x, fill_method);
    y = fillmissing(y, fill_method);
    X_all(:,k) = x; Y_all(:,k) = y;
end

% Bodypart names from CSV header
fid = fopen(dlc_csv,'r'); fgetl(fid); bp_line = fgetl(fid); fclose(fid);
bp_parts = strsplit(bp_line, ',');
bp_names = cell(1, nPts);
for k = 1:nPts
    col_idx = 2 + 3*(point_ids(k)-1);
    if col_idx <= numel(bp_parts), bp_names{k} = bp_parts{col_idx};
    else, bp_names{k} = sprintf('pt%d', point_ids(k)); end
end

%% 2 -- Per-whisker PCA
PC1 = nan(T, nPts);
PC2 = nan(T, nPts);
varExpl = nan(nPts, 2);

for k = 1:nPts
    xy = [detrend(X_all(:,k)), detrend(Y_all(:,k))];
    [coeff, sc, ~, ~, expl] = pca(xy);
    if coeff(2,1) < 0, sc(:,1) = -sc(:,1); end
    PC1(:,k) = sc(:,1) / std(sc(:,1));
    PC2(:,k) = sc(:,2) / std(sc(:,2));
    varExpl(k,:) = expl(1:2)';
end

%% 3 -- dPC1/dt (velocity)
dPC1 = diff(PC1) * fs;
t_d  = t(1:end-1) + 0.5/fs;

for k = 1:nPts
    dPC1(:,k) = (dPC1(:,k) - mean(dPC1(:,k))) / std(dPC1(:,k));
end

%% 4 -- Multitaper PSD
params_mt.Fs     = fs;
params_mt.tapers = [TW, 2*TW-1];
params_mt.pad    = 0;
params_mt.fpass  = fpass;
params_mt.err    = [2 0.05];

S_pc1 = nan([], nPts);
S_vel = nan([], nPts);
f_psd = []; f_vel = [];

for k = 1:nPts
    [Sk, fk] = mtspectrumc(PC1(:,k), params_mt);
    if isempty(f_psd)
        S_pc1 = nan(numel(Sk), nPts);
        f_psd = fk(:);
    end
    S_pc1(:,k) = Sk(:);
end

for k = 1:nPts
    [Sk, fk] = mtspectrumc(dPC1(:,k), params_mt);
    if isempty(f_vel)
        S_vel = nan(numel(Sk), nPts);
        f_vel = fk(:);
    end
    S_vel(:,k) = Sk(:);
end

%% 5 -- Figure
col = lines(nPts);
offsets = (nPts-1:-1:0) * 4;

fig = figure('Units','centimeters','Position',[1 1 30 28],'Color','w');

% Row 1: PC1 traces
ax1 = subplot(4,2,[1 2]);
hold on;
for k = 1:nPts
    plot(t, PC1(:,k) + offsets(k), 'Color', col(k,:), 'LineWidth', 0.3);
    text(t(1)-1, offsets(k), bp_names{k}, 'FontSize', 7, 'Color', col(k,:), ...
        'HorizontalAlignment','right', 'FontWeight','bold');
end
xlim(trace_xlim); ylabel('PC1 (z)'); title('PC1 (position)','FontSize',9);
set(ax1,'FontSize',7,'YTick',[]); box on;

% Row 2: dPC1/dt traces
ax2 = subplot(4,2,[3 4]);
hold on;
for k = 1:nPts
    plot(t_d, dPC1(:,k) + offsets(k), 'Color', col(k,:), 'LineWidth', 0.3);
    text(t_d(1)-1, offsets(k), bp_names{k}, 'FontSize', 7, 'Color', col(k,:), ...
        'HorizontalAlignment','right', 'FontWeight','bold');
end
xlim(trace_xlim); ylabel('dPC1/dt (z)'); title('dPC1/dt (velocity)','FontSize',9);
set(ax2,'FontSize',7,'YTick',[]); box on;
xlabel('Time (s)');

% Row 3 left: PSD — position
ax3 = subplot(4,2,5);
hold on;
for k = 1:nPts
    plot(f_psd, 10*log10(S_pc1(:,k)), 'Color', col(k,:), 'LineWidth', 1.2);
end
set(ax3,'XScale','log'); xlim(fpass);
xlabel('Hz'); ylabel('Power (dB)');
title('PC1 PSD (position)','FontSize',9);
legend(bp_names, 'FontSize', 6, 'Location','northeast');
set(ax3,'FontSize',7); box on;

% Row 3 right: PSD — velocity
ax4 = subplot(4,2,6);
hold on;
for k = 1:nPts
    plot(f_vel, 10*log10(S_vel(:,k)), 'Color', col(k,:), 'LineWidth', 1.2);
end
set(ax4,'XScale','log'); xlim(fpass);
xlabel('Hz'); ylabel('Power (dB)');
title('dPC1/dt PSD (velocity)','FontSize',9);
legend(bp_names, 'FontSize', 6, 'Location','northeast');
set(ax4,'FontSize',7); box on;

% Row 4 left: PC1 correlation
ax5 = subplot(4,2,7);
C_pos = corrcoef(PC1);
imagesc(C_pos); colorbar; caxis([0 1]);
set(ax5,'XTick',1:nPts,'YTick',1:nPts,'XTickLabel',bp_names,'YTickLabel',bp_names,'FontSize',7);
title('PC1 correlation','FontSize',9); axis square;

% Row 4 right: dPC1/dt correlation
ax6 = subplot(4,2,8);
C_vel = corrcoef(dPC1);
imagesc(C_vel); colorbar; caxis([0 1]);
set(ax6,'XTick',1:nPts,'YTick',1:nPts,'XTickLabel',bp_names,'YTickLabel',bp_names,'FontSize',7);
title('dPC1/dt correlation','FontSize',9); axis square;

annotation(fig, 'textbox', [0.01 0.96 0.98 0.04], ...
    'String', 'PC1 position vs velocity (dPC1/dt) — per whisker', ...
    'EdgeColor','none','FontSize',9,'FontWeight','bold', ...
    'HorizontalAlignment','center');

%% 6 -- Save
if doSave
    [folderPath, csvStem] = fileparts(dlc_csv);
    shortName = char(regexprep(csvStem, 'DLC.*$', ''));
    outPdf = fullfile(folderPath, [shortName 'whisker_PC1_vs_velocity.pdf']);
    print(fig, '-dpdf', '-bestfit', char(outPdf));
    fprintf('Saved: %s\n', outPdf);
end

fprintf('\nPer-whisker PC1 variance explained:\n');
for k = 1:nPts
    fprintf('  %s: PC1=%.1f%%  PC2=%.1f%%\n', bp_names{k}, varExpl(k,1), varExpl(k,2));
end
fprintf('Done.\n');
