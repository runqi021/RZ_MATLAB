clear; close all; clc;
%% Batch_breath_summary.m  –  Breathing-rate & IBI summary across sessions
%  Loops through all subfolders of masterFolder looking for breath_peak_data.mat.
%  Produces:
%    1) KDE overlay of IBI distribution per session
%    2) Poincaré plot of IBI(n) vs IBI(n+1) per session (all individual IBIs)
%  Exports both figures as vector PDF (transparent background).
%  Saves pooled + per-session data to breath_summary_master.mat.
%
%  ── USER PARAMETERS ──────────────────────────────────────────────────────
masterFolder = 'D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\maybe breathing';   % <-- EDIT: path to master folder
skipIdx = [];   % sessions to skip. % Oxygen off to test phenotype
xl = [0 3];
fs=30;
% ──────────────────────────────────────────────────────────────────────────

%% 1 ── Discover sessions with breath_peak_data.mat
allMat = dir(fullfile(masterFolder, '**', '*breath_peak_data.mat'));

% % Filter out sessions in skipFolders
skipIdx = flip(skipIdx);
keep = true(numel(allMat),1);
keep(skipIdx) = false;
allMat = allMat(keep);

%skipFolders  = {};   % <-- EDIT: folder names to skip, e.g. {'session03','bad_run'}

% if ~isempty(skipFolders)
%     parentNames = cell(numel(allMat),1);
%     for ii = 1:numel(allMat)
%         [~, parentNames{ii}] = fileparts(allMat(ii).folder);
%     end
%     keep = ~ismember(parentNames, skipFolders);
%     skipped = parentNames(~keep);
%     allMat = allMat(keep);
%     for ii = 1:numel(skipped)
%         fprintf('  SKIPPED: %s\n', skipped{ii});
%     end
% end

nSess  = numel(allMat);
fprintf('Found %d sessions with breath_peak_data.mat\n', nSess);
assert(nSess > 0, 'No breath_peak_data.mat found under %s', masterFolder);

% Pre-allocate per-session storage
sessNames = cell(nSess,1);
rate_bpm   = cell(nSess,1);   % breaths/min (scalar per session)
ibi_all    = cell(nSess,1);   % IBI vectors (seconds)
recDur_s   = zeros(nSess,1);  % recording duration

for k = 1:nSess
    D = load(fullfile(allMat(k).folder, allMat(k).name));

    % Session label (parent folder name)
    sessNames{k} = allMat(k).folder;
    [~, sessNames{k}] = fileparts(allMat(k).folder);

    % Peak times in seconds
    peak_t = D.insp_onsets_t(:);

    % Recording duration from t_breath
    recDur_s(k) = D.t_breath(end) - D.t_breath(1);

    % Breaths per minute
    rate_bpm{k} = numel(peak_t) / (recDur_s(k) / 60);

    % Inter-breath intervals
    ibi_all{k} = diff(peak_t);

    fprintf('  [%2d] %-40s  rate=%.1f bpm  nPeaks=%d  dur=%.1fs\n', ...
        k, sessNames{k}, rate_bpm{k}, numel(peak_t), recDur_s(k));
end

% Pooled vectors
rate_bpm_vec = cell2mat(rate_bpm);         % [nSess x 1]
ibi_pooled   = vertcat(ibi_all{:});        % all IBIs concatenated

%% 2 ── Colour map (rainbow across sessions)
cmap = turbo(nSess);

%% 3 ── Figure 1: KDE of IBI distribution per session (overlaid)
fig1 = figure('Color','w','Position',[100 100 700 500]);
ax1  = axes(fig1, 'Color','w'); hold(ax1,'on');

% Common evaluation grid across all IBIs
xGrid = linspace(0, max(ibi_pooled)*2, 2000);

for k = 1:nSess
    ibi = ibi_all{k};
    ibi = ibi(isfinite(ibi) & ibi > 0);
    fprintf('[%02d] %-25s n=%d  median=%.3f  min=%.3f  max=%.3f\n', ...
    k, sessNames{k}, numel(ibi), median(ibi), min(ibi), max(ibi));
    
    %if numel(ibi) < 2, continue; end
    [fKDE, xi] = ksdensity(ibi, xGrid);
    fprintf('[%02d] n=%d  maxKDE=%.4g  minKDE=%.4g\n', ...
    k, numel(ibi), max(fKDE), min(fKDE));
    plot(ax1, xi, fKDE, '-', 'Color', cmap(k,:), 'LineWidth', 1.5, ...
        'DisplayName', sessNames{k});
end
xlim(ax1, xl);
xlabel(ax1, 'Inter-Breath Interval (s)');
ylabel(ax1, 'Density');
title(ax1, 'IBI Distribution Across Sessions');
box(ax1,'on');
axis(ax1, 'square');
colormap(ax1, cmap);
cb = colorbar(ax1);
clim(ax1, [1 nSess]);
cb.Label.String = 'Session #';
cb.Ticks = round(linspace(1, nSess, min(nSess,10)));

%% 4 ── Figure 2: Poincaré plot IBI(n) vs IBI(n+1)
fig2 = figure('Color','w','Position',[100 100 650 600]);
ax2  = axes(fig2, 'Color','w'); hold(ax2,'on');

for k = 1:nSess
    ibi = ibi_all{k};
    if numel(ibi) < 2, continue; end
    ran1 =[]; ran1 = (rand(size(ibi))-0.5).*2; ran1 = ran1./fs./2;
    ran2 =[]; ran2 = (rand(size(ibi))-0.5).*2; ran2 = ran2./fs./2;

    ibi_ran1 = ibi+ran1;
    ibi_ran2 = ibi+ran2;

    scatter(ax2, ibi_ran1(1:end-1), ibi_ran2(2:end), 12, cmap(k,:), ...
        'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', sessNames{k});
end

% Identity line
plot(ax2, xl, xl, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
xlim(ax2, xl); ylim(ax2, xl);
axis(ax2, 'square');

xlabel(ax2, 'IBI_n (s)');
ylabel(ax2, 'IBI_{n+1} (s)');
title(ax2, 'Poincaré Plot — Inter-Breath Interval');
box(ax2, 'on');

colormap(ax2, cmap);
cb2 = colorbar(ax2);
clim(ax2, [1 nSess]);
cb2.Label.String = 'Session #';
cb2.Ticks = round(linspace(1, nSess, min(nSess,10)));

%% 5 ── Export PDFs
outPDF1 = fullfile(masterFolder, 'breath_IBI_KDE_summary.pdf');
outPDF2 = fullfile(masterFolder, 'breath_IBI_poincare.pdf');

exportgraphics(fig1, outPDF1, 'ContentType','vector','BackgroundColor','none');
fprintf('Saved %s\n', outPDF1);

exportgraphics(fig2, outPDF2, 'ContentType','vector','BackgroundColor','none');
fprintf('Saved %s\n', outPDF2);

%% 6 ── Save master .mat
outMat = fullfile(masterFolder, 'breath_summary_master.mat');

% Per-session struct array
for k = nSess:-1:1
    sessions(k).name      = sessNames{k};
    sessions(k).rate_bpm  = rate_bpm{k};
    sessions(k).ibi       = ibi_all{k};
    sessions(k).recDur_s  = recDur_s(k);
    sessions(k).folder    = allMat(k).folder;
end

save(outMat, 'sessions', 'rate_bpm_vec', 'ibi_pooled', 'ibi_all', ...
             'sessNames', 'recDur_s', 'masterFolder');
fprintf('Saved %s\n', outMat);

fprintf('\n── Done. %d sessions summarised. ──\n', nSess);
