%% breathing_is_motion_analysis.m
% Single-folder analysis driven by the shift-based breathing peaks produced
% by breathing_is_motion_gui.m.
%
% Produces:
%   1) |shift| breathing trace + stacked dFF (with inspiration markers)
%   2) Inspiration-triggered dFF overlay (per ROI: trials + mean)
%   3) breath <-> dFF coherence across frequency (Chronux coherencyc),
%      plus the breath power spectrum for a breathing-rate readout
%
% Inputs auto-found in folderPath:
%   *_breath_peak_data.mat   (GUI output: breath, t_breath, insp_onset_idx, ...)
%   *cpSAM_output.mat        (F  [T x N])
%   *_meta.mat               (fps)

close all; clear; clc;

%% ---------------- user params ----------------
folderPath     = "C:\Users\Admin\Desktop\260521_ChAT_soma_G8s\roi5_7x_x-1200y200z-30_3000f_23lp_00001";
BaselineWinSec = 10;        % sliding-median baseline for dFF
win_sec        = 1;       % inspiration-triggered window (+/- s)
% Coherence / spectrum multitaper (breathing band, ~99 s record)
TW    = 4;
fpass = [0.05 5];             % Hz
% Calcium event detection (for event-train coherence)
minProm_ca   = 0.1;         % min dFF prominence for a calcium peak
minDist_ca_s = 0.2;         % min seconds between calcium peaks
% ----------------------------------------------

%% Path setup (helper. + Chronux)
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

%% Locate input files
peak_hits = dir(fullfile(folderPath, '*_breath_peak_data.mat'));   % breath train (GUI)
ca_hits   = dir(fullfile(folderPath, 'ca_spike_data.mat'));        % calcium train (GUI)
dff_hits  = dir(fullfile(folderPath, '*_dFF.mat'));                % dFF used for spikes
sam_hits  = dir(fullfile(folderPath, '*cpSAM_output.mat'));        % fallback for dFF
meta_hits = dir(fullfile(folderPath, '*_meta.mat'));

if isempty(peak_hits)
    error(['No *_breath_peak_data.mat in:\n  %s\n' ...
           'Run breathing_is_motion_gui.m first to detect inspiration peaks.'], folderPath);
end
if isempty(meta_hits)
    error('No *_meta.mat in:\n  %s', folderPath);
end
if isempty(dff_hits) && isempty(sam_hits)
    error('No *_dFF.mat or *cpSAM_output.mat (need one for dFF) in:\n  %s', folderPath);
end

peakFile = fullfile(folderPath, peak_hits(1).name);
metaFile = fullfile(folderPath, meta_hits(1).name);

fprintf('Peaks : %s\n', peak_hits(1).name);
fprintf('Meta  : %s\n', meta_hits(1).name);

%% Load
BP = load(peakFile);                 % breath, t_breath, insp_onset_idx, ...
Me = load(metaFile);
fps = Me.fps;

breath      = BP.breath(:);                  % z-scored shift breathing trace
insp_onsets = BP.insp_onset_idx(:);          % frame indices
Sb          = BP.insp_onsets_train(:);       % binary train

%% dFF  (prefer the *_dFF.mat used for calcium spike detection; else recompute)
if ~isempty(dff_hits)
    D   = load(fullfile(folderPath, dff_hits(1).name), 'dFF');
    dFF = D.dFF;
    fprintf('dFF   : %s\n', dff_hits(1).name);
else
    SAMload = load(fullfile(folderPath, sam_hits(1).name));
    dFFout  = helper.dFF_RZ(SAMload.F, 'FPS', fps, 'BaselineWinSec', BaselineWinSec);
    dFF     = dFFout.dFF;
    fprintf('dFF   : recomputed from %s\n', sam_hits(1).name);
end

%% Align breath and dFF (both post-toss; allow 1-frame trim)
T = min(numel(breath), size(dFF,1));
if numel(breath) ~= size(dFF,1)
    warning('breath (%d) and dFF (%d) length differ; trimming to %d.', ...
            numel(breath), size(dFF,1), T);
end
breath = breath(1:T);
Sb     = Sb(1:T);
dFF    = dFF(1:T, :);
insp_onsets = insp_onsets(insp_onsets >= 1 & insp_onsets <= T);

[~, N_roi] = size(dFF);
t_img = (0:T-1)' / fps;
roi_ids = arrayfun(@(k) sprintf('%02d', k), 1:N_roi, 'UniformOutput', false);

fprintf('Aligned: T=%d frames, N=%d ROIs, fps=%g, %d inspiration onsets\n', ...
        T, N_roi, fps, numel(insp_onsets));

%% Calcium event train per ROI
% Primary: use the curated spikes from calcium_spike_gui.m (ca_spike_data.mat).
% Fallback: auto-detect with islocalmax only if that file is absent.
ca_idx  = cell(N_roi,1);
ifSpike = true(1, N_roi);   % whether each ROI is a spiking cell (from GUI)

if ~isempty(ca_hits)
    CA = load(fullfile(folderPath, ca_hits(1).name), 'roi_spikes', 'ifSpike');
    rs = CA.roi_spikes;
    if numel(rs) ~= N_roi
        warning('ca_spike_data has %d ROIs but dFF has %d; matching by index.', ...
                numel(rs), N_roi);
    end
    if isfield(CA,'ifSpike') && numel(CA.ifSpike) == N_roi
        ifSpike = logical(CA.ifSpike(:)');
    end
    for i = 1:min(N_roi, numel(rs))
        si = rs(i).spike_idx(:);
        ca_idx{i} = si(si >= 1 & si <= T);
    end
    fprintf('Calcium train: curated %s (%s)\n', ca_hits(1).name, ...
            mat2str(cellfun(@numel, ca_idx)'));
else
    warning(['No ca_spike_data.mat found - falling back to islocalmax ' ...
             'auto-detection. Run calcium_spike_gui.m for curated spikes.']);
    minDist_ca = max(1, round(minDist_ca_s * fps));
    for i = 1:N_roi
        tf = islocalmax(dFF(:,i), 'MinProminence', minProm_ca, ...
                        'MinSeparation', minDist_ca);
        ca_idx{i} = find(tf);
    end
    fprintf('Calcium train: islocalmax %s\n', mat2str(cellfun(@numel, ca_idx)'));
end

%% |shift| in microns (for display) recomputed from the MC output the GUI used
fp = BP.findpeak_params;
mcFile = '';
if isfield(fp,'mc_file') && ischar(fp.mc_file) && isfile(fp.mc_file)
    mcFile = fp.mc_file;
else
    mc_hits = dir(fullfile(folderPath, '*_preproc_MC_output.mat'));
    mc_hits = mc_hits(~contains({mc_hits.name}, '_MC_MC_output'));
    if ~isempty(mc_hits), mcFile = fullfile(folderPath, mc_hits(1).name); end
end

if isempty(mcFile)
    warning('No MC output found; |shift| panel falls back to z-scored breath.');
    shift_um = breath;  shift_lab = 'Breathing (z)';
else
    if isfield(fp,'pixelSize_um') && isfinite(fp.pixelSize_um)
        px = fp.pixelSize_um;
    elseif isfield(Me,'pixelSize_um')
        px = Me.pixelSize_um;
    else
        px = 1;  warning('No pixelSize_um; |shift| shown in pixels.');
    end
    MC = load(mcFile);  sh = MC.mcOut.shifts;  nF = numel(sh);  yx = zeros(nF,2);
    for k = 1:nF, s = squeeze(sh(k).shifts); yx(k,:) = s(:)'; end
    if isfield(fp,'shift_source') && strcmpi(fp.shift_source,'MC1+MC2')
        mc2 = strrep(mcFile,'_preproc_MC_output.mat','_preproc_MC_MC_output.mat');
        if isfile(mc2)
            MC2 = load(mc2); sh2 = MC2.mcOut.shifts;
            for k = 1:min(nF,numel(sh2)), s2 = squeeze(sh2(k).shifts); yx(k,:) = yx(k,:) + s2(:)'; end
        end
    end
    shift_um  = sqrt(sum((yx*px).^2, 2));
    shift_lab = '|shift| (\mum)';
end
% length-match to T
if numel(shift_um) >= T, shift_um = shift_um(1:T); else, shift_um(end+1:T,1) = NaN; end

%% ===== FIGURE 1: |shift| (um) + stacked dFF =====
amp = max(dFF(:)) - min(dFF(:));
if amp == 0, amp = 1; end
gap = 0.8 * amp;

Y = zeros(T, N_roi);
for k = 1:N_roi
    Y(:,k) = dFF(:,k) + gap*(N_roi - k);
end
yt  = 0:gap:gap*(N_roi-1);
ytl = roi_ids(end:-1:1);

idx_b = find(Sb ~= 0);

fig1 = figure('Position',[200 200 1000 750],'Color','w');
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(1);
plot(ax1, t_img, shift_um, 'k', 'LineWidth', 1); hold(ax1,'on');
if ~isempty(idx_b)
    plot(ax1, t_img(idx_b), shift_um(idx_b), 'ro', 'MarkerFaceColor','r','MarkerSize',3);
end
xlim(ax1, [t_img(1) t_img(end)]);
xlabel(ax1, 'Time (s)'); ylabel(ax1, shift_lab);
title(ax1, 'MC |shift| with detected inspirations');
box(ax1,'off'); hold(ax1,'off');

ax2 = nexttile(2, [2 1]); hold(ax2,'on');
for kk = 1:N_roi
    plot(ax2, t_img, Y(:,kk), 'Color', 'k', 'LineWidth', 0.5);
    ci = ca_idx{kk};                     % calcium peaks (QC overlay)
    if ~isempty(ci)
        plot(ax2, t_img(ci), Y(ci,kk), 'r.', 'MarkerSize', 7);
    end
end
xlim(ax2, [t_img(1) t_img(end)]);
ylim(ax2, [-gap, gap*(N_roi-1)+gap]);
yticks(ax2, yt); yticklabels(ax2, ytl);
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'ROI#');
title(ax2, '\DeltaF/F');
box(ax2,'off');

% 20% dF/F vertical scale bar (bot-right)
dFF_scale = 0.2;
xb = t_img(end) - 0.01*(t_img(end) - t_img(1));
%yb = gap*(N_roi-1);
yb = gap/2.3;
plot(ax2, [xb xb], [yb yb+dFF_scale], 'k', 'LineWidth', 3);
text(ax2, xb, yb + dFF_scale/2, sprintf('%g%% \\DeltaF/F  ', dFF_scale*100), ...
     'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',9);
hold(ax2,'off');

outFig = fullfile(folderPath, 'shift_breath_N_dFF.png');
exportgraphics(fig1, outFig, 'Resolution',150);
fprintf('Saved: %s\n', outFig);

%% ===== FIGURE 2: inspiration-triggered dFF overlay =====
win   = round(win_sec * fps);
t_evt = (-win:win) / fps;
insp_ev = insp_onsets(insp_onsets > win & insp_onsets <= T - win);
nEv = numel(insp_ev);
fprintf('Inspiration-triggered: %d usable onsets (of %d).\n', nEv, numel(insp_onsets));

dff_evt_seg_all = cell(N_roi,1);
for i = 1:N_roi
    seg = zeros(nEv, 2*win+1);
    for e = 1:nEv
        c = insp_ev(e);
        seg(e,:) = dFF(c-win : c+win, i);
    end
    dff_evt_seg_all{i} = seg;
end

fig2 = figure('Name', sprintf('Inspiration-triggered dFF (n=%d)', nEv), 'Color','w');
nCol = min(N_roi, 5);
nRow = ceil(N_roi / nCol);
for i = 1:N_roi
    seg = dff_evt_seg_all{i};
    subplot(nRow, nCol, i); hold on;
    if ~isempty(seg)
        plot(t_evt, seg', 'Color', [0.8 0.8 0.8], 'LineWidth', 0.7);
        plot(t_evt, mean(seg,1), 'k', 'LineWidth', 2);
    end
    xline(0, '--', 'Color', [0.4 0.4 0.4]);
    yline(0, '--', 'Color', [0.6 0.6 0.6]);
    xlabel('Time from insp. (s)'); ylabel('\DeltaF/F');
    title(sprintf('ROI#%d', i));
    axis square; hold off;
end

outFig = fullfile(folderPath, 'shift_triggered_avg.png');
exportgraphics(fig2, outFig, 'Resolution',150);
fprintf('Saved: %s\n', outFig);

%% ===== FIGURE 3: breath-train <-> calcium-train coherence (point process) =====
% Event-train coherence: inspiration onsets vs per-ROI calcium peaks, both as
% point processes (spike times in seconds) via Chronux coherencypt. A common
% time grid (t_grid) is passed so every ROI shares one frequency axis.
params_coh.Fs     = fps;
params_coh.tapers = [TW, 2*TW - 1];
params_coh.pad    = 0;
params_coh.fpass  = fpass;
params_coh.err    = [2 0.005];
params_coh.trialave = 0;

breath_t = t_img(insp_onsets);          % inspiration onset times (s)
t_grid   = t_img(:)';                   % common time grid for the tapers

coh = struct('C',[],'phi',[],'f',[],'confC',[],'n_ca',[]);
fig3 = figure('Name','breath-train <-> calcium-train coherence', 'Color','w');
for i = 1:N_roi
    ca_t = t_img(ca_idx{i});            % calcium peak times (s)
    coh(i).n_ca = numel(ca_t);

    subplot(nRow, nCol, i); hold on;
    if numel(ca_t) >= 2 && numel(breath_t) >= 2
        % Chronux wants point times as structs (.times) so it treats them as
        % point processes rather than length-matched continuous arrays.
        d_br.times = breath_t(:);
        d_ca.times = ca_t(:);
        % Request 9 outputs so confC (significance level) is computed.
        [C, phi, ~, ~, ~, f, ~, confC, ~] = ...
            coherencypt(d_br, d_ca, params_coh, 0, t_grid);
        coh(i).C = C; coh(i).phi = phi; coh(i).f = f; coh(i).confC = confC;
        plot(f, C, 'k', 'LineWidth', 1.2);
        yline(confC, 'r--', 'LineWidth', 0.8);   % significance level
    else
        text(0.5, 0.5, 'too few events', 'Units','normalized', ...
             'HorizontalAlignment','center');
    end
    xlim(fpass); ylim([0 1]);
    xlabel('Frequency (Hz)'); ylabel('Coherence');
    ttl = sprintf('#%d (n_{Ca}=%d, n_{br}=%d)', i, numel(ca_t), numel(breath_t));
    %if ~ifSpike(i), ttl = [ttl ' [non-spiking]']; end
    title(ttl);
    axis square; box on; hold off;
end

outFig = fullfile(folderPath, 'shift_breath_coherence.png');
exportgraphics(fig3, outFig, 'Resolution',150);
fprintf('Saved: %s\n', outFig);

%% ===== Breath power spectrum (breathing-rate readout) =====
params_spec.Fs     = fps;
params_spec.tapers = [TW, 2*TW - 1];
params_spec.pad    = 0;
params_spec.fpass  = fpass;
params_spec.err    = [2 0.05];

[S_breath, f_breath] = mtspectrumc(breath, params_spec);

fig4 = figure('Color','w');
plot(f_breath, S_breath, 'k', 'LineWidth', 1.2);
xlabel('Frequency (Hz)'); ylabel('Power');
xlim(fpass); axis square; grid on;
title('|shift| breathing spectrum');

outFig = fullfile(folderPath, 'breath_spectrum_shift.png');
exportgraphics(fig4, outFig, 'Resolution',150);
fprintf('Saved: %s\n', outFig);

%% ===== Save master =====
params = struct('BaselineWinSec',BaselineWinSec,'win_sec',win_sec, ...
                'TW',TW,'fpass',fpass,'fps',fps, ...
                'minProm_ca',minProm_ca,'minDist_ca_s',minDist_ca_s);
spectrum = struct('S',S_breath,'f',f_breath);

outMat = fullfile(folderPath, 'shift_breath_master.mat');
save(outMat, 'folderPath','fps','breath','shift_um','t_img','dFF','Sb', ...
             'insp_onsets','ca_idx','ifSpike','dff_evt_seg_all','t_evt', ...
             'coh','spectrum','params');
fprintf('Saved master: %s\n', outMat);
