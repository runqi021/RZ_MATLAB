%% phase_precession.m  --  Ca²⁺ spike phase precession relative to breathing
%  For each spiking ROI, produces one figure with 3 subplots:
%
%    (1) Phase vs Breath Cycle Number   — drift over time
%    (2) Phase vs Firing Rate           — rate dependence
%    (3) Phase vs Breath Cycle Duration — speed dependence
%
%  Phase = (spike_frame − cycle_start) / (cycle_end − cycle_start) ∈ [0,1)
%         0 = inspiration onset, 1 = next inspiration onset
%
%  Saves .fig + .pdf per ROI into folderPath/phase_precession/
%  ── USER PARAMETERS ─────────────────────────────────────────────────────
clear; close all; clc;
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

folderPath = "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing\left_pFN_roi1_z0_3x_00001";

fps_img = 30;
nDrop   = 30;
phase_range = [-3*pi/2  3*pi/2];    % display range (points wrap to fill)
% ──────────────────────────────────────────────────────────────────────────

%% ---- Auto-find files ----
sam_hits = dir(fullfile(folderPath, '*cpSAM_output.mat'));
bp_hits  = dir(fullfile(folderPath, '*breath_peak_data.mat'));

if isempty(sam_hits) || isempty(bp_hits)
    error('Missing cpSAM_output.mat or breath_peak_data.mat in:\n  %s', folderPath);
end
if ~isfile(fullfile(folderPath, 'ca_spike_data.mat'))
    error('Missing ca_spike_data.mat in:\n  %s', folderPath);
end

SAM_file = fullfile(folderPath, sam_hits(1).name);

%% ---- Load SAM → dFF ----
SAMload = load(SAM_file);
F = SAMload.F;
F(1:nDrop, :) = [];
dFFout = helper.dFF_RZ(F);
dFF = dFFout.dFF;
[T, N_roi] = size(dFF);
fs = fps_img;

fprintf('dFF: [%d x %d]\n', T, N_roi);

%% ---- Load calcium spikes ----
in = load(fullfile(folderPath, 'ca_spike_data.mat'));
ifSpk     = in.ifSpike;
roiSpk_id = find(ifSpk == 1);
ca_spk_data = in.roi_spikes;
nSpkROI   = numel(roiSpk_id);

ca_spk_id = cell(1, nSpkROI);
for k = 1:nSpkROI
    ca_spk_id{k} = ca_spk_data(roiSpk_id(k)).spike_idx;
end

fprintf('Spiking ROIs: %d of %d\n', nSpkROI, N_roi);
assert(nSpkROI > 0, 'No spiking ROIs found.');

%% ---- Load breath peaks ----
breath_peak = load(fullfile(folderPath, bp_hits(1).name));
insp_onsets = breath_peak.insp_onset_idx;

insp_onsets(insp_onsets < nDrop) = [];
insp_onsets = insp_onsets - nDrop;
insp_onsets = sort(insp_onsets(:));

b_frames = insp_onsets;
nCycles  = numel(b_frames) - 1;
cycle_dur = diff(b_frames) / fs;   % seconds

assert(nCycles >= 2, 'Fewer than 3 inspiration onsets.');
fprintf('Breath cycles: %d  (mean %.3f s)\n', nCycles, mean(cycle_dur));

%% ---- Output directory ----
outDir = fullfile(folderPath, 'phase_precession');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ---- Loop over spiking ROIs ----
for k = 1:nSpkROI
    roi_id = roiSpk_id(k);
    ca_ev  = ca_spk_id{k}(:);

    if isempty(ca_ev), continue; end

    % Assign each Ca spike to a breath cycle
    [~, ~, bin_idx] = histcounts(ca_ev, b_frames);
    keep = (bin_idx >= 1) & (bin_idx <= nCycles);
    ca_in  = ca_ev(keep);
    cyc_id = bin_idx(keep);       % which cycle each spike belongs to

    if numel(ca_in) < 3
        fprintf('ROI#%02d: <3 spikes in breath cycles, skipping.\n', roi_id);
        continue;
    end

    % Phase within cycle: canonical [0, 2*pi), then wrap copies into phase_range
    phase_frac = (ca_in - b_frames(cyc_id)) ./ (b_frames(cyc_id + 1) - b_frames(cyc_id));
    phase_base = phase_frac * 2*pi;               % [0, 2*pi)
    % Create wrapped copies so scatter fills the full display range
    phase_all = [phase_base - 2*pi; phase_base; phase_base + 2*pi];
    in_range  = phase_all >= phase_range(1) & phase_all <= phase_range(2);
    % Expand matching x-variables too
    cycle_number_3 = repmat(cyc_id, [3 1]);
    spike_dur_3    = repmat(cycle_dur(cyc_id), [3 1]);
    % (spike_rate computed below, also tripled)

    % Firing rate per cycle
    spike_rate = NaN(size(ca_in));
    for c = 1:nCycles
        mask = (cyc_id == c);
        if any(mask)
            spike_rate(mask) = sum(mask) / cycle_dur(c);
        end
    end
    spike_rate_3 = repmat(spike_rate, [3 1]);

    % Keep only wrapped points inside phase_range
    phase_plot   = phase_all(in_range);
    cyc_plot     = cycle_number_3(in_range);
    rate_plot    = spike_rate_3(in_range);
    dur_plot     = spike_dur_3(in_range);

    nEv = numel(ca_in);

    % ---- Y-axis ticks for current phase_range ----
    all_ticks  = (-2*pi : pi/2 : 2*pi)';
    all_labels = {'-2\pi','-3\pi/2','-\pi','-\pi/2','0', ...
                  '\pi/2','\pi','3\pi/2','2\pi'};
    keep_t = all_ticks >= phase_range(1) & all_ticks <= phase_range(2);
    ytk    = all_ticks(keep_t);
    ytkl   = all_labels(keep_t);

    % ---- Figure ----
    fig = figure('Color','w', 'Position', [50 100 1400 400], ...
        'Name', sprintf('Phase precession ROI#%02d', roi_id));

    % (1) Phase vs Breath Cycle Number
    ax1 = subplot(1,3,1);
    scatter(ax1, cyc_plot, phase_plot, 15, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.4);
    xlabel('Breath cycle #');
    ylabel('Phase (rad)');
    title('Phase vs Cycle #');
    ylim(phase_range); xlim([0 nCycles+1]);
    set(ax1, 'YTick', ytk, 'YTickLabel', ytkl);
    set(ax1, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');

    % (2) Phase vs Firing Rate
    ax2 = subplot(1,3,2);
    scatter(ax2, rate_plot, phase_plot, 15, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.4);
    xlabel('Ca^{2+} spike rate (spikes/s)');
    ylabel('Phase (rad)');
    title('Phase vs Firing Rate');
    ylim(phase_range);
    set(ax2, 'YTick', ytk, 'YTickLabel', ytkl);
    set(ax2, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');

    % (3) Phase vs Breath Cycle Duration
    ax3 = subplot(1,3,3);
    scatter(ax3, dur_plot, phase_plot, 15, [0.3 0.3 0.3], 'filled', 'MarkerFaceAlpha', 0.4);
    xlabel('Breath cycle duration (s)');
    ylabel('Phase (rad)');
    title('Phase vs Cycle Duration');
    ylim(phase_range);
    set(ax3, 'YTick', ytk, 'YTickLabel', ytkl);
    set(ax3, 'FontSize', 10, 'Box', 'off', 'TickDir', 'out');

    sgtitle(sprintf('Phase Precession — ROI#%02d  (%d spikes, %d cycles)', ...
        roi_id, nEv, nCycles), 'FontSize', 12, 'FontWeight', 'bold');

    % ---- Save ----
    savefig(fig, fullfile(outDir, sprintf('phase_precession_ROI%02d.fig', roi_id)));
    exportgraphics(fig, fullfile(outDir, sprintf('phase_precession_ROI%02d.pdf', roi_id)), ...
        'ContentType', 'vector', 'BackgroundColor', 'none');

    fprintf('ROI#%02d: %d spikes in %d cycles\n', roi_id, nEv, nCycles);
end

fprintf('\nDone. Figures saved to: %s\n', outDir);
