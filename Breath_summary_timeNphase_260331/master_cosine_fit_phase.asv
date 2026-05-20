% master_cosine_fit_phase.m
%  Overlays per-ROI cosine fits of P(spike | breath phase) across experiments
%  onto a single square plot.
%
%  For each ROI, cosine: y = a + b*cos(phase - mu)
%    mu = circular mean of that ROI's spike phases
%    a  = mean of its phase-binned P(spike)
%    b  = 2 * mean( P(spike) .* cos(ctrs - mu) )
%
%  Reads per-session raw files (ca_spike_data.mat, breath_peak_data.mat)
%  from each expFolder and recomputes spike phases using the same formula
%  as breath_combined_summary_260329.m. Pairs with all_is_sig_spk in
%  coherence_data.mat (matched by iteration order) for sig filtering.

clear; close all; clc;

%% ========================= USER PARAMETERS =============================
expFolders = {
    "D:\251124_live_vglut2_soma_g8s+cy5\phys\breathing"
    "D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\maybe breathing"
    "C:\Users\Admin\Desktop\260330_sst_soma_g8s\phys\maybe_breathing"
    "D:\batch_dffQC_test_260325\260323_vgat_g8s\phys\processed"
    "D:\batch_dffQC_test_260325\260224_vglut2_soma_g8s\phys\processed\breathing"
};

doSave     = true;
outDir     = "C:\Users\Admin\Desktop\RZ_MATLAB\Breath_summary_timeNphase_260331\master_summary";
sigOnly    = false;    % overlay only coherence-sig ROIs (from coherence_data.mat)
minSpikes  = 5;      % minimum spikes for a fit
nBinsHist  = 24;      % phase histogram bins (match combined summary)
nDrop      = 30;      % frames dropped from start (must match combined summary)
fps_img    = 30;      % fallback fps if detect_session_fps fails
% ========================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

% Ensure repo root on path (for detect_session_fps, helper package, etc.)
thisFile = mfilename('fullpath');
[thisDir,~] = fileparts(thisFile);
repoRoot = fileparts(thisDir);
addpath(repoRoot);

edges_hist = linspace(0, 2*pi, nBinsHist + 1);
ctrs_hist  = (edges_hist(1:end-1) + edges_hist(2:end)) / 2;

%% 1 -- Per-expFolder: iterate sessions, compute per-ROI spike phases
nExp = numel(expFolders);
all_ph_cell = {};
all_sig_vec = logical([]);
all_exp_vec = [];
expNames    = cell(nExp, 1);

% Breath waveform pool: per-cycle, phase-normalized, each cycle scaled [0,1]
% One entry per expFolder so we can plot a separate mean per folder.
nPhaseBins = 200;
cycles_per_exp = cell(nExp, 1);

for ee = 1:nExp
    inputPath = char(expFolders{ee});
    [~, folderName] = fileparts(inputPath);
    if isempty(folderName), [~, folderName] = fileparts(fileparts(inputPath)); end
    expNames{ee} = folderName;

    % Discover sessions (matches breath_combined_summary logic)
    if exist(fullfile(inputPath, 'ca_spike_data.mat'), 'file') == 2
        allMat = dir(fullfile(inputPath, 'ca_spike_data.mat'));
    else
        allMat = dir(fullfile(inputPath, '**', 'ca_spike_data.mat'));
    end
    nSess = numel(allMat);
    if nSess == 0
        fprintf('WARNING: no ca_spike_data.mat under %s\n', inputPath);
        continue
    end

    % Load fov_map.mat: fov_map{k,1}=sessName, fov_map{k,2}=sig ROI IDs
    fovMat = fullfile(inputPath, 'breath_combined_summary', 'fov_map.mat');
    has_fov = isfile(fovMat);
    if has_fov
        Sf = load(fovMat);
        fov_map = Sf.fov_map;
        fov_sess_names = string(fov_map(:,1));
    else
        fprintf('WARNING: no fov_map.mat in %s — sig filter disabled\n', inputPath);
        fov_map = {};
        fov_sess_names = strings(0,1);
    end

    exp_ph_cell = {};
    exp_sig     = logical([]);
    roi_ctr     = 0;

    for kk = 1:nSess
        folderPath = allMat(kk).folder;
        [~, sessName] = fileparts(folderPath);

        try
            bp_file = dir(fullfile(folderPath, '*breath_peak_data.mat'));
            if isempty(bp_file)
                fprintf('  [%s] no breath_peak_data.mat, skip\n', sessName);
                continue
            end

            caLoad = load(fullfile(folderPath, 'ca_spike_data.mat'));
            ifSpk  = caLoad.ifSpike;
            caData = caLoad.roi_spikes;
            roiSpk_id = find(ifSpk == 1);
            if isempty(roiSpk_id), continue; end

            bpLoad = load(fullfile(bp_file(1).folder, bp_file(1).name));
            insp_onsets_raw = sort(bpLoad.insp_onset_idx(:));
            insp_onsets = insp_onsets_raw;
            insp_onsets(insp_onsets <= nDrop) = [];
            insp_onsets = insp_onsets - nDrop;
            if numel(insp_onsets) < 2, continue; end

            % ---- Pool per-cycle breathing waveforms (phase-normalized) ----
            if isfield(bpLoad, 'breath') && numel(insp_onsets_raw) >= 2
                breath_raw = bpLoad.breath(:);
                phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1);
                phase_1cyc = phase_1cyc(1:end-1);
                nB = numel(breath_raw);
                for c = 1:numel(insp_onsets_raw)-1
                    i1 = insp_onsets_raw(c); i2 = insp_onsets_raw(c+1);
                    if i1 >= 1 && i2 <= nB && i2 > i1
                        chunk = breath_raw(i1:i2);
                        orig_ph = linspace(0, 2*pi, numel(chunk));
                        cyc = interp1(orig_ph, chunk, phase_1cyc, 'linear');
                        mn = min(cyc); mx = max(cyc);
                        if mx > mn
                            cycles_per_exp{ee}(end+1, :) = (cyc - mn) / (mx - mn); %#ok<SAGROW>
                        end
                    end
                end
            end

            % fs — use detect_session_fps if available, else fallback
            if exist('detect_session_fps','file') == 2
                try
                    fs = detect_session_fps(folderPath, fps_img);
                catch
                    fs = fps_img;
                end
            else
                fs = fps_img;
            end

            % Look up sig ROI IDs for this session from fov_map (by sessName)
            if has_fov
                mm = find(fov_sess_names == string(sessName), 1);
                if isempty(mm)
                    sig_ids_this_sess = [];
                else
                    sig_ids_this_sess = fov_map{mm, 2};
                end
            else
                sig_ids_this_sess = [];
            end

            % Compute per-ROI spike phases (same formula as combined summary)
            for jj = 1:numel(roiSpk_id)
                roi_id = roiSpk_id(jj);
                ca_ev = caData(roi_id).spike_idx(:);

                dt_l = nan(numel(ca_ev), 1);
                dt_n = nan(numel(ca_ev), 1);
                for ss = 1:numel(ca_ev)
                    c   = ca_ev(ss);
                    bef = insp_onsets(insp_onsets <= c);
                    aft = insp_onsets(insp_onsets >  c);
                    if ~isempty(bef), dt_l(ss) = (c - bef(end)) / fs; end
                    if ~isempty(aft), dt_n(ss) = (aft(1) - c)   / fs; end
                end
                keep_sp = ~isnan(dt_l) & ~isnan(dt_n);
                if any(keep_sp)
                    ph = 2*pi * (dt_l(keep_sp) + 0.5/fs) ./ ...
                         (dt_l(keep_sp) + dt_n(keep_sp));
                    ph = mod(ph, 2*pi);
                else
                    ph = [];
                end

                roi_ctr = roi_ctr + 1;
                exp_ph_cell{roi_ctr,1} = ph; %#ok<SAGROW>
                if has_fov
                    exp_sig(roi_ctr,1) = any(sig_ids_this_sess == roi_id); %#ok<SAGROW>
                else
                    exp_sig(roi_ctr,1) = true; %#ok<SAGROW>
                end
            end
        catch ME
            fprintf('  [%s] ERROR: %s\n', sessName, ME.message);
            continue
        end
    end

    all_ph_cell = [all_ph_cell; exp_ph_cell];                %#ok<AGROW>
    all_sig_vec = [all_sig_vec; exp_sig];                    %#ok<AGROW>
    all_exp_vec = [all_exp_vec; repmat(ee, numel(exp_ph_cell), 1)]; %#ok<AGROW>

    fprintf('[%s] %d ROIs (%d sig)\n', folderName, numel(exp_ph_cell), sum(exp_sig));
end

nTotal = numel(all_ph_cell);
fprintf('\n=== %d ROIs across %d experiments ===\n', nTotal, nExp);

%% 2 -- Fit cosine per ROI
mu_fit = nan(nTotal, 1);
a_fit  = nan(nTotal, 1);
b_fit  = nan(nTotal, 1);
nSpk   = zeros(nTotal, 1);

for rr = 1:nTotal
    ph = all_ph_cell{rr};
    ph = ph(~isnan(ph));
    nSpk(rr) = numel(ph);
    if nSpk(rr) < minSpikes, continue; end

    mu_fit(rr) = angle(mean(exp(1i * ph)));
    cnt  = histcounts(ph, edges_hist);
    prob = cnt / sum(cnt);
    a_fit(rr) = mean(prob);
    b_fit(rr) = 2 * mean(prob .* cos(ctrs_hist - mu_fit(rr)));
end

fit_ok = ~isnan(mu_fit);
if sigOnly
    keep = fit_ok & all_sig_vec;
else
    keep = fit_ok;
end
fprintf('Plotting %d ROIs (sigOnly=%d, minSpikes=%d)\n', sum(keep), sigOnly, minSpikes);

%% 3 -- Overlay plot (one square axis)
cmap_exp = [
    0.00 0.30 0.00;   % 1: 251124 vglut2 breathing  -> darker green
    1.00 0.85 0.00;   % 2: 260322 sst maybe bre.    -> yellow
    1.00 0.85 0.00;   % 3: 260330 sst maybe_bre.    -> yellow
    0.85 0.15 0.15;   % 4: 260323 vgat processed    -> red
    0.35 0.70 0.30;   % 5: 260224 vglut2 breathing  -> medium-light green
];
assert(size(cmap_exp,1) >= nExp, 'cmap_exp needs >= nExp rows');
ax_w = 12; ax_h = 6;  % cm, width = 2 * height
fig = figure('Color','w', 'Units','centimeters', ...
    'Position',[3 3 ax_w+3 ax_h+4], 'PaperPositionMode','auto', ...
    'Renderer','painters');
ax = axes(fig, 'Units','centimeters', 'Position',[2 2.2 ax_w ax_h]);
hold(ax,'on');

x_fit = linspace(0, 4*pi, 800);
h_line = gobjects(nExp, 1);
for rr = find(keep(:))'
    y = a_fit(rr) + b_fit(rr) * cos(x_fit - mu_fit(rr));
    c = cmap_exp(all_exp_vec(rr), :);
    if all_sig_vec(rr)
        lw = 2;  % sig -> double width
    else
        lw = 0.5;
    end
    h = plot(ax, x_fit, y, '-', 'Color', c, 'LineWidth', lw);
    if ~isgraphics(h_line(all_exp_vec(rr)))
        h_line(all_exp_vec(rr)) = h;
    end
end
xline(ax, 2*pi, 'k--', 'LineWidth', 1);

xlim(ax, [0 4*pi]);
yl = ylim(ax); ylim(ax, [0, yl(2)]);

% ---- Overlay mean breathing waveform per folder (black) ----
phase_1cyc = linspace(0, 2*pi, nPhaseBins + 1); phase_1cyc = phase_1cyc(1:end-1);
phase_full = [phase_1cyc, phase_1cyc + 2*pi];
yl2 = ylim(ax);
for ee = 1:nExp
    cyc_mat = cycles_per_exp{ee};
    if isempty(cyc_mat), continue; end
    mean_cycle = mean(cyc_mat, 1, 'omitnan');
    mc_min = min(mean_cycle); mc_max = max(mean_cycle);
    if mc_max > mc_min
        mean_cycle = (mean_cycle - mc_min) / (mc_max - mc_min);
    end
    cycle_full = [mean_cycle, mean_cycle];
    cycle_scaled = yl2(1) + cycle_full * (yl2(2) - yl2(1));
    plot(ax, phase_full, cycle_scaled, 'k-', 'LineWidth', 1.2);
    fprintf('Breath overlay [%s]: %d cycles\n', expNames{ee}, size(cyc_mat, 1));
end
set(ax, 'XTick', [0 pi 2*pi 3*pi 4*pi], ...
        'XTickLabel', {'0','\pi','2\pi','3\pi','4\pi'});
xlabel(ax, 'Breath phase (rad)');
ylabel(ax, 'P(spike | phase)');
title(ax, sprintf('Per-ROI cosine fit (n = %d)', sum(keep)));
box(ax,'on');
set(ax, 'GridAlpha', 1, 'MinorGridAlpha', 1);  % no transparency -> vector PDF
grid(ax,'on');
pbaspect(ax, [2 1 1]);

valid_exp = find(isgraphics(h_line));
if ~isempty(valid_exp)
    legend(h_line(valid_exp), expNames(valid_exp), ...
        'Location','southoutside', 'Orientation','horizontal', ...
        'Interpreter','none', 'FontSize', 7, 'Box','off');
end
hold(ax,'off');

%% 4 -- Save
if doSave
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    tag  = ternary(sigOnly, 'sig', 'all');
    pdf_path = char(fullfile(outDir, sprintf('master_cosine_fit_phase_%s.pdf', tag)));
    mat_path = char(fullfile(outDir, sprintf('master_cosine_fit_phase_%s.mat', tag)));

    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    save(mat_path, 'mu_fit','a_fit','b_fit','nSpk', ...
        'all_sig_vec','all_exp_vec','expNames','expFolders', ...
        'sigOnly','minSpikes','nBinsHist');
    fprintf('Saved figure + mat to %s\n', outDir);
end

%% ===== helpers =====
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
