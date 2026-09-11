% Ventral_surface_cycle_dFF_insp_onset_260529.m
% -----------------------------------------------------------------------
%  Quick single-cycle dF/F heatmap aligned to INSPIRATION ONSET (foot).
%  Two panels: SIG (left) vs NON-SIG (right) from the phase-coherence
%  classification (Ventral_surface_coherence_polar_phase_260529.m).
%
%  Each row = one breath cycle from one ROI, triggered on its foot
%  (insp_start_idx). Rows sorted by dt to LAST inspiration PEAK so the
%  expiration duration sweeps from short (top) to long (bottom).
%
%  Three line overlays per row:
%    - last inspiration peak (gray)   = previous peak_idx, at t < 0
%    - this insp onset      (red)    = the trigger, at t = 0
%    - next insp onset      (yellow) = next foot, at t > 0
%
%  No spike tick overlay. Time domain only. Vglut2/1124/IO excluded.
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, '2p_breathing_coherence'));

%% ===================== USER-EDITABLE =================================
rootPath = 'D:\Ventral_surface_summary';
cohMat   = fullfile(rootPath,'coherence_polar_phase_260529','coherence_polar_phase_data.mat');

exclude_substrings = { fullfile('Vglut2','1124','IO') };

t_max   = 2.5;        % seconds forward from insp onset; covers slowest cycle
t_dt    = 0.02;       % seconds per column
prcLim  = [1 99.5];   % color clamp percentiles on raw dF/F

nDrop        = 30;
fallback_fps = 30;

doSave = false;
% =====================================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

%% ---- load coherence labels + sig flags ----
assert(isfile(cohMat), 'Missing %s', cohMat);
S = load(cohMat, 'PP','labels','confC');
PP_all     = S.PP;
labels_all = S.labels;
confC      = S.confC;
sig_all    = PP_all.r >= confC;

% map label -> (sig flag, ROI idx in PP)
lookup = containers.Map(labels_all, num2cell(1:numel(labels_all)));

%% ---- common time grid: one-sided, 0 to t_max ----
t_grid = 0 : t_dt : t_max;
nW     = numel(t_grid);

%% ---- walk folders, collect per-cycle snippets ----
% Non-overlapping: per FOV per cycle ONE row in sig (mean dF/F across sig
% ROIs of that FOV) and ONE row in non-sig (mean across non-sig ROIs).
% A FOV with zero sig (or zero non-sig) ROIs simply contributes only to
% the other panel.
scan_dirs = {'ChAT','Vglut2','Vgat','Sst'};

sig_snips = [];  sig_dtL = [];  sig_dtN = [];  sig_dtF = [];
nsg_snips = [];  nsg_dtL = [];  nsg_dtN = [];  nsg_dtF = [];

for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(rootPath, sname);
    if ~isfolder(gdir), continue; end
    allMat = dir(fullfile(gdir,'**','ca_spike_data.mat'));
    fprintf('\n[%s] %d recordings\n', sname, numel(allMat));

    for kk = 1:numel(allMat)
        folderPath = allMat(kk).folder;

        if any(cellfun(@(p) ~isempty(strfind(lower(folderPath), lower(p))), exclude_substrings)) %#ok<STREMP>
            fprintf('  EXCLUDED: %s\n', folderPath); continue;
        end

        gname = sname;
        if strcmpi(sname,'ChAT') && is_io_path(folderPath, gdir), gname = 'IO'; end

        try
            df_hit = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
            if isempty(df_hit), continue; end
            peak_file = find_data_file(folderPath,'breath_peak_data.mat');
            ip_file   = find_data_file(folderPath,'breath_insp_start_data.mat');
            if isempty(peak_file) || isempty(ip_file)
                fprintf('  skip (missing): %s\n', folder_basename(folderPath)); continue;
            end

            fps = detect_session_fps(folderPath, fallback_fps);
            BP  = load(peak_file);
            IP  = load(ip_file);
            D   = load(fullfile(df_hit(1).folder, df_hit(1).name),'dFF');
            CA  = load(fullfile(folderPath,'ca_spike_data.mat'));   % nROI

            dff  = double(D.dFF);
            nROI = size(dff,2);

            % which ROIs in this FOV are sig / non-sig?
            recName = folder_basename(folderPath);
            sig_cols = false(1, nROI);
            nsg_cols = false(1, nROI);
            for rid = 1:nROI
                lab = sprintf('%s/%s#%d', gname, recName, rid);
                if ~isKey(lookup, lab), continue; end
                if sig_all(lookup(lab)), sig_cols(rid) = true;
                else,                    nsg_cols(rid) = true;
                end
            end

            % build mean traces ONCE per FOV
            tr_sig = nan(size(dff,1),1);
            tr_nsg = nan(size(dff,1),1);
            if any(sig_cols), tr_sig = mean(dff(:, sig_cols), 2, 'omitnan'); end
            if any(nsg_cols), tr_nsg = mean(dff(:, nsg_cols), 2, 'omitnan'); end

            peak_idx = sort(double(BP.insp_onset_idx(:)) - nDrop);
            peak_idx(peak_idx<1) = [];
            foot_idx = sort(double(IP.insp_start_idx(:)) - nDrop);
            foot_idx(foot_idx<1) = [];

            T = min([size(dff,1), max(peak_idx(end), foot_idx(end))+1]);
            peak_idx(peak_idx>T) = [];
            foot_idx(foot_idx>T) = [];
            if numel(foot_idx) < 2 || numel(peak_idx) < 1, continue; end

            win_fr = round(t_max*fps);

            % iterate cycles; trigger = THIS insp onset (foot) at t=0.
            % Window forward to next onset (NaN beyond).
            for fi_ = 1:numel(foot_idx)
                f  = foot_idx(fi_);
                if f + win_fr > size(dff,1), continue; end

                % previous breath peak (sort metric)
                pk_before = peak_idx(peak_idx < f);
                if isempty(pk_before), continue; end
                dt_prev_peak = (pk_before(end) - f) / fps;     % negative

                % the breath peak inside this cycle (between this onset and next onset)
                pk_after = peak_idx(peak_idx > f);
                if isempty(pk_after), continue; end
                dt_this_peak = (pk_after(1) - f) / fps;        % positive

                % next insp onset
                ft_after = foot_idx(foot_idx > f);
                if isempty(ft_after), continue; end
                dt_next_foot = (ft_after(1) - f) / fps;        % positive (= cycle dur)

                % --- extract dF/F over one cycle, NaN beyond next onset ---
                next_foot_fr = ft_after(1) - f;
                next_foot_fr = min(next_foot_fr, win_fr);      % clamp to t_max
                if any(sig_cols)
                    snip = tr_sig(f : f + next_foot_fr);
                    rs   = interp1((0:next_foot_fr)/fps, snip, t_grid, 'linear', NaN);
                    sig_snips(end+1,:) = rs;                   %#ok<AGROW>
                    sig_dtL(end+1,1)   = dt_prev_peak;         %#ok<AGROW>
                    sig_dtN(end+1,1)   = dt_this_peak;         %#ok<AGROW>
                    sig_dtF(end+1,1)   = dt_next_foot;         %#ok<AGROW>
                end
                if any(nsg_cols)
                    snip = tr_nsg(f : f + next_foot_fr);
                    rs   = interp1((0:next_foot_fr)/fps, snip, t_grid, 'linear', NaN);
                    nsg_snips(end+1,:) = rs;                   %#ok<AGROW>
                    nsg_dtL(end+1,1)   = dt_prev_peak;         %#ok<AGROW>
                    nsg_dtN(end+1,1)   = dt_this_peak;         %#ok<AGROW>
                    nsg_dtF(end+1,1)   = dt_next_foot;         %#ok<AGROW>
                end
            end
            fprintf('  ok %s  (sig ROIs=%d, nonsig ROIs=%d)\n', recName, sum(sig_cols), sum(nsg_cols));
        catch ME
            warning('  ERROR %s: %s', folderPath, ME.message);
        end
    end
end

assert(~isempty(sig_snips) || ~isempty(nsg_snips), 'No cycles collected.');
fprintf('\nNon-overlapping cycles: sig=%d, non-sig=%d\n', size(sig_snips,1), size(nsg_snips,1));

%% ---- sort by dt_prev_peak (ascending: most negative on top = longest cycle) ----
[snips_sig, dtL_sig, dtN_sig, dtF_sig] = sort_by_dtprev(sig_snips, sig_dtL, sig_dtN, sig_dtF);
[snips_nsg, dtL_nsg, dtN_nsg, dtF_nsg] = sort_by_dtprev(nsg_snips, nsg_dtL, nsg_dtN, nsg_dtF);

pool  = [snips_sig(:); snips_nsg(:)];
gclim = prctile(pool(isfinite(pool)), prcLim);

%% ---- figure ----
fig = figure('Color','w','Name','single-cycle dF/F | trig=insp onset', ...
             'Units','centimeters','Position',[2 2 26 14]);

axS = subplot(1,2,1);
draw_panel(axS, snips_sig, t_grid, dtN_sig, dtF_sig, gclim, ...
           sprintf('SIG  %d cycles', size(snips_sig,1)));

axN = subplot(1,2,2);
draw_panel(axN, snips_nsg, t_grid, dtN_nsg, dtF_nsg, gclim, ...
           sprintf('NON-SIG  %d cycles', size(snips_nsg,1)));

sgtitle('single-cycle dF/F, trigger=insp onset, sorted by dt to previous peak');

if doSave
    outDir = fullfile(rootPath,'coherence_polar_phase_260529');
    if ~isfolder(outDir), mkdir(outDir); end
    exportgraphics(fig, fullfile(outDir,'cycle_dFF_insp_onset.png'), 'Resolution',200,'BackgroundColor','white');
    exportgraphics(fig, fullfile(outDir,'cycle_dFF_insp_onset.pdf'), 'ContentType','vector','BackgroundColor','white');
end

%% ========================= LOCAL FUNCTIONS ==========================
function [M, dtL, dtN, dtF] = sort_by_dtprev(snips, dt_prev_peak, dt_next_peak, dt_next_foot)
    if isempty(snips)
        M = snips; dtL = dt_prev_peak; dtN = dt_next_peak; dtF = dt_next_foot; return;
    end
    [dtL, ord] = sort(dt_prev_peak, 'ascend');
    M   = snips(ord, :);
    dtN = dt_next_peak(ord);
    dtF = dt_next_foot(ord);
end

function draw_panel(ax, M, t_grid, dtPeak, dtFoot, clim, ttl)
% One-sided display: trigger = this insp onset at t = 0, forward to the
% next insp onset (variable per cycle). Each row's signal is NaN past
% next_foot. Three overlays: red vertical at t=0 (this onset), yellow
% curve at this peak, blue curve at next onset.
    h = imagesc(ax, t_grid, 1:size(M,1), M, 'AlphaData', ~isnan(M)); hold(ax,'on');
    set(ax,'Color',[0.96 0.96 0.96]);   % NaN regions show as light gray
    colormap(ax, flipud(gray(256))); caxis(ax, clim);
    set(ax,'YDir','reverse'); axis(ax,'tight');
    n = size(M,1); y = (1:n)';

    xline(ax, 0, 'r-', 'LineWidth', 1.4);                              % this onset
    plot(ax, dtPeak, y, '-', 'Color',[0.95 0.80 0.10], 'LineWidth',1.0);  % this peak
    plot(ax, dtFoot, y, '-', 'Color',[0   0.55 1   ], 'LineWidth',1.0);  % next onset

    xlim(ax, [t_grid(1) t_grid(end)]);
    xlabel(ax,'time from insp onset (s)');
    ylabel(ax,'cycle (sorted by dt to previous peak)');
    title(ax, ttl);
    cb = colorbar(ax); cb.Label.String = '\DeltaF/F';
end

function fp = find_data_file(folderPath, basename)
    hits = dir(fullfile(folderPath, ['*' basename]));
    if ~isempty(hits)
        fp = fullfile(hits(1).folder, hits(1).name);
    else
        bare = fullfile(folderPath, basename);
        if isfile(bare), fp = bare; else, fp = ''; end
    end
end

function tf = is_io_path(folderPath, groupRoot)
    rel = strrep(folderPath, groupRoot, '');
    rel = regexprep(rel, '^[\\/]+', '');
    parts = regexp(rel, '[\\/]', 'split');
    tf = numel(parts) >= 2 && ~isempty(regexpi(parts{2}, 'IO', 'once'));
end

function name = folder_basename(p)
    p = char(p);
    while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end)=[]; end
    [~,n,e] = fileparts(p);
    name = [n e];
end
