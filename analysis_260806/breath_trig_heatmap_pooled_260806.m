% breath_trig_heatmap_pooled_260806.m
% -----------------------------------------------------------------------
%  The same figure as breath_trig_heatmap_260806.m, but pooled over EVERY
%  recording in one experiment folder.  Same panels, same colours, same rules:
%  no sorting, no selection, every ROI of every recording is in there.
%
%      row 1   pooled breath trace                    | breath power spectra
%      row 2   per-recording mean dF/F (5% bar)       | dF/F heatmap
%              | z-scored dF/F heatmap                | z power spectra
%
%  HOW POOLING IS DONE, and why it is not a plain average
%  -----------------------------------------------------
%  Each recording breathes at its own rate (1.74-1.98 Hz here, and one recording
%  was imaged at 42 fps rather than 30).  Averaging PETHs on a shared axis in
%  SECONDS would line up the trigger at t=0 and then progressively smear the
%  neighbouring breaths, so the +/-1 cycle peaks would wash out for no reason
%  other than rate mismatch.
%  So every recording's PETH is built on ITS OWN time axis and then resampled
%  onto a common axis measured in BREATH CYCLES.  The x axis is finally labelled
%  in seconds using the pooled median IBI, i.e. "seconds, for an average breath".
%  Only the trigger and the cycle structure are forced to align; nothing else is
%  rescaled and no amplitude is touched.
%
%  ROI ordering is recording by recording, in folder order, native ROI order
%  within each -- still no sorting of any kind.  White lines mark the boundaries.
%
%  Input : every subfolder of rootPath holding breath_peak_pc1.mat (or
%          breath_insp_start_pc1.mat) plus *_ch1_dFF.mat
%  Output: rootPath\breath_trig_heatmap_POOLED_<TRIGGER>.png / .pdf / .mat
%
%  Runqi Zhang / 2026-08-06

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));

%% ===================== USER-EDITABLE =====================
rootPath = 'C:\260806_sst-soma-g8s\phys';

TRIGGER    = 'peak';       % 'peak' | 'onset'   (same meaning as the single-FOV script)
nCyc_panel = 10/3;         % breath cycles per panel in row 2
nCyc_top   = 3*nCyc_panel; % cycles in the top breath trace
nDrop      = 30;           % calcium TossFrames
prcLim     = [0.5 99.5];   % percentile clamp for both heatmap colour scales
nGrid      = 241;          % samples across a panel window on the common cycle axis

fpass      = [0.1 14];     % Hz
W_hz       = 0.25;         % multitaper half-bandwidth (smoothing), Hz
% NOTHING is normalised.  Breath PC1 is z-scored upstream and the calcium traces
% are z-scored here, so both spectra are already in z^2/Hz.

stackGain    = 1.0;        % spacing of the per-recording dF/F traces
scaleBarDFF  = 0.05;       % 5% dF/F calibration bar
fallback_fps = 30;
doSave       = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

%% ---- find the recordings ----
dd = dir(rootPath); dd = dd([dd.isdir]);
dd = dd(~ismember({dd.name}, {'.','..'}) & ~startsWith({dd.name}, '.'));
recs = {};
for i = 1:numel(dd)
    fp = fullfile(rootPath, dd(i).name);
    if ~isempty(dir(fullfile(fp,'*_ch1_dFF.mat'))) && isfile(fullfile(fp,'breath_peak_pc1.mat'))
        recs{end+1} = dd(i).name; %#ok<SAGROW>
    end
end
nRec = numel(recs);
assert(nRec > 0, 'No usable recording under %s', rootPath);
fprintf('%s\n  %d recordings\n', rootPath, nRec);

%% ---- common axes ----
u   = linspace(-nCyc_panel/2, nCyc_panel/2, nGrid);          % panels, in CYCLES
ub  = linspace(-nCyc_top/2,   nCyc_top/2,   3*nGrid);        % breath trace, in cycles
fg  = logspace(log10(fpass(1)), log10(fpass(2)), 600).';     % common frequency grid

M_all = []; Mz_all = []; recOf = [];        % pooled ROI matrices + their recording
B_rec = nan(nRec, numel(ub));               % per-recording breath PETH
D_rec = nan(nRec, nGrid);                   % per-recording MEAN dF/F trace
Sb_rec = nan(numel(fg), nRec);              % per-recording breath spectrum
Sc_all = [];                                % per-ROI z spectra (unnormalised)
ibi_rec = nan(nRec,1); nTrig_rec = nan(nRec,1); nROI_rec = nan(nRec,1);

for k = 1:nRec
    fp = fullfile(rootPath, recs{k});
    df = dir(fullfile(fp,'*_ch1_dFF.mat'));
    D  = load(fullfile(df(1).folder, df(1).name), 'dFF');
    dFF = double(D.dFF);

    switch lower(TRIGGER)
        case 'peak'
            B = load(fullfile(fp,'breath_peak_pc1.mat'));  trig = round(B.insp_onset_idx(:));
        case 'onset'
            B = load(fullfile(fp,'breath_insp_start_pc1.mat')); trig = round(B.insp_start_idx(:));
    end

    fps    = detect_session_fps(fp, fallback_fps);
    breath = B.breath(:); breath = breath(nDrop+1:end);
    T      = min(size(dFF,1), numel(breath));
    dFF    = dFF(1:T,:); breath = breath(1:T);
    trig   = trig - nDrop;

    nROI = size(dFF,2);
    Z    = (dFF - mean(dFF,1)) ./ std(dFF,0,1);

    ibi   = median(diff(sort(trig)))/fps;
    win   = round(nCyc_panel*ibi*fps/2);
    win_b = round(nCyc_top  *ibi*fps/2);
    tl_u  = (-win:win)   /fps/ibi;                 % local axis, in cycles
    tl_ub = (-win_b:win_b)/fps/ibi;

    tg   = trig(trig-win   >= 1 & trig+win   <= T);
    tg_b = trig(trig-win_b >= 1 & trig+win_b <= T);
    if isempty(tg) || isempty(tg_b), warning('%s: no usable trigger, skipped.', recs{k}); continue; end
    idx   = tg(:).'   + (-win:win).';
    idx_b = tg_b(:).' + (-win_b:win_b).';

    % breath PETH -> common cycle axis
    B_rec(k,:) = interp1(tl_ub, mean(breath(idx_b),2), ub, 'linear', NaN);

    % per-ROI PETHs -> common cycle axis
    Mk = nan(nROI, nGrid); Mzk = nan(nROI, nGrid);
    for r = 1:nROI
        tr = dFF(:,r); zr = Z(:,r);
        Mk(r,:)  = interp1(tl_u, mean(tr(idx),2), u, 'linear', NaN);
        Mzk(r,:) = interp1(tl_u, mean(zr(idx),2), u, 'linear', NaN);
    end
    M_all  = [M_all;  Mk];   %#ok<AGROW>
    Mz_all = [Mz_all; Mzk];  %#ok<AGROW>
    recOf  = [recOf; repmat(k, nROI, 1)];  %#ok<AGROW>
    D_rec(k,:) = mean(Mk, 1);

    % spectra, each recording on its own taper set, then onto the common f grid
    TW = max(2, round(W_hz * T/fps)); Kt = 2*TW - 1;
    params = struct('Fs',fps, 'tapers',[TW Kt], 'fpass',fpass, 'pad',0, 'trialave',0, 'err',0);
    [Sb, fb] = mtspectrumc(detrend(breath), params);
    [Sc, fc] = mtspectrumc(detrend(Z),      params);
    Sb_rec(:,k) = interp1(fb, Sb, fg, 'linear', NaN);
    Sc_all = [Sc_all, interp1(fc, Sc, fg, 'linear', NaN)];  %#ok<AGROW>

    ibi_rec(k) = ibi; nTrig_rec(k) = numel(tg); nROI_rec(k) = nROI;
    fprintf('  %-38s %3d ROI | IBI %.3f s | %3d triggers | fps %.1f\n', ...
            recs{k}, nROI, ibi, numel(tg), fps);
end

ibi_pool = median(ibi_rec, 'omitnan');
tt = u  * ibi_pool;                    % common axis, in "seconds of an average breath"
tb = ub * ibi_pool;
nROIall = size(M_all,1);
fprintf('  POOLED: %d ROIs, %d recordings, %d triggers, median IBI %.3f s (%.2f Hz)\n', ...
        nROIall, nRec, sum(nTrig_rec,'omitnan'), ibi_pool, 1/ibi_pool);

%% ---- figure ----
recCol = hsv(nRec);                    % one colour per RECORDING (not per ROI)
bnd    = find(diff(recOf) ~= 0) + 0.5; % ROI-row boundaries between recordings

fh = figure('Color','w','Position',[30 30 1780 990]);
tl = tiledlayout(fh, 5, 4, 'TileSpacing','compact','Padding','compact');

% ---- row 1: 1 grid row against row 2's 4, so row 2 gets the space ----
axb = nexttile(tl, 1, [1 3]); hold(axb,'on');
for k = 1:nRec, plot(axb, tb, B_rec(k,:), '-', 'Color',[recCol(k,:) 0.35], 'LineWidth',0.8); end
plot(axb, tb, mean(B_rec,1,'omitnan'), 'k-', 'LineWidth', 1.8);
xline(axb, 0, 'r-', 'LineWidth', 1);
for e = [-1 1]*tt(end), xline(axb, e, ':', 'Color',[0.5 0.5 0.5]); end
xlim(axb, [tb(1) tb(end)]); box(axb,'off');
set(axb, 'TickDir','out', 'XTick', unique(round((-4:4)*ibi_pool, 2)));
xlabel(axb, sprintf('time from inspiratory %s (s, at the pooled median breath)', TRIGGER));
ylabel(axb, 'breath (PC1, z)');
title(axb, sprintf(['pooled %s-triggered breath, %.0f cycles  ' ...
                    '(black = mean of %d recordings; dotted = span of the panels below)'], ...
                    TRIGGER, nCyc_top, nRec));

axs1 = nexttile(tl, 4, [1 1]); hold(axs1,'on');
for k = 1:nRec, loglog(axs1, fg, Sb_rec(:,k), '-', 'Color',[recCol(k,:) 0.8], 'LineWidth',0.9); end
finish_spec(axs1, fpass, 1/ibi_pool, 'power (z^2/Hz)', 'breath (PC1) power, per recording');

% ---- row 2 ----
draw_stack(nexttile(tl, 5, [4 1]), tt, M_all, recCol(recOf,:), stackGain, scaleBarDFF, ...
           TRIGGER, bnd, recs, recOf);
draw_heatmap(nexttile(tl, 6, [4 1]), tt, M_all,  prcLim, bnd, ...
             sprintf('mean dF/F, %d ROIs pooled', nROIall), 'dF/F', TRIGGER, gray(256));
draw_heatmap(nexttile(tl, 7, [4 1]), tt, Mz_all, prcLim, bnd, ...
             'mean z-scored dF/F', 'z (SD)', TRIGGER, parula(256));

axs2 = nexttile(tl, 8, [4 1]); hold(axs2,'on');
for i = 1:nROIall
    loglog(axs2, fg, Sc_all(:,i), '-', 'Color',[recCol(recOf(i),:) 0.5], 'LineWidth',0.6);
end
finish_spec(axs2, fpass, 1/ibi_pool, 'power (SD^2/Hz)', ...
            sprintf('z-scored dF/F power, all %d ROIs (not normalised)', nROIall));

[~, expName] = fileparts(rootPath);
if strcmpi(expName,'phys'), [~, expName] = fileparts(fileparts(rootPath)); end
title(tl, sprintf('%s   |   POOLED: %d recordings, %d ROIs, no sorting', expName, nRec, nROIall), ...
      'Interpreter','none','FontWeight','bold');

%% ---- save ----
if doSave
    base = fullfile(rootPath, sprintf('breath_trig_heatmap_POOLED_%s', lower(TRIGGER)));
    exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'M_all','Mz_all','D_rec','B_rec','recOf','recs','tt','tb','u','ub', ...
                        'Sb_rec','Sc_all','fg','ibi_rec','ibi_pool','nTrig_rec','nROI_rec', ...
                        'TRIGGER','nCyc_panel','nCyc_top','nDrop','prcLim','fpass','W_hz', ...
                        'rootPath');
    fprintf('  saved %s.{png,pdf,mat}\n', base);
end

%% ---- local ----
function draw_heatmap(ax, tt, A, prcLim, bnd, ttl, cbl, TRIGGER, cmap)
c = prctile(A(:), prcLim);
if ~(c(2) > c(1)), c = [min(A(:)) max(A(:))+eps]; end
imagesc(ax, tt, 1:size(A,1), A);
set(ax, 'CLim', c, 'YDir','normal', 'TickDir','out', 'Layer','top');
colormap(ax, cmap); hold(ax,'on');
xline(ax, 0, 'r-', 'LineWidth', 1);
for b = bnd(:).'                       % recording boundaries
    plot(ax, [tt(1) tt(end)], [b b], 'w-', 'LineWidth', 0.8);
end
xlabel(ax, sprintf('time from inspiratory %s (s)', TRIGGER));
ylabel(ax, 'ROI (recordings stacked, native order within each)');
title(ax, ttl);
cb = colorbar(ax); cb.Label.String = cbl;
end

function draw_stack(ax, tt, A, colPerRow, gain, barDFF, TRIGGER, bnd, recs, recOf)
% EVERY pooled ROI, one trace each, coloured by its recording.  All on ONE common
% dF/F scale -- nothing auto-scaled, so the bar applies to every trace.
nR   = size(A,1);
step = gain * median(range(A, 2), 'omitnan');
if ~(step > 0), step = 1; end
hold(ax,'on');
for k = 1:nR
    plot(ax, tt, A(k,:) - mean(A(k,:),'omitnan') + (k-1)*step, '-', ...
         'Color', colPerRow(k,:), 'LineWidth', 0.5);
end
xline(ax, 0, 'k-', 'LineWidth', 1);
xlim(ax, [tt(1) tt(end)]); ylim(ax, [-step, (nR-1)*step + step]);
% one y tick per recording, at its first ROI
lbl  = regexprep(recs, '_00001$', '');
first = arrayfun(@(k) find(recOf==k, 1, 'first'), 1:numel(recs));
set(ax, 'YTick', (first-1)*step, 'YTickLabel', lbl, 'TickDir','out', 'Box','off', ...
        'FontSize',8, 'TickLabelInterpreter','none');   % names have _ in them
for b = bnd(:).', yline(ax, (b-1)*step, '-', 'Color',[0.8 0.8 0.8]); end
xb = tt(end) - 0.02*(tt(end)-tt(1)); yb = -0.6*step;
plot(ax, [xb xb], yb + [0 barDFF], 'k-', 'LineWidth', 2.5);
text(ax, xb - 0.02*(tt(end)-tt(1)), yb + barDFF/2, sprintf('%g%% dF/F', 100*barDFF), ...
     'HorizontalAlignment','right', 'VerticalAlignment','middle', 'FontSize', 9);
xlabel(ax, sprintf('time from inspiratory %s (s)', TRIGGER));
title(ax, sprintf('mean dF/F, all %d ROIs (common scale)', nR));
end

function finish_spec(ax, fpass, f_breath, ylb, ttl)
hold(ax,'on');
xline(ax, f_breath, 'r--', 'LineWidth', 1);
set(ax, 'XScale','log', 'YScale','log', 'TickDir','out', 'Box','off');
xlim(ax, fpass);
xt = [0.1 0.2 0.5 1 2 3 5 7 10 14]; xt = xt(xt>=fpass(1) & xt<=fpass(2));
set(ax, 'XTick', xt, 'XTickLabel', compose('%g', xt));
xlabel(ax, 'frequency (Hz)'); ylabel(ax, ylb); title(ax, ttl);
end
