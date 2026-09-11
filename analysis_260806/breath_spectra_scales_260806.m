% breath_spectra_scales_260806.m
% -----------------------------------------------------------------------
%  One figure for the calcium power spectra alone: the SAME curves four times,
%  on the four axis scalings, so the shape can be judged without committing to
%  one view.
%
%      lin-lin | log-log | log x, lin y | lin x, log y
%
%  ACTIVE = green, INACTIVE = red, alpha 0.6.  Nothing is normalised: the traces
%  were z-scored before the multitaper, so every curve is in z^2/Hz already.
%
%  READ THE LOW END WITH CARE.  The multitaper half-bandwidth is W_hz (0.25 Hz by
%  default), so the whole 0.05-0.25 Hz range sits inside a single resolution
%  element.  It is plotted because you asked for the band, but structure there is
%  smearing, not signal.
%
%  Input : breath_trig_heatmap_GENOTYPES_peak.mat
%  Output: breath_spectra_scales_<TRIGGER>.png / .pdf
%
%  Runqi Zhang / 2026-08-06

clear; clc;

%% ===================== USER-EDITABLE =====================
dataFile  = 'D:\Ventral_surface_summary\breath_trig_heatmap_260806\breath_trig_heatmap_GENOTYPES_peak.mat';
specAlpha = 0.6;
colAct    = [0    0.65 0.2 ];   % active   = green
colOff    = [0.85 0.1  0.1 ];   % inactive = red
lw        = 0.5;
doSave    = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
S = load(dataFile, 'Sc_all','fg','spkOf','rayZ','zLock','fpass','W_hz','ibi_all','TRIGGER');

isAct = S.spkOf == 1 | S.rayZ >= S.zLock;      % locked ROIs are active by definition
f     = S.fg(:);
fprintf('%d ROIs: %d active (green), %d inactive (red) | %.2f-%.1f Hz, W = %.2f Hz\n', ...
        numel(isAct), nnz(isAct), nnz(~isAct), S.fpass(1), S.fpass(2), S.W_hz);

SCALES = { 'lin-lin', 'linear','linear'
           'log-log', 'log',   'log'
           'log x, lin y', 'log','linear'
           'lin x, log y', 'linear','log' };

fBreath = 1/median(S.ibi_all);

fh = figure('Color','w','Position',[30 30 1500 950]);
tl = tiledlayout(fh, 2, 2, 'TileSpacing','compact','Padding','compact');

for k = 1:size(SCALES,1)
    ax = nexttile(tl, k); hold(ax,'on');
    for i = find(~isAct).'
        plot(ax, f, S.Sc_all(:,i), '-', 'Color',[colOff specAlpha], 'LineWidth',lw);
    end
    for i = find(isAct).'
        plot(ax, f, S.Sc_all(:,i), '-', 'Color',[colAct specAlpha], 'LineWidth',lw);
    end
    xline(ax, fBreath, 'k--', 'LineWidth', 1);
    xline(ax, S.W_hz,  'k:',  'LineWidth', 1);      % below this nothing is resolved
    set(ax, 'XScale', SCALES{k,2}, 'YScale', SCALES{k,3}, 'TickDir','out', 'Box','off');
    xlim(ax, S.fpass);
    if strcmp(SCALES{k,2}, 'log')
        xt = [0.05 0.1 0.2 0.5 1 2 3 5 7 10 14]; xt = xt(xt>=S.fpass(1) & xt<=S.fpass(2));
        set(ax, 'XTick', xt, 'XTickLabel', compose('%g', xt));
    end
    xlabel(ax, 'frequency (Hz)'); ylabel(ax, 'power (SD^2/Hz)');
    title(ax, SCALES{k,1});
end

title(tl, sprintf(['z-scored dF/F power, %d ROIs: %d active (green) vs %d inactive (red)   |   ' ...
                   'dashed = breath %.2f Hz, dotted = multitaper resolution %.2f Hz'], ...
                   numel(isAct), nnz(isAct), nnz(~isAct), fBreath, S.W_hz), ...
      'FontWeight','bold','Interpreter','none');

if doSave
    base = fullfile(fileparts(dataFile), sprintf('breath_spectra_scales_%s', lower(S.TRIGGER)));
    exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    fprintf('saved %s.{png,pdf}\n', base);
end
