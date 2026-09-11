function stitch_shift_null_plot_both_260910(varargin)
%STITCH_SHIFT_NULL_PLOT_BOTH_260910  Both channels' shift nulls in one figure.
%
%   stitch_shift_null_plot_both_260910
%   stitch_shift_null_plot_both_260910('Root', ..., 'Channels', [1 3])
%
% Re-plots the saved results of stitch_shift_null_260910 for two channels on ONE
% set of axes per panel, so the somata label and the vessel label are read
% against each other rather than across two figures.  Nothing is recomputed --
% this reads the .mat each run wrote, so the figure can be restyled without
% re-scoring a thousand shifts.
%
% Layout: one row per seam, three columns
%   1  YF null distribution + the curated value
%   2  YF against shift magnitude
%   3  Pearson NCC null + the curated value  (the cross-check on YF)
%
% COLOUR CARRIES CHANNEL IDENTITY AND NOTHING ELSE.  Two categorical hues, fixed
% order, validated for colour-vision deficiency (worst adjacent pair dE 24.7
% protan, 33.6 normal -- the target is >= 8).  The curated value is NOT a third
% series: it is drawn as a heavy solid rule in its own channel's hue, so it
% inherits the identity already established, and it is labelled with its number
% directly.  The earlier per-channel figures used red for it, which is worse on
% two counts -- red is reserved for status, and with two channels overlaid a
% single red would have had to mean two different things.
%
% Histograms within a panel SHARE their bin edges across channels.  Letting each
% channel pick its own bins would make the two distributions incomparable in the
% one place the figure exists to compare them.
%
% The csv beside each .mat is the table view of this figure.
%
% Runqi Zhang / 2026-09-10

p = inputParser;
p.addParameter('Root', 'C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer\matlab_stitch');
p.addParameter('Channels', [1 3]);
p.addParameter('OutName', 'shift_null_both');
p.parse(varargin{:});
opt = p.Results;

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(here);
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

% Categorical slots 1 and 2 of the validated palette, in fixed order.
COL  = [ 42 120 214;  235 104  52] / 255;
INK  = [0.05 0.05 0.05];
INK2 = [0.32 0.32 0.31];
GRD  = [0.88 0.88 0.87];

nc = numel(opt.Channels);
D  = cell(nc,1);
for k = 1:nc
    ch = opt.Channels(k);
    fn = fullfile(opt.Root, sprintf('shift_null_ch%d', ch), sprintf('shift_null_ch%d.mat', ch));
    assert(isfile(fn), 'missing %s -- run stitch_shift_null_260910 for ch %d first', fn, ch);
    D{k} = load(fn, 'R', 'St', 'opt');
    fprintf('loaded ch%d: %d seams, %d draws each\n', ch, numel(D{k}.R), numel(D{k}.R(1).yf));
end
nS = numel(D{1}.R);
for k = 2:nc
    assert(numel(D{k}.R) == nS, 'channel %d has %d seams, channel %d has %d', ...
           opt.Channels(k), numel(D{k}.R), opt.Channels(1), nS);
end

f = figure('Color','w','Units','pixels','Position',[40 40 1400 400*nS]);
drawnow; got = get(f,'Position'); FS = min(1, got(3)/1400);
set(f,'DefaultAxesFontSize', max(7, 10*FS));
tl = tiledlayout(f, nS, 3, 'TileSpacing','compact', 'Padding','compact');

for s = 1:nS
    lbl = sprintf('seam %d-%d', D{1}.R(s).a, D{1}.R(s).b);

    % ---- YF null ---------------------------------------------------------
    ax = nexttile(tl, (s-1)*3 + 1);  hold(ax,'on');
    allv = []; for k=1:nc, allv = [allv; D{k}.R(s).yf(:); D{k}.R(s).yf0]; end %#ok<AGROW>
    ed = linspace(min(allv)*0.98, max(allv)*1.02, 34);
    for k = 1:nc
        histogram(ax, D{k}.R(s).yf, ed, 'FaceColor', COL(k,:), 'FaceAlpha', 0.55, ...
                  'EdgeColor','none', 'DisplayName', sprintf('ch %d null', opt.Channels(k)));
    end
    yl = ylim(ax);
    for k = 1:nc
        v = D{k}.R(s).yf0;
        plot(ax, [v v], [0 yl(2)], '-', 'Color', COL(k,:), 'LineWidth', 2, ...
             'HandleVisibility','off');
        label_rule(ax, v, yl(2)*(1 - 0.10*(k-1)), INK2, max(6,8*FS));
    end
    style(ax, GRD, INK2);
    xlabel(ax, 'YF  =  \Sigma min(a,b) / \Sigma max(a,b)', 'Color', INK);
    ylabel(ax, 'draws', 'Color', INK);
    title(ax, sprintf('%s    YF null', lbl), 'FontWeight','normal', 'Color', INK);
    if s == 1
        lg = legend(ax, 'Location','northwest'); lg.Box = 'off'; lg.TextColor = INK2;
    end
    hold(ax,'off');

    % ---- YF vs displacement ---------------------------------------------
    ax = nexttile(tl, (s-1)*3 + 2);  hold(ax,'on');
    for k = 1:nc
        r = D{k}.R(s);
        scatter(ax, hypot(r.dx, r.dy), r.yf, 9, COL(k,:), 'filled', ...
                'MarkerFaceAlpha', 0.40, 'HandleVisibility','off');
    end
    for k = 1:nc     % curated last, with a surface ring so it sits above the cloud
        plot(ax, 0, D{k}.R(s).yf0, 'o', 'MarkerSize', 9, 'MarkerFaceColor', COL(k,:), ...
             'MarkerEdgeColor','w', 'LineWidth', 1.5, 'HandleVisibility','off');
    end
    style(ax, GRD, INK2);
    xlabel(ax, 'shift magnitude  |(dx,dy)|   (\mum)', 'Color', INK);
    ylabel(ax, 'YF', 'Color', INK);
    title(ax, sprintf('%s    YF vs displacement', lbl), 'FontWeight','normal', 'Color', INK);
    hold(ax,'off');

    % ---- NCC null --------------------------------------------------------
    ax = nexttile(tl, (s-1)*3 + 3);  hold(ax,'on');
    allv = []; for k=1:nc, allv = [allv; D{k}.R(s).ncc(:); D{k}.R(s).ncc0]; end %#ok<AGROW>
    ed = linspace(min(allv)-0.02, max(allv)+0.02, 34);
    for k = 1:nc
        histogram(ax, D{k}.R(s).ncc, ed, 'FaceColor', COL(k,:), 'FaceAlpha', 0.55, ...
                  'EdgeColor','none', 'HandleVisibility','off');
    end
    yl = ylim(ax);
    for k = 1:nc
        v = D{k}.R(s).ncc0;
        plot(ax, [v v], [0 yl(2)], '-', 'Color', COL(k,:), 'LineWidth', 2, ...
             'HandleVisibility','off');
        label_rule(ax, v, yl(2)*(1 - 0.10*(k-1)), INK2, max(6,8*FS));
    end
    style(ax, GRD, INK2);
    xlabel(ax, 'Pearson NCC', 'Color', INK);
    ylabel(ax, 'draws', 'Color', INK);
    title(ax, sprintf('%s    NCC null  (cross-check)', lbl), 'FontWeight','normal', 'Color', INK);
    hold(ax,'off');
end

o = D{1}.opt;
% The session NAME, not the full path: at this font the absolute path runs off
% both edges of a 1400 px canvas, and the path is already in the file's location.
% St.root IS the session folder, so one fileparts -- taking two walks up to the
% parent and labels every figure "Desktop".
[~, sess] = fileparts(D{1}.St.root);
title(tl, {sprintf(['stitch shift null, ch %s  --  %d shifts per seam, |dx| and |dy| ' ...
        'uniform integer in [%d, %d] um with random sign, z held at the curated value'], ...
        strjoin(arrayfun(@(c) sprintf('%d',c), opt.Channels, 'uni',0), ' vs ch '), ...
        o.NShift, o.ShiftMin, o.ShiftMax), ...
        sprintf(['%s, curated %s   |   solid rule = the curated placement, in its own ' ...
        'channel colour   |   window eroded by %d px, tile A fixed, background removed per patch'], ...
        sess, D{1}.St.savedOn, o.ShiftMax)}, ...
      'FontWeight','bold', 'Interpreter','none', 'FontSize', max(7, 9.5*FS));

base = fullfile(opt.Root, opt.OutName);
exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
exportgraphics(f, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
d = dir([base '.png']);
fprintf('saved %s.png (%.0f kB) / .pdf\n', base, d.bytes/1024);
if d.bytes < 20000
    warning('%s.png is only %d bytes -- likely an OpenGL export stub, re-run.', base, d.bytes);
end
end

function label_rule(ax, v, y, ink2, fs)
% Put the number on whichever side of the rule has room.  A label always drawn to
% the right runs off the axes whenever the curated value sits near the upper
% limit -- which, for a good alignment, is exactly where it always sits.
xl = xlim(ax);
if v > mean(xl)
    text(ax, v, y, sprintf('%.3f ', v), 'Color', ink2, 'FontSize', fs, ...
         'VerticalAlignment','top', 'HorizontalAlignment','right');
else
    text(ax, v, y, sprintf(' %.3f', v), 'Color', ink2, 'FontSize', fs, ...
         'VerticalAlignment','top', 'HorizontalAlignment','left');
end
end

function style(ax, grd, ink2)
% Recessive furniture: the marks carry the message, the axes stay quiet.
set(ax, 'Box','off', 'TickDir','out', 'XColor', ink2, 'YColor', ink2, ...
        'GridColor', grd, 'GridAlpha', 1, 'Layer','bottom');
grid(ax, 'on');
end
