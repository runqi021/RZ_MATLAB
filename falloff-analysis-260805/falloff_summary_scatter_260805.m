function falloff_summary_scatter_260805(rootDir)
% FALLOFF_SUMMARY_SCATTER_260805  Attenuation length across sessions and groups.
%
% Scans <session>\falloff_<tag>.mat -- one per FOV/ROI -- and scatters the fitted
% attenuation length l.  Both writers are read:
%   falloff_260804.m        vessel sessions, several stacks per FOV  -> Lp, r2p
%   falloff_gcamp_260730.m  one ramped soma stack per FOV            -> ell, r2
% A session folder is taken relative to ROOT unless it is an ABSOLUTE path, so a
% session that does not live under C:\fall-off can still be listed.
%
%   left   l per FOV, one column per session
%   right  the same points pooled by GROUP (wildtype vs shiverer)
%
% COLOUR IS THE FLUORESCENT LABEL, NOT THE GROUP: green = GFP (soma GCaMP),
% magenta = vessel dye.  The group is already carried by the x axis -- panel
% order on the left, the two columns on the right -- so colour is free to show
% the other thing that differs between FOVs, and a green point sitting in the
% shiverer column reads immediately as "that one is not vessel dye".
%
% Every point is one FOV.  Mean +- SD is drawn per session and per group, but
% the individual points are always shown: n is small enough that a bar alone
% would hide the spread, and on this data the two shiverer sessions do NOT
% overlap each other.
%
% GROUP ASSIGNMENT is a property of the SESSION FOLDER and is declared in the
% table below -- it is never guessed from the folder name.  Add a row to add a
% session or a new treatment group.
%
% l here is the decay constant of (top5% - bot5%)/mW^2 vs depth.  The
% literature "effective attenuation length" is 2*l, because two-photon signal
% goes as P^2 exp(-2z/l_e); both are reported.
%
% USAGE
%   falloff_summary_scatter_260805
%   falloff_summary_scatter_260805('C:\fall-off')
%
% OUTPUT (into rootDir)
%   falloff_summary_scatter.png / .pdf / .mat

%% --------------------------- USER PARAMETERS -------------------------------
ROOT = 'C:\fall-off';

%  folder                      group        display label                        label
SESSIONS = { ...
  '251104_wt_fitc',            'wildtype',  'WT 251104 (fitc)'        , 'vessel'
  '260728_vglut2_vessel',      'wildtype',  'WT 260728 (vglut2)'      , 'vessel'
  '260114_shiver_vessel',      'shiverer',  'Shiverer 260114'         , 'vessel'
  '260804_shiver_dbh_vessel',  'shiverer',  'Shiverer 260804 (dbh)'   , 'vessel'
  '260806_sst_vessel',         'wildtype',  'WT 260806 (sst)'         , 'vessel'
  'E:\260810_vglut2-soma-g8s\roi1_stack_z-5-70_30f_00001', 'wildtype',  'WT 260810 (vglut2 soma)', 'gfp'
  'C:\260730_shiver_vgat_g8m', 'shiverer',  'Shiverer 260730 (vgat-G8m, soma)', 'gfp' };
% NOTE 260730 is a SOMA GCaMP label where every other row is vessel dye, and it
% lives outside ROOT.  It pools into the shiverer group like any other FOV.
% The 4th column is DECLARED, not inferred: '260806_sst_vessel' saves an SST
% channel AND a vessel channel, and falloff_run_all reads channel 3 (the
% vessels), so the folder name alone would label it wrong.
% NOTE 260810 is a 70 um stack, so its lambda comes from a 10-65 um fit --
% a third of the depth range behind every other WT FOV.

GROUP_ORDER = {'wildtype', 'shiverer'};
LABEL_ORDER = {'gfp', 'vessel'};
LABEL_COL   = [0.00 0.62 0.30;    % gfp     -- green
               0.85 0.10 0.62];   % vessel  -- magenta
LABEL_LEG   = {'GFP (soma GCaMP)', 'vessel dye'};
BAR_COL     = [0.35 0.35 0.35];   % mean +- SD bars: neutral, because a
                                  % group can now contain both labels
HLINE_Y  = 140;              % reference line on the lambda axis (um); [] = none
SAVE_FIG = true;

% BOXED FITS OVERRIDE. falloff_metrics_box_260901 refits l inside the 200x200 um
% sub-volume that the 3D renders show, rather than over the whole 456 um field.
% Where it has a value for a (session, FOV) it REPLACES the whole-field one here,
% so this figure and the render figures quote the same number.
%
% THE POOL IS THEN MIXED, and that is the thing to be careful about: two FOVs are
% boxed and every other FOV is whole-field. Boxed l runs slightly SHORTER because
% the box excludes large surface vessels that keep `top` high at shallow depths
% and flatten the early decay (WT 23.5 vs 24.6, shiverer 60.1 vs 61.4 -- both
% within ~1.3 um). Overridden points are drawn with a black edge so the mixed
% provenance is visible rather than silent. Set BOX_MAT = '' to use whole-field
% l for everything.
BOX_MAT = 'C:\fall-off\falloff_metrics_box_260901.mat';
%% ---------------------------------------------------------------------------

if nargin >= 1 && ~isempty(rootDir), ROOT = rootDir; end

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fileparts(thisDir));

%% --- collect -----------------------------------------------------------------
T = struct('session',{},'group',{},'label',{},'ltype',{},'tag',{},'ell',{},'r2',{},'fit',{});
for s = 1:size(SESSIONS,1)
    d = SESSIONS{s,1};
    if ~(startsWith(d, filesep) || ~isempty(regexp(d, '^[A-Za-z]:', 'once')))
        d = fullfile(ROOT, d);       % a plain name is relative to ROOT
    end
    if ~isfolder(d)
        warning('missing session folder: %s', d); continue
    end
    % falloff_<tag>.mat, excluding the *_metrics.mat cache and the proj volumes
    hits = dir(fullfile(d, 'falloff_*.mat'));
    hits = hits(~contains({hits.name}, {'_metrics.mat','_vol.mat'}) & ...
                ~startsWith({hits.name}, 'falloff_proj_'));
    for h = 1:numel(hits)
        S = load(fullfile(d, hits(h).name));
        % falloff_260804 writes the POOLED fit as Lp/r2p; falloff_gcamp_260730 has
        % only one stack so it writes it as ell/r2.  Same quantity either way.
        if     isfield(S,'Lp'),  L = S.Lp;   R2 = S.r2p;
        elseif isfield(S,'ell'), L = S.ell;  R2 = S.r2;
        else,  continue
        end
        if ~isfinite(L), continue; end
        n = numel(T) + 1;
        T(n).session = SESSIONS{s,1};
        T(n).group   = SESSIONS{s,2};
        T(n).label   = SESSIONS{s,3};
        T(n).ltype   = SESSIONS{s,4};
        T(n).tag     = regexprep(hits(h).name, '^falloff_|\.mat$', '');
        T(n).ell     = L;
        T(n).r2      = R2;
        T(n).fit     = S.FIT_DEPTH;
        T(n).boxed   = false;
    end
end

% --- apply the boxed fits, where there is one for this (session, FOV) ---
if ~isempty(BOX_MAT) && isfile(BOX_MAT)
    Bx = load(BOX_MAT, 'BOX');
    for b = 1:numel(Bx.BOX)
        j = find(strcmp({T.session}, Bx.BOX(b).session) & ...
                 strcmp({T.tag},     Bx.BOX(b).tag), 1);
        if isempty(j)
            fprintf(2, 'boxed fit for %s/%s has no matching FOV in the scatter\n', ...
                    Bx.BOX(b).session, Bx.BOX(b).tag);
            continue;
        end
        fprintf('boxed override %-26s %-6s  l %.1f -> %.1f  (fit %g-%g)\n', ...
                T(j).session, T(j).tag, T(j).ell, Bx.BOX(b).ell, Bx.BOX(b).fit);
        T(j).ell   = Bx.BOX(b).ell;
        T(j).r2    = Bx.BOX(b).r2;
        T(j).fit   = Bx.BOX(b).fit;
        T(j).boxed = true;
    end
elseif ~isempty(BOX_MAT)
    fprintf(2, 'BOX_MAT not found (%s) -- using whole-field l for every FOV\n', BOX_MAT);
end
assert(~isempty(T), 'No falloff_<tag>.mat found under %s', ROOT);

fprintf('\n%-26s %-10s %-16s %8s %8s %14s\n', ...
    'session', 'group', 'FOV', 'lambda (um)', 'R^2', 'fit range');
fprintf('%s\n', repmat('-', 1, 88));
for i = 1:numel(T)
    fprintf('%-26s %-10s %-16s %8.1f %8.3f %6.0f-%-6.0f\n', T(i).session, ...
        T(i).group, T(i).tag, T(i).ell, T(i).r2, T(i).fit(1), T(i).fit(2));
end

%% --- figure -------------------------------------------------------------------
f = figure('Color','w','Position',[80 80 1250 560]);
tl = tiledlayout(f, 1, 2, 'TileSpacing','compact','Padding','compact');

% --- per session ---
[uSess, ia] = unique({T.session}, 'stable');
uGrp  = {T(ia).group};
uLbl  = {T(ia).label};
uLtp  = {T(ia).ltype};
[~, ord] = sort(cellfun(@(g) find(strcmp(GROUP_ORDER, g)), uGrp));
uSess = uSess(ord);  uGrp = uGrp(ord);  uLbl = uLbl(ord);  uLtp = uLtp(ord);

ax = nexttile; hold(ax,'on');
hline(ax, HLINE_Y);
for k = 1:numel(uSess)
    isS = strcmp({T.session}, uSess{k});
    v   = [T(isS).ell];
    bx  = [T(isS).boxed];
    c = LABEL_COL(strcmp(LABEL_ORDER, uLtp{k}), :);
    jit = (rand(size(v)) - 0.5) * 0.28;      % jitter so overlapping FOVs show
    plot(ax, k + jit(~bx), v(~bx), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor','none', ...
        'MarkerSize', 7);
    % Boxed refits get a black edge: same quantity, different spatial support,
    % and a reader comparing columns should be able to see which is which.
    plot(ax, k + jit(bx), v(bx), 'o', 'MarkerFaceColor', c, 'MarkerEdgeColor','k', ...
        'LineWidth', 1.1, 'MarkerSize', 8);
    plot(ax, k + [-0.3 0.3], [1 1]*mean(v), '-', 'Color', c*0.6, 'LineWidth', 2.5);
    if numel(v) > 1
        plot(ax, [k k], mean(v) + [-1 1]*std(v), '-', 'Color', c*0.6, 'LineWidth', 1.4);
    end
    text(ax, k, max(v) + 4, sprintf('n=%d', numel(v)), ...
        'HorizontalAlignment','center', 'FontSize', 9, 'Color', c*0.6);
end
set(ax, 'XTick', 1:numel(uSess), 'XTickLabel', uLbl, 'XTickLabelRotation', 20);
xlim(ax, [0.4 numel(uSess)+0.6]);
grid(ax,'on'); box(ax,'on');
ylabel(ax, 'length constant  \lambda  (\mum)');
title(ax, 'per session -- each dot is one FOV', 'FontWeight','normal');

% --- pooled by group ---
ax = nexttile; hold(ax,'on');
hline(ax, HLINE_Y);
for g = 1:numel(GROUP_ORDER)
    ig = find(strcmp({T.group}, GROUP_ORDER{g}));
    v  = [T(ig).ell];
    if isempty(v), continue; end
    jit = (rand(size(v)) - 0.5) * 0.30;
    % One plot call per point: a group can mix label types, so the colour has
    % to come from the POINT, not from the group.
    for i = 1:numel(ig)
        c = LABEL_COL(strcmp(LABEL_ORDER, T(ig(i)).ltype), :);
        if T(ig(i)).boxed
            plot(ax, g + jit(i), v(i), 'o', 'MarkerFaceColor', c, ...
                'MarkerEdgeColor','k', 'LineWidth',1.1, 'MarkerSize', 9);
        else
            plot(ax, g + jit(i), v(i), 'o', 'MarkerFaceColor', c, ...
                'MarkerEdgeColor','none', 'MarkerSize', 8);
        end
    end
    plot(ax, g + [-0.32 0.32], [1 1]*mean(v), '-', 'Color', BAR_COL, 'LineWidth', 3);
    plot(ax, [g g], mean(v) + [-1 1]*std(v), '-', 'Color', BAR_COL, 'LineWidth', 1.6);
    fprintf('\n%-10s n=%2d   lambda = %.1f +- %.1f um  (range %.1f-%.1f)   EAL = 2*lambda = %.1f um\n', ...
        GROUP_ORDER{g}, numel(v), mean(v), std(v), min(v), max(v), 2*mean(v));
end
% Legend off dummy handles, so the label types are listed once instead of once
% per point.  NaN data draws nothing but still carries the marker style.
hL = gobjects(numel(LABEL_ORDER), 1);
for i = 1:numel(LABEL_ORDER)
    hL(i) = plot(ax, nan, nan, 'o', 'MarkerFaceColor', LABEL_COL(i,:), ...
        'MarkerEdgeColor','none', 'MarkerSize', 8);
end
legend(ax, hL, LABEL_LEG, 'Location','northeast', 'Box','off', 'FontSize', 9);
set(ax, 'XTick', 1:numel(GROUP_ORDER), 'XTickLabel', GROUP_ORDER);
xlim(ax, [0.4 numel(GROUP_ORDER)+0.6]);
grid(ax,'on'); box(ax,'on');
ylabel(ax, 'length constant  \lambda  (\mum)');
title(ax, 'pooled by group -- bar = mean \pm SD', 'FontWeight','normal');

% one shared y range so the two panels are directly comparable
% the reference line is INSIDE the range: with lambda ~25-86 um a 140 um line
% would otherwise be clipped off the top and silently not drawn.
yl = [0, max([[T.ell], HLINE_Y])*1.15];
arrayfun(@(a) ylim(a, yl), findobj(f,'Type','axes'));

title(tl, {sprintf('%s  --  length constant by session and group', ROOT), ...
    ['lambda = decay constant of (top 5% - bot 5%)/mW^2;  literature EAL = 2*lambda' ...
     '   |   colour = label: green GFP, magenta vessel dye']}, ...
    'FontWeight','bold','Interpreter','none','FontSize',10);

if SAVE_FIG
    base = fullfile(ROOT, 'falloff_summary_scatter');
    exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
    exportgraphics(f, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'T', 'SESSIONS');
    fprintf('\nsaved %s.png/.pdf/.mat\n', base);
end
end

% --------------------------------------------------------------------------
function hline(ax, y)
% Horizontal reference line on the lambda axis, drawn behind the data.
if isempty(y), return; end
% [num2str ' \mum'], not sprintf: sprintf would try to read \m as an escape.
yline(ax, y, '--', [num2str(y) ' \mum'], 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left', ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 9);
end
