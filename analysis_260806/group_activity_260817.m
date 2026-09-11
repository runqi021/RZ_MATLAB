% group_activity_260817.m
% -----------------------------------------------------------------------
%  Population activity per (genotype x class), consistent with the ventral map.
%
%      colour     = GENOTYPE   (genotype_colors_260817, the same table the map
%                               tints its anatomy patches and cells with)
%      line style = CLASS      solid I | dash-dot E (post-I) | dashed pre-I |
%                              dotted tonic/rhythmic null
%
%  WHY NOT REUSE polar_selected's overlay. That one takes its groups from
%  group_display_colors_260816, a list of eight (genotype, class) pairs. Any pair
%  missing from that list is filtered out with no message, and after the classes
%  were extended on 2026-08-17 that silently dropped FOURTEEN cells -- ChAT null,
%  Sert null, Vgat null and Vglut2 null, none of which had an entry. Splitting
%  colour (genotype) from line style (class) removes the need for a per-pair list
%  entirely: any class the data contains gets drawn, and a class with no style is
%  reported rather than dropped.
%
%  VENTRAL SURFACE ONLY, to match the map: cells deeper than 120 um are excluded.
%  Depth is |motorPosition(3)|, never the z in the filename.
%
%  MEAN +/- SEM ACROSS CELLS, each cell counting once regardless of how many
%  events it fired -- otherwise one 315-event Vglut2 cell outweighs twenty
%  8-event IO cells. Raw AND z-scored dF/F are both shown: raw is what was
%  recorded but is set by the brightest cells, z-scored gives every cell equal
%  weight and shows the SHAPE of the modulation.
%
%  Reads the per-cell curves already computed by polar_selected_260816.m
%  (R.dffNo / R.dffNoz / R.histNo on a fixed +/-1.5 s grid), so nothing is
%  recomputed here and the numbers cannot drift from that run.
%
%  Runqi Zhang / 2026-08-17
% -----------------------------------------------------------------------

clear; clc;
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

%% ===================== USER-EDITABLE =====================
sumRoot  = 'D:\Ventral_surface_summary';
selMat   = fullfile(sumRoot,'polar_selected_260816','polar_selected_260816.mat');
regFile  = fullfile(sumRoot,'event_latency_260811','event_latency_data.mat');
outDir   = fullfile(sumRoot,'group_activity_260817');
WIN_S    = 1.5;          % must match polar_selected_260816 (its tauN grid)
DEEP_UM  = 120;          % ventral-surface bar, same as the map
showIO   = false;        % IO is 107 of 157 cells and swamps the rest
% line style per class. A class not listed here is REPORTED, never dropped.
STY = { 'I',                    '-'      % solid
        'post-I',               '-.'     % dash-dot  (drawn as E on the map)
        'pre-I',                '--'     % dashed
        'tonic/rhythmic null',  ':' };   % dotted
LW  = 1.8;
% =========================================================

assert(isfile(selMat), 'run polar_selected_260816.m first: %s missing', selMat);
if ~isfolder(outDir), mkdir(outDir); end
S = load(selMat,'R');  R = S.R;
GC = genotype_colors_260817();

tauN = linspace(-WIN_S, WIN_S, 201);
eN   = linspace(-WIN_S, WIN_S, 42);  ctrN = (eN(1:end-1)+eN(2:end))/2;
assert(numel(R(1).dffNo)==numel(tauN), 'grid mismatch: WIN_S here disagrees with the .mat');

%% ---- classes read LIVE, not from the .mat ----
% The .mat carries each cell's class as it stood when polar_selected last ran.
% Reading that would mean a class edit silently plots the OLD grouping until the
% polar is re-run -- a wrong figure that looks right. The CURVES in the .mat do
% not depend on class, so they are still used; only the labels come from
% cell_classes_260816, which is the authority.
clear cell_classes_260816
CT = cell_classes_260816();
selCells = [R.cell].';
[inTbl, loc] = ismember(selCells, CT.cell);
% A classified cell with no curves in the .mat cannot be averaged. That happens
% when a cell is newly classified and polar_selected has not been re-run since,
% and it MUST be loud: silently averaging 7 of 8 cells and printing n=7 is
% indistinguishable from a group that really has 7.
absent = setdiff(CT.cell, selCells);
if ~isempty(absent)
    fprintf(2, ['%d classified cell(s) have NO curves in the .mat and are NOT ' ...
                'in this average: %s\n' ...
                '   re-run polar_selected_260816.m to add them.\n'], ...
            numel(absent), mat2str(absent(:).'));
end
R = R(inTbl);  loc = loc(inTbl);
grp = string(CT.group(loc));
cls = string(CT.class(loc));
cls(ismember(cls, ["tonic pre-I suppressed","null"])) = "tonic/rhythmic null";

D = load(regFile,'CELL','OBS','REC');
obsOf = pooled_obs_260814(D.CELL, D.OBS);
isDeep = false(numel(R),1);
for k = 1:numel(R)
    c = R(k).cell;
    if c > numel(obsOf) || isempty(obsOf{c}), continue; end
    dp = [];
    for o = obsOf{c}(:)'
        fp = D.REC(D.OBS(o).rec).folder;
        mh = dir(fullfile(fp,'*_ch1_meta.mat'));
        if isempty(mh), mh = dir(fullfile(fp,'*_meta.mat')); end
        if isempty(mh), continue; end
        M = load(fullfile(mh(1).folder,mh(1).name));
        if isfield(M,'motorPosition') && numel(M.motorPosition)>=3
            dp(end+1) = abs(M.motorPosition(3)); %#ok<SAGROW>
        end
    end
    isDeep(k) = ~isempty(dp) && all(dp > DEEP_UM);
end
if any(isDeep)
    fprintf('ventral surface only: dropped %d deep cell(s): %s\n', ...
            nnz(isDeep), mat2str([R(isDeep).cell]));
end
keep = ~isDeep;
if ~showIO, keep = keep & (grp ~= "IO"); end
R = R(keep);  grp = grp(keep);  cls = cls(keep);
fprintf('%d cells in the overlay\n', numel(R));

%% ---- one line per (genotype, class) present ----
key = grp + "|" + cls;
uk  = unique(key);
G = struct('lab',{},'col',{},'sty',{},'n',{},'idx',{});
for k = 1:numel(uk)
    pr = split(uk(k),"|");  g = pr(1);  c = pr(2);
    idx = find(key == uk(k));
    if ~isfield(GC, g)
        fprintf(2,'NO COLOUR for genotype "%s" -- %d cell(s) not drawn\n', g, numel(idx));
        continue;
    end
    si = find(strcmp(STY(:,1), c), 1);
    if isempty(si)
        fprintf(2,'NO LINE STYLE for class "%s" -- %d cell(s) not drawn\n', c, numel(idx));
        continue;
    end
    G(end+1) = struct('lab', sprintf('%s %s (n=%d)', g, c, numel(idx)), ...
                      'col', GC.(g), 'sty', STY{si,2}, 'n', numel(idx), ...
                      'idx', idx); %#ok<SAGROW>
end
[~,ord] = sort(string({G.lab}));  G = G(ord);

%% ---- figure ----
PAN = { 'dffNo',  tauN, 'mean \DeltaF/F \pm SEM',          '\DeltaF/F'
        'dffNoz', tauN, 'mean z-scored \DeltaF/F \pm SEM', 'z-scored \DeltaF/F'
        'histNo', ctrN, 'event histogram',                 'spk/cyc %' };
fig = figure('Color','w','Units','normalized','Position',[0.03 0.25 0.94 0.5], ...
             'Name','group activity');
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
dOnPk = median([R.dOnPk], 'omitnan');

for pi = 1:size(PAN,1)
    ax = nexttile(tl); hold(ax,'on');
    xv = PAN{pi,2};
    for q = 1:numel(G)
        Mx = cell2mat({R(G(q).idx).(PAN{pi,1})}.');
        mu = mean(Mx,1,'omitnan');
        se = std(Mx,0,1,'omitnan') ./ sqrt(max(sum(~isnan(Mx),1),1));
        ok = ~isnan(mu);
        fill(ax, [xv(ok) fliplr(xv(ok))], [mu(ok)+se(ok) fliplr(mu(ok)-se(ok))], ...
             G(q).col, 'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');
        plot(ax, xv, mu, G(q).sty, 'Color',G(q).col, 'LineWidth',LW, ...
             'DisplayName',G(q).lab);
    end
    % red = trigger (inspiration onset), blue = mean breath peak after it
    xline(ax, 0, '-', 'Color',[1.00 0.40 0.75], 'LineWidth',1.2, 'HandleVisibility','off');
    if isfinite(dOnPk)
        xline(ax, dOnPk, '-', 'Color',[0.35 0.55 1.00], 'LineWidth',1.2, 'HandleVisibility','off');
    end
    grid(ax,'on'); xlim(ax,[-WIN_S WIN_S]);
    xlabel(ax,'time from insp onset (s)'); ylabel(ax, PAN{pi,4});
    title(ax, PAN{pi,3});
    if pi == 1, legend(ax,'Location','northwest','Box','off','Interpreter','none'); end
end
title(tl, sprintf(['inspiration-onset-triggered   |   colour = genotype, ' ...
      'line style = class   |   %d cells, ventral surface only%s'], ...
      numel(R), ternary_local(showIO,'','  (IO not shown)')), ...
      'FontWeight','bold');

exportgraphics(fig, fullfile(outDir,'group_activity.png'), ...
               'Resolution',300,'BackgroundColor','white');
exportgraphics(fig, fullfile(outDir,'group_activity.pdf'), ...
               'ContentType','vector','BackgroundColor','white');
T = table(string({G.lab}).', [G.n].', 'VariableNames',{'group','n'});
writetable(T, fullfile(outDir,'group_activity_groups.csv'));
fprintf('saved -> %s\n', outDir);
disp(T);

function out = ternary_local(c,a,b), if c, out=a; else, out=b; end, end
