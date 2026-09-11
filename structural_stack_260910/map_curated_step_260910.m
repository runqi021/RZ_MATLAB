function T = map_curated_step_260910(varargin)
%MAP_CURATED_STEP_260910  Empirical tile step and overlap from a curated stitch.
%
%   T = map_curated_step_260910
%   T = map_curated_step_260910('Root', '...\map')
%
% Reads `stitch_manual_coords.mat` from stitch_manual_gui_fast_260818 and asks
% what the CURATED placement says the tile step actually is -- in pixels, in
% microns of overlap, and, by dividing into the stage travel the filenames
% record, in microns per pixel.
%
% WHY THIS IS WORTH DOING.  The stage step and the image step are two
% independent rulers for the same distance.  The stage says how far it was told
% to travel; the curated placement says how far the tissue actually moved in the
% image.  Their ratio is the lateral scale calibration, measured on real data,
% and it needs no test target -- an 81-tile curated mosaic has already measured
% it 144 times.
%
% ONLY ADJACENT PAIRS IN THE SAME ROW (or column) ARE USED, and each pair
% contributes once.  Comparing tile 1 to tile 5 would give the same answer with
% less noise, but it would also hide any DRIFT in the step, which is exactly what
% distinguishes a scale error (constant step, wrong size) from a stage fault
% (step varies with position or direction).
%
% THE STEP IS SIGNED BY THE STAGE, not assumed.  Image column runs along +stage x
% and image row along -stage y on this rig, so the per-pair micron travel is read
% from the filenames' own x/y tokens rather than taken as a constant: a skipped
% or repeated position would otherwise be silently averaged in.
%
% OUTPUT  -> <Root>\matlab_stitch\curated_step\
%   curated_step.png/.pdf   overlap and um/px distributions, plus a drift check
%   curated_step.csv        every adjacent pair
%
% Runqi Zhang / 2026-09-10

p = inputParser;
p.addParameter('Root', 'C:\Users\Admin\Desktop\260909_ChAT_g8m_Shiverer\map');
p.addParameter('NominalUmPerPx', []);   % [] = take umPerPx from the coords file
p.parse(varargin{:});
opt = p.Results;

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); addpath(here);
try, opengl('software'); catch, end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

COL  = [ 42 120 214;  235 104  52] / 255;   % categorical slots 1, 2 (validated)
INK  = [0.05 0.05 0.05];
INK2 = [0.32 0.32 0.31];
GRD  = [0.88 0.88 0.87];

cf = fullfile(opt.Root, 'matlab_stitch', 'stitch_manual_coords.mat');
assert(isfile(cf), 'no curated coords at %s', cf);
S = load(cf);

umpp = opt.NominalUmPerPx;
if isempty(umpp), umpp = S.umPerPx; end

% S.x IS ALREADY THE CURATED POSITION.  stitch_manual_gui_fast_260818 computes
% x = (colIdx-1)*stepXpx + manX and saves that, keeping manX alongside only as a
% record of the nudge -- so adding manX here counts every nudge twice.  It does
% not error, it just inflates the step (432.5 px instead of 418.75) and with it
% the scale the whole script exists to measure.  The stored canvas size is the
% guard: the placement must span exactly outW x outH.
%
% NOTE the clash with stack_stitch_gui_260910, which stores base and nudge
% SEPARATELY and there x + manX IS correct.  Two stitchers, two conventions.
X = S.x(:);
Y = S.y(:);
assert(abs((max(X)-min(X)+S.W0) - S.outW) < 1e-6 && ...
       abs((max(Y)-min(Y)+S.H0) - S.outH) < 1e-6, ...
    ['placement spans %g x %g but the file records a %g x %g canvas -- the ' ...
     'coords convention is not what this script assumes'], ...
     max(X)-min(X)+S.W0, max(Y)-min(Y)+S.H0, S.outW, S.outH);
ci = S.colIdx(:);  ri = S.rowIdx(:);
W0 = S.W0;  H0 = S.H0;

fprintf('\n============== map_curated_step_260910 ==============\n');
fprintf('coords     : %s\n', cf);
fprintf('tiles      : %d   grid %d x %d   tile %d x %d px\n', numel(X), S.nRows, S.nCols, H0, W0);
fprintf('placeMode  : %s   (%s)\n', S.placeMode, S.placeSrc);
fprintf('umPerPx    : %.4f   nominal grid step %g px = %.1f um\n', umpp, S.stepXpx, S.stepXpx*umpp);
fprintf('curated    : %d tiles nudged in x, %d in y  (already folded into x/y)\n', ...
        nnz(S.manX ~= 0), nnz(S.manY ~= 0));
fprintf(['nudge grain: %g px -- the GUI nudges in whole steps, so the per-pair spread\n' ...
         '             below is largely QUANTISATION, not measurement scatter\n'], ...
        nudge_grain([S.manX(:); S.manY(:)]));

% Stage position per tile, from the FILENAME tokens the acquisition wrote.
sx = S.stageX(:);  sy = S.stageY(:);

%% ---------------------- adjacent pairs, both axes ---------------------------
P = struct('axis',{},'i',{},'j',{},'dpx',{},'dum',{},'umpp',{},'pos',{});
for k = 1:numel(X)
    % + 1 column: the neighbour to the right
    j = find(ri == ri(k) & ci == ci(k)+1, 1);
    if ~isempty(j)
        P(end+1) = pair('x', k, j, X(j)-X(k), sx(j)-sx(k), ci(k)); %#ok<AGROW>
    end
    j = find(ci == ci(k) & ri == ri(k)+1, 1);
    if ~isempty(j)
        P(end+1) = pair('y', k, j, Y(j)-Y(k), sy(k)-sy(j), ri(k)); %#ok<AGROW>
    end
end
assert(~isempty(P), 'no adjacent pairs found -- are colIdx/rowIdx populated?');

isx = strcmp({P.axis}, 'x');
dpxX = [P(isx).dpx].';  dumX = [P(isx).dum].';  posX = [P(isx).pos].';
dpxY = [P(~isx).dpx].'; dumY = [P(~isx).dum].'; posY = [P(~isx).pos].';

ovX_px = W0 - dpxX;  ovY_px = H0 - dpxY;
uX = dumX ./ dpxX;   uY = dumY ./ dpxY;      % microns per pixel, per pair

report('x (columns)', dpxX, ovX_px, dumX, uX, umpp);
report('y (rows)',    dpxY, ovY_px, dumY, uY, umpp);

uAll = [uX; uY];
fprintf('\nPOOLED over %d adjacent pairs\n', numel(uAll));
fprintf('  empirical  %.4f +/- %.4f um/px   (median %.4f)\n', mean(uAll), std(uAll), median(uAll));
fprintf('  header     %.4f um/px\n', umpp);
fprintf('  ratio      %.4f   ->  header understates the field by %.1f%%\n', ...
        mean(uAll)/umpp, 100*(mean(uAll)/umpp - 1));
fprintf('  implied objectiveResolution %.2f um/deg  (header says 38.90)\n', 38.9*mean(uAll)/umpp);

%% -------------------------------- figure -----------------------------------
f = figure('Color','w','Units','pixels','Position',[40 40 1400 780]);
drawnow; got = get(f,'Position'); FS = min(1, got(3)/1400);
set(f,'DefaultAxesFontSize', max(7,10*FS));
tl = tiledlayout(f, 2, 3, 'TileSpacing','compact','Padding','compact');

ax = nexttile(tl,1); hold(ax,'on');
ed = binedges([ovX_px; ovY_px], 26);
histogram(ax, ovX_px, ed, 'FaceColor',COL(1,:),'FaceAlpha',0.55,'EdgeColor','none','DisplayName','x (columns)');
histogram(ax, ovY_px, ed, 'FaceColor',COL(2,:),'FaceAlpha',0.55,'EdgeColor','none','DisplayName','y (rows)');
style(ax,GRD,INK2); xlabel(ax,'overlap width (px)','Color',INK); ylabel(ax,'pairs','Color',INK);
title(ax,'overlap width, curated','FontWeight','normal','Color',INK);
lg = legend(ax,'Location','northwest'); lg.Box='off'; lg.TextColor=INK2; hold(ax,'off');

ax = nexttile(tl,2); hold(ax,'on');
ed = binedges([dpxX; dpxY], 26);
histogram(ax, dpxX, ed, 'FaceColor',COL(1,:),'FaceAlpha',0.55,'EdgeColor','none');
histogram(ax, dpxY, ed, 'FaceColor',COL(2,:),'FaceAlpha',0.55,'EdgeColor','none');
yl = ylim(ax);
plot(ax,[S.stepXpx S.stepXpx],[0 yl(2)],'-','Color',INK2,'LineWidth',1.5);
text(ax,S.stepXpx,yl(2),sprintf(' grid %g',S.stepXpx),'Color',INK2,'FontSize',max(6,8*FS), ...
     'VerticalAlignment','top');
style(ax,GRD,INK2); xlabel(ax,'step between adjacent tiles (px)','Color',INK); ylabel(ax,'pairs','Color',INK);
title(ax,'empirical step','FontWeight','normal','Color',INK); hold(ax,'off');

ax = nexttile(tl,3); hold(ax,'on');
ed = binedges(uAll, 26);
histogram(ax, uX, ed, 'FaceColor',COL(1,:),'FaceAlpha',0.55,'EdgeColor','none');
histogram(ax, uY, ed, 'FaceColor',COL(2,:),'FaceAlpha',0.55,'EdgeColor','none');
yl = ylim(ax);
plot(ax,[umpp umpp],[0 yl(2)],'-','Color',INK2,'LineWidth',1.5);
text(ax,umpp,yl(2),sprintf(' header %.4f',umpp),'Color',INK2,'FontSize',max(6,8*FS), ...
     'VerticalAlignment','top');
style(ax,GRD,INK2); xlabel(ax,'stage travel / image step   (\mum per px)','Color',INK);
ylabel(ax,'pairs','Color',INK);
title(ax,'lateral scale, measured','FontWeight','normal','Color',INK); hold(ax,'off');

ax = nexttile(tl,4); hold(ax,'on');
scatter(ax, posX, dpxX, 16, COL(1,:), 'filled', 'MarkerFaceAlpha',0.6);
scatter(ax, posY, dpxY, 16, COL(2,:), 'filled', 'MarkerFaceAlpha',0.6);
style(ax,GRD,INK2); xlabel(ax,'index of the left / upper tile','Color',INK);
ylabel(ax,'step (px)','Color',INK);
title(ax,'drift check: step vs position','FontWeight','normal','Color',INK); hold(ax,'off');

ax = nexttile(tl,5); hold(ax,'on');
scatter(ax, dumX, dpxX, 16, COL(1,:), 'filled', 'MarkerFaceAlpha',0.6);
scatter(ax, dumY, dpxY, 16, COL(2,:), 'filled', 'MarkerFaceAlpha',0.6);
xl = [min([dumX;dumY])-5, max([dumX;dumY])+5];
plot(ax, xl, xl/umpp, '-',  'Color', INK2, 'LineWidth', 1.2);
plot(ax, xl, xl/mean(uAll), '--', 'Color', INK2, 'LineWidth', 1.2);
text(ax, xl(2), xl(2)/umpp, sprintf('header %.4f ', umpp), 'Color',INK2, ...
     'FontSize',max(6,8*FS),'HorizontalAlignment','right','VerticalAlignment','bottom');
text(ax, xl(2), xl(2)/mean(uAll), sprintf('measured %.4f ', mean(uAll)), 'Color',INK2, ...
     'FontSize',max(6,8*FS),'HorizontalAlignment','right','VerticalAlignment','top');
style(ax,GRD,INK2); xlabel(ax,'stage travel (\mum)','Color',INK); ylabel(ax,'image step (px)','Color',INK);
title(ax,'two rulers for one distance','FontWeight','normal','Color',INK); hold(ax,'off');

ax = nexttile(tl,6); hold(ax,'on');
ed = binedges([ovX_px*mean(uAll); ovY_px*mean(uAll)], 26);
histogram(ax, ovX_px*mean(uAll), ed, 'FaceColor',COL(1,:),'FaceAlpha',0.55,'EdgeColor','none');
histogram(ax, ovY_px*mean(uAll), ed, 'FaceColor',COL(2,:),'FaceAlpha',0.55,'EdgeColor','none');
style(ax,GRD,INK2); xlabel(ax,'overlap width (\mum, at the MEASURED scale)','Color',INK);
ylabel(ax,'pairs','Color',INK);
title(ax,'overlap in real microns','FontWeight','normal','Color',INK); hold(ax,'off');

[~, sess] = fileparts(fileparts(opt.Root));
title(tl, {sprintf(['curated tile step and overlap  --  %s\\map  --  %d tiles, %d adjacent pairs ' ...
        '(%d in x, %d in y)'], sess, numel(X), numel(P), nnz(isx), nnz(~isx)), ...
        sprintf(['image step %.1f +/- %.1f px for %.1f +/- %.1f um of stage travel  ->  %.4f um/px ' ...
        'measured vs %.4f in the header  (%+.1f%%)'], ...
        mean([dpxX;dpxY]), std([dpxX;dpxY]), mean([dumX;dumY]), std([dumX;dumY]), ...
        mean(uAll), umpp, 100*(mean(uAll)/umpp - 1))}, ...
      'FontWeight','bold','Interpreter','none','FontSize',max(7,9.5*FS));

outDir = fullfile(opt.Root, 'matlab_stitch', 'curated_step');
if ~isfolder(outDir), mkdir(outDir); end
base = fullfile(outDir, 'curated_step');
exportgraphics(f, [base '.png'], 'Resolution', 200, 'BackgroundColor','white');
exportgraphics(f, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');

T = table({P.axis}.', [P.i].', [P.j].', [P.dpx].', [P.dum].', [P.umpp].', [P.pos].', ...
    'VariableNames', {'axis','tileA','tileB','step_px','stage_um','um_per_px','index'});
writetable(T, [base '.csv']);
fprintf('\nsaved %s.png / .pdf / .csv\n', base);
end

% =============================================================================
function s = pair(ax, i, j, dpx, dum, pos)
% MAGNITUDES.  This dataset carries flipY = true, so image rows run opposite to
% stage y and the signed row step comes out negative -- which silently turned
% `tile - step` into an overlap of 960 px on a 512 px tile.  Only the distance
% between neighbours is being measured here; its sign is a property of the
% display convention, not of the step.
dpx = abs(dpx);  dum = abs(dum);
s = struct('axis', ax, 'i', i, 'j', j, 'dpx', dpx, 'dum', dum, ...
           'umpp', dum/max(eps,dpx), 'pos', pos);
end

function report(name, dpx, ovpx, dum, u, umpp)
fprintf('\n%s   %d pairs\n', name, numel(dpx));
fprintf('  step        %.2f +/- %.2f px   (range %g .. %g)\n', mean(dpx), std(dpx), min(dpx), max(dpx));
fprintf('  overlap     %.2f +/- %.2f px   = %.1f um at the header scale\n', ...
        mean(ovpx), std(ovpx), mean(ovpx)*umpp);
fprintf('  stage       %.2f +/- %.2f um\n', mean(dum), std(dum));
fprintf('  -> scale    %.4f +/- %.4f um/px\n', mean(u), std(u));
end

function g = nudge_grain(m)
%NUDGE_GRAIN  The step the curation was actually done at, read off the nudges.
% Their GCD, because every nudge is a whole number of the GUI's move step.
m = unique(abs(m(m ~= 0)));
if isempty(m), g = NaN; return; end
g = m(1);
for k = 2:numel(m), g = gcd(g, m(k)); end
end

function ed = binedges(v, n)
lo = min(v); hi = max(v);
if hi <= lo, hi = lo + 1; end
ed = linspace(lo - 0.02*(hi-lo), hi + 0.02*(hi-lo), n);
end

function style(ax, grd, ink2)
set(ax,'Box','off','TickDir','out','XColor',ink2,'YColor',ink2, ...
       'GridColor',grd,'GridAlpha',1,'Layer','bottom');
grid(ax,'on');
end
