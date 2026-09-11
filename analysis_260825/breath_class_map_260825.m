% breath_class_map_260825.m
% -----------------------------------------------------------------------
%  Draw the average chest-motion maps per breath-cycle class, from the .mat
%  written by breath_class_map_260825.py.
%
%  WHAT IS PLOTTED
%    row 1  MEAN DISPLACEMENT MAP, peak frame minus foot frame, averaged over
%           that class's cycles. Signed, so it shows WHICH WAY the chest moved
%           and where. Same quantity as the pipeline's diff_imgs / diffImg, just
%           split by class, so it can be read against u1_img directly.
%    row 2  MAGNITUDE MAP, mean |peak - foot|. Direction-free, for the case
%           where a class's displacements partly cancel in the signed average --
%           if the magnitude is large where the signed map is near zero, the
%           motion is there but inconsistent in direction across cycles.
%
%  SCALING: normal and gasp share one symmetric colour limit inside each row, so
%  the two panels are directly comparable by eye. The difference panel gets its
%  own limit -- forcing it onto the same scale would flatten it to nothing.
%
%  POOLING ACROSS RECORDINGS is count-weighted: a recording contributing 82
%  normal cycles counts more than one contributing 21. Every recording is on the
%  same anatomical grid already, because the extractor applied each run's own
%  registration shift before cropping.
%
%  Runqi Zhang / 2026-08-25
% -----------------------------------------------------------------------

clear; close all; clc;

%% ===================== USER-EDITABLE =====================
sessionDir = 'C:\260824_Vglut2-soma-g8s_vagotomized\phys';
mapFile    = fullfile(sessionDir, 'breath_class_maps.mat');
outDir     = fullfile(sessionDir, 'breath_class_maps');
showClasses = {'normal','gasp'};     % the two to compare
doSave     = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end
assert(isfile(mapFile), 'not found: %s\nRun breath_class_map_260825.py first.', mapFile);
M = load(mapFile);
if doSave && ~isfolder(outDir), mkdir(outDir); end

names = cellstr(string(M.class_names(:)));
ia = find(strcmp(names, showClasses{1}));
ib = find(strcmp(names, showClasses{2}));
assert(~isempty(ia) && ~isempty(ib), 'class names not found in the map file');

cnt = double(M.counts);                       % [nRun x nClass]
fprintf('\n=========== breath_class_map_260825 ===========\n');
fprintf('%s: %d cycles over %d recordings\n', names{ia}, sum(cnt(:,ia)), nnz(cnt(:,ia)));
fprintf('%s: %d cycles over %d recordings\n', names{ib}, sum(cnt(:,ib)), nnz(cnt(:,ib)));

Da = pool(M.map_diff,    cnt, ia);    Db = pool(M.map_diff,    cnt, ib);
Aa = pool(M.map_absdiff, cnt, ia);    Ab = pool(M.map_absdiff, cnt, ib);
Fa = pool(M.map_foot,    cnt, ia);

% scalar summary over the ROI, so the figure is not the only evidence
fprintf('\nmean |peak-foot| over the ROI:  %s %.3f   %s %.3f   ratio %.2fx\n', ...
        names{ia}, mean(Aa(:),'omitnan'), names{ib}, mean(Ab(:),'omitnan'), ...
        mean(Ab(:),'omitnan')/max(mean(Aa(:),'omitnan'),eps));
fprintf('peak signed displacement:       %s %+.3f / %+.3f   %s %+.3f / %+.3f  (min/max)\n', ...
        names{ia}, min(Da(:)), max(Da(:)), names{ib}, min(Db(:)), max(Db(:)));

%% ---------------------------- FIGURE ----------------------------------
f = figure('Color','w','Units','centimeters','Position',[2 2 30 15]);
set(f,'DefaultAxesFontSize',8);
cmDiv = diverging_map();

lim1 = symlim([Da(:); Db(:)]);        % shared, so normal vs gasp is comparable
lim2 = [0, max([Aa(:); Ab(:)], [], 'omitnan')];

panel(f,1, Da, cmDiv, lim1, sprintf('%s   mean(peak - foot)', names{ia}), sum(cnt(:,ia)));
panel(f,2, Db, cmDiv, lim1, sprintf('%s   mean(peak - foot)', names{ib}), sum(cnt(:,ib)));
panel(f,3, Db-Da, cmDiv, symlim(Db(:)-Da(:)), sprintf('%s - %s', names{ib}, names{ia}), NaN);
panel(f,4, M.u1_img, parula, [], 'u1\_img  (PC1 basis, reference)', NaN);

panel(f,5, Aa, parula, lim2, sprintf('%s   mean |peak - foot|', names{ia}), sum(cnt(:,ia)));
panel(f,6, Ab, parula, lim2, sprintf('%s   mean |peak - foot|', names{ib}), sum(cnt(:,ib)));
panel(f,7, Ab-Aa, cmDiv, symlim(Ab(:)-Aa(:)), sprintf('|%s| - |%s|', names{ib}, names{ia}), NaN);
panel(f,8, Fa, gray, [], sprintf('mean frame at foot (%s)', names{ia}), NaN);

sgtitle(sprintf(['chest motion by breath-cycle class   |   %s   |   ' ...
   'ROI [%d %d %d %d], bin %d'], strrep(M.session,'\','/'), M.roi_xywh(1), ...
   M.roi_xywh(2), M.roi_xywh(3), M.roi_xywh(4), M.bin_factor), 'FontSize',9);

if doSave
    exportgraphics(f, fullfile(outDir,'breath_class_maps.png'), 'Resolution',300, ...
                   'BackgroundColor','white');
    exportgraphics(f, fullfile(outDir,'breath_class_maps.pdf'), 'ContentType','vector', ...
                   'BackgroundColor','white');
    fprintf('\nSaved to %s\n', outDir);
end
fprintf('Done.\n');

%% ========================= LOCAL FUNCTIONS ==========================
function P = pool(A, cnt, ci)
%POOL  Count-weighted mean across recordings of A(:,:,run,class).
%  Weighted, not plain: recordings contributed 21 to 82 cycles of a class, and
%  an unweighted mean would let the thinnest recording speak as loudly as the
%  richest.
X = squeeze(A(:,:,:,ci));                 % [H W nRun]
w = reshape(cnt(:,ci), 1, 1, []);
w(isnan(X(1,1,:))) = 0;
num = sum(X .* w, 3, 'omitnan');
den = sum(w .* ~isnan(X), 3);
P = num ./ max(den, eps);
P(den == 0) = NaN;
end

function L = symlim(v)
%SYMLIM  Symmetric limits about 0, so a diverging map is honest about sign.
v = v(~isnan(v));
if isempty(v), L = [-1 1]; return; end
m = max(abs(v));
if m <= 0, m = 1; end
L = [-m m];
end

function panel(f, k, img, cm, lim, ttl, n)
ax = subplot(2,4,k,'Parent',f);
imagesc(ax, img);
axis(ax,'image'); set(ax,'XTick',[],'YTick',[]);
colormap(ax, cm);
if ~isempty(lim) && all(isfinite(lim)) && lim(2) > lim(1), caxis(ax, lim); end
colorbar(ax);
if isnan(n), title(ax, ttl, 'FontWeight','normal','Interpreter','tex');
else,        title(ax, sprintf('%s  (n=%d)', ttl, n), 'FontWeight','normal','Interpreter','tex');
end
end

function cm = diverging_map()
%DIVERGING_MAP  blue-white-red, for signed displacement.
n = 256; h = floor(n/2);
b = [linspace(0.15,1,h)', linspace(0.35,1,h)', linspace(0.75,1,h)'];
r = [linspace(1,0.75,h)', linspace(1,0.15,h)', linspace(1,0.15,h)'];
cm = [b; r];
end
