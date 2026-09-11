function roi_ncc_plot_260801()
%% roi_ncc_plot_260801  Overlaid NEAR vs FAR NCC distributions + shift distribution.
% -----------------------------------------------------------------------
% Reads <outDir>\ncc_pooled.mat and produces:
%
%   1  NEAR vs FAR NCC, overlaid, pooled across sessions   <- the main figure
%   2  the same split by session, to check it is not one session carrying it
%   3  NCC vs separation, showing where the two populations cross
%   4  distribution of stitch shifts across all sessions
%
% The threshold is read off the FAR distribution at cfg.fpr_target, so it carries
% a stated false-positive rate rather than being chosen by eye. Separability is
% also reported as AUC, which needs no threshold at all.

cfg = ncc_cfg_260801();
assert(isfile(cfg.pooledFile), 'Run roi_ncc_collect_260801 first (%s missing).', cfg.pooledFile);
P = load(cfg.pooledFile);
NEAR = P.NEAR; FAR = P.FAR; SH = P.SH;

nv = [NEAR.ncc]'; fv = [FAR.ncc]';
nd = [NEAR.d]';   fd = [FAR.d]';
ns = string({NEAR.sess})'; fs = string({FAR.sess})';
assert(~isempty(nv) && ~isempty(fv), 'No pairs collected.');

thr = prctile(fv, 100*(1-cfg.fpr_target));
tpr = mean(nv >= thr);
auc = rocauc(nv, fv);
rt  = cfg.ncc_report_thr;

fprintf('=== roi_ncc_plot_260801 ===\n');
fprintf('NEAR (<=%g um xy, <=%g um z): n=%d  median %.2f\n', ...
    cfg.near_xy_um, cfg.near_z_um, numel(nv), median(nv));
fprintf('FAR  (%g-%g um xy, <=%g um z): n=%d  median %.2f\n', ...
    cfg.near_xy_um, cfg.far_xy_um, cfg.far_z_um, numel(fv), median(fv));
fprintf('AUC = %.3f  (0.5 = useless, 1.0 = perfect)\n\n', auc);

% Two thresholds, two questions. The fixed cut is what was asked for; the
% calibrated one says what false-positive rate that cut actually buys, which the
% fixed number alone cannot tell you.
fprintf('FIXED cut NCC > %.2f\n', rt);
fprintf('   NEAR above: %d / %d (%.0f%%)\n', nnz(nv>rt), numel(nv), 100*mean(nv>rt));
fprintf('   FAR  above: %d / %d (%.1f%%)  <-- false-positive rate at this cut\n', ...
    nnz(fv>rt), numel(fv), 100*mean(fv>rt));
fprintf('CALIBRATED cut at %.0f%% FPR = %.2f -> captures %.0f%% of near pairs\n\n', ...
    100*cfg.fpr_target, thr, 100*tpr);

% every pair above the fixed cut, both populations, most similar first
P2 = [NEAR(:); FAR(:)];
lab = [repmat("NEAR",numel(NEAR),1); repmat("FAR",numel(FAR),1)];
keep = [P2.ncc]' > rt;
Tp = table(lab(keep), string({P2(keep).sess})', [P2(keep).ncc]', ...
           [P2(keep).d]', [P2(keep).dz]', ...
           string({P2(keep).fovA})', [P2(keep).ridA]', ...
           string({P2(keep).fovB})', [P2(keep).ridB]', ...
    'VariableNames',{'population','session','ncc','d_um','dz_um', ...
                     'fovA','roiA','fovB','roiB'});
Tp = sortrows(Tp,'ncc','descend');
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
writetable(Tp, fullfile(cfg.outDir,'ncc_pairs_above_thr.csv'));
fprintf('%d pair(s) above %.2f written to ncc_pairs_above_thr.csv\n\n', height(Tp), rt);

f = figure('Color','w','Units','inches','Position',[.4 .4 17 9], ...
           'Name','roi_ncc_260801','NumberTitle','off');
t = tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');

% -- (1) pooled overlay --------------------------------------------------------
ax = nexttile(t,1); hold(ax,'on');
edges = -0.6:0.02:1;
histogram(ax, fv, edges, 'Normalization','probability', ...
    'FaceColor',[.45 .45 .45],'EdgeColor','none','FaceAlpha',.85);
histogram(ax, nv, edges, 'Normalization','probability', ...
    'FaceColor',[.85 .15 .15],'EdgeColor','none','FaceAlpha',.7);
xline(ax, rt, 'b-','LineWidth',1.8);
xline(ax, thr,'k--','LineWidth',1.4);
legend(ax, { sprintf('FAR  %g-%g um xy, <=%g um z  (n=%d)', cfg.near_xy_um, cfg.far_xy_um, cfg.far_z_um, numel(fv)), ...
             sprintf('NEAR <=%g um xy, <=%g um z  (n=%d)', cfg.near_xy_um, cfg.near_z_um, numel(nv)), ...
             sprintf('fixed %.2f: %.0f%% near, %.1f%% far', rt, 100*mean(nv>rt), 100*mean(fv>rt)), ...
             sprintf('calibrated %.2f (%.0f%% FPR)', thr, 100*cfg.fpr_target)}, ...
    'Location','northwest','Box','off');
xlabel(ax,'patch NCC'); ylabel(ax,'fraction of pairs');
title(ax, sprintf('pooled across %d sessions   AUC %.3f', numel(P.Sess), auc));
box(ax,'on'); xlim(ax,[-0.6 1]);

% -- (2) per session -----------------------------------------------------------
% One session with many pairs could carry a pooled result on its own; this is the
% check that the separation is a property of the data, not of one recording day.
ax = nexttile(t,2); hold(ax,'on');
u = unique([ns;fs]); yl = strings(numel(u),1);
for k = 1:numel(u)
    a = nv(ns==u(k)); b = fv(fs==u(k));
    if ~isempty(b), boxrow(ax, b, k-0.18, [.45 .45 .45]); end
    if ~isempty(a), boxrow(ax, a, k+0.18, [.85 .15 .15]); end
    yl(k) = sprintf('%s  (n=%d/%d)', u(k), numel(a), numel(b));
end
xline(ax, thr,'k--');
set(ax,'YTick',1:numel(u),'YTickLabel',yl,'FontSize',8,'YDir','reverse');
xlabel(ax,'patch NCC'); xlim(ax,[-0.6 1]); grid(ax,'on'); box(ax,'on');
title(ax,'per session   (grey = FAR, red = NEAR;  box = quartiles, line = median)');

% -- (3) NCC vs separation -----------------------------------------------------
ax = nexttile(t,3); hold(ax,'on');
scatter(ax, fd, fv, 6, [.55 .55 .55], 'filled', 'MarkerFaceAlpha',.25);
scatter(ax, nd, nv, 14, [.85 .15 .15], 'filled', 'MarkerFaceAlpha',.6);
[bc, bm] = binmedian([nd;fd], [nv;fv], 0:5:cfg.far_xy_um);
plot(ax, bc, bm, 'k-', 'LineWidth',1.8);
xline(ax, cfg.near_xy_um, 'b--','LineWidth',1.2);
yline(ax, thr,'k--');
xlabel(ax,'lateral separation between the two ROIs (um)'); ylabel(ax,'patch NCC');
title(ax,'appearance decays with separation; blue line = NEAR cut');
box(ax,'on'); grid(ax,'on'); xlim(ax,[0 cfg.far_xy_um]); ylim(ax,[-0.6 1]);

% -- (4) stitch shifts ---------------------------------------------------------
ax = nexttile(t,4); hold(ax,'on');
sv = [SH.shift]'; ss = string({SH.sess})'; sp = [SH.placed]';
u2 = unique(ss);
for k = 1:numel(u2)
    v = sv(ss==u2(k) & sp);
    if isempty(v), continue; end
    jitter = (rand(numel(v),1)-.5)*0.3;
    scatter(ax, v, k+jitter, 26, [.2 .45 .8],'filled','MarkerFaceAlpha',.75);
    plot(ax, median(v), k, 'k|','MarkerSize',16,'LineWidth',2);
end
set(ax,'YTick',1:numel(u2),'YTickLabel',u2,'FontSize',8,'YDir','reverse');
xlabel(ax,'stitch correction |d| (um), placed FOVs only');
title(ax, sprintf('shift distribution across sessions (median %.1f um, max %.1f um)', ...
    median(sv(sp)), max([sv(sp);0])));
grid(ax,'on'); box(ax,'on');

title(t, sprintf('Cross-session NCC separability -- %s', cfg.scanRoot), ...
    'Interpreter','none','FontWeight','bold');

if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
exportgraphics(f, fullfile(cfg.outDir,'ncc_separability.png'),'Resolution',200);
exportgraphics(f, fullfile(cfg.outDir,'ncc_separability.pdf'),'ContentType','vector');

T = table(cfg.near_xy_um, cfg.near_z_um, numel(nv), numel(fv), ...
    median(nv), median(fv), rt, mean(nv>rt), mean(fv>rt), thr, tpr, auc, ...
    'VariableNames',{'near_xy_um','near_z_um','n_near','n_far', ...
    'median_near','median_far','fixed_thr','frac_near_above','frac_far_above', ...
    'calibrated_thr','sensitivity','AUC'});
writetable(T, fullfile(cfg.outDir,'ncc_summary.csv'));
fprintf('Saved ncc_separability.png/.pdf + ncc_summary.csv in %s\n', cfg.outDir);
end

%% ---------------------------------------------------------------------------
function a = rocauc(pos, neg)
% Mann-Whitney form: the probability a random NEAR pair scores above a random FAR
% pair. Threshold-free, so it cannot be flattered by a lucky cut point.
x = [pos(:); neg(:)]; y = [true(numel(pos),1); false(numel(neg),1)];
[~,o] = sort(x); r = zeros(size(x)); r(o) = 1:numel(x);
a = (sum(r(y)) - numel(pos)*(numel(pos)+1)/2) / (numel(pos)*numel(neg));
end

function boxrow(ax, v, y, col)
q = prctile(v,[25 50 75]); w = 0.14;
patch(ax, [q(1) q(3) q(3) q(1)], y+[-w -w w w], col, 'EdgeColor','none','FaceAlpha',.65);
plot(ax, [q(2) q(2)], y+[-w w], 'k-','LineWidth',1.5);
lo = prctile(v,5); hi = prctile(v,95);
plot(ax, [lo hi], [y y], '-', 'Color',col*0.6, 'LineWidth',.8);
end

function [c, m] = binmedian(x, y, edges)
c = nan(numel(edges)-1,1); m = c;
for k = 1:numel(edges)-1
    s = x >= edges(k) & x < edges(k+1);
    if nnz(s) >= 5, c(k) = mean(edges(k:k+1)); m(k) = median(y(s)); end
end
ok = isfinite(c); c = c(ok); m = m(ok);
end
