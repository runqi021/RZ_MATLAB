function roi_stitch_verify_260801()
%% roi_stitch_verify_260801  Confirm position-matched cells by appearance.
% -----------------------------------------------------------------------
% Step 4, run AFTER stitching and matching. A cell that is both co-located and
% looks the same in every recording is about as certain as this data allows.
%
% WHY APPEARANCE ONLY WORKS AT THIS POINT
% The appearance matcher failed because it cut each patch at its OWN ROI's
% centroid. When Cellpose gives one neuron a different extent in two recordings
% the centroid moves, the patches misregister, and a true pair scored 0.57
% against a false pair's 0.45 -- the metric was measuring segmentation
% disagreement, not identity.
%
% After stitching that is fixable: both patches are cut at the CELL'S CONSENSUS
% POSITION in the stitched frame, so they are registered to each other by
% construction and segmentation extent no longer enters. On Sert/0721 the stitch
% holds cells to 0.6 um median spread, far below a soma, so the consensus centre
% is meaningful.
%
% THE THRESHOLD IS MEASURED, NOT CHOSEN
% "NCC > 0.5" would be taste. Instead a NULL is built from patches of DIFFERENT
% cells in the same pair of FOVs -- exactly the confusion we care about, two
% distinct neurons that happen to look alike. The cut is placed at a percentile
% of that null, so it carries a false-positive rate rather than a vibe.
%
% Position remains the matcher. This does not create or delete matches; it labels
% each existing cell CONFIRMED or FLAGGED so you know which to look at.
%
% Input : <outDir>\fov_table.mat, fov_offsets.mat, roi_stitch_results.mat
% Output: <outDir>\roi_stitch_verified.mat, roi_stitch_verified.csv,
%         roi_stitch_verify_qc.png/.pdf

%% ===================== USER-EDITABLE =====================
PATCH_UM   = 40;      % full width of the comparison patch, um
PATCH_PX   = 41;      % common resample grid (zoom-independent)
NULL_N     = 4000;    % random different-cell pairs for the null
NULL_PCTL  = 99;      % confirm above this percentile of the null
NULL_MIN_UM= 15;      % null pairs must be at least this far apart (not the same cell)
SAVE_PDF   = true;
%% =========================================================

cfg = stitch_cfg_260801();
tf  = fullfile(cfg.outDir,'fov_table.mat');
assert(isfile(tf), 'Run roi_stitch_scan_260801 first.');
assert(isfile(cfg.resultsFile), 'Run roi_stitch_match_260801 first.');
S = load(tf); F = S.F; nF = numel(F);
R = load(cfg.resultsFile);

X = R.roi_x_um; Y = R.roi_y_um; FV = R.roi_fov; RID = R.roi_id;
grp = R.grp; sz = R.cell_size;
assert(numel(X) == sum(arrayfun(@(f) numel(f.roi_lab), F)), ...
    'results do not match fov_table -- re-run roi_stitch_match_260801');

fprintf('=== roi_stitch_verify_260801 ===\n%s\n', cfg.datasetPath);
fprintf('%d ROIs, %d cells, %d seen >1x\n\n', numel(X), max(grp), nnz(sz>1));

% ---- average projections, in the stitch frame -------------------------------
A = cell(nF,1);
for i = 1:nF
    if isempty(F(i).avgPath) || ~isfile(F(i).avgPath), continue; end
    a = double(imread(F(i).avgPath));
    a = a - prctile(a(:),1); a = a / max(prctile(a(:),99.5),eps);
    A{i} = roi_stitch_orient_260801(min(max(a,0),1), cfg);
end

% ---- patch at an arbitrary stitched position, from a given FOV ---------------
% Cutting at a POSITION rather than at an ROI centroid is the whole point.
getp = @(fi, xq, yq) patch_at(A{fi}, F(fi), R.tx(fi), R.ty(fi), xq, yq, PATCH_UM, PATCH_PX);

% ---- observed: all-pairs NCC within each multi-observation cell --------------
multi = find(sz > 1);
obsMin = nan(max(grp),1); obsMean = nan(max(grp),1); obsN = zeros(max(grp),1);
for c = multi(:)'
    m  = find(grp == c);
    xc = median(X(m)); yc = median(Y(m));      % consensus position
    P = cell(numel(m),1);
    for k = 1:numel(m), P{k} = getp(FV(m(k)), xc, yc); end
    v = [];
    for a1 = 1:numel(m)
        for b1 = a1+1:numel(m)
            r = ncc(P{a1}, P{b1});
            if isfinite(r), v(end+1) = r; end %#ok<AGROW>
        end
    end
    if ~isempty(v), obsMin(c) = min(v); obsMean(c) = mean(v); obsN(c) = numel(v); end
end
fprintf('observed: %d cells scored (%d pairwise comparisons)\n', ...
    nnz(isfinite(obsMin)), sum(obsN));

% ---- null: DIFFERENT cells, same pair of FOVs -------------------------------
% The question this answers is "how alike do two genuinely different neurons look
% when compared the same way", which is the error mode that matters. A null built
% by displacing within one cell would instead measure something else entirely.
rng(0);
nullv = nan(NULL_N,1); k = 0; tries = 0;
while k < NULL_N && tries < NULL_N*40
    tries = tries + 1;
    p = randi(numel(X)); q = randi(numel(X));
    if FV(p) == FV(q) || grp(p) == grp(q), continue; end
    if hypot(X(p)-X(q), Y(p)-Y(q)) < NULL_MIN_UM, continue; end
    r = ncc(getp(FV(p), X(p), Y(p)), getp(FV(q), X(q), Y(q)));
    if isfinite(r), k = k + 1; nullv(k) = r; end
end
nullv = nullv(1:k);
thr = prctile(nullv, NULL_PCTL);
fprintf('null    : %d different-cell pairs | median %.2f | p%g = %.2f  <-- threshold\n', ...
    k, median(nullv), NULL_PCTL, thr);

% ---- classify ---------------------------------------------------------------
verdict = strings(max(grp),1); verdict(:) = "single";
verdict(multi) = "FLAGGED";
conf = multi(obsMin(multi) >= thr);
verdict(conf) = "CONFIRMED";
nMulti = numel(multi); nConf = numel(conf);
fprintf('\nCONFIRMED %d / %d multi-observation cells (%.0f%%)  at a %.0f%% false-positive rate\n', ...
    nConf, nMulti, 100*nConf/max(nMulti,1), 100-NULL_PCTL);
fprintf('FLAGGED   %d  -- co-located but appearance disagrees; review these\n', nMulti-nConf);
fprintf('observed median min-NCC %.2f vs null median %.2f\n', ...
    median(obsMin(multi),'omitnan'), median(nullv));

% ---- save -------------------------------------------------------------------
V = struct('cfg',cfg,'thr',thr,'nullv',nullv,'nullPctl',NULL_PCTL, ...
           'obsMin',obsMin,'obsMean',obsMean,'obsN',obsN, ...
           'verdict',verdict,'grp',grp,'cell_size',sz, ...
           'roi_fov',FV,'roi_id',RID,'roi_x_um',X,'roi_y_um',Y);
save(fullfile(cfg.outDir,'roi_stitch_verified.mat'),'-struct','V');

T = table((1:max(grp))', sz, obsMin, obsMean, obsN, verdict, ...
    'VariableNames',{'cell_id','cell_size','minNCC','meanNCC','nPairs','verdict'});
T = T(sz>1,:); T = sortrows(T,'minNCC');
writetable(T, fullfile(cfg.outDir,'roi_stitch_verified.csv'));

%% ---- QC figure --------------------------------------------------------------
f = figure('Color','w','Units','inches','Position',[.5 .5 15 5], ...
           'Name','roi_stitch_verify_qc','NumberTitle','off');
t = tiledlayout(f,1,3,'Padding','compact','TileSpacing','compact');

ax = nexttile(t,1); hold(ax,'on');
edges = -0.4:0.02:1;
histogram(ax, nullv, edges, 'Normalization','probability', ...
    'FaceColor',[.6 .6 .6],'EdgeColor','none');
histogram(ax, obsMin(multi), edges, 'Normalization','probability', ...
    'FaceColor',[.85 .2 .2],'EdgeColor','none','FaceAlpha',.75);
xline(ax, thr, 'k--','LineWidth',1.5);
legend(ax, {'null: different cells', 'observed: matched cells', ...
    sprintf('threshold p%g = %.2f',NULL_PCTL,thr)}, 'Location','northwest','Box','off');
xlabel(ax,'patch NCC at consensus position'); ylabel(ax,'fraction');
title(ax,'appearance separates matched from unmatched'); box(ax,'on');

ax = nexttile(t,2); hold(ax,'on');
sc = scatter(ax, R.cell_spread_um(multi), obsMin(multi), 26, sz(multi), 'filled');
sc.MarkerEdgeColor = 'k'; sc.LineWidth = .2;
yline(ax, thr, 'k--');
colormap(ax, parula); cb = colorbar(ax); cb.Label.String = 'times seen';
xlabel(ax,'positional spread within cell (um)'); ylabel(ax,'min pairwise NCC');
title(ax,'position vs appearance agreement'); box(ax,'on'); grid(ax,'on');

ax = nexttile(t,3); axis(ax,'off');
L = { sprintf('dataset: %s', cfg.datasetPath), '', ...
      sprintf('cells                    %d', max(grp)), ...
      sprintf('seen >1x                 %d', nMulti), '', ...
      sprintf('null pairs               %d', k), ...
      sprintf('null median NCC          %.2f', median(nullv)), ...
      sprintf('threshold (p%g)          %.2f', NULL_PCTL, thr), '', ...
      sprintf('CONFIRMED                %d  (%.0f%%)', nConf, 100*nConf/max(nMulti,1)), ...
      sprintf('FLAGGED                  %d', nMulti-nConf), '', ...
      sprintf('observed median min-NCC  %.2f', median(obsMin(multi),'omitnan')), ...
      sprintf('positional spread median %.1f um', median(R.cell_spread_um(multi),'omitnan')), ...
      '', 'position matches; appearance only labels.', ...
      'no match is created or deleted here.'};
text(ax,0,1,L,'VerticalAlignment','top','FontName','Consolas','FontSize',9,'Interpreter','none');
title(ax,'summary');
title(t, sprintf('Appearance verification of position-matched cells -- %s', cfg.datasetPath), ...
    'Interpreter','none','FontWeight','bold');

exportgraphics(f, fullfile(cfg.outDir,'roi_stitch_verify_qc.png'),'Resolution',200);
if SAVE_PDF
    exportgraphics(f, fullfile(cfg.outDir,'roi_stitch_verify_qc.pdf'),'ContentType','vector');
end
fprintf('\nSaved roi_stitch_verified.mat/.csv + roi_stitch_verify_qc.png\n  in %s\n', cfg.outDir);
end

%% ---------------------------------------------------------------------------
function P = patch_at(A, f, tx, ty, xq, yq, patch_um, npx)
% Cut a patch centred on stitched position (xq,yq) out of FOV f's ORIENTED image.
%
% No axis inversion is needed. roi_stitch_orient_260801 guarantees its output
% satisfies  B(r,c) at  X = Xmin + (c-0.5)*px,  Y = Ymin + (r-0.5)*px  -- column
% runs along X and row along Y regardless of rig sign or mount rotation. Inverting
% that contract is one subtraction per axis, and it stays correct if the frame
% convention ever changes, which a hand-rolled inverse rotation would not.
P = [];
if isempty(A), return; end
[Ho, Wo] = size(A);
Xmin = f.motor(1) + tx - f.halfW_um;
Ymin = f.motor(2) + ty - f.halfH_um;
cx = (xq - Xmin)/f.px_um + 0.5;
cy = (yq - Ymin)/f.px_um + 0.5;

half = (patch_um/2) / f.px_um;
r0 = round(cy-half); r1 = round(cy+half);
c0 = round(cx-half); c1 = round(cx+half);
if r0 < 1 || c0 < 1 || r1 > Ho || c1 > Wo, return; end
if r1-r0 < 4 || c1-c0 < 4, return; end
P = imresize(A(r0:r1, c0:c1), [npx npx]);
end

function r = ncc(P, Q)
r = NaN;
if isempty(P) || isempty(Q), return; end
a = P(:) - mean(P(:)); b = Q(:) - mean(Q(:));
na = norm(a); nb = norm(b);
if na < eps || nb < eps, return; end
r = (a'*b) / (na*nb);
end
