function Ventral_surface_tiepoint_register_260731()
%% Ventral_surface_tiepoint_register_260731  Refine FOV placement from matched cells.
% -----------------------------------------------------------------------
% A neuron matched across two FOVs has ONE true position, so the disagreement
% between its observations on the cartoon is pure placement error. The
% cross-FOV matcher has already solved the hard part -- WHICH ROIs are the same
% cell -- so those matches are tie-points and the correction is a small
% least-squares, not an image-registration problem.
%
% NO IMAGE STITCHING IS INVOLVED, deliberately. These FOVs span 1.2x-7x zoom and
% Z from +6 to -320 um, so cells visible in one plane are not the ones in
% another; pixel registration only works within a plane, while the matcher
% already handles the cross-Z case by appearance.
%
% WHAT IT SOLVES
%   one TRANSLATION per FOV. Not rotation or scale: with a handful of tie-points
%   those overfit, and they would fight the global stage->canvas affine, which is
%   a property of the vessel calibration and not something a few cells should be
%   allowed to renegotiate.
%
% GAUGE. The solution is only defined up to a constant per connected component
% (FOVs linked by shared cells), so within each component the MEAN correction is
% forced to zero. That fixes RELATIVE placement while leaving the component's
% absolute position exactly where the stage coordinates put it -- absolute
% anchoring belongs to the vessel calibration, not to this.
%
% ROBUSTNESS. A false merge would drag the fit, and we know they exist (one
% 236 um merge on Vglut2/1124). So the fit is iterated with tie-point rejection:
% any cell whose post-fit spread exceeds REJECT_UM is dropped and the fit redone.
%
% NON-DESTRUCTIVE. Writes x_um_reg / y_um_reg ALONGSIDE the originals into a new
% file. Nothing existing is modified, so the two maps can be compared before
% anything is committed.
%
% Input : <Desktop>\ROI_on_cartoon_data.mat   (written by Ventral_surface_ROI_on_cartoon_260729)
% Output: <Desktop>\ROI_on_cartoon_registered.mat  + a per-FOV shift report

%% ----------------------------- USER PARAMETERS -----------------------------
IN_MAT     = 'C:\Users\Admin\Desktop\ROI_on_cartoon_data.mat';
OUT_MAT    = 'C:\Users\Admin\Desktop\ROI_on_cartoon_registered.mat';
REJECT_UM  = 150;    % drop a tie-point cell whose members still disagree by more
                     %   than this after fitting (false merge, not a placement error)
MAX_ITER   = 500;
SAVE_FIG   = true;
%% ---------------------------------------------------------------------------

here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fileparts(here));

assert(isfile(IN_MAT), 'Not found: %s -- run Ventral_surface_ROI_on_cartoon_260729 first.', IN_MAT);
S = load(IN_MAT);
R = S.R;
X = [R.x_um]';  Y = [R.y_um]';
key  = string({R.cell_key})';
fovN = string({R.group})' + "/" + string({R.date})' + "/" + string({R.fov})';
sess = string({R.group})' + "/" + string({R.date})';

[uF, ~, fIdx] = unique(fovN);           % FOV index per ROI
nF = numel(uF);
fprintf('%d ROI observations across %d FOVs\n', numel(R), nF);

% ---- tie-points: cells seen in >1 FOV --------------------------------------
good = key ~= "" & ~ismissing(key);
[uK, ~, kIdx] = unique(key(good));
gi = find(good);
ties = {};
for c = 1:numel(uK)
    m = gi(kIdx == c);
    if numel(unique(fIdx(m))) < 2, continue; end
    ties{end+1,1} = m; %#ok<AGROW>
end
fprintf('tie-point cells (seen in >1 FOV): %d\n', numel(ties));

% ---- connected components of FOVs linked by shared cells --------------------
Adj = false(nF);
for t = 1:numel(ties)
    f = unique(fIdx(ties{t}));
    for a = 1:numel(f), for b = a+1:numel(f), Adj(f(a),f(b)) = true; Adj(f(b),f(a)) = true; end, end
end
comp = zeros(nF,1); nc = 0;
for f = 1:nF
    if comp(f), continue; end
    nc = nc + 1; stack = f;
    while ~isempty(stack)
        v = stack(end); stack(end) = [];
        if comp(v), continue; end
        comp(v) = nc;
        stack = [stack, find(Adj(v,:) & comp'==0)]; %#ok<AGROW>
    end
end
linked = accumarray(comp, 1) > 1;
fprintf('FOV components: %d  (%d FOVs are linked to at least one other, %d isolated)\n', ...
    nc, sum(linked(comp)), sum(~linked(comp)));

% ---- alternating least squares, with tie-point rejection --------------------
tx = zeros(nF,1); ty = zeros(nF,1);
active = true(numel(ties),1);
for pass = 1:5
    for it = 1:MAX_ITER
        ax = zeros(nF,1); ay = zeros(nF,1); cnt = zeros(nF,1);
        for t = find(active)'
            m = ties{t};
            cx = mean(X(m) + tx(fIdx(m)));  cy = mean(Y(m) + ty(fIdx(m)));
            for i = m'
                f = fIdx(i);
                ax(f) = ax(f) + (cx - (X(i)+tx(f)));
                ay(f) = ay(f) + (cy - (Y(i)+ty(f)));
                cnt(f) = cnt(f) + 1;
            end
        end
        k = cnt > 0;
        dx = zeros(nF,1); dy = zeros(nF,1);
        dx(k) = ax(k)./cnt(k);  dy(k) = ay(k)./cnt(k);
        tx = tx + dx;  ty = ty + dy;
        if max(hypot(dx,dy)) < 1e-4, break; end
    end
    % gauge: zero-mean correction inside each component (keep absolute position)
    for c = 1:nc
        f = comp == c;
        if ~any(f), continue; end
        tx(f) = tx(f) - mean(tx(f));  ty(f) = ty(f) - mean(ty(f));
    end
    % reject cells that still disagree -- these are match errors, not placement
    spread = nan(numel(ties),1);
    for t = 1:numel(ties)
        m = ties{t};
        xs = X(m)+tx(fIdx(m)); ys = Y(m)+ty(fIdx(m));
        spread(t) = max(hypot(xs-mean(xs), ys-mean(ys)));
    end
    bad = spread > REJECT_UM & active;
    if ~any(bad), break; end
    fprintf('  pass %d: rejected %d tie-point cell(s) over %g um\n', pass, sum(bad), REJECT_UM);
    active(bad) = false;
end
fprintf('tie-points used: %d of %d\n', sum(active), numel(ties));

% ---- before / after ---------------------------------------------------------
sp0 = nan(numel(ties),1); sp1 = nan(numel(ties),1);
for t = 1:numel(ties)
    m = ties{t};
    sp0(t) = max(hypot(X(m)-mean(X(m)), Y(m)-mean(Y(m))));
    xs = X(m)+tx(fIdx(m)); ys = Y(m)+ty(fIdx(m));
    sp1(t) = max(hypot(xs-mean(xs), ys-mean(ys)));
end
u = active;
fprintf('\n---- same-cell disagreement, tie-points used in the fit (n=%d) ----\n', sum(u));
fprintf('  %-8s %8s %8s %8s %8s\n','', 'median','75th','90th','max');
fprintf('  %-8s %8.0f %8.0f %8.0f %8.0f\n','before', median(sp0(u)), prctile(sp0(u),75), prctile(sp0(u),90), max(sp0(u)));
fprintf('  %-8s %8.0f %8.0f %8.0f %8.0f\n','after ', median(sp1(u)), prctile(sp1(u),75), prctile(sp1(u),90), max(sp1(u)));
if any(~u)
    fprintf('  rejected cells (NOT corrected, likely bad matches): n=%d, spread %.0f-%.0f um\n', ...
        sum(~u), min(sp0(~u)), max(sp0(~u)));
end

% ---- per-FOV shift ----------------------------------------------------------
shift = hypot(tx,ty);
fprintf('\n---- per-FOV correction ----\n');
fprintf('  linked FOVs: %d | median shift %.0f um | 90th %.0f | max %.0f\n', ...
    sum(linked(comp)), median(shift(linked(comp))), prctile(shift(linked(comp)),90), max(shift));
fprintf('  isolated FOVs (no tie-point, correction forced to 0): %d\n', sum(~linked(comp)));
[~,o] = sort(shift,'descend');
fprintf('\n  %-52s %9s %8s\n','FOV','shift_um','linked');
for i = o(1:min(12,nF))'
    fprintf('  %-52s %9.0f %8d\n', extractBefore(uF(i)+blanks(60),53), shift(i), linked(comp(i)));
end

% ---- write, non-destructively ----------------------------------------------
Xr = X + tx(fIdx);  Yr = Y + ty(fIdx);
for i = 1:numel(R)
    R(i).x_um_reg = Xr(i);
    R(i).y_um_reg = Yr(i);
    R(i).reg_shift_um = shift(fIdx(i));
    R(i).reg_linked   = linked(comp(fIdx(i)));
end
fovTable = table(uF, tx, ty, shift, linked(comp), comp, ...
    'VariableNames', {'fov','tx_um','ty_um','shift_um','linked','component'});
save(OUT_MAT, 'R', 'fovTable', 'sp0', 'sp1', 'active', 'REJECT_UM');
writetable(fovTable, strrep(OUT_MAT,'.mat','_fov_shifts.csv'));
fprintf('\nSaved %s\n  + _fov_shifts.csv\n', OUT_MAT);

if SAVE_FIG
    f = figure('Color','w','Position',[80 80 1100 440]);
    tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
    nexttile; hold on;
    plot(sp0(u), sp1(u), 'o', 'MarkerSize',4, 'Color',[.2 .4 .8], 'MarkerFaceColor',[.2 .4 .8]);
    if any(~u), plot(sp0(~u), sp1(~u), 'o','MarkerSize',5,'Color',[.85 .2 .2],'LineWidth',1); end
    L = max([sp0;sp1])*1.05; plot([0 L],[0 L],'k--');
    axis equal; xlim([0 L]); ylim([0 L]); grid on; box on;
    xlabel('same-cell disagreement BEFORE (\mum)'); ylabel('AFTER (\mum)');
    title(sprintf('below the line = improved  (%d used, %d rejected)', sum(u), sum(~u)), ...
        'FontWeight','normal');
    nexttile; hold on;
    histogram(shift(linked(comp)), 20, 'FaceColor',[.2 .4 .8], 'EdgeColor','none');
    grid on; box on; xlabel('per-FOV correction (\mum)'); ylabel('FOVs');
    title(sprintf('median %.0f \\mum, max %.0f \\mum', ...
        median(shift(linked(comp))), max(shift)), 'FontWeight','normal');
    title(tl,'tie-point registration from matched cells','FontWeight','bold');
    exportgraphics(f, strrep(OUT_MAT,'.mat','.png'), 'Resolution',180);
    fprintf('  + %s\n', strrep(OUT_MAT,'.mat','.png'));
end
end
