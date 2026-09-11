function roi_stitch_match_260801()
%% roi_stitch_match_260801  Match ROIs by POSITION in the stitched frame.
% -----------------------------------------------------------------------
% Two observations of one neuron must land in the same place once the FOVs are
% stitched. So: stitch first, then match on distance. Appearance is computed and
% reported, but it never accepts or rejects a match.
%
% WHY THIS WAY ROUND. The appearance matcher scores a zero-shift NCC on patches
% cut at each ROI's centroid. That is fooled by segmentation disagreement: when
% Cellpose gives one neuron a different extent in two recordings the centroid
% moves, the patches misregister, and a true pair scores 0.57 while a false pair
% of similar-looking somata scores 0.45. Position does not have that failure
% mode -- measured on this rig, true pairs agree to ~2 um after registration
% while a false merge sat at 236 um.
%
% THREE GATES, ALL NECESSARY
%   lateral  <= cfg.maxMatchDist_um   the actual identity test
%   |dz|     <= cfg.maxDz_match_um    Z spans 330 um across 21 planes here, so a
%                                     lateral coincidence between distant planes
%                                     is two DIFFERENT cells stacked in depth
%   unambiguous                       if the runner-up candidate is nearly as
%                                     close as the winner, no match is made --
%                                     in a dense field the nearest neighbour is
%                                     not evidence of identity
%
% STANDALONE: reads only this pipeline's own fov_table.mat / fov_offsets.mat and
% writes only roi_stitch_results.mat. It does not touch roi_match_out_*.
%
% Input : <outDir>\fov_table.mat   (roi_stitch_scan_260801)
%         <outDir>\fov_offsets.mat (optional; zero offsets if absent)
% Output: <outDir>\roi_stitch_results.mat + roi_stitch_cells.csv

cfg = stitch_cfg_260801();
tf = fullfile(cfg.outDir,'fov_table.mat');
assert(isfile(tf), 'Run roi_stitch_scan_260801 first (%s missing).', tf);
S = load(tf); F = S.F;
nF = numel(F);

% ---- FOV offsets from the stitch --------------------------------------------
tx = zeros(nF,1); ty = zeros(nF,1); haveOff = false;
if isfile(cfg.offsetFile)
    O = load(cfg.offsetFile);
    if isfield(O,'fovTable')
        nm = string({F.name});
        for i = 1:nF
            j = find(O.fovTable.fov == nm(i), 1);
            if ~isempty(j), tx(i) = O.fovTable.tx_um(j); ty(i) = O.fovTable.ty_um(j); end
        end
        haveOff = true;
    end
end
fprintf('=== roi_stitch_match_260801 ===\n');
fprintf('dataset : %s\n', cfg.datasetPath);
if haveOff
    fprintf('offsets : %s (%d FOVs moved, max %.1f um)\n', cfg.offsetFile, ...
        nnz(hypot(tx,ty)>1e-9), max(hypot(tx,ty)));
else
    % plain fprintf, not fprintf(2,..): stderr is buffered separately, so a
    % warning printed there arrives out of order at the end of the log and reads
    % like a crash. This is an expected first-run state, not an error.
    fprintf('offsets : !! NONE FOUND -- matching on RAW STAGE COORDINATES.\n');
    fprintf('          Stitch first; otherwise every result here inherits the stage error.\n');
end

% ---- flatten to one ROI list ------------------------------------------------
X=[]; Y=[]; Z=[]; FV=[]; RID=[]; AR=[];
for i = 1:nF
    n = numel(F(i).roi_lab);
    X  = [X;  F(i).roi_x_um + tx(i)];  Y = [Y; F(i).roi_y_um + ty(i)]; %#ok<AGROW>
    Z  = [Z;  repmat(F(i).motor(3), n, 1)];   FV = [FV; repmat(i,n,1)]; %#ok<AGROW>
    RID= [RID; F(i).roi_lab];  AR = [AR; F(i).roi_area]; %#ok<AGROW>
end
nR = numel(X);
fprintf('%d FOVs, %d ROIs\n\n', nF, nR);

% ---- candidate pairs ---------------------------------------------------------
fprintf('gates: lateral <= %g um | |dz| <= %g um | unambiguous (runner-up > %.0f%% further)\n', ...
    cfg.maxMatchDist_um, cfg.maxDz_match_um, 100*(1/cfg.ambiguousRatio - 1));
pi_=[]; pj_=[]; pd_=[];
for a = 1:nR-1
    sel = find(FV(a+1:end) ~= FV(a)) + a;
    if isempty(sel), continue; end
    d  = hypot(X(sel)-X(a), Y(sel)-Y(a));
    dz = abs(Z(sel)-Z(a));
    ok = d <= cfg.maxMatchDist_um & dz <= cfg.maxDz_match_um;
    if ~any(ok), continue; end
    pi_ = [pi_; repmat(a,nnz(ok),1)]; pj_ = [pj_; sel(ok)]; pd_ = [pd_; d(ok)]; %#ok<AGROW>
end
fprintf('candidate pairs passing distance + dz: %d\n', numel(pd_));

% ---- ambiguity rejection ------------------------------------------------------
% A match survives only if, for BOTH partners, the next-best candidate in the
% other ROI's FOV is meaningfully further away. Without this, dense fields
% produce confident nearest-neighbour matches that are just the closest of many.
keep = true(size(pd_));
for k = 1:numel(pd_)
    for side = 1:2
        if side==1, me = pi_(k); other = pj_(k); else, me = pj_(k); other = pi_(k); end
        sib = find((pi_==me & FV(pj_)==FV(other)) | (pj_==me & FV(pi_)==FV(other)));
        dd  = sort(pd_(sib));
        if numel(dd) > 1 && dd(1) > cfg.ambiguousRatio*dd(2), keep(k) = false; end
    end
end
fprintf('rejected as ambiguous (runner-up nearly as close): %d\n', nnz(~keep));
pi_=pi_(keep); pj_=pj_(keep); pd_=pd_(keep);

% ---- greedy one-per-FOV matching, closest first --------------------------------
% Closest-first, and a merge is refused if the two groups already share a FOV.
% That constraint is what stops a chain forming: one neuron can be seen at most
% once per recording, so a group holding two ROIs from the same FOV is wrong by
% construction. (This is the failure the appearance matcher had -- single-link
% chaining walked from small round blobs to large elongated somata.)
[~,o] = sort(pd_); pi_=pi_(o); pj_=pj_(o); pd_=pd_(o);
parent = (1:nR)';
fovOf  = arrayfun(@(i) FV(i), (1:nR)', 'UniformOutput', false);
nMerge = 0;
for k = 1:numel(pd_)
    ra = uf_root(pi_(k)); rb = uf_root(pj_(k));
    if ra == rb, continue; end
    if ~isempty(intersect(fovOf{ra}, fovOf{rb})), continue; end   % one per FOV
    parent(rb) = ra;  fovOf{ra} = [fovOf{ra}; fovOf{rb}];  nMerge = nMerge + 1;
end
grp = zeros(nR,1);
for i = 1:nR, grp(i) = uf_root(i); end
[~,~,grp] = unique(grp);
nCells = max(grp);
sz = accumarray(grp,1);
fprintf('\n%d merges -> %d cells (%d seen >1x, max %dx)\n', nMerge, nCells, nnz(sz>1), max(sz));

% ---- appearance, REPORTED only -------------------------------------------------
fprintf('\ncomputing appearance NCC for accepted matches (reported, not used to gate)\n');
ncc = nan(nCells,1); spread = nan(nCells,1);
for c = find(sz>1)'
    m = find(grp==c);
    spread(c) = max(hypot(X(m)-mean(X(m)), Y(m)-mean(Y(m))));
    v = [];
    for a = 1:numel(m)
        for b = a+1:numel(m)
            v(end+1) = patch_ncc(F, FV(m(a)), RID(m(a)), FV(m(b)), RID(m(b))); %#ok<AGROW>
        end
    end
    if ~isempty(v), ncc(c) = min(v); end
end
lo = ncc < cfg.minNCCReport & sz > 1;
fprintf('  multi-obs cells: %d | median min-NCC %.2f | below %.2f (review): %d\n', ...
    nnz(sz>1), median(ncc(sz>1),'omitnan'), cfg.minNCCReport, nnz(lo));
fprintf('  positional spread within a cell: median %.1f um, max %.1f um\n', ...
    median(spread(sz>1),'omitnan'), max(spread));

% ---- save ---------------------------------------------------------------------
res = struct('cfg',cfg,'F',F,'tx',tx,'ty',ty,'haveOffsets',haveOff, ...
             'roi_fov',FV,'roi_id',RID,'roi_x_um',X,'roi_y_um',Y,'roi_z_um',Z, ...
             'roi_area',AR,'grp',grp,'cell_size',sz,'cell_minNCC',ncc, ...
             'cell_spread_um',spread);
save(cfg.resultsFile,'-struct','res');

T = table(grp, string({F(FV).name})', RID, X, Y, Z, sz(grp), ncc(grp), spread(grp), ...
    'VariableNames', {'cell_id','fov','roi_index','x_um','y_um','z_um', ...
                      'cell_size','cell_minNCC','cell_spread_um'});
T = sortrows(T, {'cell_size','cell_id'}, {'descend','ascend'});
writetable(T, fullfile(cfg.outDir,'roi_stitch_cells.csv'));
fprintf('\nSaved %s\n      %s\n', cfg.resultsFile, fullfile(cfg.outDir,'roi_stitch_cells.csv'));

    % nested so it shares `parent` (union-find with path compression)
    function r = uf_root(x)
        r = x; while parent(r) ~= r, r = parent(r); end
        while parent(x) ~= x, nx = parent(x); parent(x) = r; x = nx; end
    end
end

%% ---------------------------------------------------------------------------
function r = patch_ncc(F, fa, ra, fb, rb)
% Appearance similarity, for REPORTING. Patches are cut on the mask bounding box
% and resampled to a common grid, so different zooms are comparable.
persistent cache
if isempty(cache), cache = containers.Map('KeyType','char','ValueType','any'); end
pa = getpatch(F, fa, ra, cache);  pb = getpatch(F, fb, rb, cache);
if isempty(pa) || isempty(pb), r = NaN; return; end
a = pa(:) - mean(pa(:)); b = pb(:) - mean(pb(:));
r = (a'*b) / max(norm(a)*norm(b), eps);
end

function p = getpatch(F, fi, rid, cache)
k = sprintf('%d_%d', fi, rid);
if isKey(cache,k), p = cache(k); return; end
p = [];
if isempty(F(fi).avgPath) || ~isfile(F(fi).avgPath), cache(k)=p; return; end
persistent imcache
if isempty(imcache), imcache = containers.Map('KeyType','double','ValueType','any'); end
if isKey(imcache, fi), A = imcache(fi); else, A = double(imread(F(fi).avgPath)); imcache(fi)=A; end
j = find(F(fi).roi_lab == rid, 1); if isempty(j), cache(k)=p; return; end
% 40 um half-window, resampled to a fixed 41x41 grid -> zoom-independent
half = round(20 / F(fi).px_um);
cx = round(F(fi).roi_cx_px(j)); cy = round(F(fi).roi_cy_px(j));
r0 = max(1,cy-half); r1 = min(size(A,1),cy+half);
c0 = max(1,cx-half); c1 = min(size(A,2),cx+half);
if r1-r0 < 5 || c1-c0 < 5, cache(k)=p; return; end
p = imresize(A(r0:r1, c0:c1), [41 41]);
cache(k) = p;
end
