function roi_ncc_collect_260801()
%% roi_ncc_collect_260801  Pool NCC pairs and stitch shifts across every session.
% -----------------------------------------------------------------------
% Walks cfg.scanRoot for stitched sessions, and for each one collects:
%
%   NEAR pairs  cross-recording ROI pairs within cfg.near_xy_um / near_z_um
%               -> probably the same neuron seen twice
%   FAR  pairs  cross-recording ROI pairs out to cfg.far_xy_um / far_z_um
%               -> definitely different neurons, but same tissue and optics
%   shifts      the per-FOV stitch corrections, so their distribution across
%               sessions can be inspected
%
% Each session's own cfg is read back from its fov_table.mat, so sessions with
% different mount rotations are each oriented correctly -- a single global setting
% would silently mis-orient every session that disagrees with it.
%
% Output: <outDir>\ncc_pooled.mat

cfg = ncc_cfg_260801();
if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end

tabs = dir(fullfile(cfg.scanRoot,'**',cfg.stitchDir,'fov_table.mat'));
assert(~isempty(tabs), 'No %s\\fov_table.mat found under %s', cfg.stitchDir, cfg.scanRoot);

fprintf('=== roi_ncc_collect_260801 ===\nroot: %s\n%d stitched session(s)\n\n', ...
    cfg.scanRoot, numel(tabs));

NEAR = struct('ncc',{},'d',{},'dz',{},'sess',{},'fovA',{},'ridA',{},'fovB',{},'ridB',{});
FAR  = struct('ncc',{},'d',{},'dz',{},'sess',{},'fovA',{},'ridA',{},'fovB',{},'ridB',{});
SH   = struct('sess',{},'fov',{},'shift',{},'placed',{},'z',{});
Sess = strings(0,1);

for s = 1:numel(tabs)
    outDir = tabs(s).folder;
    sessName = sessionname(outDir, cfg);
    T = load(fullfile(outDir,'fov_table.mat'));
    F = T.F; scfg = T.cfg; nF = numel(F);

    % offsets from the stitch (absent = raw stage coordinates)
    tx = zeros(nF,1); ty = zeros(nF,1); placed = false(nF,1); hasOff = false;
    ofile = fullfile(outDir,'fov_offsets.mat');
    if isfile(ofile)
        O = load(ofile);
        if isfield(O,'fovTable')
            nm = string({F.name});
            for i = 1:nF
                j = find(O.fovTable.fov == nm(i), 1);
                if ~isempty(j)
                    tx(i)=O.fovTable.tx_um(j); ty(i)=O.fovTable.ty_um(j);
                    placed(i)=O.fovTable.placed(j);
                end
            end
            hasOff = true;
        end
    end

    for i = 1:nF
        SH(end+1) = struct('sess',sessName,'fov',string(F(i).name), ...
            'shift',hypot(tx(i),ty(i)),'placed',placed(i),'z',F(i).motor(3)); %#ok<AGROW>
    end

    % ---- flatten ROIs in the stitched frame ---------------------------------
    X=[];Y=[];Z=[];FV=[];RID=[];
    for i=1:nF
        n = numel(F(i).roi_lab);
        X=[X;F(i).roi_x_um+tx(i)]; Y=[Y;F(i).roi_y_um+ty(i)]; %#ok<AGROW>
        Z=[Z;repmat(F(i).motor(3),n,1)]; FV=[FV;repmat(i,n,1)]; %#ok<AGROW>
        RID=[RID;F(i).roi_lab(:)]; %#ok<AGROW>
    end
    if isempty(X)
        fprintf('  %-28s no ROIs, skipped\n', sessName); continue;
    end

    % ---- oriented average projections ---------------------------------------
    A = cell(nF,1);
    for i = 1:nF
        if isempty(F(i).avgPath) || ~isfile(F(i).avgPath), continue; end
        a = double(imread(F(i).avgPath));
        a = a - prctile(a(:),1); a = a / max(prctile(a(:),99.5),eps);
        A{i} = roi_stitch_orient_260801(min(max(a,0),1), scfg);
    end

    % ---- candidate pairs, cross-recording only ------------------------------
    nearI=[];nearJ=[]; farI=[];farJ=[];
    for p = 1:numel(X)-1
        q = find(FV(p+1:end) ~= FV(p)) + p;
        if isempty(q), continue; end
        d  = hypot(X(q)-X(p), Y(q)-Y(p));
        dz = abs(Z(q)-Z(p));
        kn = d <= cfg.near_xy_um & dz <= cfg.near_z_um;
        kf = d >  cfg.near_xy_um & d <= cfg.far_xy_um & dz <= cfg.far_z_um;
        nearI=[nearI;repmat(p,nnz(kn),1)]; nearJ=[nearJ;q(kn)]; %#ok<AGROW>
        farI =[farI; repmat(p,nnz(kf),1)]; farJ =[farJ; q(kf)]; %#ok<AGROW>
    end
    % cap FAR so one dense session cannot dominate the pooled null
    nFarAll = numel(farI);
    if nFarAll > cfg.max_far_per_session
        rng(0); k = randperm(nFarAll, cfg.max_far_per_session);
        farI = farI(k); farJ = farJ(k);
    end

    nearv = pair_ncc(A, F, tx, ty, X, Y, FV, nearI, nearJ, cfg, true);
    farv  = pair_ncc(A, F, tx, ty, X, Y, FV, farI,  farJ,  cfg, false);

    % Pair IDENTITY is stored, not just the score, so a pair above threshold can
    % actually be looked at rather than merely counted.
    for k = 1:numel(nearI)
        if ~isfinite(nearv(k)), continue; end
        NEAR(end+1) = mkpair(nearv(k), X,Y,Z, nearI(k), nearJ(k), FV, RID, F, sessName); %#ok<AGROW>
    end
    for k = 1:numel(farI)
        if ~isfinite(farv(k)), continue; end
        FAR(end+1) = mkpair(farv(k), X,Y,Z, farI(k), farJ(k), FV, RID, F, sessName); %#ok<AGROW>
    end
    Sess(end+1) = sessName; %#ok<AGROW>

    fprintf('  %-28s %2d FOV %4d ROI | near %4d | far %5d%s | shift med %.1f max %.1f um%s\n', ...
        sessName, nF, numel(X), nnz(isfinite(nearv)), nnz(isfinite(farv)), ...
        tern(nFarAll>cfg.max_far_per_session, sprintf(' (of %d)',nFarAll), ''), ...
        median(hypot(tx,ty)), max(hypot(tx,ty)), tern(hasOff,'',' [NOT STITCHED]'));
end

save(cfg.pooledFile, 'NEAR','FAR','SH','Sess','cfg');
fprintf('\npooled: %d near pairs, %d far pairs, %d sessions\n', ...
    numel(NEAR), numel(FAR), numel(Sess));
fprintf('Saved %s\n\nNext: roi_ncc_plot_260801\n', cfg.pooledFile);
end

%% ---------------------------------------------------------------------------
function p = mkpair(v, X, Y, Z, i, j, FV, RID, F, sessName)
p = struct('ncc',v, ...
    'd',hypot(X(i)-X(j), Y(i)-Y(j)), 'dz',abs(Z(i)-Z(j)), 'sess',sessName, ...
    'fovA',string(F(FV(i)).name), 'ridA',RID(i), ...
    'fovB',string(F(FV(j)).name), 'ridB',RID(j));
end

function v = pair_ncc(A, F, tx, ty, X, Y, FV, I, J, cfg, isNear)
v = nan(numel(I),1);
useMid = isNear && strcmpi(cfg.center_mode,'midpoint');
for k = 1:numel(I)
    p = I(k); q = J(k);
    if useMid
        xp=(X(p)+X(q))/2; yp=(Y(p)+Y(q))/2; xq=xp; yq=yp;
    else
        xp=X(p); yp=Y(p); xq=X(q); yq=Y(q);
    end
    P = patch_at(A{FV(p)}, F(FV(p)), tx(FV(p)), ty(FV(p)), xp, yp, cfg);
    Q = patch_at(A{FV(q)}, F(FV(q)), tx(FV(q)), ty(FV(q)), xq, yq, cfg);
    v(k) = ncc(P,Q);
end
end

function P = patch_at(A, f, tx, ty, xq, yq, cfg)
% roi_stitch_orient_260801 guarantees B(r,c) sits at X = Xmin+(c-0.5)*px,
% Y = Ymin+(r-0.5)*px, so inverting it is one subtraction per axis and stays
% correct if the frame convention ever changes.
P = [];
if isempty(A), return; end
[Ho,Wo] = size(A);
cx = (xq - (f.motor(1)+tx-f.halfW_um))/f.px_um + 0.5;
cy = (yq - (f.motor(2)+ty-f.halfH_um))/f.px_um + 0.5;
half = (cfg.patch_um/2)/f.px_um;
r0=round(cy-half); r1=round(cy+half); c0=round(cx-half); c1=round(cx+half);
if r0<1||c0<1||r1>Ho||c1>Wo||r1-r0<4||c1-c0<4, return; end
P = imresize(A(r0:r1,c0:c1), [cfg.patch_px cfg.patch_px]);
end

function r = ncc(P,Q)
r = NaN;
if isempty(P)||isempty(Q), return; end
a=P(:)-mean(P(:)); b=Q(:)-mean(Q(:));
na=norm(a); nb=norm(b);
if na<eps||nb<eps, return; end
r = (a'*b)/(na*nb);
end

function s = sessionname(outDir, cfg)
p = erase(string(outDir), string(cfg.scanRoot));
p = erase(p, filesep + string(cfg.stitchDir));
p = regexprep(char(p), ['^\' filesep], '');
s = string(strrep(p, filesep, '/'));
if strlength(s)==0, s = "root"; end
end

function v = tern(c,a,b), if c, v=a; else, v=b; end, end
