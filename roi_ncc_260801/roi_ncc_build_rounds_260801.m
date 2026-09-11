function roi_ncc_build_rounds_260801()
%% roi_ncc_build_rounds_260801  Seed the round 1-4 curation from the MAP STITCH.
% -----------------------------------------------------------------------
% Same four-round curation workflow as before (toss -> re-home -> merge ->
% place ungrouped), but the candidate groups are now gated by STITCHED POSITION
% instead of by an appearance search over an 80 um radius.
%
% WHAT CHANGED, AND WHY IT MATTERS
%   old: candidates = any ROI pair within 80 um / 30 um dz, decided by NCC
%   new: candidates = any ROI pair within  5 um / 15 um dz in the STITCHED frame,
%        NCC then scored on those
% The manual stitch holds repeat observations of one neuron to ~2 um, so a 5 um
% gate is far outside the registration error yet well inside a soma. It removes
% the failure the old radius allowed -- a neighbouring soma 30 um away that
% happens to look similar can no longer become a candidate at all.
%
% This writes roi_match_results.mat in EXACTLY the format the existing round 1-4
% GUIs read, so those tools are reused unchanged rather than reimplemented. Only
% the seeding differs.
%
% STANDALONE: writes into its own <dataset>\roi_ncc_out_260801, never into
% roi_match_out_* or roi_stitch_out_*. Existing curation is untouched.
%
% Requires: roi_stitch_260801 on the path, and the dataset already stitched.
%
% Then run, in order:
%   roi_ncc_round1_toss_260801
%   roi_ncc_round2_assign_260801
%   roi_ncc_round3_merge_260801
%   roi_ncc_round4_ungrouped_260801

%% ===================== USER-EDITABLE =====================
% Dataset comes from ncc_round_cfg_260801 so the builder and all four rounds are
% retargeted in one place and cannot drift apart.
rcfg          = ncc_round_cfg_260801();
datasetPath   = rcfg.datasetPath;
stitchDirName = rcfg.stitchDirName;
outDirName    = rcfg.outDirName;

maxPairDist_um = 5;     % stitched-frame lateral gate  <-- the new limit
maxPairDz_um   = 15;    % stitched-frame depth gate

% Patch geometry kept identical to the morph matcher, so the review GUIs'
% patch extraction and display look exactly as they always have.
patchHalf_um   = 35;
targetUmPerPx  = 0.5;
viewHalf_um    = 35;
imgClip        = [1 99.5];
dispClip       = [1 99.9];
dispGamma      = 0.6;
nccMergeThr    = 0.5;   % recorded for the GUIs; grouping here is POSITION-gated
%% =========================================================

outDir    = fullfile(datasetPath, outDirName);
stitchDir = fullfile(datasetPath, stitchDirName);
if ~isfolder(outDir), mkdir(outDir); end
tf = fullfile(stitchDir,'fov_table.mat');
assert(isfile(tf), ['No stitch found at %s\n' ...
    'Run roi_stitch_scan_260801 / roi_stitch_map_260801 on this dataset first.'], tf);

% Only F is needed: patches are cut from the RAW FOV image at native centroid
% pixels (what the review GUIs display), while the 5 um gate uses F's already
% stitch-framed roi_x_um/roi_y_um. No image reorientation is involved here.
S = load(tf); F = S.F; nF = numel(F);
fprintf('=== roi_ncc_build_rounds_260801 ===\n%s\n%d FOVs\n', datasetPath, nF);

% ---- stitch offsets ---------------------------------------------------------
tx = zeros(nF,1); ty = zeros(nF,1); nMoved = 0;
ofile = fullfile(stitchDir,'fov_offsets.mat');
if isfile(ofile)
    O = load(ofile); nm = string({F.name});
    for i = 1:nF
        j = find(O.fovTable.fov == nm(i), 1);
        if ~isempty(j) && O.fovTable.placed(j)
            tx(i)=O.fovTable.tx_um(j); ty(i)=O.fovTable.ty_um(j);
            if hypot(tx(i),ty(i))>0, nMoved = nMoved + 1; end
        end
    end
    fprintf('stitch offsets: %d FOV(s) moved, max %.1f um\n', nMoved, max(hypot(tx,ty)));
else
    fprintf(['!! no fov_offsets.mat -- using RAW stage coordinates.\n' ...
             '   A %g um gate on unstitched coordinates will miss most repeats.\n'], maxPairDist_um);
end

% ---- flatten ROIs, with the geometry the GUIs expect ------------------------
x=[];y=[];z=[];fovv=[];roii=[];labv=[];area=[];ecc=[];sol=[];cxpx=[];cypx=[];
for i = 1:nF
    Q = load(F(i).maskPath,'maskL');
    maskL = Q.maskL;
    pr = regionprops(maskL,'Centroid','Area','Eccentricity','Solidity');
    lab = F(i).roi_lab(:);
    n = numel(lab);
    x=[x; F(i).roi_x_um+tx(i)]; y=[y; F(i).roi_y_um+ty(i)]; %#ok<AGROW>
    z=[z; repmat(F(i).motor(3),n,1)]; %#ok<AGROW>
    fovv=[fovv; repmat(i,n,1)]; roii=[roii; (1:n)']; labv=[labv; lab]; %#ok<AGROW>
    area=[area; arrayfun(@(k) pr(k).Area, lab)]; %#ok<AGROW>
    ecc =[ecc;  arrayfun(@(k) pr(k).Eccentricity, lab)]; %#ok<AGROW>
    sol =[sol;  arrayfun(@(k) pr(k).Solidity, lab)]; %#ok<AGROW>
    cxpx=[cxpx; arrayfun(@(k) pr(k).Centroid(1), lab)]; %#ok<AGROW>
    cypx=[cypx; arrayfun(@(k) pr(k).Centroid(2), lab)]; %#ok<AGROW>
end
nR = numel(x);
fprintf('%d ROIs\n', nR);

% ---- patches, oriented in the stitch frame ----------------------------------
% Extracted here once and cached, so the GUIs load instantly and score against
% the same pixels this grouping used.
fprintf('extracting patches...\n');
half = round(patchHalf_um/targetUmPerPx);  nPix = (2*half+1)^2;
vhalf = round(viewHalf_um/targetUmPerPx);
pAll = cell(nR,1); mAll = cell(nR,1); Umat = nan(nR,nPix);
for i = 1:nF
    if isempty(F(i).avgPath) || ~isfile(F(i).avgPath), continue; end
    Araw = double(imread(F(i).avgPath));
    Q = load(F(i).maskPath,'maskL');
    sel = find(fovv == i);
    Ilin = norm_img(Araw, imgClip, 1);
    Idsp = norm_img(Araw, dispClip, dispGamma);
    for k = sel(:)'
        cx = cxpx(k); cy = cypx(k);
        pAll{k} = extract_patch(Idsp, cx, cy, viewHalf_um, F(i).px_um, vhalf);
        mAll{k} = extract_patch(double(Q.maskL==labv(k)), cx, cy, viewHalf_um, F(i).px_um, vhalf) > 0.5;
        t = extract_patch(Ilin, cx, cy, patchHalf_um, F(i).px_um, half);
        if ~isempty(t)
            v = t(:)-mean(t(:)); nv = norm(v);
            if nv>0, Umat(k,:) = (v/nv).'; end
        end
    end
end

% ---- candidate pairs: STITCHED POSITION gate --------------------------------
pi_=[];pj_=[];pn_=[];pd_=[];
for a = 1:nR-1
    q = find(fovv(a+1:end) ~= fovv(a)) + a;
    if isempty(q), continue; end
    d  = hypot(x(q)-x(a), y(q)-y(a));
    dz = abs(z(q)-z(a));
    ok = d <= maxPairDist_um & dz <= maxPairDz_um;
    if ~any(ok), continue; end
    qq = q(ok);
    r = nan(numel(qq),1);
    for k = 1:numel(qq)
        if all(isfinite(Umat(a,:))) && all(isfinite(Umat(qq(k),:)))
            r(k) = Umat(a,:) * Umat(qq(k),:).';
        end
    end
    pi_=[pi_;repmat(a,numel(qq),1)]; pj_=[pj_;qq(:)]; %#ok<AGROW>
    pn_=[pn_;r]; pd_=[pd_;d(ok)]; %#ok<AGROW>
end
fprintf('candidate pairs within %g um / %g um dz: %d\n', maxPairDist_um, maxPairDz_um, numel(pd_));

% ---- group: closest first, never two ROIs from the same recording -----------
% One neuron is seen at most once per recording, so a group holding two ROIs from
% one FOV is wrong by construction. That constraint is what stopped the old
% single-link chaining from walking across neighbouring somata.
[~,o] = sort(pd_);
parent = (1:nR)';
fovOf = arrayfun(@(i) fovv(i), (1:nR)', 'UniformOutput', false);
for k = o(:)'
    ra = uf(pi_(k)); rb = uf(pj_(k));
    if ra==rb, continue; end
    if ~isempty(intersect(fovOf{ra},fovOf{rb})), continue; end
    parent(rb)=ra; fovOf{ra}=[fovOf{ra};fovOf{rb}];
end
grp = zeros(nR,1); for i=1:nR, grp(i)=uf(i); end
[~,~,grp] = unique(grp);
nCells = max(grp);
groupSize = accumarray(grp,1);
multi = find(groupSize>1);
cellMembers = arrayfun(@(c) find(grp==c), (1:nCells)', 'UniformOutput', false);

memberNCC = cell(nCells,1); cellMinNCC = nan(nCells,1); cellMedNCC = nan(nCells,1);
for c = 1:nCells
    m = cellMembers{c};
    if numel(m)<2, memberNCC{c}=[]; continue; end
    v=[];
    for a1=1:numel(m), for b1=a1+1:numel(m)
        if all(isfinite(Umat(m(a1),:))) && all(isfinite(Umat(m(b1),:)))
            v(end+1) = Umat(m(a1),:)*Umat(m(b1),:).'; %#ok<AGROW>
        end
    end, end
    memberNCC{c}=v; if ~isempty(v), cellMinNCC(c)=min(v); cellMedNCC(c)=median(v); end
end
fprintf('%d ROIs -> %d cells (%d seen >1x, max %dx)\n', nR, nCells, numel(multi), max(groupSize));
fprintf('multi-obs cells: median min-NCC %.2f\n', median(cellMinNCC(multi),'omitnan'));

% ---- write in the format the review GUIs read -------------------------------
params = struct('datasetPath',datasetPath, ...
    'maxPairDist_um',maxPairDist_um,'maxPairDz_um',maxPairDz_um, ...
    'patchHalf_um',patchHalf_um,'targetUmPerPx',targetUmPerPx, ...
    'nccMergeThr',nccMergeThr,'tossFolders',{{}}, ...
    'viewHalf_um',viewHalf_um,'imgClip',imgClip, ...
    'dispClip',dispClip,'dispGamma',dispGamma, ...
    'seededFrom','roi_stitch_260801 (map stitch)','stitchDir',stitchDir);

roiT = struct('x',x,'y',y,'z',z,'fov',fovv,'roi',roii,'lab',labv, ...
              'area',area,'ecc',ecc,'sol',sol,'cx_px',cxpx,'cy_px',cypx);

match = struct('params',params, ...
    'fov_name',{string({F.name})'}, 'fov_folder',{string({F.folder})'}, ...
    'fov_pxum',arrayfun(@(f) f.px_um, F)', ...
    'roi',roiT, 'pairs',struct('i',pi_,'j',pj_,'ncc',pn_,'d',pd_), ...
    'grp',grp, 'groupSize',groupSize, 'multi',multi, ...
    'cellMembers',{cellMembers}, 'memberNCC',{memberNCC}, ...
    'cellMinNCC',cellMinNCC, 'cellMedNCC',cellMedNCC);
save(fullfile(outDir,'roi_match_results.mat'),'match','-v7.3');

patchCache = struct('pAll',{pAll},'mAll',{mAll},'Umat',Umat);
save(fullfile(outDir,'roi_match_patches.mat'),'patchCache','-v7.3');

fprintf('\nSaved %s\n      %s\n', fullfile(outDir,'roi_match_results.mat'), ...
    fullfile(outDir,'roi_match_patches.mat'));
fprintf('\nNext: roi_ncc_round1_toss_260801  (then 2, 3, 4)\n');

    function r = uf(a)
        r=a; while parent(r)~=r, r=parent(r); end
        while parent(a)~=a, na=parent(a); parent(a)=r; a=na; end
    end
end

%% ---------------------------------------------------------------------------
function I = norm_img(A, clip, gamma)
lo = prctile(A(:),clip(1)); hi = prctile(A(:),clip(2)); if hi<=lo, hi=lo+1; end
I = min(max((A-lo)/(hi-lo),0),1); if gamma~=1, I = I.^gamma; end
end

function P = extract_patch(A, cx, cy, half_um, px_um, halfOut)
% Cut a half_um box around (cx,cy) in native pixels and resample to the common
% target pixel size, so patches from different zooms are directly comparable.
P = [];
h = half_um/px_um;
r0 = round(cy-h); r1 = round(cy+h); c0 = round(cx-h); c1 = round(cx+h);
if r0<1||c0<1||r1>size(A,1)||c1>size(A,2)||r1-r0<3||c1-c0<3, return; end
P = imresize(A(r0:r1,c0:c1), [2*halfOut+1, 2*halfOut+1]);
end
