function roi_curation_port_260727(sourceDir)
%% roi_curation_port_260727  Carry a hand curation across a folder reorganisation.
% -----------------------------------------------------------------------
%   roi_curation_port_260727()            % port from the default SOURCE below
%   roi_curation_port_260727(sourceDir)   % port from any other matcher output dir
%
% THE INVARIANT THIS RELIES ON
%   A curation decision is really a statement about
%       (recording folder NAME, ROI index)  ->  which cell
%   The ROI index comes from that recording's own cpSAM maskL, and the folder name
%   travels with the folder. NEITHER depends on where the folder sits in the
%   directory tree. So as long as folder names are unchanged, a curation survives
%   ANY amount of moving recordings between parent folders -- splitting one phys\
%   into baseline\ / roomAir\ / reoxygenated\, regrouping by condition, anything.
%   Re-curating after a reorganisation is never necessary.
%
%   What the curation does NOT survive is re-running the SEGMENTATION, because
%   that renumbers the ROIs. This script checks for exactly that (see below).
%
% WHAT IT DOES
%   1. Loads SOURCE roi_match_results.mat + roi_match_curated.mat -> the old
%      (recName, roiIdx) -> grpOf map.
%   2. Loads TARGET roi_match_results.mat (the matcher re-run on the new folder).
%   3. Matches on (recName, roiIdx) and carries grpOf and homeOf across.
%   4. INTEGRITY CHECK: for every matched pair the ROI centroid (cx_px, cy_px) must
%      agree to within centTol_px. If it does not, the ROI numbering means
%      something different on the two sides -- the segmentation was re-run -- and
%      porting would silently assign decisions to the wrong neurons. The script
%      refuses rather than producing a plausible-looking wrong answer.
%   5. Relabels group ids compactly over the surviving subset.
%   6. Writes TARGET roi_match_curated.mat + .csv, in the exact format the review
%      GUIs and cell_link_260727 expect, so you can carry on curating from there.
%
% WHAT HAPPENS TO PARTIAL GROUPS
%   A cell whose members are split across condition folders keeps only the members
%   present in the target. A group reduced to ONE surviving member stays grpOf > 0
%   rather than reverting to 0: it was a decision you already made, and reverting
%   would push it back into round 4's triage queue to be re-answered. Downstream it
%   is identical either way, since a 1-member group and an ungrouped ROI both
%   become a 1x cell.
%
% Target comes from cell_cfg_260727. Source defaults to the pre-split full-set
% curation. Nothing in SOURCE is ever written to.
%
% Runqi Zhang / 2026-07-27

%% ---- params ----
repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot); addpath(fileparts(mfilename('fullpath')));

if nargin < 1 || isempty(sourceDir)
    sourceDir = 'D:\260721_Sert_soma_G8s\phys\roi_match_out_260724';   % <-- the pre-split curation
end
centTol_px = 1e-6;      % ROI centroids must be identical, not merely close
doSave     = true;

cfg = cell_cfg_260727();
fprintf('\n============ roi_curation_port_260727 ============\n');
fprintf('SOURCE (read only): %s\n', sourceDir);
fprintf('TARGET            : %s\n', cfg.outDir);

%% ---- load ----
srcRes = fullfile(sourceDir,'roi_match_results.mat');
srcCur = fullfile(sourceDir,'roi_match_curated.mat');
assert(isfile(srcRes), 'SOURCE roi_match_results.mat not found: %s', srcRes);
assert(isfile(srcCur), 'SOURCE roi_match_curated.mat not found: %s', srcCur);
assert(isfile(cfg.resultsFile), ['TARGET roi_match_results.mat not found:\n  %s\n' ...
    'Run roi_pair_morph_match_260727.m on the new folder first.'], cfg.resultsFile);

S  = load(srcRes,'match');  sm = S.match;
SC = load(srcCur,'curated'); sc = SC.curated;
T  = load(cfg.resultsFile,'match'); tm = T.match;

sGrp = sc.grpOf(:);
if isfield(sc,'homeOf') && numel(sc.homeOf)==numel(sGrp), sHome = sc.homeOf(:); else, sHome = sGrp; end
assert(numel(sGrp) == numel(sm.roi.fov), ...
    'SOURCE curation (%d ROIs) does not match SOURCE results (%d ROIs).', numel(sGrp), numel(sm.roi.fov));

%% ---- build the (recName, roiIdx) key on both sides ----
sKey = make_key(sm);
tKey = make_key(tm);
nT   = numel(tKey);
fprintf('\nSOURCE: %d observations across %d recordings\n', numel(sKey), numel(sm.fov_folder));
fprintf('TARGET: %d observations across %d recordings\n', nT, numel(tm.fov_folder));

[tf, loc] = ismember(tKey, sKey);
fprintf('matched by (folder name, ROI index): %d of %d target observations\n', nnz(tf), nT);

%% ---- INTEGRITY: the ROI numbering must mean the same thing on both sides ----
mi = find(tf);
dc = abs(tm.roi.cx_px(mi) - sm.roi.cx_px(loc(mi)));
dr = abs(tm.roi.cy_px(mi) - sm.roi.cy_px(loc(mi)));
bad = dc > centTol_px | dr > centTol_px;
if any(bad)
    b = mi(bad);
    fprintf(2,'\n*** CENTROID MISMATCH on %d of %d matched observations ***\n', nnz(bad), numel(mi));
    for q = b(1:min(8,numel(b)))'
        fprintf(2,'   %s  roi %d :  target (%.3f, %.3f)  vs  source (%.3f, %.3f)\n', ...
            tKey(q), tm.roi.roi(q), tm.roi.cx_px(q), tm.roi.cy_px(q), ...
            sm.roi.cx_px(loc(q)), sm.roi.cy_px(loc(q)));
    end
    error(['ROI numbering differs between SOURCE and TARGET for the same recording, ' ...
           'which means the segmentation was re-run. Porting would assign your ' ...
           'decisions to the wrong neurons. Aborting -- re-curate instead.']);
end
fprintf('integrity OK: every matched ROI has an identical centroid on both sides\n');

%% ---- carry the decisions across ----
grpRaw  = zeros(nT,1);  homeRaw = zeros(nT,1);
grpRaw(mi)  = sGrp(loc(mi));
homeRaw(mi) = sHome(loc(mi));

% relabel group ids compactly over what survives
posOld = sort(unique(grpRaw(grpRaw>0)));
grpOf  = zeros(nT,1);  homeOf = zeros(nT,1);
grpOf(grpRaw < 0) = -1;                       % tossed ports directly
sizes  = zeros(numel(posOld),1);
for g = 1:numel(posOld)
    sel = grpRaw == posOld(g);
    sizes(g) = nnz(sel);
    grpOf(sel) = g;
end
for g = 1:numel(posOld)
    homeOf(homeRaw == posOld(g)) = g;
end
homeOf(grpOf > 0) = grpOf(grpOf > 0);
homeOf(homeOf == 0 & grpOf == 0) = 0;

%% ---- INHERIT THE SOURCE CELL NUMBERING ------------------------------------
% cell_id is otherwise POSITIONAL -- the rank of a group in this subset -- so it
% renumbers whenever the ROI subset changes. Splitting phys\ into condition folders
% made baseline cell 87 the same neuron as full-set cell 159, and across the 447
% shared observations the two numberings agreed only twice.
%
% The port already knows which SOURCE cell every target observation came from, so
% it can simply carry that id across. Downstream then numbers cells by the SOURCE
% (full-set) ids, leaving gaps where a cell has no member in this subset. The
% numbering is sparse but it is the SAME in baseline, roomAir and reoxygenated, so
% cell 159 means one neuron everywhere.
srcCellId = nan(nT,1);
sPos = sort(unique(sGrp(sGrp>0)));
sCid = nan(numel(sGrp),1);  k = 0;
for i = 1:numel(sPos), k = k+1; sCid(sGrp==sPos(i)) = k; end
sUng = find(sGrp==0);
for i = 1:numel(sUng), k = k+1; sCid(sUng(i)) = k; end   % exactly cell_link's rule
srcCellId(mi) = sCid(loc(mi));
srcCellId(grpOf < 0) = NaN;                               % tossed have no cell
fprintf('inherited SOURCE cell ids for %d of %d observations (max id %d)\n', ...
        nnz(~isnan(srcCellId)), nT, max(srcCellId));

%% ---- audit ----
srcSizes = arrayfun(@(i) nnz(sGrp==i), sort(unique(sGrp(sGrp>0))));
fprintf('\n---- what came across ----\n');
fprintf('  ported   : %d observations (%d grouped | %d ungrouped | %d tossed)\n', ...
        nnz(tf), nnz(grpOf>0), nnz(grpOf==0 & tf), nnz(grpOf<0));
fprintf('  NOT in the source (fresh, never reviewed): %d\n', nnz(~tf));
if any(~tf)
    u = unique(regexprep(tKey(~tf), '\|\d+$',''));
    for q = 1:numel(u), fprintf('        %s (%d ROIs)\n', u(q), nnz(strncmp(tKey(~tf), u(q), strlength(u(q))))); end
end
fprintf('  cells    : %d in source (%d multi) -> %d here (%d multi)\n', ...
        numel(srcSizes), nnz(srcSizes>1), numel(posOld), nnz(sizes>1));
fprintf('  groups reduced to a single surviving member: %d (kept as settled 1x cells)\n', nnz(sizes==1));
fprintf('  surviving size distribution: ');
for s = 1:max([sizes;0]), if nnz(sizes==s), fprintf('%dx:%d  ', s, nnz(sizes==s)); end, end
fprintf('\n');

%% ---- save in the format the review GUIs + cell_link expect ----
if doSave
    if ~isfolder(cfg.outDir), mkdir(cfg.outDir); end
    roi = tm.roi;  fov_name = tm.fov_name;
    posIds = sort(unique(grpOf(grpOf>0)));
    groups = cell(0,1);
    for g = 1:numel(posIds)
        m = find(grpOf==posIds(g)); [~,o] = sort(roi.fov(m)); groups{end+1,1} = m(o); %#ok<AGROW>
    end
    ungr = find(grpOf==0);
    for k = 1:numel(ungr), groups{end+1,1} = ungr(k); end %#ok<AGROW>

    curated = struct('grpOf',grpOf,'homeOf',homeOf, ...
        'curatedCellMembers',{groups},'original',{tm.cellMembers}, ...
        'roi',roi,'fov_name',fov_name, ...
        'srcCellId',srcCellId, ...   % <- inherited full-set numbering; cell_link uses it
        'ported_from',string(sourceDir),'ported_n',nnz(tf)); %#ok<NASGU>
    save(fullfile(cfg.outDir,'roi_match_curated.mat'),'curated','-v7.3');

    % Coordinates travel with the curation. Two frames, both needed downstream:
    %   x/y/z_um   stage coordinates -- comparable ACROSS recordings, so this is
    %              what a spatial map of the whole preparation is drawn in
    %   cx/cy_px   pixel centroid WITHIN its own FOV -- what you need to draw an
    %              ROI back onto that recording's own image
    rows = {};
    for g = 1:numel(groups)
        m = groups{g};
        for c = 1:numel(m)
            rr = m(c);
            rows(end+1,:) = {g, numel(m), char(fov_name(roi.fov(rr))), roi.roi(rr), ...
                             roi.x(rr), roi.y(rr), roi.z(rr), ...
                             roi.cx_px(rr), roi.cy_px(rr), ...
                             mean(roi.x(m)), mean(roi.y(m)), mean(roi.z(m))}; %#ok<AGROW>
        end
    end
    if ~isempty(rows)
        Tt = cell2table(rows,'VariableNames', ...
            {'cell_id','cell_size','fov_name','roi_index','x_um','y_um','z_um', ...
             'cx_px','cy_px','cell_x_um','cell_y_um','cell_z_um'});
        writetable(Tt, fullfile(cfg.outDir,'roi_match_curated.csv'));
    end
    % the round-4 reviewed marks do not transfer: the ROI set changed
    rf = fullfile(cfg.outDir,'roi_match_r4_reviewed.mat');
    if isfile(rf), delete(rf); end
    fprintf('\nSaved ported curation to\n  %s\n', fullfile(cfg.outDir,'roi_match_curated.mat'));
    fprintf('The review GUIs will now resume from it. SOURCE was not modified.\n');
end
end

%% ========================= helpers =========================
function k = make_key(m)
% "<recording folder name>|<roi index>" -- stable under any directory move.
names = strings(numel(m.fov_folder),1);
for f = 1:numel(m.fov_folder)
    p = char(m.fov_folder(f));
    while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end) = []; end
    [~,n,e] = fileparts(p);
    names(f) = string([n e]);
end
k = names(m.roi.fov(:)) + "|" + string(m.roi.roi(:));
end
