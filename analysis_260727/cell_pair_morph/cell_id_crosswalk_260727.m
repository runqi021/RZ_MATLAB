function cell_id_crosswalk_260727(dirA, dirB)
%% cell_id_crosswalk_260727  Map cell ids between two analyses of the same animal.
% -----------------------------------------------------------------------
%   cell_id_crosswalk_260727(dirA, dirB)   two matcher output folders
%   cell_id_crosswalk_260727()             defaults: full phys set  <->  baseline
%
% WHY THIS IS NEEDED
% cell_id is POSITIONAL -- the index into the sorted group list -- so it renumbers
% whenever the ROI subset changes. After phys\ was split into baseline\ / roomAir\
% / reoxygenated\, baseline cell 87 and full-set cell 159 are the SAME neuron, and
% across the 447 shared observations the two numberings agreed only twice (0.4%).
% Figures named cell_087.png in one analysis and cell_159.png in the other are the
% same cell. This builds the map.
%
% HOW
% Cells are matched by SHARED MEMBERS: the (recording folder name, ROI index)
% pairs, which are invariant to any directory reorganisation. Two cells correspond
% if they share at least one member. Reported per pair:
%   n_shared     how many observations they have in common
%   n_only_A / n_only_B   members present in one analysis and not the other, which
%                is normal when one analysis is a subset of the other
%   relation     'exact'  same member set
%                'subset' A's members are all in B (or vice versa) -- the usual
%                         case when B is the full set and A a condition folder
%                'split'  they overlap but each has members the other lacks, which
%                         means the two curations DISAGREE about this cell and is
%                         worth looking at
%
% Output: cell_id_crosswalk.csv, written into dirA.
%
% Runqi Zhang / 2026-07-28

if nargin < 1 || isempty(dirA)
    dirA = 'D:\260721_Sert_soma_G8s\phys\baseline\roi_match_out_260727';
end
if nargin < 2 || isempty(dirB)
    dirB = 'D:\260721_Sert_soma_G8s\phys\roi_match_out_260724';
end

A = load_side(dirA);  B = load_side(dirB);
fprintf('\n=========== cell_id_crosswalk_260727 ===========\n');
fprintf('A: %s\n   %d observations, %d cells\n', dirA, numel(A.key), numel(A.cells));
fprintf('B: %s\n   %d observations, %d cells\n', dirB, numel(B.key), numel(B.cells));

% member key -> cell id, on each side
mapB = containers.Map(cellstr(B.key), num2cell(B.cid));

rows = {};
nExact = 0; nSubset = 0; nSplit = 0; nOrphan = 0;
for a = 1:numel(A.cells)
    ma = A.cells{a};
    ka = A.key(ma);
    hit = false(numel(ka),1);  bids = nan(numel(ka),1);
    for i = 1:numel(ka)
        if isKey(mapB, char(ka(i)))
            v = mapB(char(ka(i)));
            if ~isnan(v), hit(i) = true; bids(i) = v; end
        end
    end
    if ~any(hit)
        nOrphan = nOrphan + 1;
        rows(end+1,:) = {A.cid_of(a), A.keyOfCell(a), NaN, "", numel(ka), 0, numel(ka), NaN, "no-match"}; %#ok<AGROW>
        continue;
    end
    ub = unique(bids(hit));
    for b = ub(:)'
        shared  = nnz(bids == b);
        onlyA   = numel(ka) - shared;
        onlyB   = numel(B.cells{b}) - shared;
        if onlyA == 0 && onlyB == 0, rel = "exact";  nExact = nExact + 1;
        elseif onlyA == 0 || onlyB == 0, rel = "subset"; nSubset = nSubset + 1;
        else, rel = "split"; nSplit = nSplit + 1;
        end
        rows(end+1,:) = {A.cid_of(a), A.keyOfCell(a), b, B.keyOfCell(b), ...
                         numel(ka), shared, onlyA, onlyB, rel}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'A_cell_id','A_cell_key','B_cell_id','B_cell_key','A_size','n_shared', ...
     'n_only_A','n_only_B','relation'});
out = fullfile(dirA,'cell_id_crosswalk.csv');
writetable(T, out);

fprintf('\n  exact  %d   subset %d   split %d   no-match %d   (%d rows)\n', ...
        nExact, nSubset, nSplit, nOrphan, height(T));
if nSplit > 0
    fprintf(['  NOTE: %d "split" pairs overlap but each side has members the other\n' ...
             '  lacks -- the two curations genuinely disagree about those cells.\n'], nSplit);
end
k = find(T.A_cell_id == 87 & T.relation ~= "no-match", 1);
if ~isempty(k)
    fprintf('\n  example: A cell %d  <->  B cell %d   (%s, %d shared)\n', ...
            T.A_cell_id(k), T.B_cell_id(k), T.relation(k), T.n_shared(k));
end
fprintf('\nSaved %s\n', out);
end

%% ========================= helpers =========================
function S = load_side(d)
fr = fullfile(d,'roi_match_results.mat');
fc = fullfile(d,'roi_match_curated.mat');
assert(isfile(fr), 'roi_match_results.mat not found in %s', d);
M = load(fr,'match'); m = M.match;
n = numel(m.roi.fov);
S.key = strings(n,1);
for r = 1:n
    p = char(m.fov_folder(m.roi.fov(r)));
    while ~isempty(p) && (p(end)=='/'||p(end)=='\'), p(end)=[]; end
    [~,a,b] = fileparts(p);
    S.key(r) = string([a b]) + "|" + string(m.roi.roi(r));
end
if isfile(fc)
    C = load(fc,'curated');  g = C.curated.grpOf(:);
else
    g = m.grp(:);
    warning('No curation in %s -- using the automatic grouping.', d);
end
S.cells = cell(0,1);  S.cid = nan(n,1);
pos = sort(unique(g(g>0)));
for i = 1:numel(pos)
    mm = find(g==pos(i));  S.cells{end+1,1} = mm;  S.cid(mm) = numel(S.cells);
end
u = find(g==0);
for i = 1:numel(u), S.cells{end+1,1} = u(i);  S.cid(u(i)) = numel(S.cells); end
S.cid_of = @(c) c;
S.keyOfCell = @(c) anchor_key(S.key, S.cells{c});
end

function k = anchor_key(keys, members)
kk = keys(members);
parts = split(kk, "|");
if isvector(parts), parts = reshape(parts,1,[]); end
[~,ord] = sortrows([parts(:,1), parts(:,2)]);
k = replace(kk(ord(1)), "|", "#");
end
