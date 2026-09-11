function [obsOfCell, mergeInfo, origCells] = apply_cell_merges_260814(CELL, OBS, G)
%APPLY_CELL_MERGES_260814  Fold hand-curated merges into the cell -> obs mapping.
%
%   [obsOfCell, mergeInfo, origCells] = apply_cell_merges_260814(CELL, OBS, G)
%
% origCells : cell array parallel to obsOfCell, listing the ORIGINAL CELL indices
%             that went into each surviving cell. Callers should label a merged
%             cell with min(origCells{i}) so that every unmerged cell keeps the
%             number -- and therefore the output filename -- it already had.
%
% CELL, OBS : as loaded from event_latency_data.mat
% G         : groups from cell_merge_overrides_260814()
%
% obsOfCell : cell array, one entry per SURVIVING cell, holding its observation
%             indices. Merged cells appear once, with their observations
%             concatenated; every other cell is passed through unchanged.
% mergeInfo : struct array describing what was merged, for logging.
%
% The archive is never modified. Matching is on (recName, roi) parsed from
% OBS.label, because cell numbering is regenerated when the registry is rebuilt
% while the recording name and ROI index are stable.
%
% A group whose members all already belong to the SAME cell is reported as
% "already merged" and does nothing -- so re-running after cell_link is fixed
% upstream is harmless rather than double-counting.

nC = numel(CELL);
cellOf = zeros(numel(OBS),1);          % obs -> cell
for c = 1:nC
    cellOf(CELL(c).obs) = c;
end

% (recName, roi) -> obs index
key = strings(numel(OBS),1);
for o = 1:numel(OBS)
    p = regexp(OBS(o).label, '/', 'split');
    if numel(p) < 3, continue; end
    key(o) = string(strjoin(p(3:end-1),'/')) + "|" + string(p{end});
end

parent = 1:nC;                          % union-find over cells
    function r = findroot(a)
        while parent(a) ~= a, parent(a) = parent(parent(a)); a = parent(a); end
        r = a;
    end

mergeInfo = struct('group',{},'cells',{},'status',{});
for gi = 1:numel(G)
    grp = G{gi};
    cs = [];  miss = {};
    for m = 1:numel(grp)
        k = string(grp{m}{1}) + "|" + string(grp{m}{2});
        o = find(key == k, 1);
        if isempty(o), miss{end+1} = char(k); continue; end %#ok<AGROW>
        cs(end+1) = cellOf(o); %#ok<AGROW>
    end
    if ~isempty(miss)
        mergeInfo(end+1) = struct('group',gi,'cells',cs, ...
            'status',['NOT FOUND: ' strjoin(miss,', ')]); %#ok<AGROW>
        continue;
    end
    cs = unique(cs);
    if numel(cs) < 2
        mergeInfo(end+1) = struct('group',gi,'cells',cs,'status','already merged'); %#ok<AGROW>
        continue;
    end
    r = findroot(cs(1));
    for m = 2:numel(cs), parent(findroot(cs(m))) = r; end
    mergeInfo(end+1) = struct('group',gi,'cells',cs,'status','merged'); %#ok<AGROW>
end

roots = arrayfun(@findroot, 1:nC);
uR = unique(roots, 'stable');
obsOfCell = cell(numel(uR),1);
origCells = cell(numel(uR),1);
for i = 1:numel(uR)
    members = find(roots == uR(i));
    o = [];
    for m = members, o = [o; CELL(m).obs(:)]; end %#ok<AGROW>
    obsOfCell{i} = unique(o(:))';
    origCells{i} = members(:)';
end
end
