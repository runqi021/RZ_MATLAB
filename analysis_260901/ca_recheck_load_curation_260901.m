function CU = ca_recheck_load_curation_260901(bundleDir)
%CA_RECHECK_LOAD_CURATION_260901  Read the re-check curation for downstream figures.
%
%   CU = ca_recheck_load_curation_260901()
%   CU = ca_recheck_load_curation_260901(bundleDir)
%
%   THE ONE PLACE downstream analysis asks "was this recording of this cell
%   thrown away, and which events does it have now". Any script that builds the
%   per-cell dF/F heatmap, the pooled PETH or the polar plots should go through
%   here rather than reading ca_spike_data.mat, so that a toss decision made in
%   the GUI is honoured everywhere by construction.
%
%   CU.T          table, one row per (cell, recording): cell_idx, stem, group,
%                 folder, roi, n_events, toss, cell_toss, kept, source, note,
%                 reviewed. `toss` is that recording's own flag, `cell_toss` is
%                 set when the WHOLE cell was thrown away, and `kept` is the one
%                 to filter on: it is false if either is set.
%   CU.get(f,roi) -> struct with .found .cell_idx .stem .spike_idx .toss .note
%                 Looked up by (folder, ROI), which is what a figure script has.
%                 .found = false means this observation was never curated -- the
%                 caller should then fall back to the archive ca_spike_data.mat
%                 and treat the recording as KEPT, never as tossed.
%   CU.kept(f,roi)-> logical: false if the recording was tossed OR its cell was.
%   CU.nCells, CU.nTossed, CU.nTossedCells, CU.dir
%
%   Runqi Zhang / 2026-09-01

if nargin < 1 || isempty(bundleDir)
    bundleDir = ['D:\Ventral_surface_summary\per-cell-summary_active_260812\' ...
                 'spike_recheck_260901'];
end
curDir = fullfile(bundleDir,'curated');
assert(isfolder(curDir), 'No curated folder in %s', bundleDir);

d = dir(fullfile(curDir,'*_curated.mat'));
rows = cell(0,13);
M    = containers.Map('KeyType','char','ValueType','any');
for k = 1:numel(d)
    q = load(fullfile(d(k).folder, d(k).name),'CUR');
    if ~isfield(q,'CUR'), continue; end
    C = q.CUR;
    ct = false;
    if isfield(C,'cell_toss'), ct = logical(C.cell_toss); end
    for i = 1:numel(C.rec)
        R = C.rec(i);
        kp = ~logical(R.toss) && ~ct;
        rows(end+1,:) = {C.cell, string(C.stem), string(C.group), ...
            string(R.folder), R.roi, R.n_events, logical(R.toss), ct, kp, ...
            string(R.source), string(R.note), logical(C.reviewed), i}; %#ok<AGROW>
        key = local_key(R.folder, R.roi);
        M(key) = struct('found',true, 'cell_idx',C.cell, 'stem',C.stem, ...
                        'spike_idx',R.spike_idx(:), 'toss',logical(R.toss), ...
                        'cell_toss',ct, 'kept',kp, ...
                        'note',R.note, 'source',R.source);
    end
end

if isempty(rows)
    T = table('Size',[0 13], ...
        'VariableTypes',{'double','string','string','string','double','double', ...
                         'logical','logical','logical','string','string', ...
                         'logical','double'}, ...
        'VariableNames',{'cell_idx','stem','group','folder','roi','n_events', ...
                         'toss','cell_toss','kept','source','note','reviewed', ...
                         'rec_in_cell'});
else
    T = cell2table(rows, 'VariableNames', {'cell_idx','stem','group','folder', ...
        'roi','n_events','toss','cell_toss','kept','source','note','reviewed', ...
        'rec_in_cell'});
    T = sortrows(T,{'cell_idx','rec_in_cell'});
end

CU = struct();
CU.dir     = bundleDir;
CU.T       = T;
CU.nCells  = numel(unique(T.cell_idx));
CU.nTossed = nnz(~T.kept);
CU.nTossedCells = numel(unique(T.cell_idx(T.cell_toss)));
CU.get     = @(f,roi) local_get(M, f, roi);
CU.kept    = @(f,roi) local_get(M, f, roi).kept;
end

% -------------------------------------------------------------------------
function s = local_get(M, f, roi)
key = local_key(f, roi);
if isKey(M, key)
    s = M(key);
else
    % Never curated: the caller falls back to the archive and treats it as KEPT.
    s = struct('found',false, 'cell_idx',NaN, 'stem','', ...
               'spike_idx',zeros(0,1), 'toss',false, 'cell_toss',false, ...
               'kept',true, 'note','', 'source','');
end
end

function k = local_key(f, roi)
k = sprintf('%s|%d', lower(strrep(char(f),'/','\')), roi);
end
