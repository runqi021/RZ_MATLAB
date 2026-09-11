function cell_link_260727()
%% cell_link_260727  Join the cross-FOV cell curation to the coherence recordings.
% -----------------------------------------------------------------------
% STEP 1 of the cell-aware coherence pipeline. This script computes NOTHING about
% breathing or calcium. Its only job is to build -- and then audit -- the map
%
%     matcher ROI observation  <->  (recording folder, ROI index)  <->  cell id
%
% and to prove that the map is lossless: every ROI the matcher knows about is
% accounted for exactly once, every recording is accounted for, and the ROI
% indices really do address the same neurons in both worlds.
%
% WHY A SEPARATE SCRIPT: the join is the only place where the two pipelines can
% silently disagree (an ROI index that means something different on each side, a
% recording present in one and missing from the other). Isolating it means the
% audit is inspectable BEFORE any coherence number is computed on top of it.
%
% -------------------------- the join key --------------------------------
%   cell side (roi_match_results/_curated.mat)
%       roi.fov(r)   -> fov_folder(f)  -> folder name = the recording name
%       roi.roi(r)   -> within-FOV ROI index = column of F/dFF = index into roi_spikes
%   coherence side (coherence_polar_general_260727.m)
%       labels{i} = 'group/date/recName/roiIdx'
%   so (recName, roiIdx) identifies the same observation on both sides.
%
% ------------------------ what "cell id" means --------------------------
% Identical to the numbering in roi_match_curated.csv, so a cell id here is the
% same cell id you saw while curating:
%       grpOf > 0  -> the curated multi-recording cells, in ascending grpOf order
%       grpOf == 0 -> ungrouped-untossed, each its OWN 1x cell, appended after
%       grpOf < 0  -> TOSSED: given cell_id = NaN, kept in the table, excluded
%                     from the cells list. Never silently deleted.
%
% Output (into cfg.cellDir):
%   cell_link.mat  - obs table + cells + rec table + audit struct
%   cell_link.csv  - the obs table, one row per ROI observation
%
% Requires: the matcher + your curation (analysis_260727\cell_pair_morph\).
% Runqi Zhang / 2026-07-27

%% ---- path setup ----
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(fileparts(scriptDir));
addpath(repoRoot); addpath(scriptDir);

cfg = coh_cfg_260727();
if ~isfolder(cfg.cellDir), mkdir(cfg.cellDir); end

fallback_fps = 30;

fprintf('\n=============== cell_link_260727 ===============\n');
fprintf('experiment : %s\n', cfg.rootPath);
fprintf('matcher out: %s\n', cfg.cellOutDir);

%% ---- 0. no cell matching for this dataset? fall back to identity grouping ----
% Cross-FOV matching is optional. Without it, every (recording, ROI) is simply its
% own cell -- which is what any single-recording analysis assumes anyway. This
% keeps the whole downstream chain (cell_pool -> breath_time) runnable on a
% dataset that has never been through the matcher, instead of hard-failing.
if ~isfile(cfg.resultsFile)
    fprintf(['\nNo roi_match_results.mat under\n  %s\n' ...
             'Falling back to IDENTITY grouping: each ROI in each recording is its own cell.\n' ...
             'Cross-recording pooling is therefore OFF for this dataset.\n'], cfg.cellOutDir);
    identity_link(cfg, fallback_fps);
    return;
end
M = load(cfg.resultsFile, 'match');  match = M.match;
roi = match.roi;  fov_name = match.fov_name;  fov_folder = match.fov_folder;
nObs = numel(roi.fov);
nFov = numel(fov_folder);

%% ---- 2. grouping: curated if present, else the matcher's auto grouping ----
if isfile(cfg.curatedFile)
    C = load(cfg.curatedFile, 'curated');
    grpOf = C.curated.grpOf(:);
    assert(numel(grpOf) == nObs, ...
        ['roi_match_curated.mat has %d ROIs but roi_match_results.mat has %d.\n' ...
         'The matcher was re-run after you curated. Re-curate, or restore the ' ...
         'matching results file -- the two cannot be joined as they stand.'], numel(grpOf), nObs);
    grpSrc = 'curated (your hand curation)';
else
    grpOf = match.grp(:);
    grpSrc = 'AUTO (roi_match_results.mat) -- no curation found';
    warning('No roi_match_curated.mat; falling back to the automatic grouping.');
end

% Cell ids. If the curation was PORTED from a larger set it carries srcCellId, the
% cell id in that source. Use it, so every condition folder (baseline / roomAir /
% reoxygenated) numbers the same neuron identically. The numbering is then SPARSE --
% ids exist for cells with no member in this subset -- which every downstream script
% already tolerates because they skip empty cells. Without srcCellId, fall back to
% the positional numbering, which is only valid within one dataset.
cellOf = nan(nObs,1);
if isfile(cfg.curatedFile) && isfield(C.curated,'srcCellId') && ...
        numel(C.curated.srcCellId) == nObs
    cellOf  = C.curated.srcCellId(:);
    nCellId = max(cellOf(~isnan(cellOf)));
    cells   = cell(nCellId,1);
    for cc = 1:nCellId
        m = find(cellOf == cc);
        if ~isempty(m), [~,o] = sort(roi.fov(m)); cells{cc} = m(o); end
    end
    grpSrc = [grpSrc ' + INHERITED source cell ids'];
    fprintf('cell ids inherited from the source curation (sparse, max id %d)\n', nCellId);
else
    posIds  = sort(unique(grpOf(grpOf > 0)));
    cells   = cell(0,1);
    for g = 1:numel(posIds)
        m = find(grpOf == posIds(g));  [~,o] = sort(roi.fov(m));  m = m(o);
        cells{end+1,1} = m;  cellOf(m) = numel(cells); %#ok<AGROW>
    end
    ungr = find(grpOf == 0);
    for k = 1:numel(ungr)
        cells{end+1,1} = ungr(k);  cellOf(ungr(k)) = numel(cells); %#ok<AGROW>
    end
end
% tossed keep cellOf = NaN by construction
cellSize = cellfun(@numel, cells);

% ---- STABLE cell key ----------------------------------------------------
% cell_id is POSITIONAL: it is the index into the sorted group list, so it
% renumbers whenever the ROI subset changes. Splitting phys\ into baseline\ /
% roomAir\ / reoxygenated\ made baseline cell 87 the same neuron as full-set cell
% 159, and across 447 shared observations the numbers agreed only 2 times. A
% positional id cannot be used to refer to a cell across analyses.
%
% cell_key is derived from CONTENT instead: the anchor member's
% "<recording folder name>#<ROI index>", the member with the lexicographically
% smallest recording name and then the smallest ROI index. Both terms are the same
% invariants the curation port relies on, so the key survives directory moves and
% survives subsetting as long as the anchor member is in the subset. Use it, not
% cell_id, to refer to a cell across datasets; use cell_id_crosswalk_260727.m when
% the anchor itself was dropped.
fovBase = strings(nFov,1);
for f = 1:nFov, fovBase(f) = folder_basename(fov_folder(f)); end
cellKey = strings(numel(cells),1);
for cc = 1:numel(cells)
    m  = cells{cc};
    if isempty(m), continue; end     % inherited ids are SPARSE: a gap is a cell whose
                                     % members all live in another condition folder
    nm = fovBase(roi.fov(m));  ri = roi.roi(m);
    [~, ord] = sortrows([string(nm(:)), string(ri(:))]);
    a = m(ord(1));
    cellKey(cc) = fovBase(roi.fov(a)) + "#" + string(roi.roi(a));
end
obs_key = strings(nObs,1);
ok2 = ~isnan(cellOf);
obs_key(ok2) = cellKey(cellOf(ok2));

status = strings(nObs,1);
status(grpOf > 0)  = "grouped";
status(grpOf == 0) = "ungrouped";
status(grpOf < 0)  = "tossed";

%% ---- 3. the recordings the coherence pipeline will scan ----
allMat  = dir(fullfile(cfg.rootPath, '**', 'ca_spike_data.mat'));
recPath = string({allMat.folder}');
recName = strings(numel(recPath),1);
for k = 1:numel(recPath), recName(k) = folder_basename(recPath(k)); end

% per-FOV -> recording folder (the matcher stores the full path, so match on that)
fovRec   = strings(nFov,1);
for f = 1:nFov, fovRec(f) = folder_basename(fov_folder(f)); end

%% ---- 4. per-recording cross-check: do the ROI indices mean the same thing? ----
% This is the guard against a distorted join. For each recording we compare
%   nROI in the matcher (from cpSAM maskL)  vs  numel(roi_spikes)  vs  size(dFF,2)
% and check the matcher's labels are contiguous 1..N (if they are not, roi.roi is
% a RANK, not a label, and would address the wrong column downstream).
nRec = numel(recPath);
rec_nSpikeROI = nan(nRec,1);  rec_nDffROI = nan(nRec,1);  rec_nCa = nan(nRec,1);
rec_fps = nan(nRec,1);        rec_nMatchROI = zeros(nRec,1);
rec_hasBreath = false(nRec,1);  rec_inMatcher = false(nRec,1);
rec_maxRoiIdx = zeros(nRec,1);  rec_labContig = true(nRec,1);

for k = 1:nRec
    f = find(fov_folder == recPath(k), 1);
    if ~isempty(f)
        rec_inMatcher(k) = true;
        m = find(roi.fov == f);
        rec_nMatchROI(k) = numel(m);
        if ~isempty(m)
            rec_maxRoiIdx(k) = max(roi.roi(m));
            % labels contiguous 1..N  =>  roi.roi == maskL label == F column
            rec_labContig(k) = isequal(sort(roi.lab(m))', 1:numel(m));
        end
    end
    try
        CA = load(fullfile(recPath(k), 'ca_spike_data.mat'), 'roi_spikes');
        if isfield(CA,'roi_spikes')
            rec_nSpikeROI(k) = numel(CA.roi_spikes);
            if rec_nSpikeROI(k) > 0
                rec_nCa(k) = numel(CA.roi_spikes(1).spike_train);
            end
        end
    catch ME
        warning('  could not read ca_spike_data.mat in %s: %s', recName(k), ME.message);
    end
    dd = dir(fullfile(recPath(k), '*_dFF.mat'));
    if ~isempty(dd)
        try
            D = load(fullfile(dd(1).folder, dd(1).name), 'dFF');
            if isfield(D,'dFF'), rec_nDffROI(k) = size(D.dFF, 2); end
        catch
        end
    end
    rec_hasBreath(k) = isfile(fullfile(recPath(k),'breath_peak_pc1.mat')) && ...
                       isfile(fullfile(recPath(k),'breath_insp_start_pc1.mat'));
    rec_fps(k) = detect_session_fps(char(recPath(k)), fallback_fps);
end

recT = table(recName, recPath, rec_inMatcher, rec_nMatchROI, rec_nSpikeROI, rec_nDffROI, ...
             rec_maxRoiIdx, rec_labContig, rec_nCa, rec_fps, rec_hasBreath, ...
    'VariableNames', {'rec_name','rec_path','in_matcher','n_roi_matcher','n_roi_spikes', ...
                      'n_roi_dff','max_roi_index','labels_contiguous','n_frames_ca','fps','has_breath'});

%% ---- 5. the observation table ----
obs        = (1:nObs)';
obs_fov    = roi.fov(:);
obs_rec    = fovRec(obs_fov);
obs_recPath= fov_folder(obs_fov);
obs_roiIdx = roi.roi(:);
obs_cell   = cellOf;
obs_size   = nan(nObs,1);
ok = ~isnan(obs_cell);  obs_size(ok) = cellSize(obs_cell(ok));

% has this observation got a spike train to contribute?
obs_hasSpikes = false(nObs,1);  obs_nSpikes = nan(nObs,1);
for k = 1:nRec
    if isnan(rec_nSpikeROI(k)), continue; end
    sel = obs_recPath == recPath(k);
    if ~any(sel), continue; end
    CA = load(fullfile(recPath(k), 'ca_spike_data.mat'), 'roi_spikes');
    idx = find(sel);
    for q = idx(:)'
        ri = obs_roiIdx(q);
        if ri >= 1 && ri <= numel(CA.roi_spikes)
            obs_hasSpikes(q) = true;
            obs_nSpikes(q)   = sum(double(CA.roi_spikes(ri).spike_train(:)) > 0);
        end
    end
end

% Per-CELL coordinates (mean over the cell's observations) alongside the per-OBSERVATION
% ones, so a downstream figure can place a cell in the preparation without re-joining.
cell_x = nan(nObs,1); cell_y = nan(nObs,1); cell_z = nan(nObs,1);
for cc = 1:numel(cells)
    m = cells{cc};
    cell_x(m) = mean(roi.x(m));  cell_y(m) = mean(roi.y(m));  cell_z(m) = mean(roi.z(m));
end
obsT = table(obs, obs_cell, obs_key, obs_size, status, obs_rec, obs_roiIdx, roi.lab(:), ...
             roi.x(:), roi.y(:), roi.z(:), roi.cx_px(:), roi.cy_px(:), ...
             cell_x, cell_y, cell_z, obs_hasSpikes, obs_nSpikes, obs_recPath, ...
    'VariableNames', {'obs','cell_id','cell_key','cell_size','status','rec_name','roi_index','maskL_label', ...
                      'x_um','y_um','z_um','cx_px','cy_px', ...
                      'cell_x_um','cell_y_um','cell_z_um','has_spikes','n_spikes','rec_path'});

%% ---- 6. AUDIT: prove nothing was lost ----
a = struct();
a.grouping_source = grpSrc;
a.n_obs_total     = nObs;
a.n_obs_grouped   = nnz(grpOf > 0);
a.n_obs_ungrouped = nnz(grpOf == 0);
a.n_obs_tossed    = nnz(grpOf < 0);
a.n_cells         = nnz(cellSize > 0);   % non-empty only; inherited ids are sparse
a.max_cell_id     = numel(cells);
a.n_cells_multi   = nnz(cellSize > 1);
a.n_cells_single  = nnz(cellSize == 1);
a.n_obs_in_cells  = sum(cellSize);
a.n_rec_scanned   = nRec;
a.rec_not_in_matcher = recName(~rec_inMatcher);
a.fov_not_a_recording = setdiff(fovRec, recName);
a.rec_missing_breath  = recName(~rec_hasBreath);
a.rec_roi_count_mismatch = recName(rec_inMatcher & ~isnan(rec_nSpikeROI) & ...
                                   rec_nMatchROI ~= rec_nSpikeROI);
a.rec_roi_index_overflow = recName(rec_inMatcher & ~isnan(rec_nSpikeROI) & ...
                                   rec_maxRoiIdx > rec_nSpikeROI);
a.rec_labels_not_contiguous = recName(rec_inMatcher & ~rec_labContig);
a.n_obs_with_spikes = nnz(obs_hasSpikes);

fprintf('\n---- grouping (%s) ----\n', grpSrc);
fprintf('  observations : %d   (grouped %d | ungrouped %d | tossed %d)\n', ...
        a.n_obs_total, a.n_obs_grouped, a.n_obs_ungrouped, a.n_obs_tossed);
fprintf('  cells        : %d   (%d multi-recording, %d singleton)\n', ...
        a.n_cells, a.n_cells_multi, a.n_cells_single);
fprintf('  accounted for: %d in cells + %d tossed = %d of %d\n', ...
        a.n_obs_in_cells, a.n_obs_tossed, a.n_obs_in_cells + a.n_obs_tossed, a.n_obs_total);
assert(a.n_obs_in_cells + a.n_obs_tossed == a.n_obs_total, ...
    'LOSSLESS CHECK FAILED: %d of %d observations are unaccounted for.', ...
    a.n_obs_total - a.n_obs_in_cells - a.n_obs_tossed, a.n_obs_total);
fprintf('  cell size distribution: ');
for s = 1:max(cellSize), if nnz(cellSize==s), fprintf('%dx:%d  ', s, nnz(cellSize==s)); end, end
fprintf('\n');

fprintf('\n---- recordings ----\n');
fprintf('  with ca_spike_data.mat under rootPath : %d\n', nRec);
fprintf('  observations with a usable spike train: %d of %d\n', a.n_obs_with_spikes, nObs);
report_list('recordings NOT in the matcher (their ROIs have no cell identity)', a.rec_not_in_matcher);
report_list('matcher FOVs with no ca_spike_data.mat (contribute no coherence)', a.fov_not_a_recording);
report_list('recordings missing breath_peak/insp_start_pc1.mat', a.rec_missing_breath);
report_list('!! ROI COUNT MISMATCH matcher vs ca_spike_data', a.rec_roi_count_mismatch);
report_list('!! ROI INDEX OVERFLOW (matcher index > numel(roi_spikes))', a.rec_roi_index_overflow);
report_list('!! maskL labels NOT contiguous (roi.roi is a rank, not a column)', a.rec_labels_not_contiguous);

if isempty(a.rec_roi_count_mismatch) && isempty(a.rec_roi_index_overflow) && ...
   isempty(a.rec_labels_not_contiguous)
    fprintf('\n  JOIN IS CLEAN: ROI indices address the same neurons on both sides.\n');
else
    fprintf(['\n  JOIN IS SUSPECT for the recordings marked !! above. Their ROI indices\n' ...
             '  may not address the same neurons in the matcher and in ca_spike_data.\n' ...
             '  cell_pool_260727.m will refuse to pool them.\n']);
end

%% ---- 7. save ----
link = struct('obsT',obsT, 'recT',recT, 'cells',{cells}, 'cellOf',cellOf, 'cellKey',cellKey, 'cellSize',cellSize, ...
              'grpOf',grpOf, 'grouping_source',grpSrc, 'roi',roi, 'fov_name',fov_name, ...
              'fov_folder',fov_folder, 'audit',a, 'cfg',cfg);
save(cfg.linkFile, 'link', '-v7.3');
writetable(obsT, fullfile(cfg.cellDir, 'cell_link.csv'));
writetable(recT, fullfile(cfg.cellDir, 'cell_link_recordings.csv'));
fprintf('\nSaved:\n  %s\n  %s\n  %s\n', cfg.linkFile, ...
    fullfile(cfg.cellDir,'cell_link.csv'), fullfile(cfg.cellDir,'cell_link_recordings.csv'));
fprintf('Next: cell_pool_260727.m\n');
end

%% ========================= helpers =========================
function identity_link(cfg, fallback_fps)
% Build the same link structures the matcher path builds, but with every
% (recording, ROI) as its own 1x cell. Field-for-field compatible with
% cell_pool_260727.m, so nothing downstream needs to know which path was taken.
allMat  = dir(fullfile(cfg.rootPath, '**', 'ca_spike_data.mat'));
assert(~isempty(allMat), 'No ca_spike_data.mat found under %s', cfg.rootPath);
recPath = string({allMat.folder}');
nRec    = numel(recPath);
recName = strings(nRec,1);
for k = 1:nRec, recName(k) = folder_basename(recPath(k)); end

rec_nSpikeROI = nan(nRec,1); rec_nDffROI = nan(nRec,1); rec_nCa = nan(nRec,1);
rec_fps = nan(nRec,1); rec_hasBreath = false(nRec,1);
obs_rec = strings(0,1); obs_path = strings(0,1); obs_ri = []; obs_ns = []; obs_has = [];
for k = 1:nRec
    fp = char(recPath(k));
    try
        CA = load(fullfile(fp,'ca_spike_data.mat'),'roi_spikes');
        if isfield(CA,'roi_spikes')
            n = numel(CA.roi_spikes);
            rec_nSpikeROI(k) = n;
            if n > 0, rec_nCa(k) = numel(CA.roi_spikes(1).spike_train); end
            for r = 1:n
                obs_rec(end+1,1)  = recName(k); %#ok<AGROW>
                obs_path(end+1,1) = recPath(k); %#ok<AGROW>
                obs_ri(end+1,1)   = r;          %#ok<AGROW>
                obs_ns(end+1,1)   = sum(double(CA.roi_spikes(r).spike_train(:)) > 0); %#ok<AGROW>
                obs_has(end+1,1)  = true;       %#ok<AGROW>
            end
        end
    catch ME
        warning('  could not read ca_spike_data.mat in %s: %s', recName(k), ME.message);
    end
    dd = dir(fullfile(fp,'*_dFF.mat'));
    if ~isempty(dd)
        try
            D = load(fullfile(dd(1).folder,dd(1).name),'dFF');
            if isfield(D,'dFF'), rec_nDffROI(k) = size(D.dFF,2); end
        catch
        end
    end
    rec_hasBreath(k) = isfile(fullfile(fp,'breath_peak_pc1.mat')) && ...
                       isfile(fullfile(fp,'breath_insp_start_pc1.mat'));
    rec_fps(k) = detect_session_fps(fp, fallback_fps);
end

nObs   = numel(obs_ri);
cellOf = (1:nObs)';
cells  = num2cell(cellOf);
obs_key = obs_rec + "#" + string(obs_ri);   % identity grouping: the ROI IS the cell
obsT = table((1:nObs)', cellOf, obs_key, ones(nObs,1), repmat("ungrouped",nObs,1), obs_rec, ...
             obs_ri, obs_ri, nan(nObs,1), nan(nObs,1), nan(nObs,1), ...
             logical(obs_has), obs_ns, obs_path, ...
    'VariableNames',{'obs','cell_id','cell_key','cell_size','status','rec_name','roi_index','maskL_label', ...
                     'x_um','y_um','z_um','has_spikes','n_spikes','rec_path'});
recT = table(recName, recPath, false(nRec,1), zeros(nRec,1), rec_nSpikeROI, rec_nDffROI, ...
             zeros(nRec,1), true(nRec,1), rec_nCa, rec_fps, rec_hasBreath, ...
    'VariableNames',{'rec_name','rec_path','in_matcher','n_roi_matcher','n_roi_spikes', ...
                     'n_roi_dff','max_roi_index','labels_contiguous','n_frames_ca','fps','has_breath'});

a = struct('grouping_source','IDENTITY (no cell matching for this dataset)', ...
           'n_obs_total',nObs,'n_obs_grouped',0,'n_obs_ungrouped',nObs,'n_obs_tossed',0, ...
           'n_cells',nObs,'n_cells_multi',0,'n_cells_single',nObs,'n_obs_in_cells',nObs, ...
           'n_rec_scanned',nRec,'rec_not_in_matcher',recName,'fov_not_a_recording',strings(0,1), ...
           'rec_missing_breath',recName(~rec_hasBreath), ...
           'rec_roi_count_mismatch',strings(0,1),'rec_roi_index_overflow',strings(0,1), ...
           'rec_labels_not_contiguous',strings(0,1),'n_obs_with_spikes',nnz(obs_has));

fprintf('\n---- identity grouping ----\n');
fprintf('  recordings   : %d\n', nRec);
fprintf('  observations : %d (each its own 1x cell)\n', nObs);
fprintf('  with a spike train: %d\n', nnz(obs_has));
report_list('recordings missing breath_peak/insp_start_pc1.mat', a.rec_missing_breath);

link = struct('obsT',obsT,'recT',recT,'cells',{cells},'cellOf',cellOf,'cellKey',obs_key, ...
              'cellSize',ones(nObs,1),'grpOf',zeros(nObs,1), ...
              'grouping_source',a.grouping_source,'roi',struct(), ...
              'fov_name',recName,'fov_folder',recPath,'audit',a,'cfg',cfg);
save(cfg.linkFile,'link','-v7.3');
writetable(obsT, fullfile(cfg.cellDir,'cell_link.csv'));
writetable(recT, fullfile(cfg.cellDir,'cell_link_recordings.csv'));
fprintf('\nSaved (identity) cell_link.mat/.csv to\n  %s\nNext: cell_pool_260727.m\n', cfg.cellDir);
end

function report_list(ttl, lst)
if isempty(lst)
    fprintf('  %-62s none\n', [ttl ':']);
else
    fprintf('  %-62s %d\n', [ttl ':'], numel(lst));
    for k = 1:numel(lst), fprintf('        %s\n', lst(k)); end
end
end

function name = folder_basename(p)
% Full last path segment (fileparts mistreats folders with dots in the name).
p = char(p);
while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end) = []; end
[~,n,e] = fileparts(p);
name = string([n e]);
end
