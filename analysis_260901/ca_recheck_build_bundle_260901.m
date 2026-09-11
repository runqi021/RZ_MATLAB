% ca_recheck_build_bundle_260901.m
% -----------------------------------------------------------------------
%  Pull every ACTIVE archive cell's dF/F, its existing event train and its
%  breath landmarks out of the archive and write them, one self-contained .mat
%  PER CELL, into a bundle folder that the re-detection GUI reads.
%
%  WHY A BUNDLE. The archive stores data per RECORDING (a folder holding one
%  dF/F matrix for all its ROIs). The re-check is per CELL: a cell is one
%  neuron, which may have been imaged in several recordings, and the question
%  being asked of it -- are these events real -- is asked of the cell, with the
%  option of throwing away individual recordings of it. Walking folders and ROIs,
%  which is what calcium_spike_gui.m does, cannot express that. So the mapping
%  from cells to (folder, ROI) observations is resolved ONCE here, and the GUI
%  never has to know the archive layout.
%
%  THE CELL NUMBERING IS THE FIGURES' NUMBERING. The cell list is rebuilt with
%  exactly the steps per_cell_summary_260812.m uses -- the same registry, the
%  same cell_merge_overrides_260814 merges, the same "a merged cell inherits the
%  LOWEST of its constituent numbers" rule -- so bundle file <group>_<date>_cellNNN
%  is the same cell as per-cell-summary_active_260812\<group>_<date>_cellNNN.png.
%  Cell IDs are positional and NOT stable across registry rebuilds, so the run
%  cross-checks its stems against the PNGs already on disk and prints both
%  directions of the mismatch rather than assuming they agree.
%
%  ACTIVE is re-measured here, not read from the stale index CSV: pooled
%  nnz(spike_train>0) / pooled duration >= 2 events/min, the same quantity and
%  the same gate per_cell_summary_260812.m applied. It is free here because the
%  traces are being loaded anyway.
%
%  ONLY dF/F. No raw F, no masks, no movie -- the re-check is a judgement about
%  traces. The breath waveform and its onset/peak trains come along because a
%  tossed recording has to drop out of the breath-triggered heatmap and the polar
%  statistics downstream, and those need the breath side to be in the bundle too.
%
%  OUTPUT  <archive>\per-cell-summary_active_260812\spike_recheck_260901\
%      cells\<group>_<date>_cellNNN.mat   one cell: cellInfo + REC struct array
%      bundle_index.mat                   the whole cell -> observation mapping
%      bundle_cells.csv                   one row per cell  (all cells, gate flag)
%      bundle_recordings.csv              one row per (cell, recording)
%      build_log.txt                      what this run did
%
%  Curation NEVER lands in these files. The GUI writes to curated\ beside them,
%  so re-running this script can only ever refresh the source snapshot.
%
%  Runqi Zhang / 2026-09-01
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(scriptDir);
addpath(fullfile(repoRoot, 'analysis_260806'));
addpath(fullfile(repoRoot, 'analysis_260727', 'coh_ca_breath'));

%% ===================== USER-EDITABLE =====================
rootPath  = 'D:\Ventral_surface_summary';
inFile    = fullfile(rootPath, 'event_latency_260811', 'event_latency_data.mat');
figDir    = fullfile(rootPath, 'per-cell-summary_active_260812');
outDir    = fullfile(figDir,  'spike_recheck_260901');

% ACTIVE gate -- pooled events/min over all of the cell's recordings. Must match
% activeMinRateHz in per_cell_summary_260812.m or the bundle and the figure
% folder describe different populations.
activeMinRateHz = 2/60;

% Bundle every cell in the registry, gate flag and all, instead of only the
% active ones. Off by default: the re-check is of the active population.
bundleAll   = false;

% Restrict to these cell NUMBERS ([] = all). For repairing individual cells.
onlyCells   = [];

overwrite   = true;    % false = leave existing cell .mat files alone

P = struct();
P.nDrop        = 30;   % breath-side frames dropped -- MUST match the figures
P.fallback_fps = 30;
% =========================================================================

if ~isfolder(outDir),                     mkdir(outDir);                     end
if ~isfolder(fullfile(outDir,'cells')),   mkdir(fullfile(outDir,'cells'));   end
if ~isfolder(fullfile(outDir,'curated')), mkdir(fullfile(outDir,'curated')); end

logFile = fullfile(outDir, 'build_log.txt');
fid = fopen(logFile, 'w');
say = @(varargin) local_say(fid, varargin{:});

say('ca_recheck_build_bundle_260901   %s\n', datestr(now)); %#ok<TNOW1,DATST>
say('archive : %s\n', rootPath);
say('bundle  : %s\n\n', outDir);

%% ===================== BUILD THE CELL LIST =====================
% Same three steps as per_cell_summary_260812.m, in the same order.
assert(isfile(inFile), 'Not found:\n  %s', inFile);
Din  = load(inFile, 'CELL', 'OBS', 'REC');
CELL = Din.CELL;  OBS = Din.OBS;  REC = Din.REC;
assert(isfield(REC,'folder'), 'REC has no .folder -- rerun the event-latency build');
say('registry: %d cells, %d observations, %d recordings\n', ...
    numel(CELL), numel(OBS), numel(REC));

[obsMerged, mergeInfo, origCells] = ...
        apply_cell_merges_260814(CELL, OBS, cell_merge_overrides_260814());
for q = 1:numel(mergeInfo)
    if ~strcmp(mergeInfo(q).status,'merged')
        say('  merge group %d: %s\n', mergeInfo(q).group, mergeInfo(q).status);
    end
end
CELL2 = struct('key',{},'obs',{});
cellNum = zeros(numel(obsMerged),1);
for i = 1:numel(obsMerged)
    cellNum(i)   = min(origCells{i});      % merged cell keeps the lowest number
    CELL2(i).obs = obsMerged{i};
    CELL2(i).key = CELL(cellNum(i)).key;
end
say('cell merges: %d groups applied, %d cells -> %d\n', ...
    nnz(strcmp({mergeInfo.status},'merged')), numel(CELL), numel(CELL2));
CELL = CELL2;

% One entry per CELL, carrying every recording it appears in.
CJ = struct('cell',{},'key',{},'group',{},'obs',{});
for c = 1:numel(CELL)
    O = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
    for o = CELL(c).obs(:)'
        parts = regexp(OBS(o).label, '/', 'split');
        if numel(parts) < 3, continue; end
        roi = str2double(parts{end});
        fp  = REC(OBS(o).rec).folder;
        if ~isfinite(roi) || roi < 1 || ~isfolder(fp), continue; end
        O(end+1) = struct('folder',fp, 'roi',roi, ...
                          'recName',strjoin(parts(3:end-1),'/'), ...
                          'group',parts{1}, 'recDate',parts{2}); %#ok<SAGROW>
    end
    if isempty(O), continue; end
    CJ(end+1) = struct('cell',cellNum(c), 'key',CELL(c).key, ...
                       'group',O(1).group, 'obs',O); %#ok<SAGROW>
end
nObsTot = sum(arrayfun(@(x) numel(x.obs), CJ));
say('%d cells with usable recordings (%d observations)\n\n', numel(CJ), nObsTot);

if ~isempty(onlyCells)
    missing = setdiff(onlyCells, [CJ.cell]);
    if ~isempty(missing)
        say('onlyCells: %s not in the cell list\n', mat2str(missing));
    end
    CJ = CJ(ismember([CJ.cell], onlyCells));
    assert(~isempty(CJ), 'onlyCells matched no cells');
    say('onlyCells: restricted to %d cell(s)\n\n', numel(CJ));
end

%% ============ FLATTEN TO OBSERVATIONS, ORDERED BY FOLDER ============
% The loader caches one folder at a time, so visiting observations in folder
% order turns 2840 file reads into one per folder.
J = []; K = [];
for j = 1:numel(CJ)
    J = [J; repmat(j, numel(CJ(j).obs), 1)];             %#ok<AGROW>
    K = [K; (1:numel(CJ(j).obs))'];                      %#ok<AGROW>
end
oFolder = strings(numel(J),1);  oRecName = strings(numel(J),1);
oGroup  = strings(numel(J),1);  oDate    = strings(numel(J),1);
oRoi    = zeros(numel(J),1);
for n = 1:numel(J)
    oo = CJ(J(n)).obs(K(n));
    oFolder(n)  = string(oo.folder);   oRecName(n) = string(oo.recName);
    oGroup(n)   = string(oo.group);    oDate(n)    = string(oo.recDate);
    oRoi(n)     = oo.roi;
end
[~, ord] = sort(oFolder);
J = J(ord); K = K(ord);
oFolder = oFolder(ord); oRecName = oRecName(ord);
oGroup  = oGroup(ord);  oDate    = oDate(ord);  oRoi = oRoi(ord);
say('%d observations over %d distinct folders\n', numel(J), numel(unique(oFolder)));

%% ===================== PASS 1: MEASURE THE ACTIVE GATE =====================
% Cheap: only nEv, T and fps are kept. Traces are not stored.
obsEv   = nan(numel(J),1);
obsDur  = nan(numel(J),1);
obsFps  = nan(numel(J),1);
obsHas  = false(numel(J),1);
obsFail = strings(numel(J),1);
ca_recheck_load_obs_260901('clearcache');
tP = tic;
for n = 1:numel(J)
    o = CJ(J(n)).obs(K(n));
    try
        r = ca_recheck_load_obs_260901(o, P);
        obsEv(n)  = r.nEv;  obsDur(n) = r.dur_s;
        obsFps(n) = r.fps;  obsHas(n) = r.has_spike;
    catch ME
        obsFail(n) = string(ME.message);
    end
    if mod(n,200)==0, fprintf('  pass 1: %d/%d\n', n, numel(J)); end
end
say('pass 1 measured %d observations in %.1f s\n', numel(J), toc(tP));
if any(obsFail ~= "")
    say('  %d observations FAILED to load:\n', nnz(obsFail ~= ""));
    for n = find(obsFail ~= "")'
        say('    cell %3d  %s  roi %d : %s\n', CJ(J(n)).cell, ...
            CJ(J(n)).obs(K(n)).recName, CJ(J(n)).obs(K(n)).roi, obsFail(n));
    end
end
say('  %d of %d observations have no ca_spike_data.mat (0 events, full duration)\n', ...
    nnz(~obsHas & obsFail==""), nnz(obsFail==""));

nEvCell  = zeros(numel(CJ),1);
durCell  = zeros(numel(CJ),1);
for n = 1:numel(J)
    if obsFail(n) ~= "", continue; end
    nEvCell(J(n)) = nEvCell(J(n)) + obsEv(n);
    durCell(J(n)) = durCell(J(n)) + obsDur(n);
end
rateCell = nEvCell ./ max(durCell, eps);
isActive = rateCell >= activeMinRateHz;
say('\nACTIVE (pooled rate >= %.3g ev/min): %d of %d cells\n', ...
    activeMinRateHz*60, nnz(isActive), numel(CJ));

if bundleAll, wanted = true(numel(CJ),1); else, wanted = isActive; end

%% ===================== PASS 2: WRITE THE BUNDLE =====================
stems  = strings(numel(CJ),1);
nWrote = 0; nSkip = 0;
bundleFile = strings(numel(CJ),1);
for j = 1:numel(CJ)
    stems(j) = string(regexprep(sprintf('%s_%s_cell%03d', CJ(j).group, ...
                     CJ(j).obs(1).recDate, CJ(j).cell), '[\\/:*?"<>|]', '_'));
end

% Which observations pass 2 needs, again in folder order.
keepObs = wanted(J) & (obsFail == "");
J2 = J(keepObs); K2 = K(keepObs);
BUF = cell(numel(CJ),1);          % accumulate each cell's recordings

ca_recheck_load_obs_260901('clearcache');
tW = tic;
for n = 1:numel(J2)
    j = J2(n);
    o = CJ(j).obs(K2(n));
    r = ca_recheck_load_obs_260901(o, P);
    if isempty(BUF{j}), BUF{j} = r; else, BUF{j}(end+1) = r; end
    if mod(n,200)==0, fprintf('  pass 2: %d/%d\n', n, numel(J2)); end
end
say('pass 2 loaded %d observations in %.1f s\n', numel(J2), toc(tW));

for j = find(wanted)'
    if isempty(BUF{j}), say('  cell %d: no loadable recording, skipped\n', CJ(j).cell); continue; end
    stem  = char(stems(j));
    fOut  = fullfile(outDir, 'cells', [stem '.mat']);
    bundleFile(j) = string(fOut);
    if ~overwrite && isfile(fOut), nSkip = nSkip + 1; continue; end

    % Put the cell's recordings back in registry order -- pass 2 walked folders.
    RECc = BUF{j};
    [~, back] = sort(arrayfun(@(x) find(strcmp({CJ(j).obs.folder}, x.folder) & ...
                                        [CJ(j).obs.roi] == x.roi, 1), RECc));
    RECc = RECc(back);

    cellInfo = struct( ...
        'cell',         CJ(j).cell, ...
        'key',          CJ(j).key, ...
        'group',        CJ(j).group, ...
        'recDate',      CJ(j).obs(1).recDate, ...
        'stem',         stem, ...
        'nRec',         numel(RECc), ...
        'pooled_nEv',   nEvCell(j), ...
        'pooled_dur_s', durCell(j), ...
        'rate_min',     rateCell(j)*60, ...
        'isActive',     isActive(j), ...
        'summary_png',  fullfile(figDir, [stem '.png']), ...
        'built',        datestr(now), ...      %#ok<TNOW1,DATST>
        'params',       P);
    save(fOut, 'cellInfo', 'RECc', '-v7.3');
    nWrote = nWrote + 1;
end
say('wrote %d cell files (%d skipped, already present)\n\n', nWrote, nSkip);

%% ===================== INDEX + CROSS-CHECK =====================
nRecCell = arrayfun(@(x) numel(x.obs), CJ)';
Tc = table([CJ.cell]', stems, string({CJ.key})', string({CJ.group})', ...
           nRecCell, nEvCell, durCell, rateCell*60, ...
           isActive, wanted, bundleFile, ...
    'VariableNames', {'cell_idx','stem','cell_key','group','n_recordings', ...
                      'pooled_n_events','pooled_dur_s','rate_per_min', ...
                      'is_active','bundled','bundle_file'});
writetable(Tc, fullfile(outDir,'bundle_cells.csv'));

Tr = table([CJ(J).cell]', stems(J), oGroup, oDate, ...
           K, oRecName, oFolder, oRoi, ...
           obsFps, obsDur, obsEv, obsEv./max(obsDur,eps)*60, obsHas, obsFail, ...
    'VariableNames', {'cell_idx','stem','group','recDate','rec_in_cell', ...
                      'recName','folder','roi','fps','dur_s','n_events', ...
                      'rate_per_min','has_ca_spike_data','load_error'});
Tr = sortrows(Tr, {'cell_idx','rec_in_cell'});
writetable(Tr, fullfile(outDir,'bundle_recordings.csv'));

BUNDLE = struct('CJ',CJ, 'stems',stems, 'isActive',isActive, 'wanted',wanted, ...
                'rate_min',rateCell*60, 'nEv',nEvCell, 'dur_s',durCell, ...
                'outDir',outDir, 'figDir',figDir, 'params',P, ...
                'activeMinRateHz',activeMinRateHz, 'built',datestr(now)); %#ok<TNOW1,DATST>
save(fullfile(outDir,'bundle_index.mat'), 'BUNDLE', '-v7.3');

% ---- do the bundle stems match the figures already on disk? --------------
pngs = dir(fullfile(figDir, '*.png'));
pngStem = string(regexprep({pngs.name}, '\.png$', ''));
pngStem = pngStem(~endsWith(pngStem, '_avgproj'));
onlyPng = setdiff(pngStem, stems(wanted));
onlyBun = setdiff(stems(wanted), pngStem);
say('cross-check against %s\n', figDir);
say('  %d summary PNGs on disk, %d cells bundled\n', numel(pngStem), nnz(wanted));
say('  %d PNG(s) with no bundled cell, %d bundled cell(s) with no PNG\n', ...
    numel(onlyPng), numel(onlyBun));
for s = onlyPng(:)', say('    PNG only   : %s\n', s); end
for s = onlyBun(:)', say('    bundle only: %s\n', s); end

say('\n---- done ----\n  %s\n', outDir);
fclose(fid);
fprintf('\nlog: %s\n', logFile);

% -------------------------------------------------------------------------
function local_say(fid, varargin)
fprintf(varargin{:});
fprintf(fid, varargin{:});
end
