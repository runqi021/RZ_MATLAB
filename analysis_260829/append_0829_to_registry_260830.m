% append_0829_to_registry_260830.m
% -----------------------------------------------------------------------
%  APPEND the 260824 vagotomised Vglut2 session to the pooled registry
%      D:\Ventral_surface_summary\event_latency_260811\event_latency_data.mat
%  so its cells get POOLED ids continuing after the existing 281, and can then
%  be rendered into per-cell-summary_active_260812 alongside every other session.
%
%  WHY APPEND RATHER THAN RE-RUN THE BUILDER
%    Ventral_surface_event_latency_260811.m scans scan_dirs in order and, within
%    a genotype, in dir() order. Vglut2\0824 sorts BETWEEN 0810 and 1124, so a
%    re-run would insert its cells mid-list and renumber everything from Vgat
%    onward -- invalidating all 817 figures already in per-cell-summary_*.
%    Appending at the end leaves cells 1..281 at exactly their current indices.
%
%  NOTHING EXISTING IS MODIFIED
%    - CELL/OBS/REC/S entries 1..end keep their indices and their contents.
%    - New REC/OBS/CELL/S are concatenated AFTER them; new cell ids start at
%      numel(CELL)+1 = 282.
%    - OBS.rec and CELL.obs of the NEW entries are offset so they index into the
%      concatenated arrays. Old entries need no adjustment because they point at
%      indices that have not moved.
%    - The .mat is rewritten with old content + new content. A timestamped backup
%      is made first and the script refuses to run without one.
%    - The script refuses to run at all if the session is already present, so a
%      second run cannot double-append.
%
%  IDENTITY COMES FROM THE SESSION'S cell_link
%    Vglut2\0824\cell_pooled\cell_link.mat (built by make_cell_link_260824.m from
%    your roi_match curation). Without it every ROI would enter the registry as
%    its own cell and the cross-FOV curation would be silently discarded, so the
%    script asserts the file exists rather than falling back.
%
%  THE STATISTICS ARE NOT RE-DERIVED. cell_latency() and draw_shifts() at the
%  bottom of this file are copied VERBATIM from lines 505-642 of
%  Ventral_surface_event_latency_260811.m. They are local functions of that
%  function file and so cannot be called from outside it; copying is the only way
%  to guarantee the appended S entries were computed by identical code. If that
%  file's statistics ever change, this copy must be updated with it.
%
%  Every analysis parameter below is copied from the builder's user block. They
%  MUST match or the appended cells are scored on a different footing from the
%  existing ones.
%
%  OUTPUT
%    - the registry, appended in place (backup alongside)
%    - <outDir>\cell_id_map_0824.csv : local cell id  <->  pooled cell id,
%      keyed on (rec_name, roi_index), including the local cells that get NO
%      pooled id because they fall below the activeMinEv gate.
%
%  Runqi Zhang / 2026-08-25
% -----------------------------------------------------------------------

clear; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir); addpath(repoRoot);
addpath(fullfile(repoRoot, 'analysis_260806'));
addpath(fullfile(repoRoot, 'analysis_260727', 'coh_ca_breath'));

%% ===================== USER-EDITABLE =====================
rootPath    = 'D:\Ventral_surface_summary';
regFile     = fullfile(rootPath, 'event_latency_260811', 'event_latency_data.mat');
sessionRel  = fullfile('Sert', '0829');         % the session to append
scanName    = 'Sert';                            % genotype for labels/keys
recDate     = '0829';
% linkFile = '' when the session was NOT cross-FOV matched. Every ROI then
% enters as its own cell, which is correct only when there really are no
% duplicates across recordings -- with a matched session, leaving this empty
% would silently split each cell into one per FOV.
%
% THIS IS THE STEP THAT GETS SKIPPED. Sert/0828 has a roi_match_curated.mat in
% BOTH its match folders and still went into the registry unlinked -- all 30 of
% its cells carry a 'roi:' key with exactly one observation, because no
% cell_pooled\cell_link.mat was ever built and this field was left empty. Same
% failure as Vglut2/0224 and 1124. 0829 has 108 cells spanning >1 recording, so
% leaving this blank would split every one of them.
linkFile    = fullfile(rootPath, sessionRel, 'cell_pooled', 'cell_link.mat');
mapCsv      = fullfile(rootPath, sessionRel, sprintf('cell_id_map_%s.csv', recDate));
dryRun      = false;      % true = compute and report, write nothing

% ---- copied verbatim from Ventral_surface_event_latency_260811.m ----
nDrop        = 30;
fallback_fps = 30;
activeMinEv  = 5;        % ACTIVE = pooled nnz(spike_train>0) > this
minTestEv    = 20;
caLagSec     = 0;
ampFrac      = 0.20;
nPerm        = 10000;
alphaPerm    = 0.001;
minShiftIBI  = 2;
smoothDiv    = 25;
rngSeed      = 260811;
% =========================================================

rng(rngSeed);
sessionDir = fullfile(rootPath, sessionRel);
assert(isfolder(sessionDir), 'session not in the archive: %s', sessionDir);
usingLink = ~isempty(linkFile);
if usingLink
    assert(isfile(linkFile), ['no cell_link for this session:%s  %s%s' ...
           'Run make_cell_link first -- without it every ROI would become ' ...
           'its own cell.'], newline, linkFile, newline);
end
assert(isfile(regFile), 'registry not found: %s', regFile);

fprintf('\n=========== append_0829_to_registry_260830 ===========\n');
fprintf('registry : %s\n', regFile);
fprintf('session  : %s\n', sessionDir);

%% ---- load the registry and refuse to double-append --------------------
R0 = load(regFile);                       % CELL, OBS, REC, S, groups, prm
req = {'CELL','OBS','REC','S','groups'};
for k = 1:numel(req)
    assert(isfield(R0, req{k}), 'registry has no %s', req{k});
end
CELL0 = R0.CELL;  OBS0 = R0.OBS;  REC0 = R0.REC;  S0 = R0.S;  groups = R0.groups;

% ---- KNOWN EXCEPTION: the S gap at cells 277-281 -----------------------
% The registry holds 281 CELL but only 276 S. Cells 277-281 are the deep Vgat
% cells added on 2026-08-15 by register_deep_recordings_260815, which appended
% CELL/OBS/REC without extending S. RZ 2026-08-25: accepted as-is, not to be
% recomputed here.
%
% BUT THE GAP CANNOT SIMPLY BE IGNORED. S is addressed positionally -- S(c)
% belongs to CELL(c). If S were left 5 short, the first appended S entry would
% land at index 277 and silently become the statistics of Vgat cell 277, and
% every 0824 cell after it would be shifted onto the wrong neuron. So the gap is
% filled with BLANK placeholders: tested = false, nEv = 0, everything else NaN --
% the same state cell_latency returns for a cell it declined to test. That keeps
% S(c) <-> CELL(c) true, leaves those five uncomputed exactly as they are now,
% and makes them visibly untested rather than quietly absent.
%
% To compute them properly later: rerun with the placeholders replaced by
% cell_latency(CELL0(c), OBS0, REC0, nPerm, minShiftIBI, minTestEv, alphaPerm,
% smoothDiv), which is the call the builder itself uses.
nSgap = numel(CELL0) - numel(S0);
if nSgap > 0
    gap = numel(S0)+1 : numel(CELL0);
    fprintf(2, 'KNOWN EXCEPTION: %d cell(s) have no S entry (%d..%d)\n', ...
            nSgap, gap(1), gap(end));
    fprintf(2, '   filling with blank placeholders to keep S(c) aligned to CELL(c);\n');
    fprintf(2, '   they stay untested, and no existing S entry is touched.\n');
    for c = gap
        S0(c) = blank_S(S0(1));
        fprintf('   placeholder S(%d)  %s\n', c, CELL0(c).label);
    end
end
assert(numel(S0) == numel(CELL0), ...
    'registry is still inconsistent: numel(S)=%d but numel(CELL)=%d', numel(S0), numel(CELL0));

already = find(contains({REC0.folder}, sessionDir, 'IgnoreCase', true));
if ~isempty(already)
    error('append_0829:alreadyPresent', ...
        ['%s is ALREADY in the registry (%d recordings, first at REC %d).\n' ...
         'Appending again would duplicate every cell. Restore a backup if you ' ...
         'need to redo it.'], sessionRel, numel(already), already(1));
end
fprintf('before   : %d cells, %d observations, %d recordings\n', ...
        numel(CELL0), numel(OBS0), numel(REC0));

%% ---- backup -----------------------------------------------------------
bkDir = fileparts(regFile);
bk    = fullfile(bkDir, sprintf('event_latency_data_backup_pre%s_%s.mat', ...
                 recDate, strrep(sessionRel, filesep, '-')));
if ~dryRun
    if ~isfile(bk), copyfile(regFile, bk); end
    d0 = dir(regFile); db = dir(bk);
    assert(~isempty(db) && db.bytes == d0.bytes, ...
        'backup missing or wrong size -- refusing to touch the registry');
    fprintf('backup   : %s (%.1f MB)\n', bk, db.bytes/1e6);
end

%% ---- cell identity from this session's cell_link ----------------------
% Same key construction as the builder: '<Genotype>/<MMDD>/<recName>/<maskL>'.
cell_map   = containers.Map('KeyType','char','ValueType','char');
tossed_set = containers.Map('KeyType','char','ValueType','logical');
pre = sprintf('%s/%s', scanName, recDate);
if ~usingLink
    Tk = table();
    fprintf(2,['no cell_link: every ROI enters as its own cell. Correct ONLY if ' ...
               'this session has no duplicate cells across recordings.%s'], newline);
else
Lk  = load(linkFile,'link');  Tk = Lk.link.obsT;
for i = 1:height(Tk)
    kk = sprintf('%s/%s/%d', pre, Tk.rec_name(i), Tk.maskL_label(i));
    if isnan(Tk.cell_id(i)), tossed_set(kk) = true;
    else,                    cell_map(kk) = sprintf('%s#c%d', pre, Tk.cell_id(i));
    end
end
fprintf('identity : %d masks -> %d cells (%d tossed)\n', height(Tk), ...
        numel(unique(Tk.cell_id(~isnan(Tk.cell_id)))), nnz(isnan(Tk.cell_id)));
end   % usingLink

gi = find(strcmp(groups, scanName), 1);
assert(~isempty(gi), 'no group "%s" in the registry''s groups list', scanName);

%% ---- SCAN the session (loop copied from the builder) ------------------
RECn = struct('fps',{},'T',{},'Trec',{},'trig',{},'gi',{},'folder',{});
OBSn = struct('rec',{},'ev',{},'gi',{},'label',{},'cellKey',{},'nEv',{});
allMat = dir(fullfile(sessionDir, '**', 'ca_spike_data.mat'));
fprintf('\n=== scanning %d recordings with spikes ===\n', numel(allMat));
for kk = 1:numel(allMat)
    folderPath = allMat(kk).folder;
    recName    = folder_basename(folderPath);
    try
        bpFile = fullfile(folderPath,'breath_peak_pc1.mat');
        ipFile = fullfile(folderPath,'breath_insp_start_pc1.mat');
        if ~isfile(bpFile) || ~isfile(ipFile)
            fprintf('  skip (no breath): %s\n', recName); continue;
        end
        fps = detect_session_fps(folderPath, fallback_fps);
        CA  = load(fullfile(folderPath,'ca_spike_data.mat'),'roi_spikes');
        nCa = numel(CA.roi_spikes(1).spike_train);

        BP = load(bpFile);  IP = load(ipFile);
        bw = detrend(double(BP.breath(:)));
        bw(1:min(nDrop,numel(bw))) = [];  bw = bw - mean(bw);
        peak_idx = round(BP.insp_onset_idx(:)) - nDrop;   % insp PEAK
        foot_idx = round(IP.insp_start_idx(:)) - nDrop;   % insp ONSET

        T = min([numel(bw), nCa]);
        peak_idx = peak_idx(peak_idx>=1 & peak_idx<=T);
        foot_idx = foot_idx(foot_idx>=1 & foot_idx<=T);
        if numel(peak_idx) < 3 || numel(foot_idx) < 3
            fprintf('  skip (too few landmarks): %s\n', recName); continue;
        end
        bw = bw(1:T);

        % breath-cycle amplitude QC, same rule as the builder
        bwz = (bw - median(bw)) / max(mad(bw,1)*1.4826, eps);
        ft  = sort(foot_idx);
        amp = nan(numel(ft)-1,1);  pkf = nan(numel(ft)-1,1);
        for i = 1:numel(ft)-1
            q = peak_idx(peak_idx>ft(i) & peak_idx<ft(i+1));
            if ~isempty(q), amp(i) = bwz(q(1)) - bwz(ft(i)); pkf(i) = q(1); end
        end
        good   = amp > ampFrac*median(amp,'omitnan');
        onTrig = ft(good);
        pkTrig = pkf(good);  pkTrig = pkTrig(~isnan(pkTrig));
        if numel(onTrig) < 3 || numel(pkTrig) < 3
            fprintf('  skip (QC left too few triggers): %s\n', recName); continue;
        end

        RECn(end+1) = struct('fps',fps,'T',T,'Trec',T/fps, ...
                             'trig',{{sort(onTrig(:))/fps, sort(pkTrig(:))/fps}}, ...
                             'gi',gi,'folder',folderPath); %#ok<SAGROW>
        recIdx = numel(RECn);

        lag  = round(caLagSec*fps);
        nInc = 0;
        for rid = 1:numel(CA.roi_spikes)
            st = double(CA.roi_spikes(rid).spike_train(:));
            st = st(1:min(T,numel(st)));
            if numel(st) < T, st(end+1:T,1) = 0; end %#ok<AGROW>
            ckin = sprintf('%s/%s/%s/%d', scanName, recDate, recName, rid);
            if isKey(tossed_set, ckin), continue; end
            if lag > 0, st = [st(1+lag:end); zeros(lag,1)]; end %#ok<AGROW>
            lab = sprintf('%s/%s/%s/%d', scanName, recDate, recName, rid);
            if isKey(cell_map, ckin), ck = cell_map(ckin); else, ck = ['roi:' lab]; end
            OBSn(end+1) = struct('rec',recIdx,'ev',find(st>0),'gi',gi, ...
                                 'label',lab,'cellKey',ck,'nEv',nnz(st>0)); %#ok<SAGROW>
            nInc = nInc + 1;
        end
        fprintf('  [%2d] %-46s fps %.2f  %3d ROI  %3d/%3d cycles pass QC\n', ...
                kk, recName(1:min(46,end)), fps, nInc, numel(onTrig), numel(ft)-1);
    catch ME
        warning('  ERROR %s: %s', recName, ME.message);
    end
end
assert(~isempty(OBSn), 'no observations collected from %s', sessionDir);

nOrphan = nnz(startsWith({OBSn.cellKey}, 'roi:'));
if nOrphan > 0
    fprintf(2, ['note: %d observation(s) had no cell_link entry and would enter as ' ...
                'their own cell\n'], nOrphan);
end

%% ---- CELLS (same gate as the builder) ---------------------------------
[uCell, ~, obsOfCell] = unique({OBSn.cellKey}, 'stable');
CELLn = struct('key',{},'gi',{},'obs',{},'nEv',{},'nObs',{},'label',{});
localOfNew = [];                      % local cell id behind each new registry cell
for c = 1:numel(uCell)
    m  = find(obsOfCell == c);
    ne = sum([OBSn(m).nEv]);
    if ne <= activeMinEv, continue; end
    if numel(m) == 1, lab = OBSn(m(1)).label;
    else,             lab = sprintf('%s [+%d rec]', OBSn(m(1)).label, numel(m)-1);
    end
    CELLn(end+1) = struct('key',uCell{c},'gi',OBSn(m(1)).gi,'obs',m(:)', ...
                          'nEv',ne,'nObs',numel(m),'label',lab); %#ok<SAGROW>
    tok = regexp(uCell{c}, '#c(\d+)$', 'tokens', 'once');
    if isempty(tok), localOfNew(end+1,1) = NaN; else, localOfNew(end+1,1) = str2double(tok{1}); end %#ok<SAGROW>
end
fprintf('\n%d observations -> %d keys -> %d ACTIVE cells (>%d events)\n', ...
        numel(OBSn), numel(uCell), numel(CELLn), activeMinEv);
assert(~isempty(CELLn), 'no cell cleared the activeMinEv gate -- nothing to append');

%% ---- offset the new indices into the concatenated arrays --------------
nREC0 = numel(REC0);  nOBS0 = numel(OBS0);  nCELL0 = numel(CELL0);
for j = 1:numel(OBSn),  OBSn(j).rec = OBSn(j).rec + nREC0;  end
for c = 1:numel(CELLn), CELLn(c).obs = CELLn(c).obs + nOBS0; end

RECc  = [REC0(:);  RECn(:)].';
OBSc  = [OBS0(:);  OBSn(:)].';
CELLc = [CELL0(:); CELLn(:)].';

%% ---- latency + permutation test for the NEW cells only ----------------
fprintf('\ntesting %d new cells x 2 triggers, %d shifts each...\n', numel(CELLn), nPerm);
tStart = tic;
Sn = S0([]);                                    % same struct type, empty
for c = 1:numel(CELLn)
    Sn(c) = cell_latency(CELLc(nCELL0+c), OBSc, RECc, nPerm, minShiftIBI, ...
                         minTestEv, alphaPerm, smoothDiv); %#ok<SAGROW>
    if mod(c,5)==0, fprintf('  %d/%d  (%.0f s)\n', c, numel(CELLn), toc(tStart)); end
end
fprintf('  done in %.0f s\n', toc(tStart));
Sc = [S0(:); Sn(:)].';

testedN = vertcat(Sn.tested);  sigN = vertcat(Sn.sig);
fprintf('new cells: onset %d tested / %d sig,  peak %d tested / %d sig\n', ...
        nnz(testedN(:,1)), nnz(sigN(:,1)), nnz(testedN(:,2)), nnz(sigN(:,2)));

%% ---- final consistency checks BEFORE writing --------------------------
% isequalN, not isequal: these structs are full of NaN, and plain isequal calls
% NaN ~= NaN, so it reports "changed" for entries that are in fact identical.
% S0 here is the POST-REPAIR baseline, so a repaired tail entry is expected to
% differ from what was on disk -- everything else must not.
bad = {};
if numel(Sc) ~= numel(CELLc), bad{end+1} = 'numel(S) ~= numel(CELL)'; end
if ~isequaln(CELLc(1:nCELL0), CELL0), bad{end+1} = 'existing CELL entries changed'; end
if ~isequaln(OBSc(1:nOBS0),  OBS0),   bad{end+1} = 'existing OBS entries changed';  end
if ~isequaln(RECc(1:nREC0),  REC0),   bad{end+1} = 'existing REC entries changed';  end
if ~isequaln(Sc(1:nCELL0),   S0),     bad{end+1} = 'existing S entries changed';    end
for c = nCELL0+1:numel(CELLc)
    o = CELLc(c).obs;
    if any(o < 1 | o > numel(OBSc)), bad{end+1} = sprintf('cell %d: obs out of range', c); end %#ok<SAGROW>
    if any([OBSc(o).rec] < 1 | [OBSc(o).rec] > numel(RECc))
        bad{end+1} = sprintf('cell %d: rec out of range', c); %#ok<SAGROW>
    end
end
if ~isempty(bad)
    fprintf(2,'CONSISTENCY CHECK FAILED:\n'); fprintf(2,'   %s\n', bad{:});
    error('append_0829:consistency','refusing to write the registry');
end
fprintf('checks   : OK -- entries 1..%d byte-identical, new indices in range\n', nCELL0);

fprintf('\nafter    : %d cells (ids %d..%d are new), %d observations, %d recordings\n', ...
        numel(CELLc), nCELL0+1, numel(CELLc), numel(OBSc), numel(RECc));

%% ---- write ------------------------------------------------------------
if dryRun
    fprintf('\nDRY RUN -- nothing written. Set dryRun = false to append.\n');
    return;
end

out = R0;                                  % preserve every original variable
out.CELL = CELLc;  out.OBS = OBSc;  out.REC = RECc;  out.S = Sc;
out.append_log = struct('session', sessionRel, 'appendedOn', '2026-08-25', ...
    'by', 'append_0829_to_registry_260830.m', ...
    'linkFile', linkFile, 'backup', bk, ...
    'firstNewCell', nCELL0+1, 'lastNewCell', numel(CELLc), ...
    'nNewRec', numel(RECn), 'nNewObs', numel(OBSn), 'nNewCell', numel(CELLn));
save(regFile, '-struct', 'out');
fprintf('registry appended: %s\n', regFile);

%% ---- the local <-> pooled map ----------------------------------------
% Only meaningful when the session HAS a local numbering to map from. Without a
% cell_link there is one cell per ROI and the pooled id is the only id, so there
% is nothing to translate.
if usingLink
localIds = unique(Tk.cell_id(~isnan(Tk.cell_id)));
fid = fopen(mapCsv,'w');
fprintf(fid,'local_cell_id,pooled_cell_id,n_recordings,rec_name,roi_index,n_events,reason\n');
nMapped = 0;
for i = 1:numel(localIds)
    lid  = localIds(i);
    rows = find(Tk.cell_id == lid);
    k    = find(localOfNew == lid, 1);
    if isempty(k)
        pid = NaN;  ne = NaN;  why = sprintf('below activeMinEv (>%d events)', activeMinEv);
    else
        pid = nCELL0 + k;  ne = CELLn(k).nEv;  why = 'pooled';  nMapped = nMapped + 1;
    end
    for r = rows(:)'
        fprintf(fid,'%d,%s,%d,%s,%d,%s,%s\n', lid, num2str_nan(pid), numel(rows), ...
                Tk.rec_name(r), Tk.roi_index(r), num2str_nan(ne), why);
    end
end
fclose(fid);
fprintf('map      : %s\n', mapCsv);
fprintf('           %d of %d local cells got a pooled id\n', nMapped, numel(localIds));
else
    fprintf('map      : skipped (no cell_link -- pooled id is the only id)%s', newline);
end

fprintf('\nNext: run per_cell_summary_260812.m with onlyCells = %d:%d and overwrite = false\n', ...
        nCELL0+1, numel(CELLc));
fprintf('Done.\n');

% =======================================================================
% ======================== LOCAL FUNCTIONS ==============================
% =======================================================================
function s = num2str_nan(v)
    if isnan(v), s = ''; else, s = sprintf('%g', v); end
end

function s = blank_S(template)
% An untested S entry, in the field ORDER of the registry's existing S (built
% from a template so concatenation cannot fail on field mismatch). Same state
% cell_latency returns for a cell it declined to test.
    s = template;
    s.ibi=[NaN NaN]; s.fps=NaN; s.L={[],[]};
    s.mode=[NaN NaN]; s.fwhm=[NaN NaN]; s.frac=[NaN NaN];
    s.dens={[],[]};   s.dctr={[],[]};
    s.med=[NaN NaN];  s.q1=[NaN NaN];  s.q3=[NaN NaN];  s.iqr=[NaN NaN];
    s.p=[NaN NaN];    s.z=[NaN NaN];
    s.nEv=[0 0];      s.nCyc=[0 0];
    s.tested=[false false]; s.sig=[false false];
end

function name = folder_basename(p)
% fileparts treats "...dir.x" as filename + ".x"; rebuild the last segment.
    p = char(p);
    while ~isempty(p) && (p(end)=='/' || p(end)=='\'), p(end)=[]; end
    [~,n,e] = fileparts(p);
    name = [n e];
end

% -----------------------------------------------------------------------
%  cell_latency and draw_shifts below are COPIED VERBATIM from lines 505-642 of
%  analysis_260806\Ventral_surface_event_latency_260811.m. Do not "improve" them
%  here: their only job is to be identical to the code that produced S for cells
%  1..281, so the appended cells are scored the same way.
% -----------------------------------------------------------------------
function S = cell_latency(C, OBS, REC, nPerm, minShiftIBI, minTestEv, alphaPerm, smoothDiv)
% Latencies and the permutation test for ONE cell, pooled over its recordings.
S = struct('ibi',[NaN NaN],'fps',NaN,'L',{{[],[]}}, ...
           'mode',[NaN NaN],'fwhm',[NaN NaN],'frac',[NaN NaN], ...
           'dens',{{[],[]}},'dctr',{{[],[]}}, ...
           'med',[NaN NaN],'q1',[NaN NaN],'q3',[NaN NaN],'iqr',[NaN NaN], ...
           'p',[NaN NaN],'z',[NaN NaN],'nEv',[0 0],'nCyc',[0 0], ...
           'tested',[false false],'sig',[false false]);

S.fps = mode(arrayfun(@(j) REC(OBS(j).rec).fps, C.obs));

for q = 1:2
    dt = []; nCyc = 0;
    for j = C.obs
        ts = REC(OBS(j).rec).trig{q};
        if numel(ts) > 1, dt = [dt; diff(ts)]; nCyc = nCyc + numel(ts) - 1; end %#ok<AGROW>
    end
    if isempty(dt), continue; end
    ibi = mean(dt);  half = ibi/2;
    S.ibi(q) = ibi;  S.nCyc(q) = nCyc;

    dtm  = 1000/S.fps;
    nf   = max(2, floor(1000*half/dtm));
    ctr  = (-nf:nf)*dtm;
    eg   = [ctr - dtm/2, ctr(end) + dtm/2];
    sg   = max(1000*ibi/smoothDiv, dtm);
    gx   = -3*sg : dtm : 3*sg;
    g    = exp(-0.5*(gx/sg).^2);  g = g/sum(g);

    Lobs = [];  cnt = zeros(numel(ctr), nPerm);
    for j = C.obs
        r  = REC(OBS(j).rec);
        ts = r.trig{q};  Trec = r.Trec;
        tE = OBS(j).ev(:)/r.fps;
        tE = tE(tE >= 0 & tE <= Trec);
        if isempty(tE), continue; end

        te   = [ts(:).'-Trec, ts(:).', ts(:).'+Trec];
        edg  = [-inf, (te(1:end-1)+te(2:end))/2, inf];

        L0 = tE - te(discretize(tE, edg)).';
        L0(abs(L0) > half) = NaN;
        Lobs = [Lobs; 1000*L0]; %#ok<AGROW>

        sh = draw_shifts(nPerm, Trec, minShiftIBI*ibi);
        M  = mod(tE - sh(:).', Trec);
        Ls = M - te(discretize(M, edg));
        Ls(abs(Ls) > half) = NaN;
        b  = discretize(1000*Ls, eg);
        ok = ~isnan(b);
        if any(ok(:))
            [~, cc] = ind2sub(size(b), find(ok));
            % DEVIATION FROM THE VERBATIM COPY (RZ, 2026-08-25), and a bug that
            % is still live in Ventral_surface_event_latency_260811.m line 578.
            % When an observation holds exactly ONE event, b is 1 x nPerm, so
            % b(ok) and cc come back as ROW vectors and [b(ok), cc] is 1 x 2K
            % instead of K x 2 -- accumarray then rejects the size argument.
            % Forcing both to columns is a no-op for every multi-event
            % observation (they are already columns) and fixes the single-event
            % case. The archive never hit this because no cell there had a
            % one-event recording; this sparse vagotomised session does.
            bb = b(ok);
            cnt = cnt + accumarray([bb(:), cc(:)], 1, [numel(ctr) nPerm]);
        end
    end

    L = Lobs(~isnan(Lobs));
    S.L{q} = L;  S.nEv(q) = numel(L);
    if numel(L) < minTestEv, continue; end

    hO = conv(histcounts(L, eg), g, 'same');
    TO = max(hO) - mean(hO);
    HS = conv2(cnt, g(:), 'same');
    TN = max(HS,[],1) - mean(HS,1);

    S.p(q)      = (1 + nnz(TN >= TO)) / (1 + nPerm);
    S.z(q)      = (TO - mean(TN)) / max(std(TN), eps);
    S.tested(q) = true;
    S.sig(q)    = S.p(q) < alphaPerm;
    S.dens{q}   = hO / max(nCyc,1);
    S.dctr{q}   = ctr;

    [pk, im] = max(hO);
    base = median(hO);  hm = base + (pk - base)/2;
    lo = im; while lo > 1         && hO(lo-1) >= hm, lo = lo - 1; end
    hi = im; while hi < numel(hO) && hO(hi+1) >= hm, hi = hi + 1; end
    S.mode(q) = ctr(im);
    S.fwhm(q) = ctr(hi) - ctr(lo);
    S.frac(q) = nnz(L >= ctr(lo) & L <= ctr(hi)) / numel(L);

    W  = 1000*ibi;
    d  = mod(L - S.mode(q) + W/2, W) - W/2;
    Qd = prctile(d, [25 50 75]);
    S.med(q) = S.mode(q) + Qd(2);
    S.q1(q)  = S.mode(q) + Qd(1);
    S.q3(q)  = S.mode(q) + Qd(3);
    S.iqr(q) = Qd(3) - Qd(1);
end
end

function sh = draw_shifts(nPerm, Trec, shMin)
% Uniform circular shifts, excluding a band of +/- shMin around no-shift.
if 2*shMin >= Trec, sh = Trec*rand(nPerm,1); return; end
sh = shMin + (Trec - 2*shMin)*rand(nPerm,1);
end
