% register_deep_recordings_260815.m
% -----------------------------------------------------------------------
%  Append the three DEEP Vgat/0730 recordings (z 250, 265, 320 um) to the event
%  latency registry so their cells get IDs and can be rendered by
%  per_cell_summary_260812.m and browsed in temporal_phase_cell_gui_260812.m.
%
%  APPEND-ONLY, BY DESIGN
%  Existing CELL entries keep their index, so every cell id already on a figure
%  filename, in pop_features.mat and in the per-cell CSVs stays valid. New cells
%  are added after the last existing one. Nothing is renumbered and nothing is
%  removed. A timestamped backup of the registry is written first.
%
%  CELL IDENTITY comes from the session's own curated matcher output, which DOES
%  cover all three deep FOVs (23 FOVs, 871 observations, 93 groups + 609
%  ungrouped, 101 tossed). So deep ROIs are grouped the same way surface ones are
%  rather than assumed one-cell-per-ROI, and curation-tossed masks are excluded.
%
%  ACTIVITY GATE = prm.activeMinEv (> 5 events), the archive's own definition, so
%  the new CELL entries mean the same thing as the existing ones. The stricter
%  rate gate (>= 2 ev/min) is applied later by per_cell_summary / pop_features,
%  exactly as it is for the surface cells.
%
%  WARNING: these depths are where the breath-locked optical artifact was
%  characterised. Registering them makes them renderable; it does NOT make their
%  phase-locking trustworthy. See deep_vgat_cells_260815.m.
%
%  Runqi Zhang / 2026-08-15
% -----------------------------------------------------------------------

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));

%% ===================== USER-EDITABLE =====================
sumRoot  = 'D:\Ventral_surface_summary';
regFile  = fullfile(sumRoot,'event_latency_260811','event_latency_data.mat');
matchDir = fullfile(sumRoot,'Vgat','0730','roi_match_out_260727');
GROUPNAME = 'Vgat';   DATESTR = '0730';

DEEP = { fullfile(sumRoot,'Vgat','0730','deep','roi1_2.4x_x1300y900_z250_3000f_30lp_00001')
         fullfile(sumRoot,'Vgat','0730','deep','roi1_3x_x1350y850_z265_3000f_31lp_00001')
         fullfile(sumRoot,'Vgat','0730','deep','roi1_2.4x_x1000y1000_z320_3000f_32lp_00001') };

dryRun = false;         % true = report what would change, write nothing
% =========================================================

assert(isfile(regFile), 'registry not found: %s', regFile);
D = load(regFile);
CELL = D.CELL; OBS = D.OBS; REC = D.REC; groups = D.groups; prm = D.prm;
nCell0 = numel(CELL); nObs0 = numel(OBS); nRec0 = numel(REC);
fprintf('registry before: %d cells, %d observations, %d recordings\n', nCell0, nObs0, nRec0);

gi = find(strcmpi(groups, GROUPNAME), 1);
assert(~isempty(gi), 'group "%s" not in registry groups', GROUPNAME);
activeMinEv = prm.activeMinEv;

%% ---- curation: (recName, roi) -> group label ----
Mm = load(fullfile(matchDir,'roi_match_results.mat'),'match');
Cc = load(fullfile(matchDir,'roi_match_curated.mat'));
if isfield(Cc,'curated') && isfield(Cc.curated,'grpOf'), gof = Cc.curated.grpOf(:);
else,                                                    gof = Mm.match.grp(:); end
fov_i = Mm.match.roi.fov(:);  roi_i = Mm.match.roi.roi(:);
curKey = containers.Map('KeyType','char','ValueType','any');
for k = 1:numel(gof)
    p = regexp(Mm.match.fov_folder{fov_i(k)}, '[\\/]', 'split');
    curKey([p{end} '|' num2str(roi_i(k))]) = gof(k);
end
fprintf('curation covers %d observations of %s/%s\n', numel(gof), GROUPNAME, DATESTR);

%% ---- build the new entries ----
existLeaf = strings(nRec0,1);
for i = 1:nRec0
    p = regexp(REC(i).folder,'[\\/]','split');  existLeaf(i) = p{end};
end
existLab = string({OBS.label}');

newCellKeys = containers.Map('KeyType','char','ValueType','double');  % key -> new CELL idx
nAddRec = 0; nAddObs = 0; nAddCell = 0; nToss = 0; nSkip = 0;

for r = 1:numel(DEEP)
    fp = DEEP{r};
    % NOT fileparts: these folder names contain dots ("roi1_2.4x_..."), and
    % fileparts eats everything after the first one as an extension. That turns
    % both the z250 and the z320 recording into "roi1_2" -- a collision that makes
    % them indistinguishable and breaks the curation lookup. Take the path leaf.
    pp = regexp(regexprep(fp,'[\\/]+$',''), '[\\/]', 'split');
    recName = pp{end};
    if any(existLeaf == string(recName))
        fprintf('SKIP (already registered): %s\n', recName);  nSkip = nSkip + 1;  continue;
    end
    sam = dir(fullfile(fp,'*_cpSAM_output.mat'));
    dff = dir(fullfile(fp,'*_ch1_dFF.mat'));
    if isempty(sam) || isempty(dff) || ~isfile(fullfile(fp,'breath_peak_pc1.mat'))
        fprintf(2,'SKIP (missing inputs): %s\n', recName);  continue;
    end

    Q  = load(fullfile(sam(1).folder,sam(1).name),'maskL');
    Fd = load(fullfile(dff(1).folder,dff(1).name),'dFF');
    [fps, ~] = detect_session_fps(fp, 30);
    Tfr  = size(Fd.dFF,1);

    % breath triggers, in SECONDS, same convention as the existing REC.trig:
    %   trig{1} = inspiration STARTS (feet), trig{2} = inspiration PEAKS
    BP = load(fullfile(fp,'breath_peak_pc1.mat'));
    pk_t = BP.insp_onsets_t(:);
    st_t = [];
    if isfile(fullfile(fp,'breath_insp_start_pc1.mat'))
        IP = load(fullfile(fp,'breath_insp_start_pc1.mat'));
        st_t = IP.insp_starts_t(:);
    end

    REC(end+1) = struct('fps',fps,'T',Tfr,'Trec',Tfr/fps,'trig',{{st_t, pk_t}}, ...
                        'gi',gi,'folder',fp); %#ok<SAGROW>
    recIdx = numel(REC);  nAddRec = nAddRec + 1;

    CA = struct('roi_spikes',[]);
    if isfile(fullfile(fp,'ca_spike_data.mat'))
        CA = load(fullfile(fp,'ca_spike_data.mat'),'roi_spikes');
    end

    rois = setdiff(unique(Q.maskL(:)),0).';
    nCellHere = 0;
    for q = rois
        lab = sprintf('%s/%s/%s/%d', GROUPNAME, DATESTR, recName, q);
        if any(existLab == string(lab)), continue; end

        ck = [recName '|' num2str(q)];
        gval = 0;                        % default: ungrouped = its own cell
        if isKey(curKey, ck), gval = curKey(ck); end
        if gval < 0 || isnan(gval)
            nToss = nToss + 1;  continue;              % rejected in curation
        end

        ev = [];
        if isfield(CA,'roi_spikes') && ~isempty(CA.roi_spikes) && q <= numel(CA.roi_spikes)
            ev = find(CA.roi_spikes(q).spike_train(:) > 0);
        end
        nEv = numel(ev);

        OBS(end+1) = struct('rec',recIdx,'ev',ev,'gi',gi,'label',lab, ...
                            'cellKey','','nEv',nEv); %#ok<SAGROW>
        obsIdx = numel(OBS);  nAddObs = nAddObs + 1;

        if nEv <= activeMinEv, continue; end           % not a cell by the archive rule

        if gval > 0, key = sprintf('%s/%s#g%d', GROUPNAME, DATESTR, gval);
        else,        key = sprintf('roi:%s', lab); end
        OBS(obsIdx).cellKey = key;

        if isKey(newCellKeys, key)                      % same cell, another FOV
            ci = newCellKeys(key);
            CELL(ci).obs  = [CELL(ci).obs, obsIdx];
            CELL(ci).nEv  = CELL(ci).nEv + nEv;
            CELL(ci).nObs = CELL(ci).nObs + 1;
        else
            CELL(end+1) = struct('key',key,'gi',gi,'obs',obsIdx,'nEv',nEv, ...
                                 'nObs',1,'label',lab); %#ok<SAGROW>
            newCellKeys(key) = numel(CELL);
            nAddCell = nAddCell + 1;  nCellHere = nCellHere + 1;
        end
    end
    fprintf('  %-46s fps %.1f, %d frames, %d ROIs -> %d new cell(s)\n', ...
            recName, fps, Tfr, numel(rois), nCellHere);
end

fprintf('\nadded: %d recordings, %d observations, %d cells (%d tossed in curation, %d recordings skipped)\n', ...
        nAddRec, nAddObs, nAddCell, nToss, nSkip);
fprintf('registry after: %d cells, %d observations, %d recordings\n', ...
        numel(CELL), numel(OBS), numel(REC));

%% ---- integrity checks BEFORE writing ----
assert(numel(CELL) >= nCell0, 'cells were removed');
for i = 1:nCell0
    assert(strcmp(CELL(i).key, D.CELL(i).key) && isequal(CELL(i).obs, D.CELL(i).obs), ...
        'existing cell %d changed -- append-only was violated', i);
end
for i = 1:nObs0
    assert(strcmp(OBS(i).label, D.OBS(i).label), 'existing observation %d changed', i);
end
allObs = [CELL.obs];
assert(numel(unique(allObs)) == numel(allObs), 'an observation belongs to two cells');
assert(max(allObs) <= numel(OBS), 'CELL.obs points past OBS');
assert(max([OBS.rec]) <= numel(REC), 'OBS.rec points past REC');
fprintf('integrity checks passed (existing %d cells untouched)\n', nCell0);

if nAddCell > 0
    fprintf('\nnew cell ids %d..%d:\n', nCell0+1, numel(CELL));
    for i = nCell0+1:numel(CELL)
        fprintf('   cell %3d  %-24s %3d events, %d obs  %s\n', i, CELL(i).key, ...
                CELL(i).nEv, CELL(i).nObs, CELL(i).label);
    end
end

%% ---- write ----
if dryRun
    fprintf('\nDRY RUN -- nothing written\n');  return;
end
bak = strrep(regFile,'.mat',sprintf('_backup_%s.mat', datestr(now,'yymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
copyfile(regFile, bak);
fprintf('\nbackup -> %s\n', bak);
S = D.S;  %#ok<NASGU>   % carried through untouched
save(regFile, 'CELL','OBS','REC','S','groups','prm','-v7.3');
fprintf('registry updated: %s\n', regFile);
