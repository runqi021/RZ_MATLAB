function archive_into_summary_260729(mode)
%% archive_into_summary_260729  Copy datasets into Ventral_surface_summary.
% -----------------------------------------------------------------------
%   archive_into_summary_260729('template')   write an EDITABLE mapping csv
%   archive_into_summary_260729('copy')       do the copy, reading that csv
%
% The archive is organised  <Genotype>\<MMDD>\<site>\<recording>\  where site is
% cell1..cellN or an anatomical label (IO, 7N, MAO). Nothing else about the
% recording folders changes -- their NAMES are preserved exactly, which is what
% keeps any existing curation portable (see roi_curation_port_260727).
%
% TWO STEPS ON PURPOSE. Assigning a site to 31 folders is a judgement call about
% the recordings, not something to infer from a filename prefix. So 'template'
% writes its best guess into a csv and you correct the `site` column; 'copy' then
% does exactly what the csv says. Nothing is written to the archive until you have
% seen and approved the mapping.
%
% COPY, NOT MOVE. Sources are left in place, so the analysis outputs already built
% there keep working and the archive copy can be verified before anything is
% deleted.
%
% Runqi Zhang / 2026-07-29

if nargin < 1, mode = 'template'; end
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here,'coh_ca_breath'));

ARCHIVE = 'D:\Ventral_surface_summary';
mapFile = fullfile(here, 'site_mapping_260729.csv');

% ---- what goes in: {source root, genotype, MMDD} ----
% Both sessions moved D: -> E: on 2026-07-30. The archive copies were made
% before the move and are unaffected (recording folder NAMES are what keep
% curation portable, not their absolute path), but 'template' re-reads the
% sources, so these must track the move or it fails its own isfolder assert.
SRC = { ...
    'E:\260721_Sert_soma_G8s\phys\baseline', 'Sert',   '0721'
    'E:\260728_vglut2_soma-g8s\phys',        'Vglut2', '0728'
    'D:\260730_vgat-g8m_shiverer\phys',      'Vgat',   '0730'   % added 2026-07-31
    };

switch lower(mode)
%% ==================================================================
case 'template'
    rows = {};
    for s = 1:size(SRC,1)
        root = SRC{s,1};  gen = SRC{s,2};  dt = SRC{s,3};
        assert(isfolder(root), 'source not found: %s', root);
        d = dir(root);  d = d([d.isdir]);
        d = d(~ismember({d.name},{'.','..'}));
        d = d(~startsWith({d.name}, {'analysis_','roi_match_','coherence_','_','.'}));
        for i = 1:numel(d)
            p = fullfile(root, d(i).name);
            nSpk = double(isfile(fullfile(p,'ca_spike_data.mat')));
            nTrg = double(isfile(fullfile(p,'breath_peak_pc1.mat')) && ...
                          isfile(fullfile(p,'breath_insp_start_pc1.mat')));
            % ACTIVE cells are the only ones the analysis uses, so count them here:
            % n_active = ROIs with >=1 event, n_analysis = ROIs clearing the >=20
            % event floor the PETH test applies.
            nEv = NaN; nRoi = NaN; nAct = NaN; nAna = NaN;
            if nSpk
                try
                    C = load(fullfile(p,'ca_spike_data.mat'),'roi_spikes');
                    ev   = arrayfun(@(r) nnz(r.spike_train>0), C.roi_spikes);
                    nRoi = numel(ev);  nEv = sum(ev);
                    nAct = nnz(ev > 5);   % VENTRAL's 'active' criterion, matched exactly:
                                          %   spike_trigger_dFF.m:34 / temporal_phase_perROI.m:56
                                          %   active ROI = nnz(spike_train>0) > minEvents, minEvents=5
                    nAna = nnz(ev >= 20); % the PETH test's own floor, stricter on purpose
                catch
                end
            end
            [st, note] = guess_site(d(i).name, nAct);
            if any(strcmp(d(i).name, coh_cfg_260727().excludeRecordings))
                st = "";  note = "EXCLUDED: different stage zero reference (see coh_cfg_260727)";
            end
            rows(end+1,:) = {string(gen), string(dt), st, string(d(i).name), ...
                             nSpk>0, nTrg>0, nRoi, nAct, nAna, nEv, note, string(p)}; %#ok<AGROW>
        end
    end
    T = cell2table(rows, 'VariableNames', ...
        {'genotype','date','site','recording','has_spikes','has_triggers', ...
         'n_rois','n_active','n_analysis','n_events','note','source_path'});
    writetable(T, mapFile);
    fprintf('\n=========== site mapping template ===========\n');
    fprintf('%d recording folders. PROPOSED sites (edit the `site` column):\n\n', height(T));
    disp(T(:,{'genotype','date','site','recording','has_spikes','has_triggers','n_events'}));
    fprintf('\nEdit: %s\n', mapFile);
    fprintf('Then: archive_into_summary_260729(''copy'')\n');
    fprintf('Nothing has been written to the archive.\n');

%% ==================================================================
case 'copy'
    assert(isfile(mapFile), ['No mapping file:\n  %s\n' ...
        'Run archive_into_summary_260729(''template'') first and edit it.'], mapFile);
    T = readtable(mapFile,'TextType','string');
    % readtable parses "0721" as the NUMBER 721, and char(721) would silently create
    % a folder with a garbage name. Force the date back to a 4-digit MMDD string.
    if isnumeric(T.date), T.date = string(compose('%04d', T.date)); end
    % A blank cell in the csv reads back as <missing>, not "", and char(<missing>)
    % errors. Normalise so the "no site given -> skip" test actually fires.
    T.site = string(T.site);  T.site(ismissing(T.site)) = "";
    T.genotype = string(T.genotype);
    assert(isfolder(ARCHIVE), 'archive not found: %s', ARCHIVE);
    fprintf('\n=========== copying into the archive ===========\n');
    fprintf('archive: %s\n', ARCHIVE);
    nDone = 0; nSkip = 0;
    for i = 1:height(T)
        if strlength(strtrim(T.site(i))) == 0
            fprintf('  SKIP (no site given): %s\n', T.recording(i));  nSkip = nSkip + 1;  continue;
        end
        dst = fullfile(ARCHIVE, char(T.genotype(i)), char(T.date(i)), ...
                       char(strtrim(T.site(i))), char(T.recording(i)));
        if isfolder(dst)
            fprintf('  SKIP (already there): %s\n', dst);  nSkip = nSkip + 1;  continue;
        end
        [ok,msg] = mkdir(dst);
        assert(ok, 'could not create %s: %s', dst, msg);
        [ok,msg] = copyfile(fullfile(char(T.source_path(i)),'*'), dst);
        if ~ok
            fprintf(2,'  FAILED %s: %s\n', T.recording(i), msg);
        else
            fprintf('  %s\\%s\\%s\\%s\n', T.genotype(i), T.date(i), strtrim(T.site(i)), T.recording(i));
            nDone = nDone + 1;
        end
    end
    fprintf('\ncopied %d, skipped %d. Sources untouched.\n', nDone, nSkip);
end
end

%% ========================= helpers =========================
function [s, note] = guess_site(nm, nActive)
% Best guess only -- the `site` column is meant to be edited. A BLANK site means
% the row is skipped by 'copy'.
%
% Convention read off the existing archive: in Vglut2\0224 the pFN recordings are
% filed as cell1..cell5, one imaging site each, and IO is its own label. MAO does
% not appear anywhere in the archive.
nm = char(nm);  note = "";
tok = regexp(nm, '^[A-Za-z]+', 'match', 'once');
num = regexp(nm, '^[A-Za-z]+(\d+)', 'tokens', 'once');
if strcmpi(tok,'MAO')
    s = "";  note = "EXCLUDED: MAO is part of IO, not a pFN site";
    return;
end
if any(strcmpi(tok, {'IO','7N'}))
    s = string(upper(tok));
elseif ~isempty(num)
    s = "cell" + string(num{1});
else
    s = "";  note = "no site could be guessed -- fill it in";
    return;
end
if ~isnan(nActive) && nActive == 0
    s = "";  note = "EXCLUDED: 0 active ROIs (>5 events)";
elseif isnan(nActive)
    s = "";  note = "EXCLUDED: no ca_spike_data yet, so no active cells";
end
end
