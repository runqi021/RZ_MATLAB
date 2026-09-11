% stage_0829_to_archive_260830.m
% -----------------------------------------------------------------------
%  Stage the 260829 vagotomised Sert session into the pooled archive
%      D:\Ventral_surface_summary\Sert\0829\<site>\<recording>
%  so that make_cell_link / append_to_registry / per_cell_summary can see it.
%
%  WHAT GETS COPIED
%    Only the 10 light files each downstream tool actually opens -- NOT the
%    multi-GB TIFFs. Sert\0828 came to 33 MB for 3 recordings; this session's 24
%    recordings land around 260 MB. The acquisition folder stays the master copy.
%
%  SITE FOLDERS ARE NAMED cell01..cellNN AND THAT NAME IS LOAD-BEARING
%    make_cell_link_260824.m:91 enumerates sites with dir(fullfile(session,'cell*')),
%    so a site folder that does not start with "cell" is invisible to it. The site
%    name does NOT enter the registry cell key (keys are
%    roi:<Genotype>/<MMDD>/<recording>/<roi>), so it is organisational only.
%
%  SITES COME FROM motorPosition, NOT FROM THE FOLDER NAME
%    The roi number in the folder name is unreliable in this session -- measured
%    2026-08-30:
%      * roi12_8x_z0_15lp_6000f sits at (-828,-1139), 28 um from roi9_6x, while the
%        other three roi12 recordings are at (-620,+910), over 2 mm away. 60% of its
%        ROIs match roi9_3x within 10 um. It is a roi9 site.
%      * roi7 and roi14 are ONE site: roi14_3x_z0 x roi7_3x_z0 match 96.6% of ROIs
%        at a 3.3-6.3 um median. All six recordings lie within ~100 um.
%      * roi4 SPLITS across two sites, (-1083,+771) and (-1241,+988), 265 um apart,
%        so sites cannot be named after roi numbers without colliding.
%    Single-linkage on stage XY at SITE_CUT_UM is therefore the authority.
%
%  Writes site_map_0829.csv alongside, recording which recording went where.
%
%  Runqi Zhang / 2026-08-30
% -----------------------------------------------------------------------

clear; clc;

%% ===================== USER-EDITABLE =====================
srcRoot     = 'C:\Users\Admin\Desktop\260829_Sert-soma-g8s_vagotomized';
archiveRoot = 'D:\Ventral_surface_summary';
genotype    = 'Sert';
dateStr     = '0829';
matchSrc    = 'roi_match_out_260828';   % the curation folder INSIDE srcRoot. Named
                                        % _260828 because cell_cfg_260727 still carried
                                        % last session's outDirName; the contents are
                                        % this session's curation (02:03, 2026-08-30).
matchDst    = 'roi_match_out_260829';   % renamed on the way in, so the archive is honest
SITE_CUT_UM = 150;                      % single-linkage cut on stage XY
dryRun      = false;                    % true = report the layout, copy nothing
%% =========================================================

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir); addpath(repoRoot);

sessionDir = fullfile(archiveRoot, genotype, dateStr);
assert(isfolder(srcRoot), 'source not found: %s', srcRoot);
assert(isfolder(archiveRoot), 'archive not found: %s', archiveRoot);

%% ---- enumerate recordings + read stage position ------------------------
d = dir(srcRoot);
d = d([d.isdir]);
d = d(~ismember({d.name}, {'.','..'}));
d = d(~startsWith({d.name}, '.'));              % .breathcache
d = d(~startsWith({d.name}, 'roi_match_out'));  % the curation folder is not a recording

rec = struct('name',{},'src',{},'motor',{});
for i = 1:numel(d)
    fp = fullfile(srcRoot, d(i).name);
    if isempty(dir(fullfile(fp, '*_cpSAM_output.mat'))), continue; end
    [~, sm] = detect_session_fps(fp, 30);
    assert(isfield(sm,'motorPosition') && numel(sm.motorPosition) >= 3, ...
        'no motorPosition for %s', d(i).name);
    rec(end+1) = struct('name', d(i).name, 'src', fp, ...
                        'motor', sm.motorPosition(:).'); %#ok<SAGROW>
end
n = numel(rec);
assert(n > 0, 'no recordings with cpSAM output under %s', srcRoot);

P = vertcat(rec.motor);
site = cluster(linkage(pdist(P(:,1:2)),'single'), 'cutoff', SITE_CUT_UM, ...
               'criterion', 'distance');

% Number sites by their mean stage X so the numbering is a property of the
% anatomy rather than of dir() order -- dir() order is exactly what makes cell
% ids unstable elsewhere (see feedback_cell_ids_not_stable).
mx = arrayfun(@(s) mean(P(site==s,1)), 1:max(site));
[~, ord] = sort(mx);
remap(ord) = 1:max(site); %#ok<AGROW>
site = remap(site).';

fprintf('=========== stage_0829_to_archive_260830 ===========\n');
fprintf('source  : %s\n', srcRoot);
fprintf('archive : %s\n', sessionDir);
fprintf('%d recordings -> %d sites at %g um single-linkage\n\n', n, max(site), SITE_CUT_UM);

%% ---- the 10 files each recording contributes ---------------------------
% Matched against the Sert\0828 staging, which is what the downstream tools were
% built and tested on. Patterns, because most carry the recording name as a stem.
patterns = { ...
    'breath_pc1.mat', 'breath_peak_pc1.mat', 'breath_insp_start_pc1.mat', ...
    'ca_spike_data.mat', ...
    '*_ch1_dFF.mat', '*_ch1_meta.mat', ...
    '*_ch1_preproc_MC_MC_AVG_for_CP.tif', ...
    '*_ch1_preproc_MC_MC_AVG_ROIlabel.tif', ...
    '*_ch1_preproc_MC_MC_AVG_ROImask.tif', ...
    '*_ch1_preproc_MC_MC_cpSAM_output.mat'};

rows = {};
nCopy = 0; nBytes = 0; missing = {};
for s = 1:max(site)
    siteName = sprintf('cell%02d', s);
    idx = find(site == s);
    fprintf('%s  (%d recording%s)\n', siteName, numel(idx), plural(numel(idx)));
    for k = idx(:)'
        dst = fullfile(sessionDir, siteName, rec(k).name);
        if ~dryRun && ~isfolder(dst), mkdir(dst); end
        got = 0;
        for p = 1:numel(patterns)
            f = dir(fullfile(rec(k).src, patterns{p}));
            if isempty(f)
                missing{end+1} = sprintf('%s : %s', rec(k).name, patterns{p}); %#ok<SAGROW>
                continue;
            end
            for q = 1:numel(f)
                srcF = fullfile(f(q).folder, f(q).name);
                if ~dryRun, copyfile(srcF, fullfile(dst, f(q).name)); end
                got = got + 1; nBytes = nBytes + f(q).bytes;
            end
        end
        nCopy = nCopy + 1;
        fprintf('    %-38s %2d files  (%8.1f %8.1f %6.1f)\n', ...
            rec(k).name, got, rec(k).motor(1), rec(k).motor(2), rec(k).motor(3));
        rows(end+1,:) = {siteName, rec(k).name, rec(k).motor(1), ...
                         rec(k).motor(2), rec(k).motor(3)}; %#ok<SAGROW>
    end
end

%% ---- carry the curation across --------------------------------------------
srcMatch = fullfile(srcRoot, matchSrc);
dstMatch = fullfile(sessionDir, matchDst);
if isfolder(srcMatch)
    if ~dryRun
        if ~isfolder(dstMatch), mkdir(dstMatch); end
        copyfile(fullfile(srcMatch,'*'), dstMatch);
    end
    fprintf('\ncuration: %s -> %s\n', matchSrc, matchDst);
else
    fprintf(2, '\nno curation folder at %s\n', srcMatch);
end

%% ---- site map --------------------------------------------------------------
T = cell2table(rows, 'VariableNames', {'site','recording','motorX','motorY','motorZ'});
csvOut = fullfile(sessionDir, sprintf('site_map_%s.csv', dateStr));
if ~dryRun, writetable(T, csvOut); end

fprintf('\nstaged %d recordings, %.0f MB\n', nCopy, nBytes/1e6);
if ~isempty(missing)
    fprintf(2, '\nMISSING (%d):\n', numel(missing));
    fprintf(2, '   %s\n', missing{:});
end
fprintf('site map: %s\n', csvOut);
if dryRun, fprintf('\nDRY RUN -- nothing was written.\n'); end

fprintf('\nnext:\n  make_cell_link_260829\n  append_0829_to_registry_260830\n');

function s = plural(n), if n == 1, s = ''; else, s = 's'; end, end
