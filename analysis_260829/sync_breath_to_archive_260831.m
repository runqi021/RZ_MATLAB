% sync_breath_to_archive_260831.m
% -----------------------------------------------------------------------
%  Push re-derived BREATH files from an acquisition folder into the pooled
%  archive, which is what every downstream script actually reads.
%
%  WHY THIS IS NEEDED AT ALL. The registry's REC.folder points into
%  D:\Ventral_surface_summary\<Genotype>\<MMDD>\<site>\<recording>, and
%  temporal_phase_cell_fig_260812.m opens breath_pc1 / breath_peak_pc1 /
%  breath_insp_start_pc1 from THAT folder. Re-running the breath pipeline on the
%  acquisition drive therefore changes NOTHING downstream, silently: the analysis
%  keeps using the archived copy and no error is raised. This is the same class
%  of failure as a curation that never became a cell_link.
%
%  IT REFUSES TO LEAVE THE THREE FILES INCONSISTENT. The foot detector
%  (find_feet_deriv in breathing_trough_gui_pc1.m) searches each INTER-PEAK gap,
%  so the feet are a function of the peaks, which are a function of the trace.
%  Copying a new trace and new peaks while leaving old feet behind produces a
%  set that no single run ever generated. The script reports any recording whose
%  breath_insp_start_pc1.mat is older than its breath_peak_pc1.mat and, unless
%  ALLOW_STALE_FEET is set, copies nothing for it.
%
%  WHAT IT DOES NOT DO. It does not touch anything DERIVED from the breath --
%  the registry's per-cell statistics in event_latency_data.mat, the popsel
%  caches, or the rendered per-cell figures. Those were computed from the old
%  breath and stay stale until they are rebuilt. The report at the end lists
%  which ones are affected so the staleness is at least visible.
%
%  Runqi Zhang / 2026-08-31
% -----------------------------------------------------------------------

clear; clc;

%% ===================== USER-EDITABLE =====================
srcRoot   = 'C:\260824_Vglut2-soma-g8s_vagotomized\phys';
archiveSession = 'D:\Ventral_surface_summary\Vglut2\0824';
FILES     = {'breath_pc1.mat','breath_peak_pc1.mat','breath_insp_start_pc1.mat'};
ALLOW_STALE_FEET = true;    % RZ 2026-08-31: reviewed all feet against the new peaks, they stand
dryRun    = false;
%% =========================================================

assert(isfolder(srcRoot), 'source not found: %s', srcRoot);
assert(isfolder(archiveSession), 'archive session not found: %s', archiveSession);

% Map recording name -> archive folder, by scanning the cell* site folders.
sites = dir(fullfile(archiveSession,'cell*'));
sites = sites([sites.isdir]);
archOf = containers.Map('KeyType','char','ValueType','char');
for s = 1:numel(sites)
    rr = dir(fullfile(archiveSession, sites(s).name));
    rr = rr([rr.isdir] & ~startsWith({rr.name},'.'));
    for k = 1:numel(rr)
        archOf(rr(k).name) = fullfile(archiveSession, sites(s).name, rr(k).name);
    end
end
fprintf('archive holds %d recordings under %d site folder(s)\n', archOf.Count, numel(sites));

d = dir(srcRoot); d = d([d.isdir] & ~startsWith({d.name},'.'));
nCopy = 0; nSkipStale = 0; nNoArch = 0; stale = {}; touched = {};

fprintf('\n%-34s %-10s %s\n','recording','action','note');
for i = 1:numel(d)
    nm = d(i).name;
    src = fullfile(srcRoot, nm);
    if ~isfile(fullfile(src,'breath_pc1.mat')), continue; end     % not a recording
    if ~isKey(archOf, nm)
        fprintf('%-34s %-10s not in the archive\n', nm(1:min(33,end)), 'SKIP');
        nNoArch = nNoArch + 1; continue;
    end
    dst = archOf(nm);

    % consistency gate: feet must not predate the peaks they were cut from
    pk = dir(fullfile(src,'breath_peak_pc1.mat'));
    ft = dir(fullfile(src,'breath_insp_start_pc1.mat'));
    isStale = ~isempty(pk) && ~isempty(ft) && ft.datenum < pk.datenum;
    if isStale
        stale{end+1} = nm; %#ok<SAGROW>
        if ~ALLOW_STALE_FEET
            fprintf('%-34s %-10s feet (%s) older than peaks (%s)\n', nm(1:min(33,end)), ...
                'HOLD', datestr(ft.datenum,'mm-dd HH:MM'), datestr(pk.datenum,'mm-dd HH:MM'));
            nSkipStale = nSkipStale + 1; continue;
        end
    end

    moved = {};
    for f = 1:numel(FILES)
        sf = fullfile(src, FILES{f});
        if ~isfile(sf), continue; end
        df = fullfile(dst, FILES{f});
        if isfile(df)
            a = dir(sf); b = dir(df);
            if a.bytes == b.bytes && abs(a.datenum-b.datenum) < 1e-6, continue; end
        end
        if ~dryRun, copyfile(sf, df); end
        moved{end+1} = FILES{f}; %#ok<SAGROW>
    end
    if isempty(moved)
        fprintf('%-34s %-10s already current\n', nm(1:min(33,end)), '-');
    else
        fprintf('%-34s %-10s %s%s\n', nm(1:min(33,end)), 'COPIED', ...
            strjoin(strrep(moved,'.mat',''),', '), ternary(isStale,'   [stale feet forced]',''));
        nCopy = nCopy + 1; touched{end+1} = nm; %#ok<SAGROW>
    end
end

fprintf('\ncopied %d recording(s); %d held for stale feet; %d not in the archive\n', ...
        nCopy, nSkipStale, nNoArch);
if ~isempty(stale)
    fprintf(2,'\n%d recording(s) have feet older than their peaks:\n', numel(stale));
    fprintf(2,'   %s\n', stale{:});
    fprintf(2,['\nRe-run breathing_trough_gui_pc1 on these (point it at %s),\n' ...
               'then run this script again. Set ALLOW_STALE_FEET = true only if you\n' ...
               'have decided the old feet are still correct for the new peaks.\n'], srcRoot);
end
if dryRun, fprintf('\nDRY RUN -- nothing was written.\n'); end

if nCopy > 0
    fprintf(['\nNOW STALE, because they were computed from the OLD breath:\n' ...
        '  * event_latency_data.mat  -- per-cell Rayleigh/phase/latency for this session\n' ...
        '  * popsel_cache_*.mat      -- rebuild with popsel_run_260831.m (REBUILD = true)\n' ...
        '  * per-cell-summary_active_260812\\ figures for these cells\n' ...
        'Nothing above is regenerated here.\n']);
end

function s = ternary(c,a,b), if c, s = a; else, s = b; end, end
