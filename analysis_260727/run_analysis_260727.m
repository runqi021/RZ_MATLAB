function run_analysis_260727(rootPath, genotype)
%% run_analysis_260727  Run the whole 260727 breath x calcium chain on one dataset.
% -----------------------------------------------------------------------
%   run_analysis_260727()                              % the default dataset in coh_cfg_260727.m
%   run_analysis_260727(rootPath, genotype)            % any other dataset
%
% e.g.  run_analysis_260727('D:\Ventral_surface_summary','Ventral')
%
% Everything lands in <rootPath>\analysis_260727\ and nothing outside it is
% touched, so a dataset can be re-analysed without disturbing earlier work.
%
% WHAT IT RUNS, IN ORDER
%   1  coherence_polar_general_260727   per-ROI coherence (kept as a secondary measure)
%   2  cell_link_260727                 join cell identity to recordings + AUDIT
%   3  cell_pool_260727                 lossless per-cell container
%   4  breath_time_overlay_260727       breath diagnostics: window + trigger QC
%   5  breath_time_peth_260727          PRIMARY: inspiration-triggered PETH per cell
%
% CELL MATCHING IS OPTIONAL. If the dataset has no roi_match_results.mat, step 2
% falls back to identity grouping (every ROI its own cell) and the chain still
% completes -- cross-recording pooling is simply off.
%
% ONE-TIME SETUP for a dataset that HAS been reorganised or is new (not part of
% this driver, because it is done once and one of the steps is interactive):
%   cell_cfg_260727('set', datasetPath)      point the matcher at it
%   roi_pair_morph_match_260727              auto-match ROIs across FOVs
%   roi_curation_port_260727(oldMatcherDir)  carry an EXISTING curation across, if
%                                            the recordings merely moved folders.
%                                            Keyed on (folder name, ROI index), so
%                                            it survives any directory reshuffle --
%                                            re-curating is never necessary.
%   roi_review1..4_260727                    curate by hand, if there is none to port
%
% THE PHASE / RAYLEIGH PIPELINE IS NOT RUN. It is retired in favour of the
% absolute-time PETH: cycle-interpolated phase is 12x more densely occupied in
% expiration than inspiration on this preparation, which needs an ECDF correction
% to interpret at all, while the time-domain null is flat with no correction. The
% scripts remain in analysis_260727\phase_rayleigh\ and can still be run by hand.
%
% Runqi Zhang / 2026-07-27

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here,'coh_ca_breath'));
addpath(fullfile(here,'breath_time'));
addpath(fullfile(here,'cell_pair_morph'));
addpath(fullfile(here,'breath_svd'));           % breath trace + peak/onset GUIs
addpath(fileparts(here));                       % repo root

if nargin >= 2
    coh_cfg_260727('set', rootPath, genotype);
end
cfg = coh_cfg_260727();

fprintf('\n################################################################\n');
fprintf('# run_analysis_260727\n');
fprintf('#   dataset : %s\n', cfg.rootPath);
fprintf('#   genotype: %s\n', cfg.genotype);
fprintf('#   output  : %s\n', cfg.outRoot);
fprintf('################################################################\n');
assert(isfolder(cfg.rootPath), 'rootPath not found: %s', cfg.rootPath);
if ~isfolder(cfg.outRoot), mkdir(cfg.outRoot); end

% Preflight first: on a new dataset it says what is present and what is missing,
% so a gap upstream is reported as itself rather than as a file-not-found several
% steps later naming a DERIVED file.
if ~analysis_preflight_260727(cfg.rootPath)
    fprintf(2, ['\nStopping: no recording has both a spike train and breath triggers.\n' ...
                'Fix the items listed above, then re-run.\n']);
    return;
end

steps = { ...
    'coherence_polar_general_260727', 'per-ROI coherence (secondary measure)'
    'cell_link_260727',               'cell identity <-> recordings, audited'
    'cell_pool_260727',               'lossless per-cell container'
    'breath_time_overlay_260727',     'breath diagnostics: window + trigger QC'
    'breath_time_both_triggers',      'PRIMARY: PETH + figures, ONSET-aligned AND PEAK-aligned'
    'breath_time_trigger_compare_260727', 'which landmark is each cell locked to?'
    'cell_metadata_260727',           'master metadata: identity + place + score, per ROI and per cell'
    };

status = strings(size(steps,1),1);
tAll = tic;
for s = 1:size(steps,1)
    fprintf('\n=================== [%d/%d] %s ===================\n', s, size(steps,1), steps{s,1});
    fprintf('%s\n', steps{s,2});
    t0 = tic;
    try
        evalin('base', sprintf('clear %s;', steps{s,1}));   % ensure a clean call
    catch
    end
    try
        feval_or_run(steps{s,1});
        status(s) = sprintf("OK (%.0f s)", toc(t0));
        fprintf('--- %s: %s\n', steps{s,1}, status(s));
    catch ME
        status(s) = "FAILED: " + string(ME.message);
        fprintf(2, '--- %s FAILED: %s\n', steps{s,1}, ME.message);
        for e = 1:min(numel(ME.stack),3)
            fprintf(2, '      at %s line %d\n', ME.stack(e).name, ME.stack(e).line);
        end
        if s <= 3
            fprintf(2, 'Step %d is a prerequisite for what follows -- stopping.\n', s);
            break;
        end
    end
end

fprintf('\n################################################################\n');
fprintf('# summary  (%.1f min total)\n', toc(tAll)/60);
for s = 1:size(steps,1)
    st = status(s); if st == "", st = "not reached"; end
    fprintf('#  %-34s %s\n', steps{s,1}, st);
end
fprintf('# output: %s\n', cfg.outRoot);
fprintf('################################################################\n');
end

function breath_time_both_triggers()
% Zero is either inspiration onset (foot_idx) or the inspiratory peak (peak_idx),
% 267 ms apart. Both are run: whichever gives the tighter latency distribution is
% the landmark that cell is actually locked to. Outputs to breath_time\<trigger>\.
for tg = {'onset','peak'}
    fprintf('\n--- %s-aligned ---\n', upper(tg{1}));
    breath_time_peth_260727(tg{1});
    breath_time_modulation_scatter_260727(tg{1});
    breath_time_peth_percell_260727(tg{1});
end
end

function feval_or_run(name)
% The pipeline mixes functions and plain scripts; call either.
% nargin() throws for a script, which is the cheapest way to tell them apart.
p = which(name);
assert(~isempty(p), 'Not on the MATLAB path: %s', name);
isFcn = true;
try
    nargin(name);
catch
    isFcn = false;
end
if isFcn
    feval(name);
else
    run(p);        % scripts run in this workspace; their `clear` cannot reach the caller
end
end
