function ok = analysis_preflight_260727(rootPath)
%% analysis_preflight_260727  What does this dataset have, and what is it missing?
% -----------------------------------------------------------------------
%   analysis_preflight_260727()            % the dataset in coh_cfg_260727
%   analysis_preflight_260727(rootPath)    % any other folder
%   ok = analysis_preflight_260727(...)    % true if the chain can run at all
%
% Run this FIRST on any new dataset. It scans every recording, reports what is
% present and what is missing, and says which upstream step to run next -- instead
% of the chain failing several steps later on a file-not-found that names a
% derived file rather than the real gap.
%
% WHY THIS EXISTS. On 260728_vglut2 the chain stopped with
%   "cell_pool.mat not found"
% but cell_pool.mat is DERIVED and would have built itself. The actual problem was
% that no recording had a breathing trace at all: the Basler .avi files were there
% but had never been through breath_svd_pc1 + the peak/trough GUIs, so there were
% no triggers to align anything to. The error named the symptom, not the cause.
%
% ------------------------- WHAT EACH STEP NEEDS --------------------------
% per recording folder:
%   *_meta.mat                 fps + zoom            (dffQC pipeline)
%   *_dFF.mat                  traces                (dffQC pipeline)
%   *_cpSAM_output.mat         masks; ONLY needed if you want cross-FOV cell
%                              matching. Without it the chain uses identity
%                              grouping -- every ROI is its own cell -- which is
%                              the right choice when there are few active cells.
%   ca_spike_data.mat          detected events       (calcium spike detector)
%   breath_pc1.mat             breathing trace       (breath_svd\, from cam1*.avi)
%   breath_peak_pc1.mat        inspiratory peaks     (peak GUI)     <- TRIGGER
%   breath_insp_start_pc1.mat  inspiration onsets    (trough GUI)   <- TRIGGER
%
% A recording contributes to the breath x Ca analysis only with ca_spike_data AND
% both breath trigger files. Everything else is optional or derived.
%
% Output: <rootPath>\analysis_260727\preflight.csv
%
% Runqi Zhang / 2026-07-29

here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here,'coh_ca_breath')); addpath(fileparts(here));
addpath(fullfile(here,'breath_svd'));           % breath trace + peak/onset GUIs
cfg = coh_cfg_260727();
if nargin >= 1 && ~isempty(rootPath), cfg.rootPath = char(rootPath); end

fprintf('\n============ analysis_preflight_260727 ============\n');
fprintf('dataset: %s\n', cfg.rootPath);
assert(isfolder(cfg.rootPath), 'rootPath not found: %s', cfg.rootPath);

%% ---- find recording folders: anything holding a meta / dFF / spike file ----
d = dir(cfg.rootPath);  d = d([d.isdir]);
d = d(~ismember({d.name}, {'.','..'}));
d = d(~startsWith({d.name}, {'analysis_','roi_match_','coherence_','_'}));
n = numel(d);
assert(n > 0, 'No recording folders under %s', cfg.rootPath);

name = strings(n,1);
has = false(n,7);   % meta dFF cpSAM spikes breathPC1 breathPeak breathStart
avi = zeros(n,1);
for i = 1:n
    p = fullfile(cfg.rootPath, d(i).name);
    name(i)  = string(d(i).name);
    has(i,1) = ~isempty(dir(fullfile(p,'*_meta.mat')));
    has(i,2) = ~isempty(dir(fullfile(p,'*_dFF.mat')));
    has(i,3) = ~isempty(dir(fullfile(p,'*_cpSAM_output.mat')));
    has(i,4) = isfile(fullfile(p,'ca_spike_data.mat'));
    has(i,5) = isfile(fullfile(p,'breath_pc1.mat'));
    has(i,6) = isfile(fullfile(p,'breath_peak_pc1.mat'));
    has(i,7) = isfile(fullfile(p,'breath_insp_start_pc1.mat'));
    avi(i)   = numel(dir(fullfile(p,'*.avi')));
end
ready = has(:,4) & has(:,6) & has(:,7);

%% ---- report ----
lbl = {'meta','dFF','cpSAM','ca_spikes','breath_pc1','breath_peak','breath_onset'};
fprintf('\n%d recording folders\n', n);
fprintf('%-16s %s\n', 'artifact', 'present');
for k = 1:numel(lbl)
    fprintf('  %-14s %3d / %-3d %s\n', lbl{k}, nnz(has(:,k)), n, bar(nnz(has(:,k)), n));
end
fprintf('  %-14s %3d / %-3d %s\n', 'cam .avi', nnz(avi>0), n, bar(nnz(avi>0), n));
fprintf('\n  READY for breath x Ca (spikes + both triggers): %d of %d\n', nnz(ready), n);

%% ---- what to do next ----
fprintf('\n---- what is blocking ----\n');
todo = strings(0,1);
if nnz(has(:,4)) < n
    m = nnz(~has(:,4));
    fprintf('  %d recording(s) have no ca_spike_data.mat\n', m);
    fprintf('      -> run the calcium spike detector on them\n');
    todo(end+1) = "ca_spike_data missing in " + m;
end
if nnz(has(:,6) & has(:,7)) < n
    m = nnz(~(has(:,6) & has(:,7)));
    fprintf('  %d recording(s) have no breath TRIGGERS\n', m);
    if nnz(has(:,5)) == 0 && nnz(avi>0) > 0
        fprintf('      the cam1*.avi files are present but were never processed:\n');
        fprintf('      all three live in analysis_260727\\breath_svd\\ (see its README)\n');
        fprintf('      -> 1. breath_svd_pc1              avi  -> breath_pc1.mat\n');
        fprintf('      -> 2. breathing_peak_gui_pc1      -> breath_peak_pc1.mat\n');
        fprintf('      -> 3. breathing_trough_gui_pc1    -> breath_insp_start_pc1.mat\n');
    elseif nnz(has(:,5)) > 0
        fprintf('      breath_pc1.mat exists but the peak/onset GUIs have not been run\n');
    else
        fprintf('      and no .avi either -- there is no breathing data for this dataset\n');
    end
    todo(end+1) = "breath triggers missing in " + m;
end
if nnz(has(:,3)) < n
    fprintf('  %d recording(s) have no cpSAM output -- cross-FOV cell matching is\n', nnz(~has(:,3)));
    fprintf('      unavailable for those. NOT a blocker: the chain falls back to\n');
    fprintf('      identity grouping (each ROI its own cell) automatically.\n');
end
if isempty(todo)
    fprintf('  nothing. run_analysis_260727() will run end to end.\n');
end

ok = nnz(ready) > 0;
if ~ok
    fprintf(['\n  ==> The breath x Ca chain CANNOT run yet: no recording has both a\n' ...
             '      spike train and breath triggers. Fix the items above first.\n']);
else
    fprintf('\n  ==> run_analysis_260727() will run, using %d of %d recordings.\n', nnz(ready), n);
end

%% ---- save ----
outRoot = fullfile(cfg.rootPath,'analysis_260727');
if ~isfolder(outRoot), mkdir(outRoot); end
T = table(name, has(:,1), has(:,2), has(:,3), has(:,4), has(:,5), has(:,6), has(:,7), avi, ready, ...
    'VariableNames',[{'rec_name'}, lbl, {'n_avi','ready'}]);
writetable(T, fullfile(outRoot,'preflight.csv'));
fprintf('\nSaved %s\n', fullfile(outRoot,'preflight.csv'));
end

function s = bar(k, n)
w = 24; f = round(w*k/max(n,1));
s = ['[' repmat('#',1,f) repmat('.',1,w-f) ']'];
end
