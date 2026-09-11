function pool = ensure_pool_260727()
%% ensure_pool_260727  Build cell_pool.mat on demand, instead of demanding it.
% -----------------------------------------------------------------------
% Any script needing the pooled container calls this instead of asserting the file
% exists. If cell_pool.mat is missing it builds the prerequisites itself:
%
%     cell_link_260727   ->  cell_pool_260727   ->  load
%
% and cell_link falls back to identity grouping (every ROI its own cell) when the
% dataset has no cross-FOV matching, so a dataset with few active cells needs no
% curation step at all. You can run any downstream script first and the chain
% assembles behind it.
%
% If the pool CANNOT be built, the error names the real gap rather than the derived
% file. That mattered on 260728_vglut2: the chain reported "cell_pool.mat not
% found", but cell_pool.mat is derived and would have built itself -- the actual
% problem was that no recording had a breathing trace, so there were no triggers to
% align to. analysis_preflight_260727 is run automatically to say so.
%
% Runqi Zhang / 2026-07-29

here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fileparts(here)); addpath(fileparts(fileparts(here)));
cfg = coh_cfg_260727();

if isfile(cfg.poolFile)
    pool = getfield(load(cfg.poolFile,'pool'), 'pool'); %#ok<GFLD>
    return;
end

fprintf('\n[ensure_pool] cell_pool.mat not found -- building it.\n');

% ---- is there anything to build FROM? ----
allMat = dir(fullfile(cfg.rootPath,'**','ca_spike_data.mat'));
if isempty(allMat)
    analysis_preflight_260727(cfg.rootPath);
    error(['No ca_spike_data.mat anywhere under\n  %s\n' ...
           'Run the calcium spike detector first. See the preflight above.'], cfg.rootPath);
end
nBreath = 0;
for k = 1:numel(allMat)
    if isfile(fullfile(allMat(k).folder,'breath_peak_pc1.mat')) && ...
       isfile(fullfile(allMat(k).folder,'breath_insp_start_pc1.mat'))
        nBreath = nBreath + 1;
    end
end
if nBreath == 0
    analysis_preflight_260727(cfg.rootPath);
    % NB: error() with a SINGLE argument does not expand \n -- it needs a format
    % spec, hence the trailing '%s'. Otherwise the message prints literal "\n".
    error(['No recording has breath triggers (breath_peak_pc1.mat +\n' ...
           'breath_insp_start_pc1.mat), so there is nothing to align calcium to.\n' ...
           'The breathing traces must be made first -- see the preflight above.%s'], '');
end

if ~isfile(cfg.linkFile)
    fprintf('[ensure_pool] running cell_link_260727 ...\n');
    cell_link_260727;
end
fprintf('[ensure_pool] running cell_pool_260727 ...\n');
cell_pool_260727;

assert(isfile(cfg.poolFile), 'cell_pool_260727 ran but produced no %s', cfg.poolFile);
pool = getfield(load(cfg.poolFile,'pool'), 'pool'); %#ok<GFLD>
fprintf('[ensure_pool] ready.\n\n');
end
