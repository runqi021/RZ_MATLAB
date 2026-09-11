function cfg = coh_cfg_260727(varargin)
% RETARGETING WITHOUT EDITING THIS FILE (added 260727):
%   coh_cfg_260727('set', rootPath, genotype)  point the whole pipeline elsewhere
%   coh_cfg_260727('reset')                    go back to the defaults below
%   coh_cfg_260727()                           read the current settings
% The override is persistent for the MATLAB session, so a driver script can set it
% once and every downstream script picks it up with no argument passing. Use this
% instead of editing the two lines below back and forth between datasets -- that
% is how outputs end up in the wrong folder.
persistent OV
if nargin > 0
    switch lower(string(varargin{1}))
        case "set"
            OV = struct('rootPath',char(varargin{2}), ...
                        'genotype',char(varargin{3}));
            fprintf('coh_cfg_260727: retargeted to %s (%s)\n', OV.rootPath, OV.genotype);
        case "reset"
            OV = [];  fprintf('coh_cfg_260727: back to defaults\n');
        otherwise
            error('coh_cfg_260727: unknown option "%s"', varargin{1});
    end
end
% coh_cfg_260727  Single source of truth for the breath x Ca coherence pipeline
%                 AND for its link to the cross-FOV cell matching.
% -----------------------------------------------------------------------
% Every script in analysis_260727\coh_ca_breath\ reads these settings, so to point
% the whole pipeline at a new experiment you edit ONLY this file.
%
%   cfg.rootPath        flat experiment folder (recordings sit directly inside)
%   cfg.genotype        label applied to every recording (e.g. 'Sert','ChAT')
%   cfg.genotype_color  marker color for that genotype
%   cfg.cohDir          foundation output dir (auto)
%   cfg.cohData         coherence_polar_data.mat the downstream scripts read (auto)
%
% Cell-matching link (new in 260727):
%   cfg.cellOutDir      the matcher's output dir, holding roi_match_curated.mat
%   cfg.curatedFile     your curation (auto)
%   cfg.resultsFile     the matcher's auto grouping (auto)
%   cfg.linkFile        cell_link.mat  written by cell_link_260727.m  (auto)
%   cfg.poolFile        cell_pool.mat  written by cell_pool_260727.m  (auto)
%
% RENAMED from coh_cfg to coh_cfg_260727 on purpose: the original
% coh_ca_breath_260722\coh_cfg.m still exists, and two files of the same name on
% the MATLAB path shadow each other depending on run order. Distinct names mean
% the 260727 pipeline can never silently read the 260722 config, and cohDir is
% tagged _260727 so it can never overwrite the 260722 results either.
%
% Runqi Zhang / 2026-07-27

cfg.rootPath       = 'D:\260728_vglut2_soma-g8s\phys';  % <-- edit per experiment
cfg.genotype       = 'Vglut2';                            % <-- edit per experiment
cfg.genotype_color = [0.90 0.45 0.10];                  % marker color for this genotype

if ~isempty(OV)
    cfg.rootPath = OV.rootPath;
    cfg.genotype = OV.genotype;
end

% ---- recordings excluded from EVERY analysis, by folder name ----
% Not a quality gate -- these are recordings whose data cannot be placed in the
% same coordinate frame as the rest, so including them would corrupt any spatial
% or pooled result.
%
% ro1_2x_x1018y1130_18lp_3000f_00001  260721_Sert, the FIRST recorded FOV. It was
%   acquired on a DIFFERENT stage zero reference, so its motor coordinates are not
%   comparable with the rest of the session. Confirmed three ways: acquired 15:24,
%   20 min before the next recording; the only recording in the session with no
%   Basler .avi (the camera was not running yet); and its motorPosition y (1018.3)
%   uniquely matches its own folder-name x while the later ro1_* recordings at the
%   nominally same site do not. It was already in the matcher's tossFolders; this
%   extends the exclusion to the analysis and to the archive.
%
% DEPTH BAR, added 2026-08-08. This archive is a VENTRAL SURFACE analysis, so the
% Vgat-G8m shiverer session (0730) -- the only one imaged deliberately deep -- is
% cut at a hard 120 um. Every 0730 FOV deeper than that is excluded from the
% analysis, the map and the polar plot.
%
% Depth is |motorPosition(3)|, NOT the z in the filename. They disagree once, and
% it matters: roi1_1.5x_x1050y1000_z100_3000f_26.5lp is named z100 but its motor
% reports 201.9 um. motorPosition is the authority (see reference_stage_axis /
% the 260728 vessel series), so that FOV is excluded despite its name.
%
% roi2_2x_x1050y2250_z120 sits at 121.6 um -- 1.6 um past the bar. Excluded,
% because the bar is on the measured depth; flip it back if you meant the
% nominal 120 in the filename.
%
% COST, stated plainly: two of these six carry Vgat's strongest phase-locked
% cells (z265 gives logZ 3.79 / 3.58 / 3.44, z250 gives 3.75). Applying the bar
% removes them, so Vgat's top of the polar plot changes. That is the intent of a
% surface-only criterion, not a side effect to be surprised by later.
cfg.excludeRecordings = [ ...
    "ro1_2x_x1018y1130_18lp_3000f_00001" ...
    "roi1_2.4x_x1300y900_z250_3000f_30lp_00001"    ... % 244.7 um
    "roi1_3x_x1350y850_z265_3000f_31lp_00001"      ... % 264.7 um
    "roi1_1.5x_x1050y1000_z100_3000f_26.5lp_00001" ... % 201.9 um (name says z100)
    "roi1_2.4x_x1000y1000_z320_3000f_32lp_00001"   ... % 314.7 um
    "roi2_2x_x1050y2250_z120_2000f_14lp_00001"     ... % 121.6 um (at the bar)
    "roi2_2x_x1050y2250_z160_2000f_16lp_00001"     ... % 161.6 um
    ];

% the matcher's output folder (must agree with cell_cfg_260727.m)
cfg.cellOutDirName = 'roi_match_out_260727';

% EVERY figure and data file this pipeline writes goes under one folder:
cfg.outRoot = fullfile(cfg.rootPath, 'analysis_260727');

% derived (leave as-is)
cfg.cohDir   = fullfile(cfg.outRoot, 'coherence_polar');   % foundation (per-ROI)
cfg.cellDir  = fullfile(cfg.outRoot, 'cell_pooled');       % cell-aware (link/pool/pooled)
cfg.phaseDir = fullfile(cfg.outRoot, 'phase_rayleigh');    % Rayleigh phase-modulation pipeline
cfg.cohData = fullfile(cfg.cohDir, 'coherence_polar_data.mat');

cfg.cellOutDir  = fullfile(cfg.rootPath, cfg.cellOutDirName);
cfg.curatedFile = fullfile(cfg.cellOutDir, 'roi_match_curated.mat');
cfg.resultsFile = fullfile(cfg.cellOutDir, 'roi_match_results.mat');
cfg.linkFile    = fullfile(cfg.cellDir, 'cell_link.mat');
cfg.poolFile    = fullfile(cfg.cellDir, 'cell_pool.mat');

% per-script figure folders, all under outRoot
cfg.figDir = @(name) fullfile(cfg.outRoot, name);
end
