function cfg = cell_cfg_260727(varargin)
% cell_cfg_260727  Single source of truth for the cross-FOV cell-matching pipeline.
% -----------------------------------------------------------------------
% RETARGETING WITHOUT EDITING THIS FILE:
%   cell_cfg_260727('set', datasetPath)   point the matcher at another folder
%   cell_cfg_260727('set', datasetPath, outDirName)
%   cell_cfg_260727('reset')              back to the defaults below
%   cell_cfg_260727()                     read the current settings
% The override is persistent for the MATLAB session, matching coh_cfg_260727 so a
% driver can retarget BOTH with one call and they cannot disagree about which
% dataset is being analysed.
%
%   cfg.datasetPath   experiment folder holding the recordings
%   cfg.outDirName    subfolder of datasetPath holding the matcher's products
%   cfg.outDir        full path to that subfolder (auto)
%   cfg.resultsFile   roi_match_results.mat   - auto grouping from the matcher (auto)
%   cfg.curatedFile   roi_match_curated.mat   - YOUR curation, from the review GUIs (auto)
%   cfg.patchFile     roi_match_patches.mat   - cached ROI patches (auto)
%
% NOTE. The original 26-recording curation lives at
%   D:\260721_Sert_soma_G8s\phys\roi_match_out_260724
% and is NOT touched by anything here -- the dataset was since split into
% baseline\ , roomAir\ and reoxygenated\ , so each condition gets its own matcher
% output under its own folder.
%
% Runqi Zhang / 2026-07-27
persistent OV
if nargin > 0
    switch lower(string(varargin{1}))
        case "set"
            OV = struct('datasetPath', char(varargin{2}), 'outDirName', '');
            if nargin >= 3, OV.outDirName = char(varargin{3}); end
            fprintf('cell_cfg_260727: retargeted to %s\n', OV.datasetPath);
        case "reset"
            OV = [];  fprintf('cell_cfg_260727: back to defaults\n');
        otherwise
            error('cell_cfg_260727: unknown option "%s"', varargin{1});
    end
end

cfg.datasetPath = 'C:\Users\Admin\Desktop\260829_Sert-soma-g8s_vagotomized';   % <-- edit per experiment
cfg.outDirName  = 'roi_match_out_260828';                     % fresh: the ROI set changed

if ~isempty(OV)
    cfg.datasetPath = OV.datasetPath;
    if ~isempty(OV.outDirName), cfg.outDirName = OV.outDirName; end
end

% derived (leave as-is)
cfg.outDir      = fullfile(cfg.datasetPath, cfg.outDirName);
cfg.resultsFile = fullfile(cfg.outDir, 'roi_match_results.mat');
cfg.curatedFile = fullfile(cfg.outDir, 'roi_match_curated.mat');
cfg.patchFile   = fullfile(cfg.outDir, 'roi_match_patches.mat');
end
