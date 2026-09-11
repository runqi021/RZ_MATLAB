function cfg = coh_cfg()
% coh_cfg  Single source of truth for the general breath x Ca coherence pipeline.
% -----------------------------------------------------------------------
% Every script in coh_ca_breath_260722/ reads these settings, so to point the
% whole pipeline at a new experiment you edit ONLY this file.
%
%   cfg.rootPath        flat experiment folder (recordings sit directly inside)
%   cfg.genotype        label applied to every recording (e.g. 'Sert','ChAT')
%   cfg.genotype_color  marker color for that genotype
%   cfg.cohDir          foundation output dir (auto)
%   cfg.cohData         coherence_polar_data.mat path the downstream scripts read (auto)
% -----------------------------------------------------------------------

cfg.rootPath       = 'D:\260721_Sert_soma_G8s\phys';   % <-- edit per experiment
cfg.genotype       = 'Sert';                            % <-- edit per experiment
cfg.genotype_color = [0.90 0.45 0.10];                  % marker color for this genotype

% derived (leave as-is)
cfg.cohDir  = fullfile(cfg.rootPath, sprintf('coherence_polar_%s_260722', lower(cfg.genotype)));
cfg.cohData = fullfile(cfg.cohDir, 'coherence_polar_data.mat');
end
