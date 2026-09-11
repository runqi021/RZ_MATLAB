function cfg = ncc_round_cfg_260801(varargin)
%% ncc_round_cfg_260801  Which dataset the stitch-seeded round 1-4 curation works on.
% -----------------------------------------------------------------------
% ONE place to retarget the whole four-round workflow. The copied review GUIs
% originally called cell_cfg_260727 -- the appearance matcher's config -- which
% would have meant retargeting one pipeline silently moved the other. They now
% call this instead, so the two are genuinely independent.
%
%   ncc_round_cfg_260801('set', datasetPath)
%   ncc_round_cfg_260801('reset')
%
% Run order:
%   roi_ncc_build_rounds_260801    seed groups from the map stitch
%   roi_ncc_round1_toss_260801     toss bad observations
%   roi_ncc_round2_assign_260801   re-home the tossed ones
%   roi_ncc_round3_merge_260801    merge same-cell splits
%   roi_ncc_round4_ungrouped_260801  place never-grouped observations

persistent OV
if nargin > 0
    switch lower(string(varargin{1}))
        case "set"
            OV = struct('datasetPath', char(varargin{2}));
            fprintf('ncc_round_cfg_260801: retargeted to %s\n', OV.datasetPath);
        case "reset"
            OV = [];
            fprintf('ncc_round_cfg_260801: reset to defaults\n');
    end
end

%% ===================== USER-EDITABLE DEFAULTS =====================
cfg.datasetPath   = 'D:\Ventral_surface_summary\Vglut2\0728';
cfg.outDirName    = 'roi_ncc_out_260801';    % NOT roi_match_out_* -- separate curation
cfg.stitchDirName = 'roi_stitch_out_260801'; % where the map stitch wrote its offsets
%% ==================================================================

if ~isempty(OV)
    f = fieldnames(OV);
    for i = 1:numel(f), cfg.(f{i}) = OV.(f{i}); end
end
cfg.outDir    = fullfile(cfg.datasetPath, cfg.outDirName);
cfg.stitchDir = fullfile(cfg.datasetPath, cfg.stitchDirName);
end
