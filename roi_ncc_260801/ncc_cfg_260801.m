function cfg = ncc_cfg_260801(varargin)
%% ncc_cfg_260801  Config for the CROSS-SESSION NCC characterisation pipeline.
% -----------------------------------------------------------------------
% Runs AFTER roi_stitch_260801 has stitched one or more sessions. Pools every
% stitched session it can find and asks one question:
%
%   how well does patch appearance separate a co-located pair (probably the same
%   neuron) from a nearby but non-co-located pair (definitely different neurons)?
%
% Answering it once across sessions gives a threshold with a measured
% false-positive rate, instead of a per-dataset guess.
%
% STANDALONE: its own config, its own output folder, its own file names. It READS
% the stitch pipeline's products (fov_table.mat, fov_offsets.mat) because it is
% downstream of them, and calls roi_stitch_orient_260801 so the image orientation
% cannot drift from the stitcher's. It writes nothing back into roi_stitch_out_*.
%
%   ncc_cfg_260801('set', rootPath)   point it at a folder to scan
%   ncc_cfg_260801('reset')           back to the default below

persistent OV
if nargin > 0
    switch lower(string(varargin{1}))
        case "set"
            OV = struct('scanRoot', char(varargin{2}));
            fprintf('ncc_cfg_260801: scanning under %s\n', OV.scanRoot);
        case "reset"
            OV = [];
            fprintf('ncc_cfg_260801: reset to defaults\n');
    end
end

%% ===================== USER-EDITABLE DEFAULTS =====================
% Every session under here with a roi_stitch_out_260801\fov_table.mat is pooled.
cfg.scanRoot   = 'D:\Ventral_surface_summary\Sert\0721';
cfg.stitchDir  = 'roi_stitch_out_260801';     % produced by the stitch pipeline
cfg.outDir     = 'C:\Users\Admin\Desktop\RZ_MATLAB\roi_ncc_260801\out';

% ---- the two pair populations ----------------------------------------------
% NEAR = co-located across recordings, so probably one neuron seen twice.
% Tightened 2026-08-01 from 10/20 now that the manual stitch is in: with a ~2 um
% registration residual, 5 um is still well outside the stitch error but well
% inside a soma, so it admits real repeats without reaching to neighbours.
cfg.near_xy_um = 5;
cfg.near_z_um  = 15;

% Report every pair scoring above this, in both populations. A fixed cut answers
% "which pairs actually look alike"; the FPR-derived threshold below answers
% "where would a calibrated cut sit". Both are printed -- they are different
% questions and it is worth seeing when they disagree.
cfg.ncc_report_thr = 0.5;

% FAR = near enough to share tissue type, laser power and optics, but far enough
% to be different neurons. This is a LOCAL null: a globally random null would be
% too easy, because it also differs in depth, brightness and region. The question
% that matters is whether appearance can separate a true pair from its NEIGHBOURS.
cfg.far_xy_um  = 150;
cfg.far_z_um   = 50;

% ---- patch extraction -------------------------------------------------------
cfg.patch_um   = 40;    % full width of the compared patch
cfg.patch_px   = 41;    % common resample grid, so different zooms are comparable

% Both patches are cut at their OWN ROI's stitched position, the same rule for
% near and far pairs. Cutting near pairs at their midpoint instead would force the
% far pairs to be centred on one shared location too, where a good stitch would
% make them near-identical and the null would collapse.
cfg.center_mode = 'own';        % 'own' | 'midpoint' (midpoint = near pairs only)

cfg.max_far_per_session = 4000; % cap so a dense session cannot dominate the pool
cfg.fpr_target = 0.01;          % report the threshold at this false-positive rate
%% ==================================================================

if ~isempty(OV)
    f = fieldnames(OV);
    for i = 1:numel(f), cfg.(f{i}) = OV.(f{i}); end
end
cfg.pooledFile = fullfile(cfg.outDir, 'ncc_pooled.mat');
end
