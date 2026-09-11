function cfg = stitch_cfg_260801(varargin)
%% stitch_cfg_260801  Config for the STITCH-FIRST ROI matching pipeline.
% -----------------------------------------------------------------------
% STANDALONE BY DESIGN. This pipeline shares NOTHING with the appearance-based
% matcher (roi_pair_morph_match_260727 / cell_cfg_260727 / roi_match_out_*):
%   - its own config function (this file), so retargeting one cannot move the other
%   - its own output folder name (roi_stitch_out_260801)
%   - its own file names (fov_offsets.mat, roi_stitch_results.mat, ...)
%   - it never reads roi_match_*.mat, and never writes into roi_match_out_*
% The two can be run on the same dataset and compared without either touching the
% other's products or curation.
%
%   stitch_cfg_260801('set', datasetPath)   point it at a dataset
%   stitch_cfg_260801('reset')              back to the default below
%   stitch_cfg_260801()                     read current settings
%
% THE IDEA. The appearance matcher uses coordinates only as a prefilter and
% decides identity by patch NCC. This one inverts that: you stitch the FOVs
% first, and then two ROIs are the same neuron if they land in the SAME PLACE.
% Appearance becomes a check on the result rather than the thing being trusted.
% That is worth doing here because the tie-point work showed registration good to
% ~2 um, while zero-shift NCC was demonstrably fooled by segmentation differences.

persistent OV
if nargin > 0
    switch lower(string(varargin{1}))
        case "set"
            OV = struct('datasetPath', char(varargin{2}));
            fprintf('stitch_cfg_260801: retargeted to %s\n', OV.datasetPath);
        case "reset"
            OV = [];
            fprintf('stitch_cfg_260801: reset to defaults\n');
    end
end

%% ===================== USER-EDITABLE DEFAULTS =====================
cfg.datasetPath = 'D:\Ventral_surface_summary\Sst\0322';
cfg.outDirName  = 'roi_stitch_out_260801';

% ---- FOV placement ---------------------------------------------------------
% motorPosition is the stitch frame. No anatomical axis map and no vessel
% calibration: this pipeline only needs the FOVs to be mutually consistent WITHIN
% one animal. Whatever anatomical mapping the cartoon applies is a separate,
% later, global step.
%
% ONE correction is applied, and it is measured, not assumed: image ROW runs
% along NEGATIVE stage y. Verified by phase-correlating 15 overlapping pairs and
% scoring all four axis conventions against each other -- median residual
% 44.7 um with row=+y versus 9.1 um with row=-y, and the best pairs fall to
% 1-5 um (the true registration residual). Re-measure with
% scratchpad/test_stage_convention.m if a rig or mounting change is suspected.
cfg.stageYSign = -1;    % +1 = image row along +stage y, -1 = along -stage y

% ---- prep mounting: stage -> anatomy ---------------------------------------
% SECOND, INDEPENDENT layer. The one above is about the RIG (where an image
% column lands in stage space). This one is about how the ANIMAL was placed.
%
% Sessions from 260721 onward were mounted on the new coordinate axis, 90 deg
% from the old one:  lateral = +stage y,  rostral = -stage x
% i.e. "x is y, and +y is -x". Same convention as AXMAP_ROT90 in
% Ventral_surface_ROI_on_cartoon_260729.m, so this pipeline's output and the
% ventral cartoon are read in the same frame.
%
% This is a rigid rotation, so it CANNOT change which ROIs match -- distances are
% invariant under it. It is applied so the map is read in anatomical axes, and so
% that a FOV's footprint is reported on the correct axis (a 320x512 px frame
% spans its long side ROSTRALLY once rotated, not laterally).
%
% Keyed per SESSION, not genotype: how the animal was mounted is a property of
% that day's prep.
cfg.rotatedMount = false;                % Vgat/0730 is post-260721
cfg.AXMAP_ROT90  = [ 0  1 ; -1  0 ];    % [lat; ros] = axmap * [stage x; stage y]
cfg.AXMAP_ID     = [ 1  0 ;  0  1 ];

% ---- which FOV pairs can be stitched by eye --------------------------------
cfg.minOverlapFrac = 0.10;   % footprint intersection / smaller footprint

% Two depth regimes, because they are aligned on different evidence.
%
% CELL regime (|dz| <= maxDz_align_um). Somata recur between the two planes, so
% you can align on the cells themselves and the yellow match count is meaningful.
% 15 um is set by the optics, not by taste: the two-photon axial PSF here is
% ~5-10 um FWHM and somata are ~15 um, so beyond ~15 um of defocus a cell present
% in one plane is largely absent from the other and there is nothing to pair.
cfg.maxDz_align_um = 15;

% VESSEL regime (maxDz_align_um < |dz| <= maxDz_vessel_um). Used ONLY to rescue a
% FOV that would otherwise have no partner at all. Vasculature is continuous in z
% so it still registers, but expect few or no cell pairings -- the GUI says so in
% the title rather than leaving you wondering why nothing lines up.
cfg.maxDz_vessel_um = 30;

% ---- ROI matching in the stitched frame ------------------------------------
% maxDist is the real parameter. Set it from the registration residual, not from
% soma size: two observations of one neuron should coincide to within the stitch
% error (~2 um measured on this rig), so a few um of slack is right and 30+ um
% starts absorbing genuinely different neighbouring cells.
cfg.maxMatchDist_um = 12;
cfg.maxDz_match_um  = 25;    % beyond this the fields are different planes; a
                             %   lateral coincidence is two cells stacked in depth
cfg.ambiguousRatio  = 0.7;   % if the 2nd-nearest candidate is closer than this
                             %   FRACTION of... see roi_stitch_match_260801: a
                             %   match is rejected as ambiguous when the runner-up
                             %   is nearly as close as the winner
cfg.minNCCReport    = 0.30;  % appearance is REPORTED for every accepted match,
                             %   never used to accept or reject. Low values are
                             %   flagged for review, not dropped.
%% ==================================================================

if ~isempty(OV)
    f = fieldnames(OV);
    for i = 1:numel(f), cfg.(f{i}) = OV.(f{i}); end
end
cfg.outDir      = fullfile(cfg.datasetPath, cfg.outDirName);
cfg.offsetFile  = fullfile(cfg.outDir, 'fov_offsets.mat');
cfg.resultsFile = fullfile(cfg.outDir, 'roi_stitch_results.mat');
end
