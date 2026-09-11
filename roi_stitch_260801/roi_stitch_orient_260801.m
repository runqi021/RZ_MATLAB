function B = roi_stitch_orient_260801(A, cfg)
%% roi_stitch_orient_260801  Put an image into the stitch frame.
% -----------------------------------------------------------------------
% The scan reports FOV centres and ROI positions in the stitch frame (rig y-sign
% + prep rotation). Image PIXELS have to be carried through the same transform,
% or the picture and the coordinates disagree -- tiles are drawn in the wrong
% orientation while their ROI markers sit correctly, and any registration against
% them is meaningless.
%
% Every consumer (map, overview, auto-register) calls THIS, so the transform is
% defined in exactly one place and they cannot drift apart.
%
% DERIVATION. A pixel (row, col) of the raw image sits at stage
%     stage_x = cx + (col - W/2)*px
%     stage_y = cy + sy*(row - H/2)*px          sy = cfg.stageYSign
% and the frame is [X; Y] = axmap * [stage_x; stage_y].
%
% Rotated mount, axmap = [0 1; -1 0]:
%     X =  stage_y = const - row*px     -> X DECREASES with row
%     Y = -stage_x = const - col*px     -> Y DECREASES with col
% so the output array needs column <-> reversed row and row <-> reversed col:
%     B = rot90(A.', 2)                 (W rows x H cols -- the footprint transposes)
%
% Standard mount, axmap = I:
%     X = stage_x  -> increases with col   (already correct)
%     Y = stage_y  = const - row*px        -> decreases with row
%     B = flipud(A)
%
% Output convention: B(r,c) sits at X = Xmin + (c-0.5)*px, Y = Ymin + (r-0.5)*px,
% i.e. both indices increase with the frame axes.

if nargin < 2, cfg = stitch_cfg_260801(); end

if cfg.stageYSign > 0
    % row already runs along +stage y; nothing to undo on that axis
    if cfg.rotatedMount, B = A.'; else, B = A; end
    if ~cfg.rotatedMount, return; end
    B = fliplr(B);
    return;
end

if cfg.rotatedMount
    B = rot90(A.', 2);      % footprint transposes: H x W -> W x H
else
    B = flipud(A);
end
end
