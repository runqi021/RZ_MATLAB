function recover_vessel_calib_260729()
%% recover_vessel_calib_260729  Rebuild ROI_on_vessel_data.mat from a prior run.
% -----------------------------------------------------------------------
% The cartoon map needs ROI_on_vessel_data.mat -- the stage->anatomy calibration:
%
%   ax, ay            affine, stage um -> stitched-vessel-canvas pixels
%   bx, by            the basilar artery BRANCH POINT in canvas pixels (clicked
%                     by hand; this is the anatomical origin)
%   theta_rad         rotation aligning the vessel axis to rostro-caudal
%   um_per_px_stitch  scale of the stitched canvas
%
% That file is missing from this machine. It is normally produced by
% Ventral_surface_ROI_on_vessel_260528.m, which is INTERACTIVE -- you click the
% branch. But the cartoon script SAVED all six variables into its own output,
% Desktop\ROI_on_cartoon_data.mat, when it last ran successfully. So the
% calibration can be recovered exactly rather than re-clicked.
%
% This is recovery, not re-derivation: the numbers are the ones the previous map
% was drawn with, so the new map lands in the SAME anatomical frame as before,
% which is what you want for the added Sert cells to be comparable.
%
% If you ever DO want to re-derive it (new stitch, better branch click), run
% Ventral_surface_ROI_on_vessel_260528.m and it will overwrite this.
%
% Runqi Zhang / 2026-07-29

SRC = 'C:\Users\Admin\Desktop\ROI_on_cartoon_data.mat';
DST = 'C:\Users\Admin\Desktop\ChAT_MAP\matlab_stitch\ROI_on_vessel_data.mat';

need = {'ax','ay','bx','by','theta_rad','um_per_px_stitch'};
assert(isfile(SRC), ['Cannot recover: %s is missing too.\n' ...
    'Then the only route is to run Ventral_surface_ROI_on_vessel_260528.m by hand.'], SRC);

w = whos('-file', SRC);
missing = setdiff(need, {w.name});
assert(isempty(missing), 'Prior run is missing: %s', strjoin(missing, ', '));

S = load(SRC, need{:});
fprintf('\n===== recovered vessel calibration =====\n');
fprintf('  source: %s\n\n', SRC);
fprintf('  ax               = [%.6g  %.6g]   (stage x -> canvas px)\n', S.ax(1), S.ax(2));
fprintf('  ay               = [%.6g  %.6g]   (stage y -> canvas px)\n', S.ay(1), S.ay(2));
fprintf('  bx, by           = %.2f, %.2f  px  (basilar branch = anatomical origin)\n', S.bx, S.by);
fprintf('  theta_rad        = %.6f  (%.2f deg vessel-axis rotation)\n', S.theta_rad, rad2deg(S.theta_rad));
fprintf('  um_per_px_stitch = %.6g um/px\n', S.um_per_px_stitch);

% sanity: the affine must be invertible and the scale positive
assert(abs(S.ax(1)) > 0 && abs(S.ay(1)) > 0, 'degenerate affine -- refusing to write');
assert(S.um_per_px_stitch > 0, 'non-positive scale -- refusing to write');

d = fileparts(DST);
if ~isfolder(d), mkdir(d); end
if isfile(DST)
    fprintf('\n  %s already exists -- NOT overwriting.\n', DST);
    return;
end
ax = S.ax; ay = S.ay; bx = S.bx; by = S.by;                       %#ok<NASGU>
theta_rad = S.theta_rad; um_per_px_stitch = S.um_per_px_stitch;    %#ok<NASGU>
recovered_from = string(SRC);                                     %#ok<NASGU>
recovered_on   = "2026-07-29";                                    %#ok<NASGU>
save(DST, 'ax','ay','bx','by','theta_rad','um_per_px_stitch', ...
          'recovered_from','recovered_on');
fprintf('\n  wrote %s\n', DST);
fprintf('  now run: Ventral_surface_ROI_on_cartoon_260729.m  (the copy that includes Sert)\n');
end
