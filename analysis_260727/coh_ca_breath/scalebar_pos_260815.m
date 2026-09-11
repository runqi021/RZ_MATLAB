function pos = scalebar_pos_260815(W, H, barLen, corner, thickFrac)
%SCALEBAR_POS_260815  rectangle Position for a scale bar in a chosen corner.
%
%   pos = scalebar_pos_260815(W, H, barLen, corner, thickFrac)
%
% W, H      : crop size in image pixels
% barLen    : bar length in image pixels
% corner    : 'upper left' | 'upper right' | 'lower left' | 'lower right'
%             (default 'lower left', which is where the bar has always been)
% thickFrac : bar thickness as a FRACTION of the crop height (default 0.02).
%             A fraction rather than a pixel count so it keeps its proportions if
%             crop_um or cropPxUm change; with the crop resampled to a fixed
%             pixel size it also comes out identical on every cell.
%
% Returns a 4-element Position for rectangle().
%
% The projection axes is drawn with YDir = 'reverse', so row 1 is the TOP of the
% displayed image and row H is the bottom. "upper" is therefore SMALL y and
% "lower" is LARGE y -- the opposite of the usual axes convention, which is the
% easy thing to get backwards here.
%
% Shared by temporal_phase_cell_fig_260812 (the panel) and
% save_avgproj_png_260812 (the standalone PNG) so the bar cannot end up in
% different corners in the two outputs.

if nargin < 4 || isempty(corner), corner = 'lower left'; end
if nargin < 5 || isempty(thickFrac), thickFrac = 0.02; end
mgn = max(2, round(0.05*H));
thk = max(1, round(thickFrac*H));
barLen = max(1, min(barLen, W - 2*mgn));      % never run off the crop

switch lower(strrep(strtrim(char(corner)), '-', ' '))
    case 'upper left',  x = mgn;              y = mgn;
    case 'upper right', x = W - mgn - barLen; y = mgn;
    case 'lower right', x = W - mgn - barLen; y = H - mgn - thk;
    otherwise,          x = mgn;              y = H - mgn - thk;   % lower left
end
pos = [x y barLen thk];
end
