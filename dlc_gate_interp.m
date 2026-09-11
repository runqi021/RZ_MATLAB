function M = dlc_gate_interp(csvPath, thr)
%DLC_GATE_INTERP  Read a DLC csv and linear-interpolate low-likelihood points.
%   M = dlc_gate_interp(csvPath, thr) reads a DeepLabCut analyzed csv (3 header
%   lines; columns = [frame, (x,y,likelihood) per bodypart]) and, for every
%   bodypart, replaces (x,y) on frames with likelihood < thr by NaN and then
%   LINEARLY interpolates over frame index (endpoints held = 'nearest').
%   Likelihood columns are returned unchanged, so M is a drop-in replacement for
%   `readmatrix(csvPath,'NumHeaderLines',3)` in the whisker/nose analysis scripts.
%
%   Matches thermal_nostril_breath_single.py gate_interp (LIK_THRESH = 0.6).
%   Default thr = 0.6.

if nargin < 2 || isempty(thr), thr = 0.6; end
M   = readmatrix(csvPath, 'NumHeaderLines', 3);
nbp = floor((size(M,2)-1)/3);                 % bodyparts (frame col + xyl triplets)
for j = 1:nbp
    xc = 3*j-1; yc = 3*j; lc = 3*j+1;         % x, y, likelihood columns
    bad = ~(M(:,lc) >= thr);                  % low-likelihood (NaN likelihood -> bad too)
    M(bad,xc) = NaN; M(bad,yc) = NaN;
    M(:,xc) = fillmissing(M(:,xc), 'linear', 'EndValues', 'nearest');
    M(:,yc) = fillmissing(M(:,yc), 'linear', 'EndValues', 'nearest');
end
end
