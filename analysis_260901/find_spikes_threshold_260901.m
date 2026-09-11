function locs = find_spikes_threshold_260901(trace, fps, threshold, minDistS, minWS, minProm)
%FIND_SPIKES_THRESHOLD_260901  Detect calcium spikes using MATLAB findpeaks.
%   Uses MinPeakHeight, MinPeakDistance, MinPeakWidth, MinPeakProminence.
%
%   VERBATIM COPY of the local function find_spikes_threshold() at the bottom of
%   calcium_spike_gui.m -- only the name is different, because a local function
%   in a script file cannot be called from anywhere else. The body below is
%   byte-for-byte the old GUI's, so ca_recheck_gui_260901 detects with the same
%   code and not merely the same idea. If the old GUI's version ever changes,
%   change this one to match it.

trace = trace(:);
locs  = [];
if length(trace) < 3, return; end

minDistFr = max(1, round(minDistS * fps));
minWFr    = max(0, round(minWS * fps));

pkArgs = {'MinPeakHeight',    threshold, ...
          'MinPeakDistance',  minDistFr};
if minWFr > 0
    pkArgs = [pkArgs, {'MinPeakWidth', minWFr}];
end
if minProm > 0
    pkArgs = [pkArgs, {'MinPeakProminence', minProm}];
end

[~, locs] = findpeaks(trace, pkArgs{:});
locs = locs(:);
end  % find_spikes_threshold_260901
