function [obsOf, mergeInfo] = pooled_obs_260814(CELL, OBS)
%POOLED_OBS_260814  Pooled observation list, keyed by the cell number in use.
%
%   [obsOf, mergeInfo] = pooled_obs_260814(CELL, OBS)
%
%   obsOf : cell array indexed by the ORIGINAL cell number -- the number that
%           appears in T.cell, in the figure filenames (..._cell%03d.png) and in
%           every CSV. obsOf{c} holds ALL observations of that cell after the
%           260814 merges. For a cell that was merged AWAY (the higher-numbered
%           partner) the entry is empty.
%
% WHY THIS EXISTS
%   Several scripts load pop_features.mat (post-merge) together with
%   event_latency_data.mat (pre-merge) and then do
%
%       c = T.cell(ci);  for o = CELL(c).obs(:)'   % <-- WRONG after merging
%
%   That still runs, and still returns observations, so nothing errors -- but for
%   a merged cell it returns only the lower-numbered constituent's recordings and
%   silently drops the partner's. The cell is the unit of analysis, so a pooled
%   quantity computed that way (IBI, event counts, phase histograms) would be
%   built from half the data for exactly the cells we just took the trouble to
%   merge. Use this instead:
%
%       [obsOf, ~] = pooled_obs_260814(CELL, OBS);
%       c = T.cell(ci);  for o = obsOf{c}(:)'
%
%   Callers that build their own cell list from scratch (per_cell_summary_260812,
%   pop_features_260813) apply the merges directly and do not need this.

[obsMerged, mergeInfo, origCells] = ...
        apply_cell_merges_260814(CELL, OBS, cell_merge_overrides_260814());

obsOf = cell(numel(CELL),1);
for i = 1:numel(obsMerged)
    obsOf{ min(origCells{i}) } = obsMerged{i};
end
end
