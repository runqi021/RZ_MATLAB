function G = cell_merge_overrides_260814()
%CELL_MERGE_OVERRIDES_260814  Hand-curated cell identity corrections.
%
%   G = cell_merge_overrides_260814()
%
% Returns a cell array. Each element is a GROUP: a list of (recording name, ROI
% index) pairs that are the SAME physical cell but were left unmerged by
% cell_link. Anything listed together is pooled into one cell downstream.
%
% WHY THIS FILE EXISTS, AND WHY IT IS KEYED THIS WAY
%   cell_link matches ROIs across FOVs using stage position as a prefilter and
%   morphology as the metric. The gate is deliberately tight -- a loose one has
%   already produced a 236 um false merge -- so genuine matches across FOVs that
%   differ in zoom or depth can fall outside it and never get compared.
%   Corrections therefore have to be made by eye, and this is where they live.
%
%   The key is (recName, roi), NOT a cell index: cell numbering is regenerated
%   whenever the registry is rebuilt, whereas the recording folder name plus the
%   ROI index inside that recording is stable across reshuffles.
%
%   Nothing in the archive is modified. Scripts call this and apply it after
%   loading event_latency_data.mat, so the correction is visible, versioned and
%   reversible, rather than baked into a .mat no one can audit.
%
% TO ADD A MERGE
%   append a group: { {'<recName>', <roi>}, {'<recName>', <roi>}, ... }
%   Record WHY in the comment above it -- an identity claim with no stated
%   evidence cannot be checked by anyone later.

G = {};

% 2026-08-14, RZ, by eye from the avg-projection crops: the same cell appears in
% both pFN FOVs of Vglut2/0224. Not merged automatically because the two
% recordings sit 83 um apart in stage x/y, 5 um apart in z, and were imaged at
% different zoom (6x vs 3x), which puts the pair outside cell_link's position
% prefilter. Was cell 25 (60 events, logZ 3.54) and cell 29 (109 events, logZ 3.84).
G{end+1} = { {'pFN_roi1_z0_512x512_6x_2000f_00001',  9}, ...
             {'pFN_roi2_z5_25x_512x512_3x_6000f_00001', 16} };

% 2026-08-14, RZ, by eye from the avg-projection crops: same cell in both Sst/0807
% roi1 FOVs, imaged at 2.5x and 5x at the same depth (z5). NOTE the two stage
% positions differ by 161 um -- larger than the 83 um of the Vglut2 pair above,
% and within the range where a loose gate has previously produced a FALSE merge
% (236 um). This merge rests on the visual match, not on position.
% Was cell 218 (43 events, logZ 0.91) and cell 220 (56 events, logZ 0.07).
G{end+1} = { {'roi1_z5_2.5x_3000f_lp15_00001',   2}, ...
             {'roi1_z5_5x_6000f_lp15.5_00001',   1} };

% 2026-08-14, RZ, by eye from the avg-projection crops: same cell in two Vglut2/0810
% roi1 FOVs (5x at z-36, 7x at z-30). Checked against the matcher output for this
% session: both ROIs ARE present (observations 120 and 153 of 160) and both were
% left UNGROUPED by the curation, so this is a genuine matcher miss and not a
% missing cell_link -- 0810 does have cell_link.mat.
% Was cell 48 and cell 52.
G{end+1} = { {'roi1_5x_z-36_24lp_3000f_00001',  7}, ...
             {'roi1_7x_z-30_21lp_3000f_00001',  6} };

% ---------------------------------------------------------------------------
% Recovered from curation that was never joined
%
% Vglut2/0224 and Vglut2/1124 both have a curated cross-FOV matcher output
% (roi_match_curated.mat) but no cell_link.mat -- the matching was done and the
% curation exists, it simply was never joined into the cell registry. The groups
% below are read straight out of those curations, so they carry the same
% authority as any other linked session; they are NOT fresh by-eye claims.
%
% Transcribed rather than fixed upstream because re-running cell_link renumbers
% every cell in the registry, which would invalidate every figure filename and
% the pop_features table. The override key is stable across rebuilds, so this
% reaches the same grouping without the renumber. If cell_link is ever run for
% these two sessions, these groups become no-ops ("already merged") rather than
% double-counting.
%
% Verified 2026-08-14 against roi_match_curated.mat in each session.
% ---------------------------------------------------------------------------

% Vglut2/0224, curated group 14 (was cells 24 and 27)
G{end+1} = { {'pFN_roi1_z0_512x512_6x_2000f_00001',        6}, ...
             {'pFN_roi2_z5_25x_512x512_3x_6000f_00001',    7} };

% Vglut2/1124, curated group 2 (was cells 178 and 183)
G{end+1} = { {'roi5_1400-1230-0_x4.4_15lp_6000f_00001',              5}, ...
             {'roi5_1400-1300-0_18lp_930_x4_512x512_6000f_00002',   12} };

% Vglut2/1124, curated group 3 (was cells 179 and 184)
G{end+1} = { {'roi5_1400-1230-0_x4.4_15lp_6000f_00001',              7}, ...
             {'roi5_1400-1300-0_18lp_930_x4_512x512_6000f_00002',   14} };

% Vglut2/1124, curated group 4 (was cells 176 and 182)
G{end+1} = { {'roi5_1400-1230-0_x4.4_15lp_6000f_00001',              3}, ...
             {'roi5_1400-1300-0_18lp_930_x4_512x512_6000f_00002',   10} };

% NOT INCLUDED, needs a convention decision:
%   Four further curated groups pair an ACTIVE cell with an ROI that did not
%   pass the >5-event active gate, so it sits in OBS but in no CELL:
%       Vglut2/0224 grp  9  -> cell 26     Vglut2/0224 grp 16 -> cell 28
%       Vglut2/0224 grp 17  -> cell 30     Vglut2/1124 grp  9 -> cell 177
%   Pooling those in would not change the cell count but WOULD add their events
%   to the cell, shifting its pooled event count, logZ and permutation p. That
%   is a real change to the numbers, so it is left out until decided.

end

% ---------------------------------------------------------------------------
% KNOWN HAZARD FOR ANYONE EXTENDING THIS FILE
%   The (recName, roi) key is taken from OBS.label in event_latency_data.mat,
%   which stores the FULL recording folder name -- that is unambiguous and is
%   what this file relies on.
%
%   Do NOT key anything off the fov_name field of roi_match_results.mat. Those
%   names were produced with fileparts() on folder names, so everything after
%   the first dot is dropped as if it were a file extension:
%       roi1_z5_2.5x_3000f_lp15_00001    -> roi1_z5_2
%       fov1_1.7x_y1650x1320_z30_...     -> fov1_1
%   In Vglut2/0728 this collapses three distinct z-planes (z30, z35, z40) onto
%   the single name "fov1_1", and in Vgat/0730 it collapses a cell1 recording
%   and a quiet recording onto "roi1_1". Use match.fov_folder (full path) there.
% ---------------------------------------------------------------------------
