% ventral_map_genotype_260817.m   (copy of ventral_map_selected_260816)
% -----------------------------------------------------------------------
%  SAME MAP, RECODED: colour is GENOTYPE, marker is RESPIRATORY CLASS.
%      purple Sst | orange Sert | red ChAT | blue Vgat | lemon green Vglut2
%      star pre-I | triangle I | square post-I | OPEN circle tonic
%  Splitting the two lets you read "which genotype" and "which class" off one
%  dot independently; the 260816 map folded both into a single colour per
%  (genotype, class) pair, so Vglut2 pre-I and Vglut2 I were two unrelated
%  greens and nothing said they were the same genotype.
%
%  VENTRAL SURFACE ONLY. Cells whose observations all lie in depth-barred FOVs
%  are dropped from the dots, not just from the anatomy. The 260816 map drew them
%  deliberately ("depth is not a reason to leave a selected cell off the map");
%  this one is a ventral-surface figure, so they go. In practice that is the
%  three deep Vgat I cells, which is why Vgat here is the tonic/null population.
%
%  Axes cropped to x [-2500 2500], y [-1000 3000].
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
%  SAME MAP, DIFFERENT RIGHT SIDE. The anatomy (left hemisphere: FOV patches,
%  tints, stitch offsets, site merge) is the 260808 placement code unaltered --
%  verified against its saved F: the same 61 patches, centres equal to 0 um. Not
%  pixel-identical, and cannot be: the soft blob mask is built from an unseeded
%  randn, so its edges differ between any two runs of either script. The white
%  sig-ROI outlines the 260808 map drew are GONE -- the significant cells are
%  shown as circles on the right, and marking them twice was clutter.
%
%  Right-hemisphere dots, one per CELL, filled, coloured by GROUP using the same
%  table as the population-average overlays (group_display_colors_260816):
%      IO      not shown at all -- 107 of the 147, in one dense clump
%      Sert    significant only
%      others  the full selected set
%
%  The selection is the 147 cells of polar_selected_260816.mat -- the same cells
%  the polar plot and the group overlays describe -- so all three figures now
%  talk about one population.
%
%  ONE DOT PER CELL, not per mask. A cell seen in several FOVs used to paint one
%  dot per observation, which drew the same neuron two or three times; here the
%  observations are averaged to a single position and the spread is printed as a
%  check.
%
%  DEPTH BAR: coh_cfg's 120 um exclusion still suppresses the anatomy patch for
%  the two deep Vgat recordings, so the left side is unchanged -- but their cells
%  (277/278/279, three of the four Vgat I) DO get dots, because depth is not a
%  reason to leave a selected cell off the map. Flagged at run time.
% -----------------------------------------------------------------------
%  FOV patches (per-FOV AVG image, circular soft-mask, group tint, white
%  outline on sig ROIs) overlaid on the brainstem cartoon background.
%  Coordinate convention = same as the working dot version:
%      data x = +lateral  (μm)
%      data y = +rostral  (μm)         y-up convention (YDir='normal')
%  ROIs from every ca_spike_data.mat under <sumRoot>\{groups} are computed
%  in μm relative to the branch via the vessel-derived calibration.
% -----------------------------------------------------------------------

clear; clc;   % NOT close all -- the polar figure is meant to stay open beside this one

%% ===================== USER-EDITABLE =====================
cartoonMat = 'C:\Users\Admin\Desktop\brainstem_map_cartoon_coordframe.mat';
mouseRoot  = 'C:\Users\Admin\Desktop\ChAT_MAP';
% REPOINTED. The old path (Desktop\260522_ChAT_soma_G8s) no longer exists -- that
% data moved to D:\260522_ChAT_soma_G8s, whose matlab_stitch holds only
% stitched_ch3_avgz.tif. ChAT_MAP\matlab_stitch is the only place with BOTH
% stitch_manual_coords.mat and stitched_ch3_avgproj.tif, which is what
% Ventral_surface_ROI_on_vessel_260528.m needs to regenerate the calibration.
%
% BLOCKER: ROI_on_vessel_data.mat does not exist anywhere on this machine, and it is
% what converts stage coordinates to cartoon-canvas coordinates. It is produced by
% Ventral_surface_ROI_on_vessel_260528.m, which is INTERACTIVE (you click the vessel
% branch to set the reference). So that has to be run once, by hand, before this
% script can draw anything -- including the Sert cells added below.
%
% Why the exclusion in coh_cfg_260727 matters here: stage_to_canvas is ONE affine
% applied to every animal, so the whole map assumes every session zeroed its stage
% on the same anatomical landmark. The first 260721_Sert FOV did not, which is
% exactly why it must stay excluded from anything spatial.
vesselMat  = fullfile(mouseRoot,'matlab_stitch','ROI_on_vessel_data.mat');

sumRoot    = 'D:\Ventral_surface_summary';
% coherence/significance source = SVD-motion pipeline output (was the older
% coherence_polar_260528). Same schema (PP.r, PP.th, labels, confC) and same
% label format 'group/date/recName/rid', so this is a drop-in source swap.
dataMat    = fullfile(sumRoot,'polar_coh_vs_rayleigh_260808','polar_coh_vs_rayleigh_data.mat');
% REPOINTED 2026-08-08 to the run that (a) includes Sst 0806/0807 and (b) applies
% the 120 um depth bar. Leaving it on 260729 would have drawn phase dots from a
% cell set that no longer matches the FOVs this map draws.
% REPOINTED from the coherence run to the RAYLEIGH-on-linear-phase run.
%   significance : logZ >= logZ_crit, i.e. alpha = 0.05   (was coherence r >= confC
%                  at alpha = 0.001 -- only 14 of 368 ROIs cleared it)
%   phase        : the occupancy-weighted circular MEDIAN of the cell's event
%                  phases (RY.thMed), on the same convention the colorbar already
%                  uses -- 0 = insp onset, pi = insp peak.
% Same linear peak<->trough phase, same 100 ms lead correction, so the map and the
% polar plot now say the same thing about the same cells.
%
% NOTE the Rayleigh result is per CELL while the map draws per ROI observation, so
% obsOfCell is used to hand every observation of a cell that cell's phase and
% verdict. A cell imaged in three recordings therefore paints three FOVs the same
% colour, which is the point.

% ---- cell identity: report CELLS, not raw Cellpose masks --------------------
% A mask is not a cell. The same neuron re-imaged at a second Z or a second zoom
% produces one mask per recording, so summing masks over recordings overcounts.
% The cross-FOV matcher resolves that, and its cell_link.mat carries the mapping
%     (recording name, maskL label) -> cell_id
% which is exactly the portable key the curation is built on. Registered per
% session below.
%
% Two things the matcher output encodes that must be honoured:
%   status 'tossed'  (cell_id NaN)  -> rejected in curation, NOT a cell. Dropped.
%   status 'ungrouped'              -> a real cell, just seen in one recording.
%                                      Counts as one cell.
% Sessions with no matcher run fall back to one cell per ROI, which is what
% cell_link_260727 itself does without curation -- so the count degrades
% gracefully rather than silently reverting to mask counts.
% Read from the ARCHIVE's own copies. These pointed at D:, the sessions then
% moved to E:, and E: is removable -- so the lookup silently failed and every
% group fell back to "1 cell per mask", overcounting Sert 402 instead of 296.
% The failure is a warning, not an error, which is exactly why it has to be a
% path that cannot go missing. Copied in beside each session 2026-07-31.
cell_link_sources = { ...
    'Sert',   '0721', fullfile(sumRoot, 'Sert',   '0721', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0728', fullfile(sumRoot, 'Vglut2', '0728', 'cell_pooled', 'cell_link.mat')
    'Vgat',   '0730', fullfile(sumRoot, 'Vgat',   '0730', 'cell_pooled', 'cell_link.mat')
    'Sst',    '0806', fullfile(sumRoot, 'Sst',    '0806', 'cell_pooled', 'cell_link.mat')
    'Sst',    '0807', fullfile(sumRoot, 'Sst',    '0807', 'cell_pooled', 'cell_link.mat')
    'Vglut2', '0810', fullfile(sumRoot, 'Vglut2', '0810', 'cell_pooled', 'cell_link.mat')
    };

scan_dirs   = {'ChAT','Vglut2','Vgat','Sst','Sert'};   % + Sert (260721 baseline)
groups      = {'IO','ChAT','Vglut2','Vgat','Sst','Sert'};

PixelSizeBase = 1.7778;

% ---- stage-axis -> anatomy map, PER SESSION --------------------------------
% The map's whole geometry assumes the animal was mounted the same way every
% session, so that stage x is the medio-lateral axis and stage y the
% rostro-caudal one. Two sessions were mounted ROTATED 90 deg, so their stage
% axes carry the opposite anatomical meaning:
%
%     -stage x = anterior (rostral)  ;  stage y = lateral (+y = left)
%
% Encoded as a 2x2 taking stage coordinates to anatomical ones:
%
%     [lateral ; rostral] = axmap * [stage x ; stage y]
%
% The lateral SIGN is not consequential -- everything is folded into the left
% hemisphere below with -abs() -- so the real content of the rotated row is
% "rostral comes from -stage x", and the FOV footprint transposes with it.
%
% Keyed on SESSION (scan dir + date), not genotype: how the animal was mounted
% is a property of that day's prep. Vglut2 0224 and 1124 are NOT rotated, only
% 0728 is -- so a genotype-level flag would have moved the older pFN data too.
%
% Independent corroboration for both, from the experimenter's own filenames:
% the rotated sessions write the coordinate pair y-first (`y1900x1320`), and
% their motorPosition satisfies motorX = -nameY, motorY = nameX.
AXMAP_DEFAULT = [ 1  0 ;  0  1 ];
AXMAP_ROT90   = [ 0  1 ; -1  0 ];
% Sessions acquired on the NEW mount / coordinate axis. This is a fact about the
% rig on the day, NOT something to infer from where the stage coordinates land --
% Vgat/0730's footprint happens to overlap the old-mount cloud, so guessing from
% coordinates got it wrong. Everything from 260721 onward is on the new axis.
rotated_sessions = { 'Sert','0721' ; 'Vglut2','0728' ; 'Vgat','0730' ; ...
                     'Sst','0806'  ; 'Sst','0807'   ; 'Vglut2','0810' };
% Vglut2/0810 added 2026-08-11 on the same rule: everything from 260721 onward is
% on the new mount. NOTE it is the same genotype as Vglut2 0224/1124, which are NOT
% rotated -- so this had to be keyed per SESSION, never per genotype.
% Sst 0806/0807 added 2026-08-08: acquired on the NEW mount, same as everything
% from 260721 onward, per the experimenter. NOTE these two sessions have no
% cross-FOV matcher run, so every mask counts as its own cell -- the Sst cell
% count is therefore an OVERCOUNT (the same neuron re-imaged at another z or
% zoom appears twice). Accepted deliberately: this map is for orientation.

% ---- stitch corrections for the ANATOMY (left) side -------------------------
% 2026-08-01. Overlapping FOVs were manually stitched with roi_stitch_260801, so
% their positions are now known far better than the raw stage reports them
% (measured residual ~2 um, versus stage corrections of up to ~26 um).
%
% SCOPE: these corrections are applied to the FOV CENTRE ONLY, which moves the
% left-hemisphere anatomy patches. The per-ROI positions in R -- and therefore
% the right-hemisphere phase dots -- are deliberately left untouched, so the
% right side of this figure is bit-identical to the 260729 version.
%
% The offsets live in the stitch pipeline's frame, so they are rotated back into
% stage coordinates before being added to motorPosition (see stitch_offset_um).
USE_STITCH_OFFSETS = true;
STITCH_DIR = 'roi_stitch_out_260801';

stitchOff = containers.Map('KeyType','char','ValueType','any');
if USE_STITCH_OFFSETS
    ot = dir(fullfile(sumRoot,'**',STITCH_DIR,'fov_offsets.mat'));
    for k = 1:numel(ot)
        Ok = load(fullfile(ot(k).folder, ot(k).name));
        Tk = load(fullfile(ot(k).folder, 'fov_table.mat'));   % folder path per FOV
        if ~isfield(Ok,'fovTable') || ~isfield(Tk,'F'), continue; end
        for i = 1:numel(Tk.F)
            j = find(Ok.fovTable.fov == string(Tk.F(i).name), 1);
            if isempty(j) || ~Ok.fovTable.placed(j), continue; end
            if abs(Ok.fovTable.tx_um(j)) < 1e-9 && abs(Ok.fovTable.ty_um(j)) < 1e-9, continue; end
            stitchOff(lower(Tk.F(i).folder)) = [Ok.fovTable.tx_um(j), Ok.fovTable.ty_um(j)];
        end
    end
    fprintf('stitch offsets loaded for %d FOV(s) from %d session(s)\n', ...
        stitchOff.Count, numel(ot));
end
nStitched = 0; stitchMag = [];

% ---- one patch per imaging SITE --------------------------------------------
% Several sites were re-imaged at a second Z plane, or (ro7/roi7) recorded twice
% under near-identical names. Each such recording used to stamp its own tinted
% patch at essentially the same place, which is what made a site look drawn
% "repetitively" -- 0.85-alpha patches piling up read as one dark smear.
% So the PATCH is drawn once per site, from the recording with the most ROIs.
% Every ROI of every recording is still kept in R: two Z planes are two
% different sets of cells and only the rendering was duplicated, not the data.
mergeSiteUm = 100;      % same site if FOV centres are closer than this (um)

% per-group tints (FOV patch colors)
% Vgat's tint is taken from the shared display table's "Vgat I" entry rather than
% the 260808 literal, so the Vgat anatomy on the left is the same blue as the
% Vgat I dots on the right and the Vgat I trace in the population overlays.
% The other genotypes keep their 260808 tints: each of them maps to more than one
% display group (Vglut2 splits pre-I / I, Vgat also has the tonic group), so
% there is no single group colour to inherit -- only Vgat I was asked for.
DTMP = group_display_colors_260816();
vgatI_col = DTMP{strcmp(DTMP(:,1),'Vgat I'), 2};
patch_tints = struct( ...
    'IO',     [0.5    0.5    0.5   ], ...   % black
    'ChAT',   [0.85 0.10 0.10], ...   % red
    'Vglut2', [0.10 0.65 0.20], ...   % green
    'Vgat',   vgatI_col,        ...   % light blue, = "Vgat I" everywhere else
    'Sst',    [0.55 0.20 0.75], ...   % purple
    'Sert',   [0.90 0.45 0.10]);      % orange, matching cfg.genotype_color
patch_alpha    = 0.85;
patch_clip_pct = [1 99.5];
patch_gamma    = 0.6;
% =========================================================

%% ---- load cartoon coordframe ----
C = load(cartoonMat);
assert(isfile(C.imgPath), 'Cartoon image missing: %s', C.imgPath);
Icart = im2double(imread(C.imgPath));
if size(Icart,3) == 1, Icart = repmat(Icart,1,1,3); end
% Light gaussian blur to soften pixelation (the cartoon was a raster)
bg_blur_sigma = 1.2;            % px ; bump to 2-3 if still pixelated
Icart = imgaussfilt(Icart, bg_blur_sigma);
[Hc, Wc, ~] = size(Icart);
% cartoon data extent (μm). pix_per_y carries sign of user's click.
xd = ([1 Wc] - C.origin_pix(1)) / C.pix_per_x;
yd = ([1 Hc] - C.origin_pix(2)) / C.pix_per_y;

%% ---- load vessel-derived calibration ----
V = load(vesselMat, 'bx','by','theta_rad','ax','ay','um_per_px_stitch');
bx = V.bx; by = V.by; theta_rad = V.theta_rad;
um_per_px_stitch = V.um_per_px_stitch;
ax = V.ax;  ay = V.ay;
stage_to_canvas = @(mx, my) deal(ax(1)*mx + ax(2), ay(1)*my + ay(2));
Rmat = [cos(theta_rad)  sin(theta_rad);
       -sin(theta_rad)  cos(theta_rad)];

%% ---- significance + coherence lookup ----
sig_map = containers.Map('KeyType','char','ValueType','logical');
r_map   = containers.Map('KeyType','char','ValueType','double');
th_map  = containers.Map('KeyType','char','ValueType','double');
crit_val = 0;  critName = '';
% Significance level for THIS FIGURE, set here rather than inherited from the
% polar run. The polar draws two circles (alpha 0.05 inner, 0.001 outer) and
% saves the 0.05 one as logZ_crit; the map wants the stricter threshold, so it
% derives its own. Z = n*Rbar^2 with p ~ exp(-Z), hence logZ_crit = log(-log(alpha))
% -- n-independent, so one cut is valid for every cell regardless of event count.
sig_alpha = 0.001;      % 0.05 -> logZ 1.10 | 0.001 -> logZ 1.93
if isfile(dataMat)
    Sray = load(dataMat,'RY','labels','obsOfCell','logZ_crit','alpha_ray');
    critZ   = log(-log(sig_alpha));
    sigCell = Sray.RY.logZ >= critZ;                    % per CELL
    for i = 1:numel(Sray.labels)
        c = Sray.obsOfCell(i);                          % this observation's cell
        lab = Sray.labels{i};
        sig_map(lab) = sigCell(c);
        r_map(lab)   = Sray.RY.logZ(c);
        th_map(lab)  = mod(Sray.RY.th(c), 2*pi);        % MEAN direction, 0..2pi (was thMed, 2026-08-11)
    end
    crit_val = critZ;
    critName = sprintf('logZ\\geq%.2f, \\alpha=%.3g', critZ, sig_alpha);
    fprintf('  phase/significance from %s\n', dataMat);
    fprintf('    %d observations -> %d cells, %d significant at alpha=%.3g (logZ>=%.2f)\n', ...
            numel(Sray.labels), numel(sigCell), nnz(sigCell), sig_alpha, critZ);
    fprintf('    (polar run itself used alpha=%.2g, logZ>=%.2f -> %d cells)\n', ...
            Sray.alpha_ray, Sray.logZ_crit, nnz(Sray.RY.logZ >= Sray.logZ_crit));
else
    warning('No polar_coh_vs_rayleigh_data.mat -- everyone non-sig.');
end

%% ---- THE SELECTED SET: cells, their phase, their verdict -------------------
% Source of truth for the right side. polar_selected_260816.mat carries, per
% CELL: the occupancy-weighted circular MEDIAN phase (thMed) and logZ. Median,
% not mean, because that is what the polar figure plots -- the two must agree.
%
% The registry (event_latency_data.mat) is what turns a cell number back into
% the (group/date/recName/rid) observations this map can place, and it is read
% through pooled_obs_260814 so the merged cells resolve to ALL their masks
% rather than the first one only.
addpath(fileparts(mfilename('fullpath')));
selMat = fullfile(sumRoot,'polar_selected_260816','polar_selected_260816.mat');
assert(isfile(selMat), ['Run polar_selected_260816.m first -- ' selMat ' is missing.']);
SEL = load(selMat,'R');
regD = load(fullfile(sumRoot,'event_latency_260811','event_latency_data.mat'), ...
            'CELL','OBS','REC');
obsOfSel = pooled_obs_260814(regD.CELL, regD.OBS);

% Significance for the DOTS. Same rule and same alpha as the left-side outlines
% (logZ >= log(-log(alpha)), n-independent), so "significant" means the same
% thing on both halves of the figure.
sel_alpha  = sig_alpha;                 % 0.001 -> logZ 1.93 | 0.05 -> logZ 1.10
sel_critZ  = log(-log(sel_alpha));

selCell  = containers.Map('KeyType','char','ValueType','double');  % label -> cell
selPhase = containers.Map('KeyType','double','ValueType','double');
selLogZ  = containers.Map('KeyType','double','ValueType','double');
selGroup = containers.Map('KeyType','double','ValueType','char');
selClass = containers.Map('KeyType','double','ValueType','char');
for i = 1:numel(SEL.R)
    c = SEL.R(i).cell;
    selPhase(c) = mod(SEL.R(i).thMed, 2*pi);
    selLogZ(c)  = SEL.R(i).logZ;
    selGroup(c) = char(SEL.R(i).group);
    selClass(c) = char(SEL.R(i).class);
    for o = obsOfSel{c}(:)'
        selCell(regD.OBS(o).label) = c;      % 'group/date/recName/rid'
    end
end
nSelSig = nnz(cellfun(@(c) selLogZ(c) >= sel_critZ, num2cell(cell2mat(selLogZ.keys))));
fprintf(['  selected set: %d cells, %d observations, %d significant ' ...
         '(logZ>=%.2f, alpha=%.3g)\n'], ...
        numel(SEL.R), selCell.Count, nSelSig, sel_critZ, sel_alpha);

% ---- WHICH selected cells actually get drawn on the right ------------------
% Two restrictions, both asked for, both applied here so that the dots and the
% white outlines cannot disagree:
%   IO      dropped entirely. 107 of the 147 are IO and they land in one dense
%           clump, so they dominate a figure that is about the other groups.
%   Sert    significant only.
% Everything else keeps its full selected set, significant or not.
DROP_GROUPS = "IO";
% NO extra significance gate (2026-08-17). Sert used to need Rayleigh logZ >=
% 1.93 on TOP of being in the class list, which silently dropped 6 of its 10
% cells while every other genotype drew its full selected set -- so Sert
% looked sparse for a reason that had nothing to do with anatomy. The class
% list is now the only selection: whatever cell_classes_260816 assigns, the
% map draws.
SIG_ONLY_GROUPS = strings(0,1);
selDraw = containers.Map('KeyType','double','ValueType','logical');
for i = 1:numel(SEL.R)
    c = SEL.R(i).cell;  g = string(SEL.R(i).group);
    keepC = ~any(g == DROP_GROUPS);
    if any(g == SIG_ONLY_GROUPS), keepC = keepC && (selLogZ(c) >= sel_critZ); end
    selDraw(c) = keepC;
end
nDraw = nnz(cell2mat(selDraw.values));
fprintf('  drawn on the right: %d of %d cells (IO dropped)\n', ...
        nDraw, numel(SEL.R));

%% ---- cell-identity lookup ----
% keys are '<scanDir>/<date>/<recName>/<maskL label>'; values are a cell key that
% is namespaced by session, since cell ids restart at 1 in every matcher run.
cell_map   = containers.Map('KeyType','char','ValueType','char');
tossed_set = containers.Map('KeyType','char','ValueType','logical');
for s = 1:size(cell_link_sources,1)
    pre = sprintf('%s/%s', cell_link_sources{s,1}, cell_link_sources{s,2});
    f   = cell_link_sources{s,3};
    if ~isfile(f)
        warning('cell_link missing for %s (%s) -- that session will count masks as cells.', pre, f);
        continue;
    end
    Lk = load(f,'link');  T = Lk.link.obsT;
    for i = 1:height(T)
        kk = sprintf('%s/%s/%d', pre, T.rec_name(i), T.maskL_label(i));
        if isnan(T.cell_id(i)), tossed_set(kk) = true;
        else,                   cell_map(kk) = sprintf('%s#c%d', pre, T.cell_id(i));
        end
    end
    fprintf('  cell identity %-14s %4d masks -> %4d cells (%d tossed in curation)\n', ...
            pre, height(T), numel(unique(T.cell_id(~isnan(T.cell_id)))), nnz(isnan(T.cell_id)));
end

%% ---- walk every FOV ----
R = struct('group',{},'date',{},'fov',{},'rid',{},'is_sig',{}, ...
           'x_um',{},'y_um',{},'r_coh',{},'th_coh',{},'cell_key',{},'sel_cell',{});
F = struct('fov',{},'gname',{},'avgPath',{},'maskPath',{}, ...
           'x_um',{},'y_um',{}, ...
           'H_fov',{},'W_fov',{},'px_um_fov',{},'sig_rids',{}, ...
           'sess',{},'site',{},'nroi',{},'rot90',{},'drawPatch',{},'cfgEx',{}, ...
           'depth_um',{});

nCfgSkip = 0;
for sg = 1:numel(scan_dirs)
    sname = scan_dirs{sg};
    gdir  = fullfile(sumRoot, sname);
    if ~isfolder(gdir), continue; end
    sam_hits = dir(fullfile(gdir,'**','*cpSAM_output.mat'));
    for hh = 1:numel(sam_hits)
        fp  = sam_hits(hh).folder;
        rel = strrep(fp, gdir, ''); rel = regexprep(rel,'^[\\/]+','');
        parts = regexp(rel,'[\\/]','split');
        if numel(parts) < 3, continue; end
        dateStr = parts{1}; subcell = parts{2};
        [~,n,e] = fileparts(fp); fov = [n e];

        % 120 um ventral-surface depth bar + the stage-zero exclusion, both from
        % coh_cfg_260727.  This script does not read that list on its own, so it
        % is applied here to keep the map, the polar plot and the per-cell
        % summary describing the SAME set of recordings.
        % CHANGED from 260808. coh_cfg's exclusion list mixes two unrelated
        % reasons, and only ONE of them is about this figure:
        %
        %   DEPTH BAR (6 recordings, all Vgat/0730 deeper than 120 um) -- a
        %   scope choice for a ventral-SURFACE analysis. Their coordinates are
        %   fine. Now INCLUDED: patch and dots both, on request.
        %
        %   STAGE ZERO (ro1_2x_x1018y1130_18lp_3000f_00001) -- that field did
        %   not zero its stage on the same landmark as every other session, so
        %   its coordinates are simply WRONG. stage_to_canvas is one affine
        %   applied to every animal, so drawing it would put a patch in a place
        %   the tissue never was. Stays excluded, and must.
        STAGE_ZERO_EX = "ro1_2x_x1018y1130_18lp_3000f_00001";
        cfgEx = any(strcmp(fov, STAGE_ZERO_EX));
        if cfgEx, nCfgSkip = nCfgSkip + 1; end

        gname = sname;
        if ~isempty(regexpi(subcell,'IO','once')), gname = 'IO'; end
        gi = find(strcmp(groups, gname), 1);
        if isempty(gi), continue; end

        mh = dir(fullfile(fp,'*_ch1_meta.mat'));
        if isempty(mh), mh = dir(fullfile(fp,'*_meta.mat')); end
        if isempty(mh), continue; end
        M = load(fullfile(mh(1).folder, mh(1).name));
        if ~isfield(M,'motorPosition') || numel(M.motorPosition) < 3, continue; end
        motor = M.motorPosition(:).';
        if isfield(M,'pixelSize_um') && M.pixelSize_um>0
            px_um_fov = M.pixelSize_um;
        elseif isfield(M,'zoomFactor') && M.zoomFactor>0
            px_um_fov = PixelSizeBase / M.zoomFactor;
        else, continue;
        end

        maskPath = fullfile(sam_hits(hh).folder, sam_hits(hh).name);
        S2 = load(maskPath, 'maskL');
        if ~isfield(S2,'maskL') || isempty(S2.maskL), continue; end
        maskL = S2.maskL; [Hf, Wf] = size(maskL);

        % which way was this animal mounted?
        isRot = any(strcmp(rotated_sessions(:,1), sname) & ...
                    strcmp(rotated_sessions(:,2), dateStr));
        if isRot, axmap = AXMAP_ROT90; else, axmap = AXMAP_DEFAULT; end

        % Manual stitch correction, ANATOMY (left) side only. The offset is stored
        % in the stitch pipeline's frame, so it is rotated back to stage
        % coordinates with axmap' (the map is orthogonal) before being added.
        motorA = motor;
        if isKey(stitchOff, lower(fp))
            d = stitchOff(lower(fp));
            ds = axmap.' * d(:);
            motorA(1:2) = motor(1:2) + ds(:).';
            nStitched = nStitched + 1;
            stitchMag(end+1) = hypot(d(1), d(2)); %#ok<SAGROW>
        end

        % FOV centre: stage->anatomical axes->canvas->rotated->μm; then y-up
        [cx_F, cy_F] = stage_to_canvas_ax(stage_to_canvas, axmap, motorA(1), motorA(2));
        dxy = Rmat * [cx_F - bx; cy_F - by];
        x_um = dxy(1) * um_per_px_stitch;
        y_um = -dxy(2) * um_per_px_stitch;            % negate so +y = rostral

        ah = dir(fullfile(fp,'*_AVG_for_CP.tif'));
        if ~isempty(ah), avgPath = fullfile(ah(1).folder, ah(1).name); else, avgPath = ''; end

        lbls = setdiff(unique(maskL(:)), 0);
        props = regionprops(maskL, 'Centroid'); props = props(lbls);
        sig_rids_here = [];
        nTossHere = 0;
        for j = 1:numel(lbls)
            rid = lbls(j);
            lab = sprintf('%s/%s/%s/%d', gname, dateStr, fov, rid);
            % cell-identity key uses the SCAN DIR (sname), not gname: the matcher
            % ran per session, before the IO folder override was applied.
            ckin = sprintf('%s/%s/%s/%d', sname, dateStr, fov, rid);
            if isKey(tossed_set, ckin)
                nTossHere = nTossHere + 1;  continue;   % rejected in curation
            end
            if isKey(cell_map, ckin), cell_key = cell_map(ckin);
            else,                     cell_key = ['mask:' lab];   % no matcher run
            end
            % CHANGED from 260808: the white outline now marks EXACTLY the cells
            % that get a phase-coloured dot on the right -- selected AND
            % significant -- instead of the old map's own sig_map. The two halves
            % of the figure therefore point at the same neurons, and an outline
            % on the left can be paired with a coloured dot on the right.
            is_sig = isKey(selCell, lab) && ...
                     selDraw(selCell(lab)) && ...
                     selLogZ(selCell(lab)) >= sel_critZ;
            if is_sig, sig_rids_here(end+1) = rid; end %#ok<AGROW>
            if isKey(r_map,  lab), r_coh  = r_map(lab);  else, r_coh  = NaN; end
            if isKey(th_map, lab), th_coh = th_map(lab); else, th_coh = NaN; end

            % per-ROI position: centroid -> stage um -> anatomical axes -> canvas
            % -> rotated -> y-up um.  The scan mirrors are fixed to the RIG, so a
            % FOV column is always stage x and a row always stage y no matter how
            % the animal was mounted; the mount only enters at the stage->anatomy
            % step, which is where axmap is applied.
            c = props(j).Centroid;                       % [column row]
            mx_r = motor(1) + (c(1) - Wf/2) * px_um_fov;
            my_r = motor(2) + (c(2) - Hf/2) * px_um_fov;
            [cxR, cyR] = stage_to_canvas_ax(stage_to_canvas, axmap, mx_r, my_r);
            dxyR = Rmat * [cxR - bx; cyR - by];
            x_umR =  dxyR(1) * um_per_px_stitch;
            y_umR = -dxyR(2) * um_per_px_stitch;

            % Which SELECTED cell is this mask, if any. NaN = not in the 147.
            if isKey(selCell, lab), selc = selCell(lab); else, selc = NaN; end

            R(end+1) = struct('group',gname,'date',dateStr,'fov',fov, ...
                              'rid',rid,'is_sig',is_sig, ...
                              'x_um',x_umR,'y_um',y_umR, ...
                              'r_coh',r_coh,'th_coh',th_coh, ...
                              'cell_key',cell_key,'sel_cell',selc); %#ok<SAGROW>
        end
        if nTossHere > 0
            fprintf('  %-42s dropped %d mask(s) tossed in curation\n', fov(1:min(42,end)), nTossHere);
        end

        F(end+1) = struct('fov',fov,'gname',gname,'avgPath',avgPath, ...
                          'maskPath',maskPath, ...
                          'x_um',x_um,'y_um',y_um, ...
                          'H_fov',Hf,'W_fov',Wf,'px_um_fov',px_um_fov, ...
                          'sig_rids',sig_rids_here, ...
                          'sess',[sname '/' dateStr],'site',subcell, ...
                          'nroi',numel(lbls),'rot90',isRot, ...
                          'drawPatch',~cfgEx,'cfgEx',cfgEx, ...
                          'depth_um',abs(motor(3))); %#ok<SAGROW>
    end
end
fprintf('excluded by coh_cfg (stage-zero + 120 um depth bar): %d FOV(s)\n', nCfgSkip);
assert(~isempty(F), 'No FOVs.');

%% ---- optional: tie-point registered coordinates -----------------------------
% Ventral_surface_tiepoint_register_260731 solves one translation per FOV from
% cells the matcher saw in more than one FOV. Applied HERE, before the site
% merge, deliberately: the merge groups FOVs by centre distance, so correcting
% centres afterwards regroups nothing and can pull a group apart (measured: one
% 8-member group went 80 -> 174 um when the correction was applied after).
% FOVs with no tie-point get zero shift and are left on stage coordinates, so a
% map made this way MIXES corrected and uncorrected placements -- F.reg_linked
% records which is which rather than hiding it.
USE_REGISTERED = false;      % <-- set true to draw from corrected coords
% Which offset file. Both carry the same fovTable schema:
%   ROI_on_cartoon_registered.mat        automatic, tie-points from matched cells
%   ROI_on_cartoon_manual_offsets.mat    hand-nudged in Ventral_surface_manual_align
REG_SOURCE = 'manual';       % 'auto' | 'manual'
if strcmpi(REG_SOURCE,'manual')
    regMat = fullfile(fileparts(cartoonMat), 'ROI_on_cartoon_manual_offsets.mat');
else
    regMat = fullfile(fileparts(cartoonMat), 'ROI_on_cartoon_registered.mat');
end
if USE_REGISTERED
    assert(isfile(regMat), 'USE_REGISTERED but %s is missing.', regMat);
    R2 = load(regMat, 'fovTable');  T2 = R2.fovTable;
    kF = string({F.sess}) + "/" + string({F.fov});
    nMoved = 0;
    for i = 1:numel(F)
        j = find(T2.fov == kF(i), 1);
        F(i).reg_linked = false;
        if isempty(j), continue; end
        F(i).x_um = F(i).x_um + T2.tx_um(j);
        F(i).y_um = F(i).y_um + T2.ty_um(j);
        F(i).reg_linked = T2.linked(j);
        if T2.linked(j), nMoved = nMoved + 1; end
    end
    kR = string({R.group}) + "/" + string({R.date}) + "/" + string({R.fov});
    for i = 1:numel(R)
        j = find(T2.fov == kR(i), 1);
        if isempty(j), continue; end
        R(i).x_um = R(i).x_um + T2.tx_um(j);
        R(i).y_um = R(i).y_um + T2.ty_um(j);
    end
    fprintf('  REGISTERED coordinates applied: %d of %d FOVs moved (rest have no tie-point)\n', ...
        nMoved, numel(F));
else
    for i = 1:numel(F), F(i).reg_linked = false; end
end

%% ---- collapse re-imaged sites to ONE patch each ----
% MERGE_SITES = false -> no collapsing at all: every FOV draws its own patch and
% overlapping fields simply overlay. That shows the true coverage instead of one
% representative field per site. The reason merging existed is that 0.85-alpha
% patches stack into a dark smear, so alpha is scaled down here to compensate --
% with 9 fields on one site the accumulated opacity is what needs controlling,
% not the number of fields.
MERGE_SITES = true;
if ~MERGE_SITES
    nOv = zeros(numel(F),1);
    for u = 1:numel(F)
        for v = 1:numel(F)
            if u==v || ~strcmp(F(u).sess,F(v).sess), continue; end
            if hypot(F(u).x_um-F(v).x_um, F(u).y_um-F(v).y_um) <= mergeSiteUm
                nOv(u) = nOv(u) + 1;
            end
        end
    end
    patch_alpha = patch_alpha / max(1, sqrt(median(nOv(nOv>0)+1)));
    fprintf('  MERGE_SITES=false: all %d patches drawn, alpha scaled to %.2f\n', ...
        numel(F), patch_alpha);
end
if MERGE_SITES
% Same session + same group + FOV centres within mergeSiteUm  ->  one site. The
% survivor is the recording with the most ROIs, so the patch shown is the one
% carrying the most cells. R is untouched: this suppresses duplicate DRAWING
% only, it does not drop any cell.
key = strcat({F.sess}, '|', {F.gname});
for k = unique(key)
    idx = find(strcmp(key, k{1}));
    for a = 1:numel(idx)
        ia = idx(a);
        if ~F(ia).drawPatch, continue; end
        for b = a+1:numel(idx)
            ib = idx(b);
            if ~F(ib).drawPatch, continue; end
            if hypot(F(ia).x_um-F(ib).x_um, F(ia).y_um-F(ib).y_um) > mergeSiteUm, continue; end
            if F(ib).nroi > F(ia).nroi, drop = ia; keepi = ib; else, drop = ib; keepi = ia; end
            F(drop).drawPatch = false;
            fprintf('  merged patch: %-42s -> %-42s (%s, %.0f um apart)\n', ...
                    F(drop).fov(1:min(42,end)), F(keepi).fov(1:min(42,end)), ...
                    F(drop).sess, hypot(F(ia).x_um-F(ib).x_um, F(ia).y_um-F(ib).y_um));
            if drop == ia, break; end
        end
    end
end
end   % if MERGE_SITES
nRot = nnz([F.rot90]);
fprintf('FOVs: %d (patches drawn: %d, %d on rotated-mount sessions)   ROIs: %d\n', ...
        numel(F), nnz([F.drawPatch]), nRot, numel(R));
if USE_STITCH_OFFSETS
    if nStitched > 0
        fprintf(['anatomy (left) side: %d FOV centre(s) moved by the manual stitch ' ...
                 '(median %.1f um, max %.1f um)\n'], ...
                nStitched, median(stitchMag), max(stitchMag));
    else
        fprintf('anatomy (left) side: no stitch offsets matched any FOV\n');
    end
    fprintf('right side (phase dots) unchanged: per-ROI positions were not touched\n');
end

%% ---- figure: cartoon background, y-up ----
figOut = figure('Color','w','Name','ROIs on cartoon', ...
                'Units','normalized','Position',[0.04 0.04 0.9 0.9]);
ax_p = axes(figOut);
image(ax_p, xd, yd, Icart); hold(ax_p,'on');
set(ax_p,'YDir','normal');
axis(ax_p,'image'); box(ax_p,'on'); grid(ax_p,'on');
xlabel(ax_p,'X (\mum)');
ylabel(ax_p,'Y (\mum)');

% --- overlay each FOV patch as its own image() in data coords ---
% DRAW ORDER, changed from 260808. The old script drew in scan_dirs order --
% ChAT, Vglut2, Vgat, Sst, Sert -- so Sert painted last and every one of the 19
% Sst patches ended up underneath one, which is why Sst looked absent from the
% map even though it was always being drawn. Sorting by footprint AREA and
% drawing the largest first means a small field can no longer be buried under a
% large one, and it is genotype-neutral: nothing gets priority for being a
% particular label, only for being small enough to be hidden otherwise.
areaF = arrayfun(@(f) (f.W_fov*f.px_um_fov) * (f.H_fov*f.px_um_fov), F);
[~, drawOrder] = sort(areaF, 'descend');
for u = drawOrder(:).'
    if isempty(F(u).avgPath) || ~isfile(F(u).avgPath), continue; end
    % The stage-zero field draws NOTHING -- its coordinates are wrong, so a
    % patch would land where the tissue never was.
    if F(u).cfgEx, continue; end
    gname_u = F(u).gname;
    if isfield(patch_tints, gname_u), tintU = patch_tints.(gname_u);
    else, tintU = [0.10 0.80 0.30]; end

    A = double(imread(F(u).avgPath));
    a_lo = prctile(A(:), patch_clip_pct(1));
    a_hi = prctile(A(:), patch_clip_pct(2));
    aN   = max(0, min(1, (A - a_lo)/max(a_hi-a_lo,eps))).^patch_gamma;

    % Rotated mount -> the FOV's own axes are 90 deg off the anatomy, so the
    % footprint transposes: a 320x512 px FOV that spanned 569 um of stage x now
    % spans 569 um ROSTRALLY and only 356 um laterally. rot90 carries both the
    % texture and the extent across. (The residual handedness of the texture is
    % not resolved -- these patches are decorative soft blobs, and the existing
    % non-rotated ones already draw row 1 at high y. Position and footprint are
    % what carry information here, and both are correct.)
    if F(u).rot90, aN = rot90(aN, -1); end

    % patch tinted RGB in source pixel size
    [hp, wp] = size(aN);
    [Yp, Xp]  = ndgrid(linspace(-1,1,hp), linspace(-1,1,wp));
    rho       = sqrt(Xp.^2 + Yp.^2);
    % Strong low-frequency perturbation -> clearly non-circular blob
    rho_noise = imgaussfilt(randn(hp, wp), max(hp,wp)/8) * 0.7;
    circ_mask = max(0, min(1, (1 - (rho + rho_noise))*3));
    alphaPatch = patch_alpha .* circ_mask .* aN;

    rgbPatch = zeros(hp, wp, 3);
    for ch = 1:3, rgbPatch(:,:,ch) = tintU(ch); end

    % world extent in μm  (centre at FOV's x_um, y_um; y-up so flip vertically)
    halfW_um = (wp/2) * F(u).px_um_fov;
    halfH_um = (hp/2) * F(u).px_um_fov;
    xext = [F(u).x_um - halfW_um, F(u).x_um + halfW_um];
    yext = [F(u).y_um + halfH_um, F(u).y_um - halfH_um];  % image row 1 at top (high y)

    % Only the TINTED BLOB is suppressed at a re-imaged site -- that is what piled
    % up and smeared. The sig-ROI outlines below still draw for every recording,
    % each in its OWN extent (they carry which cells were significant, and a thin
    % white outline does not smear). Masks from different recordings cannot be
    % overlaid on one patch anyway: different zoom, pixel size and Z.
    if F(u).drawPatch
        image(ax_p, xext, yext, rgbPatch, 'AlphaData', alphaPatch);
    end

    % NO ROI outlines on the anatomy. The significant cells are already shown as
    % circles on the right, and marking them a second time on the left added
    % nothing but clutter. F(u).sig_rids is still recorded in the saved .mat, so
    % which masks were significant per FOV is not lost -- only the drawing went.
end

% --- right hemisphere: ONE DOT PER SELECTED CELL ------------------------------
% Black for every cell in the selected set; the significant ones filled with
% their tuned phase. Phase convention as before: 0 = insp onset, pi = insp peak,
% HSV because it wraps.
%
% A cell can appear in several FOVs, so its observations are averaged to ONE
% position -- drawing a dot per mask would draw the same neuron twice. The
% observations of a genuine cell sit within a few um of each other; the spread is
% printed below so a bad match would show up rather than average away quietly.
phase_rot = 0;        % rotate HSV so phase 0 lands at desired hue (radians)
selIds = unique([R(~isnan([R.sel_cell])).sel_cell]);
selIds = selIds(arrayfun(@(c) selDraw(c), selIds));    % IO out, Sert sig-only
% VENTRAL SURFACE ONLY (2026-08-17). A cell is dropped when EVERY one of its
% observations sits in a depth-barred FOV. Cells with at least one surface
% observation are kept and keep all their observations in the position average --
% they are surface cells that were also imaged deep, not deep cells.
% Depth is |motorPosition(3)|, NEVER the z in the filename -- the two disagree
% (roi1_1.5x_..._z100 reads 201.9 um on the motor) and the motor is the authority.
% coh_cfg's own depth bar is NOT used: it flagged 0 FOVs in this run, so relying
% on it would have silently kept every deep cell on a ventral-surface figure.
DEEP_UM = 120;
if isfield(F,'depth_um')
    deepFovs = {F([F.depth_um] > DEEP_UM).fov};
    fprintf('depth bar %g um: %d of %d FOV(s) are deep\n', ...
            DEEP_UM, numel(deepFovs), numel(F));
    isDeep = false(size(selIds));
    for i = 1:numel(selIds)
        m = find([R.sel_cell] == selIds(i));
        isDeep(i) = ~isempty(m) && all(ismember({R(m).fov}, deepFovs));
    end
    if any(isDeep)
        fprintf('ventral-surface only: dropped %d deep cell(s): %s\n', ...
                nnz(isDeep), mat2str(selIds(isDeep)));
        selIds = selIds(~isDeep);
    end
end
nSel   = numel(selIds);
sx = nan(nSel,1); sy = nan(nSel,1); sth = nan(nSel,1);
ssig = false(nSel,1); sspread = nan(nSel,1); snobs = zeros(nSel,1);
for i = 1:nSel
    m = find([R.sel_cell] == selIds(i));
    xs = [R(m).x_um];  ys = [R(m).y_um];
    sx(i) = mean(xs);  sy(i) = mean(ys);  snobs(i) = numel(m);
    if numel(m) > 1
        sspread(i) = max(hypot(xs - sx(i), ys - sy(i)));
    end
    sth(i)  = selPhase(selIds(i));
    ssig(i) = selLogZ(selIds(i)) >= sel_critZ;
end
rx = -sx;                            % mirror LEFT->RIGHT
ry =  sy;
th = mod(sth - phase_rot, 2*pi);     % kept for the CSV only, no longer drawn

% COLOUR = GENOTYPE, MARKER = CLASS (2026-08-17).
% Two independent codes rather than one colour per (genotype, class) pair: the
% reader can then answer "which genotype" and "which class" separately, and a
% genotype keeps ONE identity across all its classes.
% COLOURS COME FROM patch_tints -- the SAME table the left-hemisphere anatomy
% patches are tinted with. Derived, not copied: a second hardcoded list is how
% the dot for a genotype ends up a different green from its own anatomy, and the
% reader then reads them as two populations. Change a colour in patch_tints and
% both sides move together.
gn_ord = {'Sst','Sert','ChAT','Vgat','Vglut2'};
GCOL = cell(numel(gn_ord),2);
for gg = 1:numel(gn_ord)
    GCOL{gg,1} = gn_ord{gg};
    assert(isfield(patch_tints, gn_ord{gg}), ...
           'patch_tints has no colour for %s', gn_ord{gg});
    GCOL{gg,2} = patch_tints.(gn_ord{gg});
end
% third column: true = filled, false = OPEN (edge only).
% NOTE 'null' is not a class the brief named a marker for. It is drawn as an open
% circle, the same as tonic, on the reading that both mean "not respiratory
% modulated". That covers 8 Sst and 3 Vgat cells -- change this one line if they
% should be a separate shape or left off entirely.
% {class, kind, glyph, filled}. kind 'txt' draws the LETTER itself at the cell's
% position instead of a marker: I = inspiratory, E = expiratory. A letter states
% the tuning directly, where a triangle-vs-square only does so via the legend.
% Roman (upright serif) so I and E read as letters, not as italic variables.
% pre-I keeps the star, at a LARGER size than the other glyphs -- a star's ink
% sits inside its bounding box, so at equal 'size' it reads smaller than a disc,
% and its dark edge eats more of what is left.
MRK  = { 'pre-I',                'mkr', 'p', true       % star, enlarged below
         'I',                    'txt', 'I', true       % letter I
         'post-I',               'txt', 'E', true       % letter E
         'tonic/rhythmic null',  'mkr', 'o', false };   % open circle
% NOTE scatter's size argument is an AREA, so a linear "x times bigger" is a
% factor of x^2 here. Written out rather than as a bare number so the intent
% survives the next adjustment.
SZ_MKR  = 52  * (2/3)^2;   % open circle: 2/3 the RADIUS  -> 23.1
SZ_STAR = 95  * 1.3^2;     % pre-I star : 1.3x the RADIUS -> 160.6
FS_TXT  = 12  * 1.2;       % I / E letters: 1.2x, and FontSize IS linear -> 14.4
gName = strings(nSel,1); cName = strings(nSel,1);
for i = 1:nSel
    gName(i) = string(selGroup(selIds(i)));
    cName(i) = string(selClass(selIds(i)));
end
% 'tonic pre-I suppressed' and 'null' are ONE display class. Both mean the cell
% is not phase-locked to the breath and both already draw the same open circle,
% so keeping the names apart split the legend without splitting the symbol.
% NOTE this is a DISPLAY merge only -- cell_classes_260816 still records which
% is which, and the CSV written below keeps the original class string.
cName(ismember(cName, ["tonic pre-I suppressed","null"])) = "tonic/rhythmic null";
combos = unique(gName + "|" + cName);
hG = gobjects(numel(combos),1); lG = strings(numel(combos),1);
for k = 1:numel(combos)
    pr = split(combos(k), "|");  g = pr(1);  c = pr(2);
    m  = find(gName == g & cName == c);
    ci = find(strcmpi(GCOL(:,1), g), 1);
    mi = find(strcmpi(MRK(:,1),  c), 1);
    % An unmapped genotype or class is REPORTED, never silently skipped -- a cell
    % missing from the map because no one gave it a colour looks identical to a
    % cell that was never selected.
    if isempty(ci) || isempty(mi)
        fprintf(2, 'NO STYLE for "%s %s" -- %d cell(s) not drawn\n', g, c, numel(m));
        continue;
    end
    col = GCOL{ci,2};  kind = MRK{mi,2};  gl = MRK{mi,3};  filled = MRK{mi,4};
    if strcmp(kind,'txt')
        % The letter IS the marker, centred both ways so it sits ON the cell's
        % coordinate rather than beside it.
        text(ax_p, rx(m), ry(m), repmat({gl}, numel(m), 1), ...
             'Color', col, 'FontName','Times New Roman', 'FontSize', FS_TXT, ...
             'FontWeight','bold', 'FontAngle','normal', ...
             'HorizontalAlignment','center', 'VerticalAlignment','middle');
        % A legend cannot use a text glyph as its icon, so the key is a filled
        % SQUARE in the genotype colour -- a colour swatch only -- and the letter
        % is named in the label. Square is free: no class draws one any more.
        hG(k) = scatter(ax_p, NaN, NaN, SZ_MKR, col, 's', 'filled', ...
                'MarkerEdgeColor','k', 'LineWidth',0.4);
    elseif filled
        sz = SZ_MKR;  if strcmp(gl,'p'), sz = SZ_STAR; end
        scatter(ax_p, rx(m), ry(m), sz, col, gl, 'filled', ...
                'MarkerEdgeColor','k', 'LineWidth',0.4);
        hG(k) = scatter(ax_p, NaN, NaN, sz, col, gl, 'filled', ...
                'MarkerEdgeColor','k', 'LineWidth',0.4);
    else
        scatter(ax_p, rx(m), ry(m), SZ_MKR, gl, ...
                'MarkerEdgeColor',col, 'LineWidth',1.3);
        hG(k) = scatter(ax_p, NaN, NaN, SZ_MKR, gl, ...
                'MarkerEdgeColor',col, 'LineWidth',1.3);
    end
    if strcmp(kind,'txt')
        lG(k) = sprintf('%s  %s = %s (n=%d)', g, gl, c, numel(m));
    else
        lG(k) = sprintf('%s %s (n=%d)', g, c, numel(m));
    end
end

% Placement check: a cell whose observations disagree by more than a soma is a
% cross-FOV match to be suspicious of, so it is named rather than averaged over.
nMulti = nnz(snobs > 1);
if nMulti > 0
    fprintf(['\nselected cells seen in >1 FOV: %d (position = mean of %d..%d obs)\n' ...
             '  spread median %.1f um, max %.1f um\n'], ...
            nMulti, min(snobs(snobs>1)), max(snobs), ...
            median(sspread(snobs>1)), max(sspread(snobs>1)));
    bad = find(sspread > 30);
    for b = bad(:)'
        fprintf('  WARNING cell %d: %d obs spread over %.0f um\n', ...
                selIds(b), snobs(b), sspread(b));
    end
end
% Only cells that were SUPPOSED to be drawn count as a placement failure. The IO
% and non-significant-Sert cells are absent on purpose and must not be reported
% as if the map had lost them.
wanted  = [SEL.R.cell];
wanted  = wanted(arrayfun(@(c) selDraw(c), wanted));
% Deep cells were removed deliberately above, so they are not placement
% failures -- reporting them as 'could NOT be placed' would read as a bug.
missing = setdiff(wanted, selIds);
if exist('deepFovs','var')
    dropped = wanted(arrayfun(@(c) ~isempty(find([R.sel_cell]==c,1)) && ...
        all(ismember({R([R.sel_cell]==c).fov}, deepFovs)), wanted));
    missing = setdiff(missing, dropped);
end
if ~isempty(missing)
    fprintf('\n%d cell(s) meant to be drawn could NOT be placed: %s\n', ...
            numel(missing), mat2str(missing));
else
    fprintf('\nall %d cells meant to be drawn were placed\n', numel(wanted));
end
% Sert is no longer gated on significance, so IO is the only withheld group.
fprintf('withheld on purpose: %d IO\n', nnz(string({SEL.R.group}) == "IO"));
ie = find([F.cfgEx]);
if ~isempty(ie)
    deepCells = unique([R(ismember({R.fov}, {F(ie).fov}) & ~isnan([R.sel_cell])).sel_cell]);
    fprintf(['\ndepth-barred fields: %d drawn as dots but NOT as anatomy ' ...
             '(cells %s)\n'], numel(ie), mat2str(deepCells));
end

% No phase colorbar: the dots carry GROUP now, not phase, and leaving an HSV bar
% beside them would invite reading a colour off it that is not on it.

% origin + title
plot(ax_p, 0, 0, 'r+', 'MarkerSize',16, 'LineWidth',2);
title(ax_p, sprintf(['left: anatomy, %d FOVs / %d sites   |   right: %d cells ' ...
      'colour = genotype, marker = class.  IO not shown; ventral surface only'], ...
      nnz(~[F.cfgEx]), nnz([F.drawPatch]), nSel));

% legend
gn = fieldnames(patch_tints);
fprintf('\n--- masks vs CELLS per group ---\n');
fprintf('%-8s %8s %8s  %s\n','group','masks','cells','note');
for gi = 1:numel(gn)
    m    = strcmp({R.group}, gn{gi});
    nAll = sum(m);
    if nAll == 0, continue; end
    nCell = numel(unique({R(m).cell_key}));
    matched = nnz(~startsWith({R(m).cell_key}, 'mask:'));
    fprintf('%-8s %8d %8d  %s\n', gn{gi}, nAll, nCell, ...
        ternary(matched==0, 'no matcher run -- 1 cell per mask', ...
                sprintf('%d masks cross-FOV matched', matched)));
end
% The legend now describes the DOTS, not the anatomy patches. The patch tints are
% one-per-genotype and the dots are one-per-group, so a single legend cannot label
% both; the dots are what the figure is about, and the patch colours are already
% obvious from the blobs themselves.
ok = isgraphics(hG);
hL = hG(ok);  lL = cellstr(lG(ok));
% inside the axes, not 'eastoutside': outside it lands on top of the colorbar and
% hides its 'preferred phase (rad)' label. The caudal-right quadrant is empty.
if ~isempty(hL)
    legend(ax_p, hL, lL, 'Location','southeast','FontSize',8, 'Color','w');
end

%% ---- save ----
% Own directory and own names. The 260808 map's files are NOT touched: this is a
% different population on the same anatomy and the two must stay comparable.
% Own folder: the 260816 map is a different encoding of a different (deep-
% inclusive) population and has to stay on disk to compare against.
outDir = fullfile(sumRoot,'ventral_map_genotype_260817');
if ~isfolder(outDir), mkdir(outDir); end
opengl('software');   % NVIDIA path silently writes blank 9,223-byte PNG stubs

% CROP, as specified: x -2500..2500, y -1000..3000. Set AFTER every plot call so
% nothing can widen it again, and used for the main export.
xlim(ax_p, [-2500 2500]); ylim(ax_p, [-1000 3000]);
exportgraphics(figOut, fullfile(outDir,'ventral_map_genotype.png'), ...
               'Resolution',300,'BackgroundColor','white');
exportgraphics(figOut, fullfile(outDir,'ventral_map_genotype.pdf'), ...
               'ContentType','vector','BackgroundColor','white');

% ZOOM. At full-brainstem scale the whole dataset occupies about a fifth of the
% canvas, so the white outlines and the smaller Sst fields are below a pixel.
% Same figure, limits pulled in to the region that actually has data, exported
% separately -- the full map stays the one that shows where that region sits.
xl0 = xlim(ax_p); yl0 = ylim(ax_p);
padZ = 150;
xz = [min([-rx(:); [F(~[F.cfgEx]).x_um].']) - padZ, max(rx(:)) + padZ];
yz = [min([ry(:); [F(~[F.cfgEx]).y_um].']) - padZ, max([ry(:); [F(~[F.cfgEx]).y_um].']) + padZ];
xlim(ax_p, xz); ylim(ax_p, yz);
title(ax_p, sprintf(['zoom: %d cells, colour = genotype, marker = class ' ...
      '(IO not shown, ventral surface only)'], nSel));
% GUARDED. exportgraphics on this machine intermittently throws a Java
% NullPointerException out of HGRasterOutputHelper (the NVIDIA/OpenGL path), and
% it did so here on 2026-08-17 at Resolution 400. Unguarded, that aborts the
% script BEFORE the CSV and the .mat below -- so a cosmetic zoom export took the
% actual data products with it. The zoom is the least important output on the
% page; it must not be able to cost the rest.
try
    exportgraphics(figOut, fullfile(outDir,'ventral_map_genotype_zoom.png'), ...
                   'Resolution',400,'BackgroundColor','white');
    exportgraphics(figOut, fullfile(outDir,'ventral_map_genotype_zoom.pdf'), ...
                   'ContentType','vector','BackgroundColor','white');
catch ME
    fprintf(2,'ZOOM EXPORT FAILED (%s) -- main map, CSV and .mat still written\n', ME.message);
end
% Guarded: closing the figure window mid-run must not cost the CSV and the .mat
% that come after it. Restoring the limits is cosmetic, saving the data is not.
if isgraphics(ax_p), xlim(ax_p, [-2500 2500]); ylim(ax_p, [-1000 3000]); end

% Per-cell placement table, so a dot can be traced back to a cell.
Tm = table(selIds(:), sx, sy, snobs, sspread, rad2deg(th(:)), ssig, ...
           'VariableNames',{'cell','x_um_left','y_um','n_obs','spread_um', ...
                            'phase_deg','significant'});
Tm.group = string(cellfun(@(c) selGroup(c), num2cell(selIds(:)), 'uni',0));
Tm.class = string(cellfun(@(c) selClass(c), num2cell(selIds(:)), 'uni',0));
Tm = sortrows(Tm,{'group','class','cell'});
writetable(Tm, fullfile(outDir,'ventral_map_genotype_cells.csv'));
save(fullfile(outDir,'ventral_map_genotype_data.mat'), 'R','F','Tm','xd','yd', ...
     'bx','by','theta_rad','um_per_px_stitch','ax','ay', ...
     'rotated_sessions','AXMAP_DEFAULT','AXMAP_ROT90','mergeSiteUm', ...
     'sel_alpha','sel_critZ');
fprintf('Saved map PNG/PDF + CSV + .mat to %s\n', outDir);

%% ---- where every rotated-mount FOV landed (so the fix is checkable) ----
ir = find([F.rot90]);
if ~isempty(ir)
    fprintf('\n--- rotated-mount sessions: placement ---\n');
    fprintf('%-8s %-7s %-40s %9s %9s %6s %s\n', ...
            'sess','site','recording','lat(um)','ros(um)','patch','footprint');
    for i = ir
        fprintf('%-8s %-7s %-40s %9.0f %9.0f %6s %4.0f x %4.0f um\n', ...
            F(i).sess, F(i).site, F(i).fov(1:min(40,end)), F(i).x_um, F(i).y_um, ...
            string(F(i).drawPatch), F(i).H_fov*F(i).px_um_fov, F(i).W_fov*F(i).px_um_fov);
    end
end

%% ========================= local helpers =========================
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

% smooth_closed_contour was here. It smoothed the vector ROI outlines that this
% map no longer draws, so it went with them rather than sitting dead in the file.

function [cx, cy] = stage_to_canvas_ax(s2c, axmap, mx, my)
% Stage coordinates -> stitched-canvas pixels, via the session's stage->anatomy
% axis map. axmap sends [stage x; stage y] to [lateral; rostral]; the lateral
% component is then folded into the LEFT hemisphere, exactly as before, so a
% non-rotated session reduces to the original stage_to_canvas(-abs(mx), my).
v = axmap * [mx; my];
[cx, cy] = s2c(-abs(v(1)), v(2));
end
