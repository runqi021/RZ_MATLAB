% vglut2_nonIO_cellcount_map_260815.m
% -----------------------------------------------------------------------
%  Non-IO Vglut2 only: how many cells were RECORDED in each experiment, how many
%  of those were ACTIVE, and where they sat on the brainstem cartoon.
%
%  TWO COUNTS, TWO DIFFERENT SOURCES -- they are not the same population:
%
%  RECORDED = every segmented ROI that survived cross-FOV curation, collapsed to
%             cells by cell_key. Read from ROI_on_cartoon_data_stitch260801.mat,
%             which is the run that already applied the 260801 stitch offsets and
%             dropped masks tossed during curation. A mask is not a cell: the same
%             neuron re-imaged at another z or zoom yields one mask per recording,
%             so unique cell_key is counted, never masks.
%
%  ACTIVE   = cells in the event-latency registry passing the rate gate
%             (>= activeMinRateHz, pooled events / pooled duration). The registry
%             only ever contained cells that cleared the archive's own event
%             floor, so ACTIVE is a subset of RECORDED by construction.
%
%  IO IS EXCLUDED by group label, not by folder. The Vglut2/1124 session imaged
%  both pFN and IO sites; the IO sites carry group 'IO' and are a different
%  structure, so they are not Vglut2 cortex-adjacent cells and do not belong in
%  this count.
%
%  THE MAP draws POSITIONS ONLY -- no phase, no coherence colouring. Every cell is
%  one marker at the mean position of its ROIs. Cells are folded onto ONE side
%  (x -> -|x|) and the cartoon is cropped to that half. In this dataset every ROI
%  is already at negative x, so the fold changes nothing here; it is there so the
%  figure stays correct if a session on the other side is ever added.
%
%  Runqi Zhang / 2026-08-15
% -----------------------------------------------------------------------

clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);  addpath(repoRoot);
addpath(fullfile(repoRoot,'analysis_260806'));
addpath(fullfile(repoRoot,'analysis_260727','coh_ca_breath'));

%% ===================== USER-EDITABLE =====================
cartoonMat = 'C:\Users\Admin\Desktop\brainstem_map_cartoon_coordframe.mat';
roiDataMat = 'C:\Users\Admin\Desktop\ROI_on_cartoon_data_stitch260801.mat';
sumRoot    = 'D:\Ventral_surface_summary';
regFile    = fullfile(sumRoot,'event_latency_260811','event_latency_data.mat');
outDir     = fullfile(sumRoot,'vglut2_nonIO_map_260815');

GROUP      = 'Vglut2';       % group label to keep; 'IO' sites are a different label
activeMinRateHz = 2/60;      % must match per_cell_summary_260812 / pop_features_260813

markRecorded = [0.72 0.72 0.72];
markActive   = [0.85 0.10 0.10];
% =========================================================

if ~isfolder(outDir), mkdir(outDir); end
set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end     % NVIDIA export path drops ~1% of PNGs

%% ===================== RECORDED CELLS =====================
assert(isfile(roiDataMat), 'ROI map data not found:\n  %s', roiDataMat);
S = load(roiDataMat, 'R','xd','yd');
R = S.R;
gAll = string({R.group}');
R    = R(gAll == GROUP);
assert(~isempty(R), 'no ROIs with group "%s"', GROUP);

dts   = string({R.date}');
fovs  = string({R.fov}');
rids  = [R.rid]';
ckey  = string({R.cell_key}');
xum   = [R.x_um]';
yum   = [R.y_um]';

% ROI-level key, to join against the registry
roiKey = GROUP + "/" + dts + "/" + fovs + "/" + string(rids);

%% ---- repair cell identity where cell_link was never run --------------------
% Vglut2/0224 and Vglut2/1124 have a CURATED cross-FOV matcher output but no
% cell_link.mat, so the map data fell back to one cell per mask ('mask:...') and
% counting unique cell_key would report 107 and 48 cells where the curation says
% 73 and fewer. The curation is read directly here.
%
% TRAP: match.fov_name is fileparts-truncated (everything after the first dot is
% eaten as an extension) and COLLIDES between recordings in some sessions, so the
% recording name is taken from match.fov_folder's path leaf instead.
curated_sources = { '0224', fullfile(sumRoot,GROUP,'0224','pFN','roi_match_out_260727')
                    '1124', fullfile(sumRoot,GROUP,'1124','pFN','roi_match_out_260727') };
for s = 1:size(curated_sources,1)
    dstr = curated_sources{s,1};  base = curated_sources{s,2};
    fRes = fullfile(base,'roi_match_results.mat');
    fCur = fullfile(base,'roi_match_curated.mat');
    if ~isfile(fRes) || ~isfile(fCur), continue; end
    M = load(fRes,'match');  Cu = load(fCur);
    if isfield(Cu,'curated') && isfield(Cu.curated,'grpOf'), gof = Cu.curated.grpOf(:);
    else,                                                    gof = M.match.grp(:); end
    fov_i = M.match.roi.fov(:);  roi_i = M.match.roi.roi(:);
    % grpOf has THREE states, and conflating them inflates the cell count:
    %   > 0   grouped   -> one cell shared by every member
    %   == 0  ungrouped -> a real cell, seen in one recording only
    %   < 0   TOSSED    -> rejected during curation, NOT a cell. Must be dropped.
    % Rejection is encoded NEGATIVE, not NaN. Counting those as ungrouped gave
    % 0224 87 cells instead of 73 (14 rejected) and 1124 36 instead of 30 (6
    % rejected), because each rejected mask became its own cell.
    % Cross-check: 0224 is 56 grouped + 37 ungrouped + 14 negative = 107 masks,
    % and 36 groups + 37 ungrouped = 73 cells.
    nFixed = 0; nToss = 0;
    for k = 1:numel(gof)
        parts   = regexp(M.match.fov_folder{fov_i(k)}, '[\\/]', 'split');
        recName = parts{end};
        key = GROUP + "/" + string(dstr) + "/" + string(recName) + "/" + string(roi_i(k));
        j = find(roiKey == key, 1);
        if isempty(j), continue; end
        if isnan(gof(k)) || gof(k) < 0
            ckey(j) = "";  nToss = nToss + 1;              % dropped below
        elseif gof(k) > 0
            ckey(j) = "cur:" + string(dstr) + "#g" + string(gof(k));
        else
            ckey(j) = "cur:" + string(dstr) + "#u" + string(k);
        end
        nFixed = nFixed + 1;
    end
    fprintf('curation applied to %s/%s : %d ROIs re-keyed, %d tossed in curation\n', ...
            GROUP, dstr, nFixed, nToss);
end

% Drop curation-tossed ROIs everywhere before anything is counted or plotted.
drop = (ckey == "");
if any(drop)
    fprintf('dropping %d curation-tossed ROI(s)\n', nnz(drop));
    keepM = ~drop;
    dts = dts(keepM); fovs = fovs(keepM); rids = rids(keepM);
    ckey = ckey(keepM); xum = xum(keepM); yum = yum(keepM);
    roiKey = roiKey(keepM);
end

%% ===================== ACTIVE CELLS =====================
Din   = load(regFile, 'CELL','OBS','REC');
obsOf = pooled_obs_260814(Din.CELL, Din.OBS);
labs  = string({Din.OBS.label}');

% Rate per registry cell, measured the same way the figures measure it. Cheap:
% statsOnly, and only the cells of this group.
Pst = struct('doCoh',false,'nDrop',30,'fallback_fps',30,'TW_spec',6,'alpha_sig',0.01, ...
    'minSpikes',2,'ca_lag_sec',0,'f_breath_search',[0.2 4],'fwhm_factor',0.6, ...
    'min_bw',0.05,'fmin',0.05,'fmax',15,'trigWin_sec',[],'trigWinIBI',2, ...
    'ylim_dff',[],'ylim_epc',[],'histBinFrames',2,'nShuffle',0,'shiftMinCyc',3, ...
    'pad_um',20,'clip_pct',[0.5 99.9],'scalebar_um',50,'rayPhaseBins',36,'featFs',30, ...
    'gamma_val',1,'PixelSizeBase',1.7778,'outlineLW',0.8,'sortMode','none', ...
    'dffColor',[0.2 0.7 0.2],'onsetCol',[0.9 0.1 0.1],'peakCol',[0.35 0.75 1], ...
    'statsOnly',true);

activeRoiKey = strings(0,1);
nCellChecked = 0;  nCellActive = 0;
for c = 1:numel(obsOf)
    if isempty(obsOf{c}), continue; end
    p1 = regexp(Din.OBS(obsOf{c}(1)).label,'/','split');
    if ~strcmp(p1{1}, GROUP), continue; end          % IO sites carry group 'IO'
    O = struct('folder',{},'roi',{},'recName',{},'group',{},'recDate',{});
    keysHere = strings(0,1);
    for o = obsOf{c}(:)'
        p = regexp(Din.OBS(o).label,'/','split');
        fp = Din.REC(Din.OBS(o).rec).folder;
        if ~isfolder(fp), continue; end
        O(end+1) = struct('folder',fp,'roi',str2double(p{end}), ...
            'recName',strjoin(p(3:end-1),'/'),'group',p{1},'recDate',p{2}); %#ok<SAGROW>
        keysHere(end+1,1) = labs(o); %#ok<SAGROW>
    end
    if isempty(O), continue; end
    nCellChecked = nCellChecked + 1;
    try
        [~,~,st] = temporal_phase_cell_fig_260812(O, Pst);
    catch ME
        fprintf(2,'  stats failed for a %s cell: %s\n', GROUP, ME.message);  continue;
    end
    if isfinite(st.rateHz) && st.rateHz >= activeMinRateHz
        nCellActive  = nCellActive + 1;
        activeRoiKey = [activeRoiKey; keysHere]; %#ok<AGROW>
    end
end
fprintf('registry: %d %s cells measured, %d active (>= %.3g ev/min)\n\n', ...
        nCellChecked, GROUP, nCellActive, activeMinRateHz*60);

isActiveRoi = ismember(roiKey, activeRoiKey);

%% ===================== COLLAPSE ROIs -> CELLS =====================
[uk, ~, ic] = unique(ckey);
nCell = numel(uk);
cx = accumarray(ic, xum, [nCell 1], @mean);
cy = accumarray(ic, yum, [nCell 1], @mean);
cAct = accumarray(ic, double(isActiveRoi), [nCell 1], @max) > 0;
cSess = strings(nCell,1);
for i = 1:nCell
    j = find(ic == i, 1);
    cSess(i) = dts(j);
end

%% ===================== TABLE =====================
uSess = unique(cSess);
T = table('Size',[numel(uSess) 5], ...
          'VariableTypes',{'string','double','double','double','double'}, ...
          'VariableNames',{'session','n_ROI','n_recorded_cells','n_active_cells','pct_active'});
for i = 1:numel(uSess)
    m  = cSess == uSess(i);
    mr = dts   == uSess(i);
    T.session(i)          = GROUP + "/" + uSess(i);
    T.n_ROI(i)            = nnz(mr);
    T.n_recorded_cells(i) = nnz(m);
    T.n_active_cells(i)   = nnz(m & cAct);
    T.pct_active(i)       = 100*nnz(m & cAct)/max(nnz(m),1);
end
T = sortrows(T,'session');
T(end+1,:) = {"TOTAL", sum(T.n_ROI), sum(T.n_recorded_cells), sum(T.n_active_cells), ...
              100*sum(T.n_active_cells)/max(sum(T.n_recorded_cells),1)};
disp(T);

% NOT every session's cell count means the same thing. Say so in the output
% rather than letting an incomparable number pass as if it were comparable.
fprintf(['\nCAVEAT on n_recorded_cells:\n' ...
    '  0810  full cross-FOV link (160 obs -> 81 cells, 28 seen in >1 recording)\n' ...
    '  0224  curated matcher output applied here (14 tossed masks excluded)\n' ...
    '  1124  curated matcher output applied here (6 tossed masks excluded)\n' ...
    '  0728  cell_link registers only 155 of 466 masks, from 6 of 16 recordings,\n' ...
    '        and merges NOTHING (every observation is its own cell). Its count is\n' ...
    '        therefore one-cell-per-mask and is an OVERCOUNT relative to the other\n' ...
    '        three. Re-running the matcher on 0728 is what would fix it.\n']);
writetable(T, fullfile(outDir,'vglut2_nonIO_cell_counts.csv'));

%% ===================== MAP =====================
C = load(cartoonMat);
assert(isfile(C.imgPath), 'cartoon image missing: %s', C.imgPath);
I = im2double(imread(C.imgPath));
if size(I,3) == 1, I = repmat(I,1,1,3); end
I = imgaussfilt(I, 1.2);
[Hc, Wc, ~] = size(I);
xd = ([1 Wc] - C.origin_pix(1)) / C.pix_per_x;
yd = ([1 Hc] - C.origin_pix(2)) / C.pix_per_y;

% Fold every cell onto the negative-x side, then crop the cartoon to that half.
cxF = -abs(cx);
nFlip = nnz(cx > 0);

fig = figure('Color','w','Units','normalized','Position',[0.05 0.05 0.42 0.86]);
ax  = axes(fig); hold(ax,'on');
image(ax, 'XData',xd, 'YData',yd, 'CData',I);
set(ax,'YDir','normal'); axis(ax,'image');
xlim(ax,[min(xd) 0]);                 % the other side is cut off
ylim(ax,[min(yd) max(yd)]);
plot(ax,[0 0], ylim(ax), '-', 'Color',[0 0 0 0.35], 'LineWidth',1);   % midline

hR = scatter(ax, cxF(~cAct), cy(~cAct), 16, 'o', ...
    'MarkerFaceColor',markRecorded,'MarkerEdgeColor','none','MarkerFaceAlpha',0.75);
hA = scatter(ax, cxF(cAct), cy(cAct), 26, 'o', ...
    'MarkerFaceColor',markActive,'MarkerEdgeColor','none','MarkerFaceAlpha',0.9);

legend(ax, [hR hA], ...
    {sprintf('recorded, not active (n=%d)', nnz(~cAct)), ...
     sprintf('active (n=%d)', nnz(cAct))}, ...
    'Location','southwest','Box','off','FontSize',9);
xlabel(ax,'lateral (\mum)'); ylabel(ax,'rostral (\mum)');
title(ax, sprintf('%s, non-IO: %d cells recorded, %d active, %d experiments', ...
      GROUP, nCell, nnz(cAct), numel(uSess)), 'FontWeight','normal');
box(ax,'on');

exportgraphics(fig, fullfile(outDir,'vglut2_nonIO_cell_map.png'), ...
               'Resolution',600,'BackgroundColor','white');
exportgraphics(fig, fullfile(outDir,'vglut2_nonIO_cell_map.pdf'), ...
               'ContentType','vector','BackgroundColor','white');

save(fullfile(outDir,'vglut2_nonIO_map_data.mat'), 'T','cx','cy','cxF','cAct','cSess','uk');
fprintf('\nfolded %d cell(s) from the positive-x side\n', nFlip);
fprintf('saved -> %s\n', outDir);
