% breath_cycle_dff_heatmap_260825.m
% -----------------------------------------------------------------------
%  TEST SCRIPT, one cell.  Per-cycle dF/F heatmap, inspiration-onset triggered,
%  shown TWICE side by side in ONE figure -- the same cycles, the same colour
%  scale, differing only in the row order:
%
%     LEFT   rows sorted by the cycle's BREATH amplitude
%     RIGHT  rows sorted by the cycle's dF/F peak
%
%  Each side carries BOTH metrics as connecting lines, so whichever ordered the
%  rows, you can see what the other one did. If the two ranked the cycles the
%  same way the two heatmaps would be identical and both lines monotonic on
%  both sides; the extent to which they are not is the answer.
%
%  BOTH METRICS ARE MEASURED OVER THE SAME SPAN: this cycle's inspiration ONSET
%  to this cycle's inspiratory PEAK.
%      breath amplitude = bwz(peak_i) - bwz(foot_i)      (robust-z breath trace)
%      dF/F peak        = max( dFF over [foot_i , peak_i] )
%  Measuring them over the same window is what makes the comparison fair -- a
%  dF/F peak taken over a wider span would pick up the response to a different
%  breath and could rank cycles by something the breath amplitude never saw.
%  dffWin below can widen it to a peak-to-peak span if you want the response
%  that FOLLOWS the inspiration rather than the one during it.
%
%  Breath amplitude is the VALUE at the two detected landmarks, the same
%  definition the labelling GUI thresholds on and the pipeline's ampFrac QC
%  uses, so a number here means the same thing as a number there. ampDef =
%  'minmax' switches to max-min within the cycle, which is robust to a
%  mis-detected landmark but also picks up noise.
%
%  Heatmap is GRAYSCALE, as in the per-cell summary.
%
%  INDEX BASE: cyc_*_idx are RAW pc1 frames; the Ca traces have nDrop frames
%  tossed off the front, so the Ca-frame onset is cyc_foot_idx - nDrop. nDrop is
%  read from the label file's params rather than assumed.
%
%  Runqi Zhang / 2026-08-26
% -----------------------------------------------------------------------

clear; close all; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir); addpath(repoRoot);
addpath(fullfile(repoRoot, 'analysis_260806'));

%% ===================== USER-EDITABLE =====================
regFile  = 'D:\Ventral_surface_summary\event_latency_260811\event_latency_data.mat';
srcRoot  = 'C:\260824_Vglut2-soma-g8s_vagotomized\phys';
outDir   = fullfile(srcRoot, 'breath_class_dff_260825');

cellList = [283 284 285 286 287 289 290 292 296 297];   % one figure per cell

trigWinIBI = 2;             % display window, total width in median IBIs
sortDir    = 'descend';     % largest at the TOP

% Span the dF/F peak is taken over:
%   'onset_peak'  [ foot_i , peak_i ]              same span as the breath amp
%   'peak_next'   [ peak_i , peak_{i+1} ]          the response that FOLLOWS it
%   'peak_prev'   [ peak_{i-1} , peak_i ]
dffWin = 'onset_peak';

% Breath amplitude definition:
%   'peak_value'  the breath trace VALUE at the inspiratory peak, on a trace
%                 normalised 0-1 within its recording. 0 = the recording's
%                 lowest point, 1 = its largest breath. Comparable across the
%                 recordings of one cell without any absolute PC1 units.
%   'peak_onset'  bwz(peak_i) - bwz(foot_i) on the robust-z trace (the rise
%                 during that inspiration; what the labelling GUI thresholds on)
%   'minmax'      max - min of bwz within the cycle
ampDef = 'peak_value';
% 'log_peak' takes the 0-1 peak value, logs it, then re-normalises the logs to
% 0-1 across that recording's cycles. The raw peak distribution is skewed --
% gasps sit near 1 and normal breaths are squashed into the bottom fifth -- so a
% linear scale spends most of its range on a handful of cycles. The log spreads
% the small breaths out; it changes the SPACING of the rows' sort key, never
% their ORDER, since log is monotonic.

fallback_fps = 30;
cmapName  = 'gray';
clipPct   = [1 99];         % colour limits from these percentiles
markClass = false;          % true = colour the line markers by class
tickFrac  = 0.42;           % raster tick half-height, in heatmap rows
doSave    = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end
if doSave && ~isfolder(outDir), mkdir(outDir); end

NORM = int8(0); GASP = int8(1);
COL_B = [0 0 0];            % breath amplitude line
% Onset/peak colours copied from temporal_phase_cell_gui_260812 (P.onsetCol /
% P.peakCol), so an onset mark means the same thing here as it does there.
COL_ONSET = [1.00 0.40 0.75];   % pink
COL_PEAK  = [0.35 0.55 1.00];   % cornflower
COL_D = [0 0 0];            % dF/F peak line

D = load(regFile,'CELL','OBS','REC');
obsOfCell = pooled_obs_260814(D.CELL, D.OBS);
fprintf('\n=========== breath_cycle_dff_heatmap_260825 ===========\n');
fprintf('dF/F peak over %s | breath amp %s | sorted %s\n', dffWin, ampDef, sortDir);

for cellId = cellList(:)'
if cellId > numel(obsOfCell) || isempty(obsOfCell{cellId})
    fprintf(2,'cell %d is not in the registry -- skipped\n', cellId); continue;
end

%% ---- one common time axis for this cell -------------------------------
% A cell imaged twice has two recordings with different median IBIs and
% slightly different frame rates, so a window sized per recording gives rows of
% different length that cannot be stacked. Size it once from the cell's own
% median IBI and interpolate every cycle onto it -- which also absorbs the
% 29.95-30.01 fps differences.
ibis = [];
for o = obsOfCell{cellId}(:)'
    q  = regexp(D.OBS(o).label,'/','split');
    cfq = fullfile(srcRoot, strjoin(q(3:end-1),'/'), 'breath_cycle_class_pc1.mat');
    if ~isfile(cfq), continue; end
    Lq = load(cfq,'cyc_ibi_s');
    ibis(end+1,1) = median(Lq.cyc_ibi_s,'omitnan'); %#ok<SAGROW>
end
if isempty(ibis)
    fprintf(2,'cell %d: no labelled recording -- skipped\n', cellId); continue;
end
winSec = trigWinIBI * median(ibis);
tAx    = linspace(-winSec/2, winSec/2, round(winSec*30)+1);

%% ---- gather every cycle of this cell ---------------------------------
M = []; ampZ = []; dffPk = []; cls = []; pkLat = []; recOf = []; recUsed = {};   % tAx is already set above
for o = obsOfCell{cellId}(:)'
    p = regexp(D.OBS(o).label,'/','split');
    recName = strjoin(p(3:end-1),'/');  roi = str2double(p{end});
    fp = fullfile(srcRoot, recName);
    cf = fullfile(fp,'breath_cycle_class_pc1.mat');
    df = dir(fullfile(fp,'*_ch1_dFF.mat'));
    bp = fullfile(fp,'breath_peak_pc1.mat');
    if ~isfile(cf) || isempty(df) || ~isfile(bp)
        fprintf(2,'  %s: missing labels / dFF / breath -- skipped\n', recName); continue;
    end
    L  = load(cf);
    Dd = load(fullfile(df(1).folder, df(1).name),'dFF');
    BP = load(bp,'breath');
    dff = double(Dd.dFF);
    if roi < 1 || roi > size(dff,2), continue; end

    fps = detect_session_fps(fp, fallback_fps);
    nDrop = 30;
    if isfield(L,'class_params') && isfield(L.class_params,'nDrop')
        nDrop = L.class_params.nDrop;
    end
    pad = ceil(winSec/2*fps) + 2;
    T = size(dff,1);  x = dff(:,roi);

    % Two normalisations of the same detrended trace, for the two amplitude
    % definitions. Detrended first in both cases: the PC1 baseline drifts over a
    % recording, so a min-max on a drifting trace would normalise the drift
    % rather than the breathing.
    bwd = detrend(double(BP.breath(:)));
    bwz = (bwd - median(bwd)) / max(mad(bwd,1)*1.4826, eps);          % robust z
    bw01 = (bwd - min(bwd)) / max(max(bwd) - min(bwd), eps);          % 0-1 per recording

    onRw = L.cyc_foot_idx(:);        % RAW pc1 frames
    pkRw = L.cyc_peak_idx(:);
    on   = onRw - nDrop;             % Ca frames
    pk   = pkRw - nDrop;
    cl   = int8(L.class_final(:));
    nC   = numel(on);

    for i = 1:nC
        % ---- dF/F peak span, in Ca frames ---------------------------
        switch lower(dffWin)
            case 'onset_peak', a = on(i);   b = pk(i);
            case 'peak_next',  if i >= nC, continue; end, a = pk(i); b = pk(i+1);
            case 'peak_prev',  if i <= 1,  continue; end, a = pk(i-1); b = pk(i);
            otherwise, error('dffWin must be onset_peak, peak_next or peak_prev');
        end
        if a < 1 || b > T || b <= a, continue; end
        if on(i)-pad < 1 || on(i)+pad > T, continue; end
        if onRw(i) < 1 || pkRw(i) > numel(bwz), continue; end

        % ---- breath amplitude, same cycle ---------------------------
        switch lower(ampDef)
            case {'peak_value','log_peak'}
                aZ = bw01(pkRw(i));    % log_peak is transformed after the scan
            case 'peak_onset'
                aZ = bwz(pkRw(i)) - bwz(onRw(i));
            case 'minmax'
                if i >= nC, continue; end
                sRw = onRw(i):onRw(i+1);
                sRw = sRw(sRw >= 1 & sRw <= numel(bwz));
                if numel(sRw) < 3, continue; end
                aZ = max(bwz(sRw)) - min(bwz(sRw));
            otherwise, error('ampDef must be log_peak, peak_value, peak_onset or minmax');
        end

        % sample at the common offsets, in this recording's own frame rate
        M(end+1,:)     = interp1((1:T)', x, double(on(i)) + tAx*fps, 'linear'); %#ok<SAGROW>
        dffPk(end+1,1) = max(x(a:b)); %#ok<SAGROW>
        ampZ(end+1,1)  = aZ; %#ok<SAGROW>
        cls(end+1,1)   = cl(i); %#ok<SAGROW>
        pkLat(end+1,1) = (pk(i) - on(i))/fps;   % inspiratory peak, s after onset %#ok<SAGROW>
        recOf(end+1,1) = numel(recUsed) + 1; %#ok<SAGROW>
    end
    recUsed{end+1} = recName; %#ok<SAGROW>
end
if isempty(M)
    fprintf(2,'cell %d: no usable cycles -- skipped\n', cellId); continue;
end
if strcmpi(ampDef,'log_peak')
    % Per RECORDING, so a cell imaged twice is not put on one linear scale by
    % the back door. floor at 1e-3 of that recording's largest peak: a cycle
    % whose peak sits essentially at the trace minimum would otherwise send
    % log() to -Inf and collapse every other cycle to the top of the range.
    for r = 1:max(recOf)
        m = recOf == r;
        if nnz(m) < 2, continue; end
        v = ampZ(m);
        v = log(max(v, max(v)*1e-3));
        ampZ(m) = (v - min(v)) / max(max(v) - min(v), eps);
    end
end
nR = size(M,1);
fprintf('%d cycles from %d recording(s): %d normal, %d gasp\n', ...
        nR, numel(recUsed), nnz(cls==NORM), nnz(cls==GASP));

rho = corr(dffPk, ampZ, 'type','Spearman','rows','complete');
fprintf('Spearman(dF/F peak, breath amplitude) = %+.3f\n', rho);

[~, ordB] = sort(ampZ,  sortDir);      % LEFT  block
[~, ordD] = sort(dffPk, sortDir);      % RIGHT block
cl2 = prctile(M(:), clipPct);

%% ---------------------------- FIGURE ----------------------------------
f = figure('Color','w','Units','centimeters','Position',[1 1 28 14]);
set(f,'DefaultAxesFontSize',8);

axA = draw_block(f, 0.060, tAx, M(ordB,:), cl2, cmapName, ...
        ampZ(ordB), dffPk(ordB), cls(ordB), markClass, ...
        sprintf('sorted by %s', ternary_lab(ampDef)), COL_B, COL_D, ampDef, pkLat(ordB), COL_ONSET, COL_PEAK, tickFrac);

draw_block(f, 0.545, tAx, M(ordD,:), cl2, cmapName, ...
        ampZ(ordD), dffPk(ordD), cls(ordD), markClass, ...
        'sorted by max(\DeltaF/F)', COL_B, COL_D, ampDef, pkLat(ordD), COL_ONSET, COL_PEAK, tickFrac);

cb = colorbar(axA,'Location','southoutside');
cb.Label.String = '\DeltaF/F';
cb.Position = [0.060 0.075 0.205 0.020];   % same width as the heatmap

if markClass
    annotation(f,'textbox',[0.58 0.015 0.40 0.045], ...
        'String','markers: black = normal breath,  red = GASP', ...
        'EdgeColor','none','FontSize',8,'HorizontalAlignment','center');
end
sgtitle(sprintf('cell %d   |   %d cycles   |   %s', cellId, nR, ...
    strjoin(recUsed,', ')), 'FontSize',8.5, 'Interpreter','none');

if doSave
    stem = sprintf('cycle_dff_heatmap_cell%03d_bothsorts_%s', cellId, dffWin);
    exportgraphics(f, fullfile(outDir,[stem '.png']),'Resolution',300,'BackgroundColor','white');
    exportgraphics(f, fullfile(outDir,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');
    fprintf('saved %s\n', fullfile(outDir,[stem '.png']));
end
close(f);        % 10 cells would otherwise leave 10 figures open
end
fprintf('Done.\n');

% =======================================================================
function axH = draw_block(f, x0, tAx, M, cl2, cmapName, ampZ, dffPk, cls, ...
                          markClass, ttl, colB, colD, ampDef, pkLat, colOn, colPk, tickFrac)
%DRAW_BLOCK  Heatmap + BOTH metric traces, sharing the heatmap's row axis.
%  Both traces are drawn on every block regardless of which one sorted it: the
%  sorted one is monotonic by construction and tells you nothing new, the other
%  one is the whole point of the panel.
% The heatmap axes is the same width as the colourbar drawn under it -- the
% image itself carries no more information for being stretched, and the two
% blocks read as a pair more easily when neither dominates.
nR = size(M,1);
wH = 0.205; wK = 0.085; y0 = 0.17; hh = 0.70;

axH = axes(f,'Position',[x0 y0 wH hh]);
imagesc(axH, tAx, 1:nR, M);
if cl2(2) > cl2(1), caxis(axH, cl2); end
colormap(axH, cmapName);
hold(axH,'on');
% ONSET and PEAK as RASTERS -- one short tick per cycle, nothing joining
% them, the same style as the GUI's per-cycle onset row
% (P.onsetRowStyle = 'raster'). A connected line would imply the peak times
% form a continuous function of row order, which they do not: the rows are
% ordered by amplitude, not by time.
raster(axH, zeros(nR,1), tickFrac, colOn, 1.2);   % inspiration ONSET
raster(axH, pkLat,       tickFrac, colPk, 1.2);   % inspiratory PEAK, per cycle
set(axH,'YDir','reverse');            % row 1 = largest, at the TOP
xlabel(axH,'time from inspiration onset (s)');
ylabel(axH,'breath cycle');
title(axH, ttl, 'FontWeight','normal','Interpreter','tex');

axB = keyAxes(f, [x0+wH+0.008 y0 wK hh], ampZ, cls, markClass, colB, nR);
xlabel(axB, ternary_lab(ampDef));

axD = keyAxes(f, [x0+wH+wK+0.016 y0 wK hh], dffPk, cls, markClass, colD, nR);
xlabel(axD, 'max(\DeltaF/F)', 'Interpreter','tex');
end

function ax = keyAxes(f, pos, v, cls, markClass, col, nR)
ax = axes(f,'Position',pos); hold(ax,'on');
plot(ax, v, 1:nR, '-','Color',col,'LineWidth',1);
if markClass
    mN = cls==0; mG = cls==1;
    plot(ax, v(mN), find(mN), 'o','MarkerFaceColor',[0 0 0], ...
         'MarkerEdgeColor','none','MarkerSize',2.5);
    plot(ax, v(mG), find(mG), 'o','MarkerFaceColor',[0.85 0.10 0.10], ...
         'MarkerEdgeColor','none','MarkerSize',2.5);
end
ylim(ax,[0.5 nR+0.5]); set(ax,'YTick',[],'YDir','reverse');
box(ax,'off');
end

function raster(ax, x, half, col, lw)
%RASTER  One short vertical tick per row, drawn as a single NaN-separated line.
%  One object instead of nR objects: with a few hundred cycles per panel that is
%  the difference between an instant redraw and a visible stall.
n = numel(x);
X = [x(:)'; x(:)'; nan(1,n)];
Y = [(1:n)-half; (1:n)+half; nan(1,n)];
plot(ax, X(:), Y(:), '-', 'Color', col, 'LineWidth', lw);
end

function lb = ternary_lab(ampDef)
if strcmpi(ampDef,'log_peak'), lb = 'log max(motion)'; else, lb = 'max(motion)'; end
end
