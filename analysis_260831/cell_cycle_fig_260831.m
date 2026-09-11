function [fig, proj, stats] = cell_cycle_fig_260831(OBS, P)
%CELL_CYCLE_FIG_260831  Per-cell PER-CYCLE figure. Drop-in for the per-cell GUI.
%
%   [fig, proj, stats] = cell_cycle_fig_260831(OBS, P)
%
% Same signature, same OBS, same P as temporal_phase_cell_fig_260812, so
% cell_cycle_gui_260831.m is that GUI with one line changed and every control --
% folder/id/Load, pooled vs single ROI, trig win, dFF ylim, crop um, scale-bar
% corner, clip, window, Prev/Next, save -- behaves identically.
%
% WHAT IS DIFFERENT: the figure. Five panels instead of the full summary.
%
%   row 1, four SQUARE plot boxes, all the same size and aligned:
%     1  avg projection, ROI outline, scale bar
%     2  breath strip (top 1/8) over the per-cycle dF/F heatmap (bottom 7/8)
%     3  onset-triggered, EVERY CYCLE as its own line
%     4  peak-triggered, likewise
%   row 2, one box the SAME height as a square, full width:
%     5  wide dF/F trace with the breath waveform over it
%
% WHY SPAGHETTI. The summary figure draws mean +/- SD. A handful of large cycles
% and a consistent small one give the same mean; drawing every cycle tells them
% apart. Nothing is smoothed and no mean is drawn on top -- a bold mean over
% faint cycles reads as the data.
%
% THE CROP AND THE STATISTICS ARE NOT REIMPLEMENTED. This calls
% temporal_phase_cell_fig_260812 once, with nShuffle forced to 0, and takes its
% `proj` (crop, ROI outline, scale bar) and `stats` (IBI, dOnPk, counts), then
% closes the figure it drew. That costs a few seconds per render and is worth it:
% the crop here is the crop the batch figures use, at the same microns per pixel
% and the same contrast, and cannot drift from them. Only the per-cycle epoch
% matrices are rebuilt here, because that function keeps them local.
%
% THE EPOCHS MIRROR ITS POOLED-EPOCH BLOCK: bw = detrend(BP.breath) with the
% first P.nDrop dropped and de-meaned; PEAK = insp_onset_idx (breath_peak_pc1),
% ONSET = insp_start_idx (breath_insp_start_pc1); T = min(dFF, breath, events)
% applied before epoching; the Vglut2/1124 rising-edge one-frame shift, gated on
% the genotype FOLDER because that session's IO sites carry the group label 'IO';
% and cycles cut wider than the window then interpolated inwards so interp1
% never extrapolates. The magenta dots are each cycle's OWN inspiration onset,
% each peak taking the onset that PRECEDES it, NaN (undrawn) if there is none.
%
% Runqi Zhang / 2026-08-31

%% ---- defaults ----
if ~isfield(P,'nDrop'),        P.nDrop = 30;          end
if ~isfield(P,'fallback_fps'), P.fallback_fps = 30;   end
if ~isfield(P,'dffColor'),     P.dffColor = [0.35 0.75 0.35]; end
% Copied from cell_cycle_gui_260831.m:243-247 so a direct call looks like a GUI
% render. These are the pink/cornflower pair RZ set on 2026-08-17; the BATCH
% (per_cell_summary_260812.m) still uses a red onset, so do not "fix" one to
% match the other without deciding which is canonical.
if ~isfield(P,'onsetCol'),     P.onsetCol = [1.00 0.40 0.75]; end
if ~isfield(P,'peakCol'),      P.peakCol  = [0.35 0.55 1.00]; end
traceObs = 1; if isfield(P,'traceObs') && ~isempty(P.traceObs), traceObs = P.traceObs; end
traceObs = max(1, min(traceObs, numel(OBS)));

%% ---- canonical pass: crop + stats ----
% This pass DRAWS the full summary figure -- that is the only way to get `proj`,
% because roi_crop_local runs after the statsOnly early return. The figure is
% closed immediately, but without the visibility guard below it still flashes on
% screen, which reads as "the GUI made two figures". DefaultFigureVisible is set
% off for the duration and restored by onCleanup, so it is restored even if the
% call errors.
Pc = P; Pc.nShuffle = 0; Pc.statsOnly = false;
proj = struct('have_proj',false); stats = struct();
vis0 = get(0,'DefaultFigureVisible');
restoreVis = onCleanup(@() set(0,'DefaultFigureVisible',vis0));
set(0,'DefaultFigureVisible','off');
f0 = [];
try
    [f0, proj, stats] = temporal_phase_cell_fig_260812(OBS, Pc);
catch ME
    fprintf(2,'cell_cycle_fig_260831: canonical pass failed (%s)\n', ME.message);
end
if ~isempty(f0) && all(ishandle(f0)), close(f0); end
clear restoreVis;                      % restore visibility before OUR figure

%% ---- trigger window ----
% P.trigWin_sec is a TOTAL width, matching temporal_phase_cell_fig_260812:63-69.
IBI = 1; if isfield(stats,'IBI') && isfinite(stats.IBI) && stats.IBI > 0, IBI = stats.IBI; end
if isfield(P,'trigWin_sec') && ~isempty(P.trigWin_sec) && P.trigWin_sec > 0
    winTot = P.trigWin_sec;
else
    nIBI = 2; if isfield(P,'trigWinIBI') && P.trigWinIBI > 0, nIBI = P.trigWinIBI; end
    winTot = nIBI * IBI;
end
winSec = winTot/2;

%% ---- per-cycle epochs ----
D = cell_epochs(OBS, P, winSec, traceObs);
if isempty(D.Ep) && isempty(D.E)
    fig = figure('Color','w'); axis off;
    text(0.5,0.5,'no complete breath cycles','Horizontal','center'); return;
end

%% ---- canvas ----
% Designed on 1700x1150 and scaled to whatever the window actually gets. It must
% be read back: a MATLAB figure is capped at the screen, and applying the design
% unscaled to a clamped canvas silently makes the "squares" rectangular and puts
% anything positioned in raw pixels tens of px out of place.
% Margins are only what the decorations need: ~65 px under row 1 for its tick
% labels and xlabel, ~55 under row 2 for the same, ~40 above row 1 for the title.
% Anything more is dead space in the export.
% desH leaves TWO bands above the row-1 boxes, not one: ax_st carries its own
% axes title ("peak-triggered, n=... cycles") which is drawn ABOVE its box, so
% the figure title has to clear that, not just clear the boxes.
desW = 1580; desH = 810;
% SCALE shrinks the whole figure -- canvas, boxes, gaps AND fonts -- by one
% factor, so the proportions the layout was designed at are preserved. Scaling
% the canvas alone would leave the text at absolute point sizes and it would
% swallow the panels.
SCALE = 0.7; if isfield(P,'figScale') && ~isempty(P.figScale), SCALE = P.figScale; end
fig = figure('Color','w','Units','pixels', ...
             'Position',[20 40 round(desW*SCALE) round(desH*SCALE)], ...
             'Name','per-cycle cell figure','NumberTitle','off');
drawnow;
q = get(fig,'Position'); figW = q(3); figH = q(4);
% One uniform factor for both axes, so a design square stays square whatever the
% window ends up as -- a MATLAB canvas is capped at the screen and does not
% always grant the size asked for.
s  = min(figW/desW, figH/desH);
px = @(x,y,w,h) [(x*s)/figW, (y*s)/figH, (w*s)/figW, (h*s)/figH];
% Fonts scale with the layout. Set BEFORE any axes is created, or the defaults
% do not apply to it.
set(fig,'DefaultAxesFontSize', max(5, 10*s), ...
        'DefaultTextFontSize', max(5, 10*s));

A    = 300;                     % square side, design px
yB   = 421;                     % bottom edge of the row-1 squares
yB2  = 55;                      % bottom edge of row 2
% Column lefts. Each gap holds the ylabel + tick labels of the panel that
% FOLLOWS it; the gap after panel 2 is wider because it has to hold the
% colourbar, the colourbar's ticks AND its label before panel 3's ylabel starts.
% Too narrow and the colourbar label lands on top of panel 3's.
colX = [62 430 860 1225];
hStrip = A/8;  hHeat = A - hStrip;

ax_pr = axes(fig,'Position',px(colX(1), yB,       A, A));
ax_hm = axes(fig,'Position',px(colX(2), yB,       A, hHeat));
ax_st = axes(fig,'Position',px(colX(2), yB+hHeat, A, hStrip));
ax_on = axes(fig,'Position',px(colX(3), yB,       A, A));
ax_pk = axes(fig,'Position',px(colX(4), yB,       A, A));
ax_tr = axes(fig,'Position',px(colX(1), yB2, colX(4)+A-colX(1), A));   % SAME height as a square
for a = [ax_pr ax_hm ax_st ax_on ax_pk ax_tr]
    try
        a.PositionConstraint = 'innerposition';
    catch
        try, a.ActivePositionProperty = 'position'; catch, end
    end
end

%% ---- 1. avg projection ----
if isstruct(proj) && isfield(proj,'have_proj') && proj.have_proj
    % Padded to a square, NOT axis image: that would resize the axes to the
    % image's aspect ratio and break the alignment. Padding leaves the pixels
    % undistorted instead of stretching a non-square crop to fill the box.
    [img, dx, dy] = pad_square(proj.crop_img);
    imagesc(ax_pr, img); colormap(ax_pr, gray(256)); caxis(ax_pr,[0 1]);
    set(ax_pr,'YDir','reverse','XTick',[],'YTick',[],'XColor','none','YColor','none');
    xlim(ax_pr,[0.5 size(img,2)+0.5]); ylim(ax_pr,[0.5 size(img,1)+0.5]);
    hold(ax_pr,'on');
    % NO ROI OUTLINE (RZ, 2026-08-31). proj.bnd_crop still carries it, so
    % restoring it is a two-line change, but the crop is already centred on the
    % ROI and the outline only obscures the cell it is pointing at.
    if isfinite(proj.barLen_pr)
        W2 = size(img,2); H2 = size(img,1);
        plot(ax_pr, [W2*0.08 W2*0.08+proj.barLen_pr], [H2*0.93 H2*0.93], ...
             '-','Color','w','LineWidth',3);
    end
else
    axis(ax_pr,'off'); text(ax_pr,0.5,0.5,'(no projection)','Horizontal','center');
end

%% ---- 2. breath strip + per-cycle heatmap ----
if ~isempty(D.Epb)
    mb = mean(D.Epb,1,'omitnan');
    plot(ax_st, D.tau, mb, '-','Color',[0.3 0.3 0.3],'LineWidth',1.2);
    xlim(ax_st,[D.tau(1) D.tau(end)]); axis(ax_st,'off');
    hold(ax_st,'on');
    yl = [min(mb) max(mb)] + [-1 1]*0.08*max(range_(mb),eps);
    plot(ax_st,[0 0],yl,'-','Color',P.peakCol,'LineWidth',1);
    if isfinite(D.dOnPk), plot(ax_st,[D.dOnPk D.dOnPk],yl,'-','Color',P.onsetCol,'LineWidth',1); end
    ylim(ax_st,yl);
    text(ax_st,D.tau(1),yl(2),'breath','FontSize',max(5,8*s),'VerticalAlignment','top');
    title(ax_st,sprintf('peak-triggered, n=%d cycles',size(D.Ep,1)), ...
          'FontWeight','normal','FontSize',max(5,9*s));
else
    axis(ax_st,'off');
end
cb = [];
if ~isempty(D.Ep)
    cl = prctile(D.Ep(:),[0.5 99.5]); if ~all(isfinite(cl)) || cl(2)<=cl(1), cl = [min(D.Ep(:)) max(D.Ep(:))]; end
    imagesc(ax_hm, D.tau, 1:size(D.Ep,1), D.Ep); colormap(ax_hm, gray(256));
    caxis(ax_hm, cl); set(ax_hm,'YDir','reverse'); xlim(ax_hm,[D.tau(1) D.tau(end)]);
    hold(ax_hm,'on');
    ok = isfinite(D.onsRow);
    plot(ax_hm, D.onsRow(ok), find(ok), '.', 'Color',P.onsetCol,'MarkerSize',6);
    xlabel(ax_hm,'time from insp peak (s)'); ylabel(ax_hm,'breath # (chronological)');
    cb = colorbar(ax_hm); cb.Label.String = '\DeltaF/F';
else
    axis(ax_hm,'off');
end

%% ---- 3 & 4. per-cycle spaghetti ----
ylD = [];
if isfield(P,'ylim_dff') && numel(P.ylim_dff) == 2, ylD = P.ylim_dff; end
spaghetti(ax_on, D.tau, D.E,  P.dffColor, P.onsetCol, 'onset-triggered','time from insp onset (s)', ylD);
spaghetti(ax_pk, D.tau, D.Ep, P.dffColor, P.peakCol,  'peak-triggered', 'time from insp peak (s)',  ylD);
if isempty(ylD)
    yl = [min([ylim(ax_on) ylim(ax_pk)]) max([ylim(ax_on) ylim(ax_pk)])];
    ylim(ax_on,yl); ylim(ax_pk,yl);
end

%% ---- 5. wide trace ----
r = D.rec;
if ~isempty(r)
    tAll = (0:numel(r.dff)-1)/r.fps;
    if isfield(P,'trace_xlim_sp') && numel(P.trace_xlim_sp) == 2
        w0 = P.trace_xlim_sp(1); w1 = P.trace_xlim_sp(2);
    else
        w0 = max(0, tAll(end)/2 - 15); w1 = min(tAll(end), w0 + 30);
    end
    m = tAll >= w0 & tAll <= w1;
    yyaxis(ax_tr,'left');
    plot(ax_tr, tAll(m), r.dff(m), '-','Color',P.dffColor,'LineWidth',1.4);
    ylabel(ax_tr,'\DeltaF/F'); ax_tr.YColor = [0.20 0.45 0.20];
    if ~isempty(ylD), ylim(ax_tr, ylD); end
    yyaxis(ax_tr,'right');
    bb = r.bw(m); bb = (bb-min(bb))/max(max(bb)-min(bb),eps);
    plot(ax_tr, tAll(m), bb, '-','Color','k','LineWidth',1.0);
    % 0-1: the breath is min-max normalised, so that is its full range and any
    % headroom above it is just dead space. It does mean the breath and the dF/F
    % overlap across the whole panel height, which is the intent -- they are read
    % against each other, not stacked.
    ylabel(ax_tr,'breath'); ax_tr.YColor = [0.2 0.2 0.2]; ylim(ax_tr,[0 1]);
    yyaxis(ax_tr,'left');
    xlim(ax_tr,[w0 w1]); xlabel(ax_tr,'Time (s)'); box(ax_tr,'off');
    grid(ax_tr,'on'); ax_tr.GridAlpha = 0.12;
end

%% ---- title ----
cid = NaN; if isfield(P,'cellId'), cid = P.cellId; end
if isnan(cid), idStr = sprintf('ROI %d', OBS(1).roi); else, idStr = sprintf('cell %d', cid); end
ttl = sprintf('%s/%s   %s   |   %d cycles over %d recording(s)   |   trace: %s', ...
    OBS(1).group, OBS(1).recDate, idStr, size(D.Ep,1), numel(OBS), OBS(traceObs).recName);
% Clears the row-1 boxes AND ax_st's axes title above them. +8 put it straight
% through that label; +38 puts it in its own band.
annotation(fig,'textbox',px(0, yB+A+38, desW, 32), ...
    'String',ttl,'EdgeColor','none','HorizontalAlignment','center', ...
    'FontSize',max(6,11*s),'Interpreter','none','VerticalAlignment','middle');

% Colourbar slot re-applied LAST, and the heatmap box WITH it. Creating a
% colourbar shrinks its peer axes to make room -- PositionConstraint does not
% prevent that -- which took the heatmap from 233 px wide to 176 and broke the
% alignment with the other three squares. Restoring the axes Position after the
% colourbar is placed is what makes the four boxes actually identical.
if ~isempty(cb) && isvalid(cb)
    cb.Units = 'normalized';
    cb.Position = px(colX(2)+A+12, yB, 11, hHeat);
    ax_hm.Position = px(colX(2), yB, A, hHeat);
end
end

% =========================================================================
function spaghetti(ax, tau, M, col, trigCol, ttl, xlab, ylD)
cla(ax); hold(ax,'on');
if isempty(M), axis(ax,'off'); return; end
% Opacity scaled to the cycle count. At a fixed alpha, 30 cycles are readable
% and 130 saturate into a solid block that shows only the envelope -- which is
% what the mean already says. This keeps ink per pixel roughly constant.
n  = size(M,1);
aL = max(0.04, min(0.45, 12/n));
for k = 1:n
    plot(ax, tau, M(k,:), '-', 'Color',[col aL], 'LineWidth',0.5);
end
if numel(ylD) == 2
    yl = ylD;
else
    yl = prctile(M(:),[0.5 99.5]);
    if ~all(isfinite(yl)) || yl(2) <= yl(1), yl = [min(M(:)) max(M(:))]; end
    yl = yl + [-1 1]*0.08*max(diff(yl),eps);
end
plot(ax,[0 0],yl,'-','Color',trigCol,'LineWidth',1.4);
ylim(ax,yl); xlim(ax,[tau(1) tau(end)]);
title(ax,ttl,'FontWeight','normal'); xlabel(ax,xlab); ylabel(ax,'\DeltaF/F');
grid(ax,'on'); ax.GridAlpha = 0.15; box(ax,'off');
end

% =========================================================================
function D = cell_epochs(OBS, P, winSec, traceObs)
D = struct('tau',[],'E',[],'Ep',[],'Epb',[],'onsRow',[],'dOnPk',NaN,'rec',[]);
R = []; fpsAll = [];
for i = 1:numel(OBS)
    rr = load_rec(OBS(i), P.nDrop, P.fallback_fps);
    if isempty(rr), continue; end
    if isempty(R), R = rr; else, R(end+1) = rr; end %#ok<AGROW>
    fpsAll(end+1) = rr.fps; %#ok<AGROW>
end
if isempty(R), return; end
fpsRef = median(fpsAll);
nH  = max(1, round(winSec*fpsRef));
tau = (-nH:nH)/fpsRef;
E = []; Ep = []; Epb = []; onsRow = [];
useB = true(1,numel(R));
if isfield(P,'breathObs') && ~isempty(P.breathObs)
    useB = false(1,numel(R));
    b = P.breathObs(P.breathObs >= 1 & P.breathObs <= numel(R));
    useB(b) = true;
    if ~any(useB), useB = true(1,numel(R)); end
end
for i = 1:numel(R)
    r = R(i);
    w  = max(1, round(winSec*r.fps));
    tk = (-w:w)/r.fps;
    on = r.foot(r.foot-w>=1 & r.foot+w<=r.T);
    pk = r.peak(r.peak-w>=1 & r.peak+w<=r.T);
    E  = [E;  regrid(cut(r.dff,on,w), tk, tau)]; %#ok<AGROW>
    Ep = [Ep; regrid(cut(r.dff,pk,w), tk, tau)]; %#ok<AGROW>
    if useB(i)
        bwz = (r.bw - mean(r.bw))/max(std(r.bw),eps);
        Epb = [Epb; regrid(cut(bwz,pk,w), tk, tau)]; %#ok<AGROW>
    end
    if ~isempty(pk)
        oo = nan(numel(pk),1);
        bb = discretize(pk,[r.foot(:); inf]); okp = ~isnan(bb);
        oo(okp) = -(pk(okp) - r.foot(bb(okp)))/r.fps;
        onsRow = [onsRow; oo]; %#ok<AGROW>
    end
end
D.tau = tau; D.E = E; D.Ep = Ep; D.Epb = Epb; D.onsRow = onsRow;
D.dOnPk = median(onsRow,'omitnan');
D.rec = R(min(traceObs,numel(R)));
end

% =========================================================================
function r = load_rec(o, nDrop, fbFps)
r = [];
df = dir(fullfile(o.folder,'*_ch1_dFF.mat'));
bp = dir(fullfile(o.folder,'breath_peak_pc1.mat'));
ip = dir(fullfile(o.folder,'breath_insp_start_pc1.mat'));
if isempty(df) || isempty(bp), return; end
fps = detect_session_fps(o.folder, fbFps);
Dd  = load(fullfile(df(1).folder, df(1).name),'dFF');
BP  = load(fullfile(bp(1).folder, bp(1).name));
dff_all = double(Dd.dFF);
if o.roi < 1 || o.roi > size(dff_all,2), return; end

bw = detrend(double(BP.breath(:))); bw(1:min(nDrop,numel(bw))) = []; bw = bw - mean(bw);
nB = numel(BP.breath);
ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;
ev(1:min(nDrop,numel(ev))) = [];
ev_foot = [];
if ~isempty(ip)
    IP = load(fullfile(ip(1).folder, ip(1).name));
    ev_foot = zeros(nB,1); fi = round(IP.insp_start_idx(:));
    ev_foot(fi(fi>=1 & fi<=nB)) = 1; ev_foot(1:min(nDrop,numel(ev_foot))) = [];
end
geno = o.group;
if strcmpi(geno,'IO') && contains(o.folder,[filesep 'Vglut2' filesep],'IgnoreCase',true)
    geno = 'Vglut2';
end
if strcmpi(geno,'Vglut2') && strcmp(o.recDate,'1124')
    bw = [bw(1); bw(1:end-1)];
    ev = [0; ev(1:end-1)];
    if ~isempty(ev_foot), ev_foot = [0; ev_foot(1:end-1)]; end
end
T = min([size(dff_all,1), numel(bw), numel(ev)]);
if ~isempty(ev_foot), T = min(T, numel(ev_foot)); end
if T < 10, return; end
dff = dff_all(1:T, o.roi); bw = bw(1:T); ev = ev(1:T);
if ~isempty(ev_foot), ev_foot = ev_foot(1:T); else, ev_foot = zeros(T,1); end
r = struct('dff',dff,'bw',bw,'fps',fps,'T',T, ...
           'peak',find(ev>0),'foot',find(ev_foot>0),'name',o.recName);
end

% =========================================================================
function Eo = cut(x, idx, w)
Eo = zeros(numel(idx), 2*w+1);
for k = 1:numel(idx), Eo(k,:) = x(idx(k)-w : idx(k)+w); end
end

function Y = regrid(Eo, tk, tau)
if isempty(Eo), Y = zeros(0,numel(tau)); return; end
Y = interp1(tk, Eo.', tau, 'linear', NaN).';
if size(Y,2) ~= numel(tau), Y = reshape(Y, [], numel(tau)); end
end

function [img, dx, dy] = pad_square(im)
[h,w] = size(im);
n = max(h,w);
img = zeros(n,n,'like',im) + min(im(:));
dy = floor((n-h)/2); dx = floor((n-w)/2);
img(dy+(1:h), dx+(1:w)) = im;
end

function y = range_(x), y = max(x(:)) - min(x(:)); end
