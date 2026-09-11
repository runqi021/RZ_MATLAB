function [fig, info] = trace_rawF_panel_260817(OBS, P)
%TRACE_RAWF_PANEL_260817  One panel: the GUI's wide trace, but RAW F not dF/F.
%
%   [fig, info] = trace_rawF_panel_260817(OBS, P)
%
% OBS : the same struct array the summary figure takes -- .folder .roi .recName
%       .group .recDate. Only ONE of them is drawn (see P.traceObs), because a
%       raw trace cannot be pooled: F is in camera counts and two recordings of
%       the same cell at different zoom, laser power and PMT gain are simply not
%       on the same scale. dF/F is what makes them comparable, and that is
%       exactly what has been taken out here.
%
% P   : the same parameter struct. Used: traceObs, trace_xlim_sp, nDrop,
%       fallback_fps, dffColor, traceBreathCol, cellId, showAcqTime.
%
% Returns the figure and an info struct (recording, roi, window, fps, counts).
%
% WHY A SEPARATE FILE. This draws nothing but the trace, and it must not be able
% to change what temporal_phase_cell_fig_260812 renders for the batch. Adding a
% raw-F branch inside that function would put a switch in the middle of the panel
% every batch figure is made from; a standalone panel cannot affect it at all.
%
% RAW F IS NOT BASELINE-CORRECTED, and that is the point of asking for it: slow
% bleaching, focus drift and any illumination change stay in the trace. Read the
% y-axis as counts, not as activity. F_roi_raw, F_roi and the cpSAM F are all the
% same array (verified identical, 2026-08-17) -- one signal under three names --
% so "raw F" is unambiguous here: the mask-mean fluorescence Cellpose extracted.
%
% ALIGNMENT is taken from the SAME file and the SAME truncation the dF/F trace
% uses (F_roi_raw sits beside dFF in *_ch1_dFF.mat with identical [T x N] shape),
% so the raw panel and the dF/F panel cannot drift apart by a frame.
%
% Runqi Zhang / 2026-08-17

%% ---- which recording ----
iT = 1;
if isfield(P,'traceObs') && ~isempty(P.traceObs)
    iT = min(max(round(P.traceObs), 1), numel(OBS));
end
o = OBS(iT);

%% ---- load, exactly as the summary figure's loader does ----
df = dir(fullfile(o.folder,'*_ch1_dFF.mat'));
bp = dir(fullfile(o.folder,'breath_peak_pc1.mat'));
assert(~isempty(df),'No *_ch1_dFF.mat in %s', o.folder);
assert(~isempty(bp),'No breath_peak_pc1.mat in %s', o.folder);

nDrop = 30;  if isfield(P,'nDrop') && ~isempty(P.nDrop), nDrop = P.nDrop; end
fbFps = 30;  if isfield(P,'fallback_fps') && ~isempty(P.fallback_fps), fbFps = P.fallback_fps; end
fps = detect_session_fps(o.folder, fbFps);

D = load(fullfile(df(1).folder, df(1).name));
assert(isfield(D,'F_roi_raw'), 'No F_roi_raw in %s', df(1).name);
F_all = double(D.F_roi_raw);
assert(o.roi>=1 && o.roi<=size(F_all,2), 'ROI %d out of range (1..%d) in %s', ...
       o.roi, size(F_all,2), o.folder);

BP = load(fullfile(bp(1).folder, bp(1).name));
bw = detrend(double(BP.breath(:)));
bw(1:min(nDrop,numel(bw))) = [];
bw = bw - mean(bw);

nB = numel(BP.breath);
if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==nB
    ev = double(BP.insp_onsets_train(:) ~= 0);
else
    ev = zeros(nB,1); oi = round(BP.insp_onset_idx(:)); ev(oi(oi>=1 & oi<=nB)) = 1;
end
ev(1:min(nDrop,numel(ev))) = [];

% Vglut2/1124 rising-edge trigger: breath and events are one frame late. Same
% gate as the summary figure -- on the GENOTYPE FOLDER, not the group label,
% because that session's IO sites carry the group name 'IO'.
gFix = '';
if isfield(o,'group'), gFix = o.group; end
if isempty(gFix) || strcmpi(gFix,'IO'), gFix = local_genotype(o.folder); end
if strcmpi(gFix,'Vglut2') && isfield(o,'recDate') && strcmp(o.recDate,'1124')
    bw = [bw(1); bw(1:end-1)];
    ev = [0; ev(1:end-1)];
end

T  = min([size(F_all,1), numel(bw), numel(ev)]);
Fr = F_all(1:T, o.roi);  bw = bw(1:T);
t  = (0:T-1)'/fps;

%% ---- window ----
if isfield(P,'trace_xlim_sp') && ~isempty(P.trace_xlim_sp)
    wsp = [max(min(P.trace_xlim_sp), t(1)), min(max(P.trace_xlim_sp), t(end))];
else
    mid = (t(1)+t(end))/2;  halfW = min(15, (t(end)-t(1))/2);
    wsp = [mid-halfW, mid+halfW];
end
m = t >= wsp(1) & t <= wsp(2);
assert(nnz(m) > 2, 'window %.1f-%.1f s contains %d samples', wsp(1), wsp(2), nnz(m));

%% ---- draw ----
dffColor = [0.2 0.7 0.2];
if isfield(P,'dffColor') && ~isempty(P.dffColor), dffColor = P.dffColor; end
bwCol = [0.6 0.6 0.6];
if isfield(P,'traceBreathCol') && ~isempty(P.traceBreathCol), bwCol = P.traceBreathCol; end

fig = figure('Color','w','Units','centimeters','Position',[2 2 24 6.5], ...
             'Name','raw F trace','NumberTitle','off');
set(fig,'DefaultAxesFontSize',9,'DefaultTextFontSize',9);
ax = axes(fig);
tRel = t - wsp(1);

% Single-axis layout, the same trick the GUI trace uses: yyaxis always paints the
% right side above the left, so a dark breath line on its own axis would cover
% the data it is context for. The breath is rescaled onto the F axis instead and
% pushed to the bottom of the child order; the right axis is kept for its LABEL.
yyaxis(ax,'right'); set(ax,'YColor',bwCol,'YTick',[]); ylabel(ax,'breath');
yyaxis(ax,'left');  hold(ax,'on');
plot(ax, tRel, Fr, '-', 'Color',dffColor, 'LineWidth',0.8);
set(ax,'YColor','k'); ylabel(ax,'F');
lo = min(Fr(m)); hi = max(Fr(m)); pad = 0.05*max(hi-lo, eps);
ylim(ax, [lo-pad hi+pad]);
dl = ylim(ax);
wbw = bw(m); if isempty(wbw), wbw = bw; end
bwS = (bw - min(wbw)) / max(max(wbw)-min(wbw), eps);
bwS = dl(1) + bwS * (dl(2) - dl(1));
hBw = plot(ax, tRel, bwS, '-', 'Color',bwCol, 'LineWidth',0.6);
uistack(hBw,'bottom');
ylim(ax, dl);                       % the breath must not rescale the F axis
hold(ax,'off');
xlim(ax, [0 wsp(2)-wsp(1)]); xlabel(ax,'Time (s)'); box(ax,'off');

%% ---- title, same identity the summary figure prints ----
cellTag = '';
if isfield(P,'cellId') && ~isempty(P.cellId) && ~isnan(P.cellId)
    cellTag = sprintf('cell %d | ', P.cellId);
end
genoStr = '';
if isfield(o,'group') && ~isempty(o.group), genoStr = o.group; end
showAcq = false;
if isfield(P,'showAcqTime') && ~isempty(P.showAcqTime), showAcq = logical(P.showAcqTime); end
if showAcq
    dt = local_acq_datetime(o.folder);
    if ~isnat(dt)
        genoStr = sprintf('%s %s %s', genoStr, datestr(dt,'yymmdd'), datestr(dt,'HH:MM'));
    elseif isfield(o,'recDate') && ~isempty(o.recDate)
        genoStr = sprintf('%s %s', genoStr, o.recDate);
    end
end
if ~isempty(genoStr), genoStr = [genoStr ' | ']; end
title(ax, {sprintf('%s%s%s  ROI%d', cellTag, genoStr, o.recName, o.roi), ...
           sprintf('raw F | %gHz | trace %.1f-%.1f s', ...
                   round(fps), wsp(1), wsp(2))}, ...
      'Interpreter','none','FontWeight','bold','FontSize',9);

info = struct('folder',o.folder,'recName',o.recName,'roi',o.roi, ...
              'traceObs',iT,'nObsAvailable',numel(OBS),'fps',fps, ...
              'window_s',wsp,'T',T,'F_min',min(Fr(m)),'F_max',max(Fr(m)), ...
              'F_median',median(Fr(m)));
end

% =======================================================================
function g = local_genotype(folderPath)
%LOCAL_GENOTYPE  Genotype from the archive path, for the Vglut2/1124 frame fix.
g = '';
pp = regexp(regexprep(folderPath,'[\\/]+$',''), '[\\/]', 'split');
known = {'ChAT','Vglut2','Vgat','Sst','Sert'};
for k = numel(pp):-1:1
    j = find(strcmpi(pp{k}, known), 1);
    if ~isempty(j), g = known{j}; return; end
end
end

function dt = local_acq_datetime(folderPath)
%LOCAL_ACQ_DATETIME  ScanImage epoch from the recording's _meta.mat.
dt = NaT;
mh = dir(fullfile(folderPath,'*_ch1_meta.mat'));
if isempty(mh), mh = dir(fullfile(folderPath,'*_meta.mat')); end
if isempty(mh), return; end
M = load(fullfile(mh(1).folder, mh(1).name));
if isfield(M,'epoch') && numel(M.epoch) >= 6
    e = M.epoch;
    try, dt = datetime(e(1),e(2),e(3),e(4),e(5),floor(e(6))); catch, dt = NaT; end
end
end
