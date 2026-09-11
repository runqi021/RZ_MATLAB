% breath_class_dff_traces_260825.m
% -----------------------------------------------------------------------
%  TWO FIGURES of the same individual-cycle dF/F, inspiration-onset triggered,
%  drawn in the per-cell summary GUI's style: every trial faint, mean bold on
%  top (P.trialOverlay = true), not a mean +/- SD ribbon.
%
%    FIGURE 1  class_dff_all_trials.png
%       coloured by CLASS.  GASP = red,  normal breath = black.
%
%    FIGURE 2  class_dff_ampcolor.png
%       the same traces coloured by that cycle's BREATH AMPLITUDE, jet.
%       Amplitude is peak-minus-foot on the breath trace after that trace has
%       been normalised 0-1 WITHIN ITS RECORDING, so cycles from different
%       recordings of one cell sit on a common scale. Traces are drawn in
%       ascending amplitude, so the big breaths end up on top rather than
%       buried under the small ones.
%
%  WINDOW: one common axis for every cell -- total width = trigWinIBI x the
%  median IBI across all labelled recordings, centred on inspiration onset.
%  Every cycle is interpolated onto that axis, which also absorbs the small fps
%  differences between recordings (29.95-30.01) and is what lets cycles from
%  different recordings of the same cell be stacked at all.
%
%  dF/F IS PLOTTED RAW, not baseline-subtracted -- it is a dF/F plot. Set
%  subtractBaseline = true to remove a short pre-onset mean instead, which takes
%  out the GCaMP carry-over offset from the previous cycle.
%
%  INDEX BASE: cyc_*_idx are RAW pc1 frames and the Ca traces have nDrop frames
%  tossed off the front, so the Ca-frame onset is cyc_foot_idx - nDrop. nDrop is
%  read from the label file's own params rather than assumed.
%
%  Runqi Zhang / 2026-08-25
% -----------------------------------------------------------------------

clear; close all; clc;

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir); addpath(repoRoot);
addpath(fullfile(repoRoot, 'analysis_260806'));

%% ===================== USER-EDITABLE =====================
regFile  = 'D:\Ventral_surface_summary\event_latency_260811\event_latency_data.mat';
srcRoot  = 'C:\260824_Vglut2-soma-g8s_vagotomized\phys';
outDir   = fullfile(srcRoot, 'breath_class_dff_260825');

cellList = [283 284 285 286 287 289 290 292 296 297];

trigWinIBI       = 2;       % total window width, in median IBIs
subtractBaseline = false;   % true = subtract mean over basePre_sec before onset
basePre_sec      = 0.10;
fallback_fps     = 30;

colNorm = [0 0 0];          % black
colGasp = [0.85 0.10 0.10]; % red
trialAlpha = 0.12;
trialLW    = 0.4;
meanLW     = 1.8;

% Figure 2 colour limits. 'data' = the observed amplitude range, which uses the
% whole colormap; 'unit' = [0 1], the full normalisation range, which is more
% literal but leaves most of jet unused because no cycle spans the whole trace.
cLimMode  = 'data';
ampAlpha  = 0.55;
ampLW     = 0.5;

nCols  = 5;
doSave = true;
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');
try, opengl('software'); catch, end
if doSave && ~isfolder(outDir), mkdir(outDir); end

NORM = int8(0); GASP = int8(1);

D = load(regFile,'CELL','OBS','REC');
obsOfCell = pooled_obs_260814(D.CELL, D.OBS);
fprintf('\n=========== breath_class_dff_traces_260825 ===========\n');

%% ---- ONE common time axis for every cell ------------------------------
ibiAll = [];  recSeen = {};
for c = cellList
    if c > numel(obsOfCell) || isempty(obsOfCell{c}), continue; end
    for o = obsOfCell{c}(:)'
        q  = regexp(D.OBS(o).label,'/','split');
        rn = strjoin(q(3:end-1),'/');
        if any(strcmp(recSeen, rn)), continue; end
        cf = fullfile(srcRoot, rn, 'breath_cycle_class_pc1.mat');
        if ~isfile(cf), continue; end
        L = load(cf,'cyc_ibi_s');
        ibiAll(end+1,1) = median(L.cyc_ibi_s,'omitnan'); %#ok<SAGROW>
        recSeen{end+1}  = rn; %#ok<SAGROW>
    end
end
assert(~isempty(ibiAll), 'no labelled recording found under %s', srcRoot);
winSec = trigWinIBI * median(ibiAll);
tAx    = linspace(-winSec/2, winSec/2, round(winSec*30)+1);
fprintf('window: %.2f s total (%g x median IBI %.3f s), %d samples\n', ...
        winSec, trigWinIBI, median(ibiAll), numel(tAx));

%% ---- collect every cycle ----------------------------------------------
C = struct('cell',{},'t',{},'A',{},'cls',{},'amp',{});
for c = cellList
    if c > numel(obsOfCell) || isempty(obsOfCell{c}), continue; end
    A = []; cls = []; amp = [];
    for o = obsOfCell{c}(:)'
        p = regexp(D.OBS(o).label,'/','split');
        recName = strjoin(p(3:end-1),'/');  roi = str2double(p{end});
        fp = fullfile(srcRoot, recName);
        cf = fullfile(fp,'breath_cycle_class_pc1.mat');
        df = dir(fullfile(fp,'*_ch1_dFF.mat'));
        bp = fullfile(fp,'breath_peak_pc1.mat');
        if ~isfile(cf) || isempty(df) || ~isfile(bp)
            fprintf(2,'  cell %d: %s missing labels / dFF / breath\n', c, recName); continue;
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

        % ---- breath trace normalised 0-1 WITHIN this recording ----------
        % Detrended first: the PC1/fb baseline drifts over a recording, and a
        % min-max on a drifting trace normalises the drift rather than the
        % breathing.
        bw = detrend(double(BP.breath(:)));
        bw = (bw - min(bw)) / max(max(bw) - min(bw), eps);

        T = size(dff,1);  x = dff(:,roi);
        b = max(1, round(basePre_sec*fps));
        on   = L.cyc_foot_idx(:) - nDrop;      % RAW pc1 frames -> Ca frames
        onRw = L.cyc_foot_idx(:);              % stay RAW for the breath trace
        pkRw = L.cyc_peak_idx(:);
        cl   = int8(L.class_final(:));
        pad  = ceil(winSec/2*fps) + 2;
        ok   = on - pad >= 1 & on + pad <= T & on - b >= 1 & ...
               onRw >= 1 & pkRw <= numel(bw);
        on = on(ok); onRw = onRw(ok); pkRw = pkRw(ok); cl = cl(ok);

        for i = 1:numel(on)
            fr  = double(on(i)) + tAx*fps;          % common offsets, own fps
            seg = interp1((1:T)', x, fr, 'linear');
            if subtractBaseline
                seg = seg - mean(x(on(i)-b : on(i)-1), 'omitnan');
            end
            A(end+1,:)   = seg; %#ok<SAGROW>
            cls(end+1,1) = cl(i); %#ok<SAGROW>
            amp(end+1,1) = bw(pkRw(i)) - bw(onRw(i)); %#ok<SAGROW>
        end
    end
    if isempty(A), continue; end
    C(end+1) = struct('cell',c,'t',tAx,'A',A,'cls',cls,'amp',amp); %#ok<SAGROW>
    fprintf('  cell %3d : %3d normal, %3d gasp  |  norm amp %.2f - %.2f\n', ...
            c, nnz(cls==NORM), nnz(cls==GASP), min(amp), max(amp));
end
assert(~isempty(C), 'nothing to plot');

%% ======================= FIGURE 1: by class ===========================
n = numel(C); nr = ceil(n/nCols);
f1 = figure('Color','w','Units','centimeters','Position',[1 1 5.2*nCols, 4.6*nr+1.2]);
set(f1,'DefaultAxesFontSize',7.5);
for i = 1:n
    ax = subplot(nr,nCols,i); hold(ax,'on');
    t = C(i).t;  N = C(i).A(C(i).cls==NORM,:);  G = C(i).A(C(i).cls==GASP,:);
    % normal first, gasp on top: whichever is drawn second wins every overlap,
    % and gasp is the set being asked about
    drawTrials(ax, t, N, colNorm, trialAlpha, trialLW);
    drawTrials(ax, t, G, colGasp, trialAlpha, trialLW);
    if ~isempty(N), plot(ax, t, mean(N,1,'omitnan'), '-','Color',colNorm,'LineWidth',meanLW); end
    if ~isempty(G), plot(ax, t, mean(G,1,'omitnan'), '-','Color',colGasp,'LineWidth',meanLW); end
    xline(ax, 0, ':','Color',[0.4 0.4 0.4]);
    xlim(ax,[t(1) t(end)]);
    title(ax, sprintf('cell %d    n=%d / %d', C(i).cell, size(N,1), size(G,1)), 'FontWeight','normal');
    if mod(i-1,nCols)==0, ylabel(ax,'\DeltaF/F'); end
    if i > n-nCols, xlabel(ax,'s from inspiration onset'); end
    box(ax,'off');
end
sgtitle('individual cycles:  black = normal breath,  red = GASP   (bold = mean)','FontSize',9);
saveFig(f1, outDir, 'class_dff_all_trials', subtractBaseline, doSave);

%% ================= FIGURE 2: coloured by breath amplitude =============
allAmp = vertcat(C.amp);
if strcmpi(cLimMode,'unit'), clim2 = [0 1];
else,                        clim2 = [min(allAmp) max(allAmp)];
end
f2 = figure('Color','w','Units','centimeters','Position',[1 1 5.2*nCols, 4.6*nr+1.6]);
set(f2,'DefaultAxesFontSize',7.5);
cm = jet(256);
for i = 1:n
    ax = subplot(nr,nCols,i); hold(ax,'on');
    t = C(i).t;  A = C(i).A;  a = C(i).amp;
    [~, ord] = sort(a, 'ascend');       % big breaths drawn last, so on top
    for k = ord(:)'
        u   = (a(k)-clim2(1)) / max(clim2(2)-clim2(1), eps);
        col = cm(1 + round(255*min(max(u,0),1)), :);
        plot(ax, t, A(k,:), '-', 'Color',[col ampAlpha], 'LineWidth',ampLW);
    end
    xline(ax, 0, ':','Color',[0.3 0.3 0.3]);
    xlim(ax,[t(1) t(end)]);
    caxis(ax, clim2); colormap(ax, cm);
    title(ax, sprintf('cell %d    n=%d', C(i).cell, size(A,1)), 'FontWeight','normal');
    if mod(i-1,nCols)==0, ylabel(ax,'\DeltaF/F'); end
    if i > n-nCols, xlabel(ax,'s from inspiration onset'); end
    box(ax,'off');
    if i == n
        cb = colorbar(ax); cb.Label.String = 'breath amplitude (0-1 per recording)';
    end
end
sgtitle(sprintf(['individual cycles coloured by breath amplitude   |   ' ...
    'peak-foot on the 0-1 normalised breath trace   |   colour range %.2f-%.2f'], ...
    clim2(1), clim2(2)), 'FontSize',9);
saveFig(f2, outDir, 'class_dff_ampcolor', subtractBaseline, doSave);

fprintf('Done.\n');

% -----------------------------------------------------------------------
function drawTrials(ax, t, M, col, a, lw)
if isempty(M), return; end
for k = 1:size(M,1)
    plot(ax, t, M(k,:), '-', 'Color',[col a], 'LineWidth',lw);
end
end

function saveFig(f, outDir, stem, subBase, doSave)
if subBase, stem = [stem '_baselinesub']; end
if ~doSave, return; end
exportgraphics(f, fullfile(outDir,[stem '.png']),'Resolution',300,'BackgroundColor','white');
exportgraphics(f, fullfile(outDir,[stem '.pdf']),'ContentType','vector','BackgroundColor','white');
fprintf('saved %s\n', fullfile(outDir,[stem '.png']));
close(f);
end
