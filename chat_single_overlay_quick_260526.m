% chat_single_overlay_quick_260526.m
% -----------------------------------------------------------------------
%  Quick single-recording overlays:
%    Fig 1: breath (top) + dF/F stacked (bottom)
%    Fig 2: dF/F (color) vs raw F (gray) overlaid, per ROI, stacked
%           (each trace normalized to [0,1] so they share a scale)
%  Just set folderPath and run.
%  Inputs: *_ch1_dFF.mat (dFF, F_roi), *DLC_*breath_peak_data.mat
%  Dependency: detect_session_fps.m
% -----------------------------------------------------------------------

clear; close all; clc;
addpath(fileparts(mfilename('fullpath')));

%% ===================== USER-EDITABLE =====================
folderPath = 'D:\batch_dffQC_test_260325\260322_sst_soma_g8s\phys\processed\maybe breathing\5x_-743-840-z30_00001';
nDrop      = 30;
fallback_fps = 30;
% =========================================================

[~, nm] = fileparts(folderPath);
df = dir(fullfile(folderPath, '*_ch1_dFF.mat'));
bp = dir(fullfile(folderPath, '*DLC*breath_peak_data.mat'));
assert(~isempty(df), 'No *_ch1_dFF.mat in %s', folderPath);

D   = load(fullfile(df(1).folder, df(1).name));
dFF = double(D.dFF);
if isfield(D,'F_roi'), Fraw = double(D.F_roi);
elseif isfield(D,'F_roi_raw'), Fraw = double(D.F_roi_raw);
else, Fraw = []; end
fps = detect_session_fps(folderPath, fallback_fps);
[Tca, N] = size(dFF);
t_ca = (0:Tca-1)'/fps;

norm01 = @(x) (x - min(x)) ./ (max(x) - min(x) + eps);

%% ===== Fig 1: breath + dF/F =====
if ~isempty(bp)
    BP = load(fullfile(bp(1).folder, bp(1).name));
    bw = detrend(double(BP.breath(:))); bw(1:nDrop) = []; bw = bw - mean(bw);
    if isfield(BP,'insp_onsets_train') && numel(BP.insp_onsets_train)==numel(BP.breath)
        ev = double(BP.insp_onsets_train(:) ~= 0);
    else
        ev = zeros(numel(BP.breath),1); oi = round(BP.insp_onset_idx(:));
        ev(oi(oi>=1 & oi<=numel(ev))) = 1;
    end
    ev(1:nDrop) = [];
    Tb = min([Tca, numel(bw), numel(ev)]);
    bw = bw(1:Tb); ev = ev(1:Tb); tb = (0:Tb-1)'/fps; ot = tb(ev>0);

    rng = max(dFF,[],1)-min(dFF,[],1); sp = max(prctile(rng,80),0.3);
    f1 = figure('Color','w','Name',[nm ' : dFF+breath'], ...
                'Units','normalized','Position',[0.05 0.1 0.9 0.8]);
    ax1 = subplot(4,1,1); hold(ax1,'on');
    plot(ax1, tb, bw, 'Color',[0 0.35 1], 'LineWidth',0.7);
    plot(ax1, ot, interp1(tb,bw,ot,'nearest'), 'r.', 'MarkerSize',8);
    ylabel(ax1,'breath'); title(ax1, sprintf('%s   (fps=%.2f, %d ROIs)', nm, fps, N), 'Interpreter','none');
    xlim(ax1,[0 tb(end)]);
    ax2 = subplot(4,1,2:4); hold(ax2,'on');
    for r=1:N, plot(ax2, t_ca, dFF(:,r)+(r-1)*sp, 'LineWidth',0.5); end
    ylim(ax2,[-sp N*sp]); set(ax2,'YTick',(0:N-1)*sp,'YTickLabel',1:N);
    xlabel(ax2,'Time (s)'); ylabel(ax2,'dF/F (stacked)'); xlim(ax2,[0 t_ca(end)]);
    linkaxes([ax1 ax2],'x');
else
    warning('No breath_peak_data.mat -> skipping Fig 1.');
end

%% ===== Fig 2: dF/F vs raw F overlay (per ROI) =====
if ~isempty(Fraw)
    f2 = figure('Color','w','Name',[nm ' : dFF+F'], ...
                'Units','normalized','Position',[0.05 0.1 0.9 0.8]);
    ax = axes(f2); hold(ax,'on');
    for r=1:N
        plot(ax, t_ca, norm01(Fraw(:,r)) + (r-1), 'Color',[0.6 0.6 0.6], 'LineWidth',0.5);
        plot(ax, t_ca, norm01(dFF(:,r))  + (r-1), 'Color',[0.85 0.1 0.1], 'LineWidth',0.5);
    end
    ylim(ax,[-0.2 N]); set(ax,'YTick',(0:N-1)+0.5,'YTickLabel',1:N);
    xlabel(ax,'Time (s)'); ylabel(ax,'ROI (each norm. 0–1)');
    title(ax, sprintf('%s   raw F (gray) vs dF/F (red)', nm), 'Interpreter','none');
    xlim(ax,[0 t_ca(end)]);
else
    warning('No F_roi in dFF.mat -> skipping Fig 2.');
end
