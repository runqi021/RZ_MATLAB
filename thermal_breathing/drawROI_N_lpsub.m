% drawROI_N_lpsub.m
% Draw L/R nostril ROIs AND baseline-remove (LP-subtraction) in one script.
%
% For each nostril (L=LEFT/L1, R=RIGHT/L2): draw an ellipse/circle on the avg
% projection, extract the raw deg-C ROI trace from the tracking-aligned stack,
% then baseline = lowpass(raw, LP_CUT), detrended = raw - baseline (inhale-up).
%
% Plots THREE traces: detrended LEFT, detrended RIGHT, and their AVERAGE.
% Saves the drawn ROI params -> <stem>_nostrilROI.mat
%   and the averaged breathing signal -> <stem>_breath.mat (canonical).
%
% RUN IN MATLAB YOURSELF (interactive: you draw the ROI). Zero-phase
% (butter+filtfilt). No flips: imagesc shows the native .ats orientation.
close all; clc; clear;

% Give ONLY the DLC csv; the _nostrilC.mat is resolved via thermal_resolve_paths.
dlcCsv    = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos\5916296_nose_n2DLC_Resnet50_260615_thermal_noseJun17shuffle2_snapshot_best-100.csv";
dataRoot  = "D:\260615_thermalNbasler";
ROI_SHAPE = 'circle';   % 'ellipse' | 'circle'
STAT      = 'mean';     % per-frame ROI readout: 'mean' | 'median' | 'max'
LP_CUT    = 1;          % Hz, low-pass cutoff = the slow baseline to subtract off
INVERT    = true;       % inhale cools the nostril -> invert so inhale reads as a rise
FINE_BP   = [1 20];     % finer band saved as breath_bp ([] = skip)

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here),'mod','bluewhitered'));
P = thermal_resolve_paths(dlcCsv, dataRoot);
S = load(P.nostrilC);
sides = {'L','R'};
fps = double(S.L_stack_fps);
[lb,la] = butter(2, LP_CUT/(fps/2), 'low');
R = struct();

% ---- draw ROI + extract + LP-subtract, per nostril ----
for i = 1:2
    sd  = sides{i};
    avg = double(S.(sd+"_avg"));
    stk = double(S.(sd+"_stack"));
    sidelbl = char(S.(sd+"_side"));
    win = size(avg,1); cc = (win+1)/2;

    figure('Color','w','Position',[120 120 620 560],'Name',sprintf('%s nostril ROI',sidelbl));
    ax = axes;
    imagesc(ax, avg); axis(ax,'image'); colormap(ax, bluewhitered(256)); colorbar(ax);
    hold(ax,'on'); plot(ax, cc, cc, 'c+','MarkerSize',10,'LineWidth',1.2);
    title(ax, sprintf('%s avg proj (\\circC) — DRAW ROI, double-click to confirm', sidelbl));
    if strcmpi(ROI_SHAPE,'circle'), h = drawcircle(ax); else, h = drawellipse(ax); end
    wait(h);
    mask = createMask(h);

    flat = reshape(stk, size(stk,1), []);
    roipix = flat(:, mask(:));
    switch STAT
        case 'mean',   tr = mean(roipix, 2, 'omitnan');
        case 'median', tr = median(roipix, 2, 'omitnan');
        case 'max',    tr = max(roipix, [], 2, 'omitnan');
        otherwise,     error('STAT must be mean | median | max');
    end

    base = filtfilt(lb, la, tr);     % zero-phase low-pass baseline
    detr = tr - base;                % baseline-removed
    if INVERT, detr = -detr; end

    R.(sd).side = sidelbl; R.(sd).mask = mask; R.(sd).trace = tr; R.(sd).detr = detr;
    R.(sd).fps = double(S.(sd+"_stack_fps")); R.(sd).npix = nnz(mask);
    R.(sd).shape = ROI_SHAPE; R.(sd).stat = STAT;
    if strcmpi(ROI_SHAPE,'circle')
        R.(sd).center = h.Center; R.(sd).radius = h.Radius;
    else
        R.(sd).center = h.Center; R.(sd).semiaxes = h.SemiAxes; R.(sd).angle = h.RotationAngle;
    end
    fprintf('%s: ROI %d px, raw mean %.2f C\n', sidelbl, nnz(mask), mean(tr,'omitnan'));
end
R.src = char(dlcCsv); R.nostrilC = P.nostrilC;
save(P.nostrilROI, '-struct', 'R');

% ---- the three traces: detrended L, R, and their average ----
detrA = (R.L.detr + R.R.detr) / 2;
t = (0:numel(detrA)-1)/fps;
figure('Color','w','Position',[200 120 1100 640]);
ax1 = subplot(3,1,1); plot(ax1, t, R.L.detr,'-'); grid(ax1,'on'); ylabel(ax1,'\circC');
title(ax1, sprintf('LEFT detrended (raw - LP %g Hz)', LP_CUT));
ax2 = subplot(3,1,2); plot(ax2, t, R.R.detr,'-'); grid(ax2,'on'); ylabel(ax2,'\circC');
title(ax2, 'RIGHT detrended');
ax3 = subplot(3,1,3); plot(ax3, t, detrA,'-'); grid(ax3,'on'); ylabel(ax3,'\circC'); xlabel(ax3,'s');
% breathing peak of the average (info only)
x = detrA - mean(detrA,'omitnan'); nf = 2^nextpow2(numel(x));
Pw = abs(fft(x,nf)).^2; fr = (0:nf-1)*(fps/nf); inb = fr>=LP_CUT & fr<=15;
[~,ip] = max(Pw(inb)); frb = fr(inb);
title(ax3, sprintf('L+R AVERAGE (breathing, peak %.2f Hz)', frb(ip)));
linkaxes([ax1 ax2 ax3],'x');
invstr = ''; if INVERT, invstr = ' (inhale up)'; end
sgtitle(sprintf('Nostril breathing — LP-subtraction (cutoff %g Hz)%s', LP_CUT, invstr));

% ---- write THE breathing signal (L+R averaged) ----
if ~isempty(FINE_BP)
    [fbb,faa] = butter(2, FINE_BP/(fps/2), 'bandpass');
    detrA_bp = filtfilt(fbb, faa, detrA);
else
    detrA_bp = [];
end
B = struct();
B.breath    = detrA(:);
B.breath_bp = detrA_bp(:);
B.fps       = fps;
B.t         = t(:);
B.method    = 'lpsub';
B.lp_cut    = LP_CUT;
B.fine_bp   = FINE_BP;
B.inverted  = INVERT;
B.roi_stat  = STAT;
B.roi_shape = ROI_SHAPE;
B.animal    = P.animal; B.run = P.k;
B.src_csv   = char(dlcCsv); B.src_ats = P.ats;
save(P.breath, '-struct', 'B');
fprintf('saved %s (ROI) and %s (breathing)\n', P.nostrilROI, P.breath);
