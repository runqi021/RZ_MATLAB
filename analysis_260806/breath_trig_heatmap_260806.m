% breath_trig_heatmap_260806.m
% -----------------------------------------------------------------------
%  Inspiration-PEAK-triggered dF/F heatmap for ONE recording folder.
%
%  Deliberately minimal:
%    - NO sorting.  ROIs are plotted in their native cpSAM order (1..N),
%      top to bottom.  Nothing is reordered, ranked, or selected.
%    - NO significance test, no event detection, no coherence.
%    - Every ROI is shown, including silent ones.
%
%  Two panels of the SAME data:
%      left   trigger-averaged dF/F                (dF/F units)
%      right  trigger-averaged z-scored dF/F       (SD units)
%             z is taken per ROI over the WHOLE trace before averaging, so
%             bright and dim ROIs become comparable.  This is the panel to
%             read for shape; the left one keeps the real amplitude.
%  Both use clim = [0.5 99.5] percentile of the matrix being displayed.
%
%  A trigger-averaged breath trace is drawn above so the calcium columns can
%  be read against the breath cycle.
%
%  ALIGNMENT.  The breath camera is 2P-frame-triggered, so breath video frame
%  k IS imaging frame k.  dFF has TossFrames removed from its front, so the
%  breath frame index is shifted by nDrop (=TossFrames) and both are then
%  truncated to the common length.
%
%  Input  (all already in the recording folder):
%      *_ch1_dFF.mat              dFF [T x N]
%      breath_peak_pc1.mat        insp_onset_idx  <- this is the PEAK
%      breath_insp_start_pc1.mat  insp_start_idx  <- this is the ONSET (foot)
%  Output (written into the same folder):
%      breath_trig_heatmap_<TRIGGER>.png / .pdf
%      breath_trig_heatmap_<TRIGGER>.mat
%
%  Runqi Zhang / 2026-08-06

clear; clc;   % NOT close all -- each run opens its own figure so FOVs can be
              % compared side by side.  Every panel is drawn through an explicit
              % axes handle, so an older figure being current cannot steal a plot.
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fileparts(scriptDir);
addpath(repoRoot);                       % detect_session_fps
addpath(genpath(fullfile(repoRoot,'chronux_2_12')));   % mtspectrumc

%% ===================== USER-EDITABLE =====================
folderPath = 'C:\260806_sst-soma-g8s\phys\roi1_6x_z5_14lp_3000f_00001';
% folderPath = 'C:\260806_sst-soma-g8s\phys\roi1_2x_z5_15lp_6000f_00001';
% folderPath = 'C:\260806_sst-soma-g8s\phys\roi1_5x_z10_16lp_3000f_00001';

TRIGGER   = 'peak';        % 'peak'  = inspiratory peak  (breath_peak_pc1.insp_onset_idx)
                           % 'onset' = inspiration foot  (breath_insp_start_pc1.insp_start_idx)

% --- window, set in BREATH CYCLES, not seconds ---
% The window is derived from the measured IBI (median trigger-to-trigger interval
% of THIS recording), so the panels always hold the same number of cycles whatever
% the animal's breath rate is.  The top breath trace is as wide as the three panels
% under it, so it gets 3x their span.
nCyc_panel = 10/3;         % cycles per panel in the row below (~3.3)
nCyc_top   = 3*nCyc_panel; % cycles in the top breath trace (= 10)
nDrop     = 30;            % frames tossed off the front of the calcium (TossFrames)
prcLim    = [0.5 99.5];    % percentile clamp for both color scales
fallback_fps = 30;
doSave    = true;

% --- power spectra (bottom row) ---
fpass     = [0.1 14];      % Hz.  Nyquist is fps/2 = 15 Hz, so the top of this
                           %   band is only ~1 Hz below it -- see note below.
W_hz      = 0.25;          % multitaper HALF-BANDWIDTH in Hz -- this is the smoothing
                           %   knob.  TW and the taper count are derived from it and
                           %   the recording length, so the spectra are smoothed by
                           %   the same amount in Hz no matter how long the file is.
                           %   Bigger W = smoother but blurrier; 0.25 Hz is still far
                           %   narrower than the ~1.75 Hz spacing of the breath
                           %   harmonics, so nothing real is merged.
% NOTHING is normalised.  Breath PC1 is detrended and z-scored upstream and the
% calcium traces are z-scored here, so both spectra are already in comparable
% units (z^2 per Hz) and dividing by a peak would only throw that away.

% --- dF/F trace stack (left panel of row 2) ---
% NOT auto-scaled per trace: every ROI is drawn on ONE common dF/F scale, so a
% big response looks big and a flat ROI looks flat.  Only the spacing between
% traces is chosen for you.
stackGain    = 1.0;        % trace spacing = stackGain x median peak-to-peak dF/F.
                           %   Lower = more overlap.
scaleBarDFF  = 0.05;       % 5% dF/F calibration bar
% =========================================================

set(0,'DefaultAxesFontName','Arial'); set(0,'DefaultTextFontName','Arial');

%% ---- load ----
df = dir(fullfile(folderPath,'*_ch1_dFF.mat'));
assert(~isempty(df), 'No *_ch1_dFF.mat in %s', folderPath);
D   = load(fullfile(df(1).folder, df(1).name), 'dFF', 'F_roi');
dFF = double(D.dFF);                                   % [T x N]
Fr  = double(D.F_roi);                                 % [T x N] RAW ROI fluorescence

switch lower(TRIGGER)
    case 'peak'
        B = load(fullfile(folderPath,'breath_peak_pc1.mat'));
        trig_raw = round(B.insp_onset_idx(:));         % field name says onset, it is the PEAK
        trigLabel = 'inspiratory peak';
    case 'onset'
        B = load(fullfile(folderPath,'breath_insp_start_pc1.mat'));
        trig_raw = round(B.insp_start_idx(:));
        trigLabel = 'inspiration onset';
    otherwise
        error('TRIGGER must be ''peak'' or ''onset''.');
end
breath_raw = B.breath(:);

fps = detect_session_fps(folderPath, fallback_fps);

%% ---- align breath frames to dFF frames ----
% breath frame k -> dFF row k-nDrop, then truncate both to the common length
breath = breath_raw(nDrop+1:end);
T      = min([size(dFF,1), size(Fr,1), numel(breath)]);
dFF    = dFF(1:T,:);
Fr     = Fr(1:T,:);
breath = breath(1:T);
trig   = trig_raw - nDrop;

nROI = size(dFF,2);

% window from the measured inter-breath interval of THIS recording
ibi_s = median(diff(sort(trig)))/fps;
win   = round(nCyc_panel*ibi_s*fps/2);   % half-window, panels
win_b = round(nCyc_top  *ibi_s*fps/2);   % half-window, top breath trace
tt    = (-win:win)/fps;
tb    = (-win_b:win_b)/fps;

keep   = trig-win   >= 1 & trig+win   <= T;
keep_b = trig-win_b >= 1 & trig+win_b <= T;
trig_b = trig(keep_b);
trig   = trig(keep);
fprintf('%s\n  fps %.3f | %d frames (%.1f s) | %d ROIs\n', folderPath, fps, T, T/fps, nROI);
fprintf('  IBI %.3f s (%.2f Hz) -> panels +/-%.2f s (%.2f cycles), breath trace +/-%.2f s (%.1f cycles)\n', ...
        ibi_s, 1/ibi_s, win/fps, nCyc_panel, win_b/fps, nCyc_top);
fprintf('  %s triggers: %d total, %d fit the panel window, %d fit the breath window\n', ...
        trigLabel, numel(keep), numel(trig), numel(trig_b));
assert(~isempty(trig) && ~isempty(trig_b), 'No trigger has a full window inside the recording.');

%% ---- trigger-average ----
idx   = trig(:).'   + (-win:win).';        % [nWin  x nTrig] frame indices
idx_b = trig_b(:).' + (-win_b:win_b).';    % wider, breath only

Z = (dFF - mean(dFF,1)) ./ std(dFF,0,1);  % per-ROI z of the WHOLE trace

M  = zeros(nROI, 2*win+1);               % mean dF/F
Mz = zeros(nROI, 2*win+1);               % mean z-scored dF/F
Mf = zeros(nROI, 2*win+1);               % mean RAW fluorescence (counts)
for r = 1:nROI
    tr = dFF(:,r); zr = Z(:,r); fr = Fr(:,r);
    M(r,:)  = mean(tr(idx), 2).';
    Mz(r,:) = mean(zr(idx), 2).';
    Mf(r,:) = mean(fr(idx), 2).';
end
bAvg = mean(breath(idx_b), 2).';         % trigger-averaged breath, 10 cycles wide

%% ---- power spectra (whole aligned trace, NOT the triggered window) ----
TW  = max(2, round(W_hz * T/fps));        % half-bandwidth in Hz -> time-bandwidth
Kt  = 2*TW - 1;                           % all the usable Slepians; no error bars
params = struct('Fs',fps, 'tapers',[TW Kt], 'fpass',fpass, 'pad',0, ...
                'trialave',0, 'err',0);
[Sb, fb_hz] = mtspectrumc(detrend(breath), params);        % breath              [F x 1]
[Sc, fc_hz] = mtspectrumc(detrend(Z),      params);        % every ROI, Z-SCORED [F x N]
[~,ib] = max(Sb); f_breath = fb_hz(ib);
fprintf('  breath fundamental %.2f Hz | spectra: W=%.2f Hz -> TW=%d, %d tapers, %.1f-%.1f Hz\n', ...
        f_breath, W_hz, TW, Kt, fpass(1), fpass(2));

%% ---- figure ----
roiCol = hsv(nROI);                      % ONE colour per ROI, reused everywhere

% TWO rows in total:
%   row 1  breath trace | breath power spectrum
%   row 2  raw F stack  | dF/F heatmap | z heatmap | z power spectrum
% Both power panels sit in the right-hand column.  The 5x4 grid only exists to
% keep row 1 short: 1 grid row against 4, so row 2 gets the space.
fh = figure('Color','w','Position',[30 30 1780 990]);
tl = tiledlayout(fh, 5, 4, 'TileSpacing','compact','Padding','compact');

% ---- row 1 ----
axb = nexttile(tl, 1, [1 3]);
plot(axb, tb, bAvg, 'k', 'LineWidth', 1.4); hold(axb,'on');
xline(axb, 0, 'r-', 'LineWidth', 1);
% the panel edges below fall every nCyc_panel cycles -- mark them
for e = [-1 1]*win/fps, xline(axb, e, ':', 'Color',[0.5 0.5 0.5]); end
xlim(axb, [tb(1) tb(end)]); box(axb,'off');
set(axb, 'TickDir','out', 'XTick', unique(round((-4:4)*ibi_s, 2)));
xlabel(axb, sprintf('time from %s (s)', trigLabel));
ylabel(axb, 'breath (PC1, z)');
title(axb, sprintf('%s-triggered mean breath, %.0f cycles  (n = %d triggers; dotted = span of each panel below)', ...
                   trigLabel, nCyc_top, numel(trig_b)));

% breath power spectrum, right-hand column of row 1
axs1 = nexttile(tl, 4, [1 1]);
loglog(axs1, fb_hz, Sb, 'k', 'LineWidth', 1.6);
finish_spec(axs1, fpass, f_breath, 'power (z^2/Hz)', 'breath (PC1) power');

% ---- row 2 (4x the height of row 1) ----
draw_stack(nexttile(tl, 5, [4 1]), tt, M, roiCol, stackGain, scaleBarDFF, trigLabel);
draw_heatmap(nexttile(tl, 6, [4 1]), tt, M,  prcLim, 'mean dF/F', 'dF/F', trigLabel, gray(256));
draw_heatmap(nexttile(tl, 7, [4 1]), tt, Mz, prcLim, 'mean z-scored dF/F', 'z (SD)', trigLabel, parula(256));

% z-scored calcium power spectrum, right-hand column of row 2
axs2 = nexttile(tl, 8, [4 1]);
hold(axs2,'on');
for r = 1:nROI
    loglog(axs2, fc_hz, Sc(:,r), '-', 'Color',[roiCol(r,:) 0.7], 'LineWidth',0.9);
end
finish_spec(axs2, fpass, f_breath, 'power (SD^2/Hz)', ...
            sprintf('z-scored dF/F power, %d ROIs (not normalised)', nROI));

[~, recName] = fileparts(folderPath);
title(tl, sprintf('%s   |   %d ROIs, no sorting', recName, nROI), ...
      'Interpreter','none','FontWeight','bold');

%% ---- save ----
if doSave
    base = fullfile(folderPath, sprintf('breath_trig_heatmap_%s', lower(TRIGGER)));
    exportgraphics(fh, [base '.png'], 'Resolution',200, 'BackgroundColor','white');
    exportgraphics(fh, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    save([base '.mat'], 'M','Mz','Mf','bAvg','tt','tb','trig','trig_b','fps','nROI','T', ...
                        'Sb','Sc','fb_hz','fc_hz','f_breath','ibi_s', ...
                        'TRIGGER','nCyc_panel','nCyc_top','nDrop','prcLim', ...
                        'fpass','W_hz','TW','Kt','folderPath');
    fprintf('  saved %s.{png,pdf,mat}\n', base);
end

%% ---- local ----
function draw_heatmap(ax, tt, A, prcLim, ttl, cbl, trigLabel, cmap)
% cmap: gray(256)   -> black at the low clim, white at the high one (dF/F panel)
%       parula(256) -> blue at the low clim, yellow at the high one (z panel)
c = prctile(A(:), prcLim);
if ~(c(2) > c(1)), c = [min(A(:)) max(A(:))+eps]; end
imagesc(ax, tt, 1:size(A,1), A);
set(ax, 'CLim', c, 'YDir','normal', 'TickDir','out', 'Layer','top');
colormap(ax, cmap);
hold(ax,'on'); xline(ax, 0, 'r-', 'LineWidth', 1);
xlabel(ax, sprintf('time from %s (s)', trigLabel));
ylabel(ax, 'ROI # (native order)');
title(ax, ttl);
cb = colorbar(ax); cb.Label.String = cbl;
end

function draw_stack(ax, tt, A, roiCol, gain, barDFF, trigLabel)
% Trigger-averaged dF/F, one trace per ROI, stacked bottom-to-top in native ROI
% order.  NOTHING is auto-scaled: all traces share one dF/F scale, so amplitudes
% are directly comparable between ROIs and the scale bar applies to all of them.
% Only the vertical SPACING is chosen, from the median peak-to-peak of the panel.
nROI = size(A,1);
step = gain * median(range(A, 2));
if ~(step > 0), step = 1; end
hold(ax,'on');
for r = 1:nROI
    plot(ax, tt, A(r,:) - mean(A(r,:)) + (r-1)*step, '-', ...
         'Color', roiCol(r,:), 'LineWidth', 0.9);
end
xline(ax, 0, 'k-', 'LineWidth', 1);
xlim(ax, [tt(1) tt(end)]);
ylim(ax, [-step, (nROI-1)*step + step]);
set(ax, 'YTick', (0:4:nROI-1)*step, 'YTickLabel', compose('%d', 1:4:nROI), ...
        'TickDir','out', 'Box','off');

% calibration bar, drawn in real dF/F units at the right edge
xb = tt(end) - 0.02*(tt(end)-tt(1));
yb = -0.6*step;
plot(ax, [xb xb], yb + [0 barDFF], 'k-', 'LineWidth', 2.5);
text(ax, xb - 0.02*(tt(end)-tt(1)), yb + barDFF/2, sprintf('%g%% dF/F', 100*barDFF), ...
     'HorizontalAlignment','right', 'VerticalAlignment','middle', 'FontSize', 9);

xlabel(ax, sprintf('time from %s (s)', trigLabel));
ylabel(ax, 'ROI # (native order)');
title(ax, 'mean dF/F (common scale)');
end

function finish_spec(ax, fpass, f_breath, ylb, ttl)
hold(ax,'on');
xline(ax, f_breath, 'r--', 'LineWidth', 1);
set(ax, 'XScale','log', 'YScale','log', 'TickDir','out', 'Box','off');
xlim(ax, fpass);
xt = [0.1 0.2 0.5 1 2 3 5 7 10 14]; xt = xt(xt>=fpass(1) & xt<=fpass(2));
set(ax, 'XTick', xt, 'XTickLabel', compose('%g', xt));
xlabel(ax, 'frequency (Hz)');
ylabel(ax, ylb);
title(ax, ttl);
end
