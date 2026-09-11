% IO_selROI_psd_overlay_260527.m
% -----------------------------------------------------------------------
%  PSD overlay of the SELECTED IO ROIs (the n_sel keep-list from
%  select_ROI_dFF_plot_260117_v2.m).
%
%  Each ROI's dF/F power spectrum is computed on the time-derivative of the
%  trace (diff(dF/F)*fps, units dF/F per s) in dB, log-x -- same convention
%  as breath_combined_summary_260329.m.  All selected ROIs overlaid, colored
%  by a colormap, legend = original ROI label.
%
%  Inputs (IO folder): *_ch1_dFF.mat (dFF [T x N]).
%  Dependencies: detect_session_fps.m, Chronux (mtspectrumc).
% -----------------------------------------------------------------------

clear; close all; clc;

%% ========================= PATH SETUP ================================
scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(genpath(fullfile(scriptDir, 'chronux_2_12')));

%% ===================== USER-EDITABLE PARAMETERS ======================
ioDir      = 'C:\Users\Admin\Desktop\ChAT_analysis\0124\IO\roi3_R_-1000_140_2x_34lp_512x256_00001';
n_sel      = [7, 8, 9, 12, 13, 14, 15, 17, 18, 19, 23, 22, 25, 26, 30, 29, 28, 31, 32, 34];

tossFrames = 30;          % drop first N calcium frames
fallback_fps = 30;
TW         = 6;           % multitaper time-bandwidth (K = 2*TW-1 tapers)
fmin       = 0.05;        % Hz
fmax       = 15;          % Hz
doSave     = true;
% =====================================================================

set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');

%% ---- load dF/F ----
df = dir(fullfile(ioDir, '*_ch1_dFF.mat'));
assert(~isempty(df), 'No *_ch1_dFF.mat in %s', ioDir);
D   = load(fullfile(df(1).folder, df(1).name), 'dFF');
dFF = double(D.dFF);
if tossFrames > 0 && size(dFF,1) > tossFrames, dFF(1:tossFrames,:) = []; end
[T, N] = size(dFF);
fps = detect_session_fps(ioDir, fallback_fps);

bad = n_sel(n_sel < 1 | n_sel > N);
assert(isempty(bad), 'n_sel has IDs outside 1..%d: %s', N, mat2str(bad));
fprintf('IO %s: fps=%.2f, %d frames, %d/%d ROIs selected.\n', ...
        df(1).name, fps, T, numel(n_sel), N);

%% ---- breath waveform (optional overlay) ----
bp = dir(fullfile(ioDir, '*DLC*breath_peak_data.mat'));
hasBreath = ~isempty(bp);
if hasBreath
    BP = load(fullfile(bp(1).folder, bp(1).name), 'breath');
    bw = detrend(double(BP.breath(:)));
    if tossFrames > 0 && numel(bw) > tossFrames, bw(1:tossFrames) = []; end
    bw = bw - mean(bw);
end

%% ---- per-ROI derivative PSD ----
p.Fs=fps; p.tapers=[TW,2*TW-1]; p.pad=0; p.fpass=[fmin,min(fmax,fps/2)]; p.err=0;
cols = turbo(numel(n_sel));

f1 = figure('Color','w','Name','IO selected-ROI PSD overlay', ...
            'Units','normalized','Position',[0.1 0.1 0.7 0.78]);
ax = axes(f1); hold(ax,'on');
h = gobjects(numel(n_sel),1);
for i = 1:numel(n_sel)
    roi = n_sel(i);
    [S,f] = mtspectrumc(diff(dFF(:,roi))*fps, p);
    h(i) = plot(ax, f(:), 10*log10(S(:)), '-', 'Color', cols(i,:), 'LineWidth', 1.0);
end

legH = h; legL = compose('ROI %d', n_sel(:));
if hasBreath                                   % breath waveform PSD, black
    [Sb,fb] = mtspectrumc(bw, p);
    hb = plot(ax, fb(:), 10*log10(Sb(:)), 'k-', 'LineWidth', 1.8);
    legH = [legH; hb]; legL = [legL; "breath"];
end

set(ax,'XScale','log'); xlim(ax,[fmin fmax]);
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'power (dB)');
title(ax, sprintf('IO selected-ROI dF/F'' PSD overlay + breath (black)  (N=%d ROIs, fps=%.2f)', ...
      numel(n_sel), fps), 'Interpreter','none');
legend(legH, legL, 'Location','eastoutside', 'NumColumns', 2, 'Box','off');
grid(ax,'on');

if doSave
    exportgraphics(f1, fullfile(ioDir,'IO_selROI_psd_overlay.png'), ...
        'Resolution',150, 'BackgroundColor','white');
    fprintf('Saved IO_selROI_psd_overlay.png to %s\n', ioDir);
end
