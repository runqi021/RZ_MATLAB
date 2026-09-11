% thermal_breath_raw_batch.m
% Build the UN-DETRENDED (raw deg-C) nostril breathing signal for every session
% that already has a manual ROI, WITHOUT re-reading the .ats: it reuses the raw
% ROI traces drawn by drawROI_N_lpsub (stored in <stem>_nostrilROI.mat as
% L.trace / R.trace, full-rate deg C). No LP-subtraction, no invert -- the literal
% temperature the ROI sees.
%
% Per session saves <stem>_breath_raw.mat (L_raw, R_raw, breath_raw=(L+R)/2, fps, t)
% and a QC figure <stem>_breath_raw.png. Run in MATLAB.
close all; clc; clear;

DATA_ROOT = "D:\260615_thermalNbasler";

rois = dir(fullfile(char(DATA_ROOT), '*', 'cam1_*', 'Rec-*_nostrilROI.mat'));
assert(~isempty(rois), 'no _nostrilROI.mat under %s', DATA_ROOT);
fprintf('RAW breath signal: %d sessions with a manual ROI\n', numel(rois));

for i = 1:numel(rois)
    roiPath = fullfile(rois(i).folder, rois(i).name);
    stem    = erase(roiPath, '_nostrilROI.mat');
    R = load(roiPath);                       % top-level fields L, R (structs)
    if ~isfield(R,'L') || ~isfield(R,'R') || ~isfield(R.L,'trace')
        fprintf('  [%2d/%2d] SKIP %s (no raw trace)\n', i, numel(rois), rois(i).name); continue;
    end

    Lraw = double(R.L.trace(:));             % raw deg C, full rate (NOT detrended)
    Rraw = double(R.R.trace(:));
    fps  = double(R.L.fps);
    n    = min(numel(Lraw), numel(Rraw));
    Lraw = Lraw(1:n); Rraw = Rraw(1:n);
    breath_raw = (Lraw + Rraw) / 2;
    t = (0:n-1)' / fps;

    B = struct('L_raw',Lraw,'R_raw',Rraw,'breath_raw',breath_raw, ...
               'fps',fps,'t',t,'method','raw_undetrended','inverted',false, ...
               'roi_shape',R.L.shape,'roi_stat',R.L.stat,'src_roi',roiPath);
    save([char(stem) '_breath_raw.mat'], '-struct', 'B');

    fh = figure('Color','w','Position',[120 120 1100 620],'Visible','off');
    a1 = subplot(3,1,1); plot(a1,t,Lraw,'-'); grid(a1,'on'); ylabel(a1,'\circC');
    title(a1, sprintf('LEFT raw (un-detrended)  [%s]', rois(i).name),'Interpreter','none');
    a2 = subplot(3,1,2); plot(a2,t,Rraw,'-'); grid(a2,'on'); ylabel(a2,'\circC'); title(a2,'RIGHT raw');
    a3 = subplot(3,1,3); plot(a3,t,breath_raw,'-'); grid(a3,'on');
    ylabel(a3,'\circC'); xlabel(a3,'s'); title(a3,'L+R AVERAGE raw (un-detrended)');
    linkaxes([a1 a2 a3],'x');
    sgtitle('Nostril breathing — RAW un-detrended deg C (no LP-subtraction)');
    exportgraphics(fh, [char(stem) '_breath_raw.png'], 'Resolution', 120);
    close(fh);
    fprintf('  [%2d/%2d] %s  fps=%.1f  n=%d  mean L %.2f / R %.2f C\n', ...
            i, numel(rois), rois(i).name, fps, n, mean(Lraw), mean(Rraw));
end
fprintf('done.\n');
