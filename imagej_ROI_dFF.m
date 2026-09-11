%% imagej_ROI_dFF.m
% Compute dF/F for a single ImageJ-ROI fluorescence trace and plot
% F (black, left axis) and dFF (blue, right axis) overlaid.
%
% Input: ImageJ "Values.csv" -> col 2 = Mean F (col 1 ignored, treated as frames).

clear; clc;

%% ---------------- user params ----------------
csvFile        = "C:\Users\Admin\Desktop\ChAT_analysis\0523\cell1\roi1_4x_x-900y700z-15_6000f_12lp_00001\ROI1_largeROI_Values.csv";
FPS            = 30;            % imaging rate (Hz)
DropFirstSec   = 0;            % seconds to drop at start
BaselineWinSec = 20;           % sliding-median baseline window (s)
% ----------------------------------------------

addpath(fileparts(mfilename('fullpath')));   % so helper. is on path

%% Load: just take F (col 2), time comes from FPS
M = readmatrix(csvFile);
F = M(:,2);

%% dF/F (project helper, sliding-median baseline)
out = helper.dFF_RZ(F, ...
    'FPS', FPS, ...
    'DropFirstSec', DropFirstSec, ...
    'BaselineWinSec', BaselineWinSec);

t_dff = out.t_dff;     % time (s) after dropping frames
dFF   = out.dFF;       % dF/F
F_dff = out.F_dff;     % matching raw F

%% Plot: F (black, left) + dFF (blue, right)
figure('Color','w','Position',[100 100 1100 420]);

yyaxis left
plot(t_dff, F_dff, 'k', 'LineWidth', 0.8);
ylabel('F (a.u.)');
ax = gca; ax.YColor = 'k';

yyaxis right
plot(t_dff, dFF, 'b', 'LineWidth', 0.8);
ylabel('\DeltaF/F');
ax.YColor = 'b';

xlabel('Time (s)');
title('F (black) and \DeltaF/F (blue)');
box off; xlim([t_dff(1) t_dff(end)]);
