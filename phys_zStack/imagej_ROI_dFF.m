%% imagej_ROI_dFF.m
% Compute dF/F for a single ImageJ-ROI fluorescence trace and plot
% F (black, left axis) and dFF (blue, right axis) overlaid.
%
% Input: a 2-column text file exported from ImageJ:
%   col 1 = time (s), col 2 = mean F   (header line "[sec]  Mean" is skipped)

clear; clc;

%% ---------------- user params ----------------
txtFile        = 'roi_F.txt';   % <-- your pasted data saved as text
FPS            = 30;            % imaging rate (Hz); 0.0333 s step => 30 Hz
DropFirstSec   = 0;            % seconds to drop at start
BaselineWinSec = 20;           % sliding-median baseline window (s)
% ----------------------------------------------

addpath(fileparts(mfilename('fullpath')));   % so helper. is on path

%% Load
M = readmatrix(txtFile);          % [N x 2]: time, F
t = M(:,1);
F = M(:,2);

%% dF/F (project helper, sliding-median baseline)
out = helper.dFF_RZ(F, ...
    'FPS', FPS, ...
    'DropFirstSec', DropFirstSec, ...
    'BaselineWinSec', BaselineWinSec);

t_dff = out.t_dff;     % time after dropping frames
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
title('ImageJ ROI: F (black) and \DeltaF/F (blue)');
box off; xlim([t_dff(1) t_dff(end)]);
