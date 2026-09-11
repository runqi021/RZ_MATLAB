% thermal_nostril_linescan.m
% Check pixel CONTINUITY around the tracked nostril center: instead of averaging
% an ROI, draw a LINE across the avg projection and read the temperature
% cross-section along it -- both on the average image (spatial profile) and over
% time (kymograph = line profile vs frame).
%
% Per nostril (L=LEFT/L1, R=RIGHT/L2):
%   FIG: (1) avg proj (deg C) with your line + tracked center
%        (2) cross-section of the AVG along the line (deg C vs position)
%        (3) kymograph: position-along-line x time, from the aligned stack
%
% RUN IN MATLAB YOURSELF (interactive: you draw the line). No flips.
close all; clc; clear;

dlcCsv   = "D:\260615_thermalNbasler\nose_thermal\260615_thermal_nose-RZ-2026-06-17\videos\5916296_nose_n4DLC_Resnet50_260615_thermal_noseJun17shuffle2_snapshot_best-100.csv";
dataRoot = "D:\260615_thermalNbasler";

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here),'mod','bluewhitered'));
P = thermal_resolve_paths(dlcCsv, dataRoot);
S = load(P.nostrilC);
sides = {'L','R'};

for i = 1:2
    sd  = sides{i};
    avg = double(S.(sd+"_avg"));
    stk = double(S.(sd+"_stack"));
    sfps = double(S.(sd+"_stack_fps"));
    sidelbl = char(S.(sd+"_side"));
    win = size(avg,1); cc = (win+1)/2;

    % --- draw the line on the avg projection ---
    figure('Color','w','Position',[120 120 600 560],'Name',sprintf('%s line',sidelbl));
    axd = axes;
    imagesc(axd, avg); axis(axd,'image'); colormap(axd, bluewhitered(256)); colorbar(axd);
    hold(axd,'on'); plot(axd, cc, cc, 'g+','MarkerSize',12,'LineWidth',1.5);
    title(axd, sprintf('%s avg proj — DRAW A LINE through the center, double-click', sidelbl));
    h = drawline(axd);
    wait(h);
    pos = h.Position;                       % [x1 y1; x2 y2] in win-px

    % --- sample points along the line (~1 per pixel) ---
    len = hypot(diff(pos(:,1)), diff(pos(:,2)));
    np  = max(2, round(len) + 1);
    xs  = linspace(pos(1,1), pos(2,1), np);
    ys  = linspace(pos(1,2), pos(2,2), np);
    d   = linspace(0, len, np);             % distance along line (native px)

    csAvg = interp2(avg, xs, ys, 'linear'); % cross-section of the average image

    % --- kymograph: sample the line in every frame of the stack ---
    T = size(stk,1);
    kymo = zeros(np, T);
    for k = 1:T
        kymo(:,k) = interp2(squeeze(stk(k,:,:)), xs, ys, 'linear');
    end
    tk = (0:T-1)/sfps;

    % --- figure ---
    figure('Color','w','Position',[180 90 1150 760],'Name',sprintf('%s linescan',sidelbl));
    subplot(2,2,1);
    imagesc(avg); axis image; colormap(gca, bluewhitered(256)); colorbar; hold on;
    plot(pos(:,1), pos(:,2), 'k-', 'LineWidth', 1.5);
    plot(pos(1,1), pos(1,2), 'ko', 'MarkerFaceColor','w');   % line start
    plot(cc, cc, 'g+', 'MarkerSize',12,'LineWidth',1.5);     % tracked center
    title(sprintf('%s avg proj + line', sidelbl));

    subplot(2,2,2);
    plot(d, csAvg, '-o', 'MarkerSize',3); grid on;
    xlabel('position along line (px)'); ylabel('\circC');
    title(sprintf('%s cross-section of AVG (start = circle)', sidelbl));

    subplot(2,2,[3 4]);
    imagesc(tk, d, kymo); axis xy; colormap(gca,'hot'); colorbar;
    xlabel('time (s)'); ylabel('position along line (px)');
    title(sprintf('%s kymograph (line profile vs time, \\circC)', sidelbl));

    sgtitle(sprintf('%s nostril line-scan — pixel continuity around center', sidelbl));
end
