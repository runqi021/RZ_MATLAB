% cartoon_label_coords_260530.m
% -----------------------------------------------------------------------
%  Define a DATA coordinate frame on a cartoon image.
%  THREE clicks, TRANSLATION + PER-AXIS SCALE only (no rotation, no shear):
%     Click 1: ORIGIN -> data (0, 0)
%     Click 2: a point you want to label x = xRef
%              (its pixel column relative to origin sets x-scale; pixel row ignored)
%     Click 3: a point you want to label y = yRef
%              (its pixel row    relative to origin sets y-scale; pixel column ignored)
%  Mapping:
%     pix_x = origin_pix_x + data_x * pix_per_x
%     pix_y = origin_pix_y + data_y * pix_per_y
%
%  Output : <stem>_coordframe.mat with:
%     origin_pix             1x2 clicked origin in pixels
%     pix_per_x, pix_per_y   scale (px per data unit)
%     anchors_pix            3x2 clicked pixel coords [origin; x_click; y_click]
%     xRef, yRef             data values used
%     A                      2x3 affine for compat: pix = A * [dx dy 1]'
%     imgPath                source cartoon path
% -----------------------------------------------------------------------

clear; close all; clc;

%% ===================== USER-EDITABLE =====================
imgPath = "C:\Users\Admin\Desktop\brainstem_map_cartoon.png";   % outline image
xRef    = 1500;          % data x-coord of click 2
yRef    = 1500;          % data y-coord of click 3
showGrid     = true;     % overlay verification grid
gridStepData = 200;      % grid step in data units
% =========================================================

imgPath = char(imgPath);
assert(isfile(imgPath), 'Image not found: %s', imgPath);
[stem, base, ~] = fileparts(imgPath);
outMat = fullfile(stem, [base '_coordframe.mat']);

I = imread(imgPath);
fig = figure('Color','w','Name','Click anchor points','Units','normalized', ...
             'Position',[0.05 0.05 0.85 0.85]);
ax = axes(fig);
imshow(I,'Parent',ax); hold(ax,'on');

labels = {'ORIGIN (0,0)', ...
          sprintf('X scale -> x = %d (pixel column sets scale)', xRef), ...
          sprintf('Y scale -> y = %d (pixel row sets scale)',    yRef)};
colors = [0 0 0; 0.85 0.10 0.10; 0.10 0.30 0.85];
anchors_pix = zeros(3,2);
for k = 1:3
    title(ax, sprintf('Click: %s', labels{k}), 'Interpreter','none');
    [xp, yp] = ginput(1);
    anchors_pix(k,:) = [xp, yp];
    plot(ax, xp, yp, 'o', 'MarkerFaceColor', colors(k,:), ...
        'MarkerEdgeColor','k', 'MarkerSize',9, 'LineWidth',1);
    text(ax, xp+18, yp, labels{k}, 'Color',colors(k,:), ...
        'FontSize',10, 'FontWeight','bold', 'Interpreter','none');
end

origin_pix = anchors_pix(1,:);
pix_per_x  = (anchors_pix(2,1) - origin_pix(1)) / xRef;
pix_per_y  = (anchors_pix(3,2) - origin_pix(2)) / yRef;
fprintf('Origin px = (%.1f, %.1f)\n', origin_pix);
fprintf('Scale: %.3f px / data-x   |  %.3f px / data-y\n', pix_per_x, pix_per_y);

% Affine for compatibility: pix = A * [dx; dy; 1]
A = [pix_per_x  0          origin_pix(1);
     0          pix_per_y  origin_pix(2)];

% --- verification grid ---
if showGrid
    dx_vals = 0 : gridStepData : xRef;
    dy_vals = 0 : gridStepData : yRef;
    for dx = dx_vals
        p1 = A * [dx; 0;    1]; p2 = A * [dx; yRef; 1];
        plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], ':', 'Color', [0.4 0.4 0.4]);
        text(ax, p1(1), p1(2)-12, sprintf('x=%d', dx), 'Color',[0.4 0.4 0.4], 'FontSize',8);
    end
    for dy = dy_vals
        p1 = A * [0;    dy; 1]; p2 = A * [xRef; dy; 1];
        plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], ':', 'Color', [0.4 0.4 0.4]);
        text(ax, p1(1)-60, p1(2), sprintf('y=%d', dy), 'Color',[0.4 0.4 0.4], 'FontSize',8);
    end
end
title(ax, sprintf('done  |  saved -> %s', outMat), 'Interpreter','none');

save(outMat, 'anchors_pix','origin_pix','pix_per_x','pix_per_y','A', ...
     'xRef','yRef','imgPath');
fprintf('Saved %s\n', outMat);

%% --- usage helper ---
% Pixel coords from a DATA point (dx, dy):
%   pix_x = origin_pix(1) + dx * pix_per_x;
%   pix_y = origin_pix(2) + dy * pix_per_y;
% Or via A:
%   pix = A * [dx; dy; 1];
function pix = data2pix(xy, A)  %#ok<DEFNU>
    n = size(xy,1);
    pix = (A * [xy ones(n,1)]')';   % n x 2
end
