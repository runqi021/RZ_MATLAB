% phase_colorwheel_260530.m
% -----------------------------------------------------------------------
%  Draw a circular HSV colorwheel as a legend for preferred phase.
%  Annulus (r in [r_in, r_out]) is colored by angle; centre is white.
%  Saves transparent PNG + vector PDF next to this script.
%
%  Convention used by the polar/cartoon scripts: foot = 0, peak = pi.
%  Tick marks at 0, pi/2, pi, 3pi/2 with labels.
% -----------------------------------------------------------------------

clear; close all; clc;

%% ===================== USER-EDITABLE =====================
outDir   = 'C:\Users\Admin\Desktop';   % where to save
outName  = 'phase_colorwheel';
N        = 800;        % image resolution (px)
r_out    = 0.95;       % outer radius (normalized to half-width)
r_in     = 0.82;       % ring 1.5x thinner (~0.13 wide)
phase_rot = 0;         % rotate HSV so phase 0 lands at desired hue (rad)
% mathematical angle convention: theta = atan2(y, x) so theta=0 points right,
% theta=pi/2 points up.  Set theta_zero_dir to flip if you want 0 at the top.
theta_zero_dir = 'right';   % 'right' | 'top'
% =========================================================

%% ---- build the image (RGBA) ----
[X, Y] = meshgrid(linspace(-1, 1, N), linspace(-1, 1, N));
R = sqrt(X.^2 + Y.^2);

% angle: take atan2(y, x), then re-orient to chosen zero direction.
if strcmpi(theta_zero_dir,'right')
    th = atan2(Y, X);                       % 0 = right, increases CCW
else   % 'top'
    th = atan2(X, -Y);                      % 0 = top,   increases CW->CCW
end
th = mod(th - phase_rot, 2*pi);             % wrap & rotate

% HSV color from theta
H = th / (2*pi);
S = ones(N);
V = ones(N);
rgb = hsv2rgb(cat(3, H, S, V));

% alpha: 1 inside annulus, 0 outside outer + (white) inside inner
alpha = double(R <= r_out);

% White centre (mix rgb toward white where r < r_in, hard cut)
ringMask = R >= r_in & R <= r_out;
centerMask = R < r_in;
for c = 1:3
    ch = rgb(:,:,c);
    ch(centerMask) = 1;        % white centre
    rgb(:,:,c) = ch;
end
% smooth alpha at outer rim for a clean edge
alpha = max(0, min(1, (r_out - R)/0.02));
alpha(R > r_out) = 0;

%% ---- write transparent PNG + show preview ----
outBase = fullfile(outDir, outName);
imwrite(uint8(rgb*255), [outBase '.png'], 'Alpha', uint8(alpha*255));
fprintf('Saved %s.png\n', outBase);

% --- vector PDF: re-render via a figure + tick labels ---
fig = figure('Color','none','Units','pixels','Position',[100 100 600 600]);
ax  = axes(fig,'Position',[0.05 0.05 0.9 0.9]);
image(ax, [-1 1], [-1 1], rgb, 'AlphaData', alpha); hold(ax,'on');
axis(ax,'image','off');
set(ax,'YDir','normal');

% tick labels at the four cardinal phases (taking phase_rot into account)
tick_phases = [0 pi/2 pi 3*pi/2];
tick_labels = {'0 (foot)','\pi/2','\pi (peak)','3\pi/2'};
for k = 1:numel(tick_phases)
    phi = mod(tick_phases(k) + phase_rot, 2*pi);
    if strcmpi(theta_zero_dir,'right')
        xt = cos(phi)*(r_out + 0.06);
        yt = sin(phi)*(r_out + 0.06);
        xi = cos(phi)*r_out;    yi = sin(phi)*r_out;
        xo = cos(phi)*(r_out + 0.04);  yo = sin(phi)*(r_out + 0.04);
    else
        xt = sin(phi)*(r_out + 0.06);
        yt = cos(phi)*(r_out + 0.06);
        xi = sin(phi)*r_out;    yi = cos(phi)*r_out;
        xo = sin(phi)*(r_out + 0.04);  yo = cos(phi)*(r_out + 0.04);
    end
    plot(ax, [xi xo], [yi yo], 'k-', 'LineWidth', 1);
    text(ax, xt, yt, tick_labels{k}, 'HorizontalAlignment','center', ...
         'FontSize',12,'FontWeight','bold');
end

xlim(ax,[-1.25 1.25]); ylim(ax,[-1.25 1.25]);
exportgraphics(fig, [outBase '.pdf'], 'ContentType','vector', 'BackgroundColor','none');
fprintf('Saved %s.pdf\n', outBase);
