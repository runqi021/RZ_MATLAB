% stitch_ch1_ch3_overlay_260528.m
% -----------------------------------------------------------------------
%  Two-channel anatomy overlay: dark blood-vessel skeleton (ch3) on a
%  light background + magenta-tinted neuron channel (ch1) on top.
%  Standalone, no dependencies beyond Image Processing Toolbox.
% -----------------------------------------------------------------------

clear; close all; clc;

%% ===================== USER-EDITABLE =================================
mouseRoot = 'C:\Users\Admin\Desktop\260522_ChAT_soma_G8s';
ch3Tif    = fullfile(mouseRoot,'matlab_stitch','stitched_ch3_avgproj.tif');   % vessels
ch1Tif    = fullfile(mouseRoot,'matlab_stitch','stitched_ch1_avgproj.tif');   % neurons
outDir    = fullfile(mouseRoot,'matlab_stitch');

% ---- channel preprocessing ----
clip_pct_ch3   = [0.5 99.5];
clip_pct_ch1   = [5 99.5];
gamma_ch3      = 0.6;
gamma_ch1      = 0.5;
clahe_tiles    = [32 32];
clahe_clip     = 0.01;
use_clahe_ch1  = false;      % CLAHE adds haze on punctate ch1; off by default
use_tophat_ch1 = true;       % subtract anything wider than a soma -> kills diffuse bg
tophat_radius_px = 25;       % px; somas ~17-22 px at 2x stitch resolution

% ---- composite look ----
invert_vessel  = true;          % dark vessels on light bg
neuron_tint    = [0.10 0.85 0.30];   % green (RGB)
neuron_alpha   = 0.85;           % opacity of neuron layer (0..1)
neuron_thresh  = 0.05;           % dim pixels below this don't show (lets bg through)
% =====================================================================

%% ---- load + preprocess vessels (ch3) ----
imV = double(imread(ch3Tif));
lo = prctile(imV(:), clip_pct_ch3(1)); hi = prctile(imV(:), clip_pct_ch3(2));
imVn = max(0, min(1, (imV - lo)/max(hi-lo,eps))).^gamma_ch3;
imVn = adapthisteq(imVn, 'NumTiles', clahe_tiles, 'ClipLimit', clahe_clip);
if invert_vessel, imVn = 1 - imVn; end

%% ---- load + preprocess neurons (ch1) ----
imN = double(imread(ch1Tif));
if use_tophat_ch1
    imN = imtophat(imN, strel('disk', tophat_radius_px));   % drop anything wider than a soma
end
lo = prctile(imN(:), clip_pct_ch1(1)); hi = prctile(imN(:), clip_pct_ch1(2));
imNn = max(0, min(1, (imN - lo)/max(hi-lo,eps))).^gamma_ch1;
if use_clahe_ch1
    imNn = adapthisteq(imNn, 'NumTiles', clahe_tiles, 'ClipLimit', clahe_clip);
end

%% ---- composite ----
% RGB background = inverted vessel grayscale
bg = repmat(imVn, [1 1 3]);

% per-pixel alpha from neuron intensity (with a soft floor to suppress background)
a = max(0, (imNn - neuron_thresh) / (1 - neuron_thresh));
a = a * neuron_alpha;

% alpha-blend magenta-tinted neuron layer onto vessel background
out = bg;
for ch = 1:3
    out(:,:,ch) = bg(:,:,ch) .* (1 - a) + neuron_tint(ch) * a;
end

%% ---- show + save ----
fig = figure('Color','w','Name','ch1 + ch3 overlay','Units','normalized','Position',[0.04 0.04 0.92 0.92]);
image(out); axis image off;
title('blood vessels (gray) + neurons (magenta)','FontSize',12);

exportgraphics(fig, fullfile(outDir,'stitch_ch1ch3_overlay.png'), 'Resolution',200,'BackgroundColor','white');
imwrite(uint8(round(out*255)), fullfile(outDir,'stitch_ch1ch3_overlay_native.png'));
fprintf('Saved:\n  %s\n  %s\n', ...
        fullfile(outDir,'stitch_ch1ch3_overlay.png'), ...
        fullfile(outDir,'stitch_ch1ch3_overlay_native.png'));
