function save_avgproj_png_260812(proj, P, outpng)
%SAVE_AVGPROJ_PNG_260812  Full-resolution PNG of the ROI avg-proj crop.
%
%   save_avgproj_png_260812(proj, P, outpng)
%
% proj   : second output of temporal_phase_cell_fig_260812 (crop_img, bnd_crop,
%          barLen_pr). Do nothing unless proj.have_proj is true.
% P      : needs outlineLW only.
% outpng : destination .png
%
% The crop embedded in the combined figure is rasterized at that figure's
% resolution; this writes the same crop at 600 dpi with the axes sized to the
% crop's own pixel dimensions, so the outline and scalebar stay crisp.
%
% Shared by temporal_phase_cell_gui_260812.m and per_cell_summary_260812.m.
[H,W] = size(proj.crop_img);
fg = figure('Color','w','Visible','off','Units','pixels','Position',[100 100 max(W,60) max(H,60)]);
axp = axes('Parent',fg,'Position',[0 0 1 1]);
imagesc(axp, proj.crop_img); colormap(axp, gray(256)); caxis(axp,[0 1]);
axis(axp,'image','off'); hold(axp,'on');
% Match the outline thickness rule used inside the combined figure. LineWidth is
% in POINTS, so the same value renders a different fraction of the image here --
% this axes is sized to the crop's own pixel dimensions -- than it does in the
% figure panel: measured 8.0 px here vs 3.0 px there for one 195x195 crop.
% P.outlineImgPx (image pixels) removes that; empty keeps the old behaviour.
wantPx = [];
if isfield(P,'outlineFrac') && ~isempty(P.outlineFrac) && P.outlineFrac > 0
    wantPx = P.outlineFrac * W;
elseif isfield(P,'outlineImgPx') && ~isempty(P.outlineImgPx)
    wantPx = P.outlineImgPx;
end
lwOut = P.outlineLW;
if ~isempty(wantPx) && W > 0
    drawnow;
    old = get(axp,'Units'); set(axp,'Units','points');
    pp  = get(axp,'Position'); set(axp,'Units',old);
    if pp(3) > 0, lwOut = max(0.15, wantPx * pp(3) / W); end
end
% Match the figure panel: a 4-element Color is [r g b alpha] on a line.
outCol = [1 1 0];
if isfield(P,'outlineAlpha') && ~isempty(P.outlineAlpha) && P.outlineAlpha < 1
    outCol = [1 1 0 max(0, min(1, P.outlineAlpha))];
end
for k = 1:numel(proj.bnd_crop)
    plot(axp, proj.bnd_crop{k}(:,2), proj.bnd_crop{k}(:,1), '-','Color',outCol,'LineWidth',lwOut);
end
sbCorner = 'lower left';
if isfield(P,'scalebarCorner') && ~isempty(P.scalebarCorner), sbCorner = P.scalebarCorner; end
sbThick = 0.02;
if isfield(P,'scalebarThickFrac') && ~isempty(P.scalebarThickFrac)
    sbThick = P.scalebarThickFrac;
end
rectangle(axp,'Position',scalebar_pos_260815(W, H, proj.barLen_pr, sbCorner, sbThick), ...
          'FaceColor','w','EdgeColor','none');
hold(axp,'off');
exportgraphics(axp, outpng, 'Resolution', 600, 'BackgroundColor','white');
close(fg);
end
