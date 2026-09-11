function breath_diff_image(folder, nFr)
% breath_diff_image  Show the camera image next to an inspiration-vs-baseline
% difference heatmap, computed from the cropped cube (breath_crop.mat) ranked by
% the breath signal (breath_pc1.mat). A DLC-free "what moves during a breath" map
% -- an interpretable alternative to the SVD eigenimage.
%
%   left  : mean camera frame (the image DLC would track on)
%   right : mean(top nFr highest-signal frames) - mean(bottom nFr lowest)
%           = inspiration - baseline; diverging heatmap, symmetric scale
%
%   breath_diff_image(folder)        % nFr = 100
%   breath_diff_image(folder, nFr)

if nargin < 2 || isempty(nFr), nFr = 100; end

C = load(fullfile(folder,'breath_crop.mat'));   mov = single(C.mov);
P = load(fullfile(folder,'breath_pc1.mat'));    bt  = double(P.breathTrace(:));
[H,W,T] = size(mov);
T2 = min(T, numel(bt)); mov = mov(:,:,1:T2); bt = bt(1:T2);

% orient so inspiration (the sharp, brief excursion) is the HIGH side
z = (bt-mean(bt))/std(bt); if mean(z.^3) < 0, bt = -bt; end
nFr = min(nFr, floor(T2/2));
[~, ord] = sort(bt, 'descend');
topIdx = ord(1:nFr);                    % peak-inspiration frames
botIdx = ord(end-nFr+1:end);            % baseline frames

camImg  = mean(mov, 3);
diffImg = mean(mov(:,:,topIdx), 3) - mean(mov(:,:,botIdx), 3);

% blue-white-red diverging colormap, symmetric about 0
m = 128;
cmap = [ [linspace(0,1,m)' linspace(0,1,m)' ones(m,1)]; ...
         [ones(m,1) linspace(1,0,m)' linspace(1,0,m)'] ];
cl = max(abs(diffImg(:))) + eps;

[~, leaf] = fileparts(folder);
f = figure('Color','w','Position',[80 120 1100 520]);
ax1 = subplot(1,2,1); imagesc(ax1, camImg); axis(ax1,'image','off'); colormap(ax1, gray);
colorbar(ax1); title(ax1, 'camera image (mean frame)','Interpreter','none');
ax2 = subplot(1,2,2); imagesc(ax2, diffImg, [-cl cl]); axis(ax2,'image','off'); colormap(ax2, cmap);
colorbar(ax2); title(ax2, sprintf('top %d - bottom %d frames  (inspiration - baseline)', nFr, nFr),'Interpreter','none');
sgtitle(strrep(leaf,'_','\_'));

exportgraphics(f, fullfile(folder,'breath_diff_image.png'), 'Resolution',130);
fprintf('Saved -> %s\n', fullfile(folder,'breath_diff_image.png'));
end
